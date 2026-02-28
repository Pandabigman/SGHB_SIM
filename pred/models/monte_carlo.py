"""Monte Carlo wrapper for SGHB population simulation.

Runs a model N times with stochastic=True to produce
statistical distributions over outcomes, including extinction
probability and confidence intervals on all genetic metrics.

Two modes:
  run_monte_carlo()  - blocking, returns full results dict (used by old POST endpoint)
  iter_monte_carlo() - generator, yields partial results for SSE streaming
"""

import numpy as np
from . import run_model

_METRICS = ['Ho', 'He', 'Na', 'F', 'population']


def _run_single_replicate(args):
    """
    Top-level worker function for ProcessPoolExecutor.
    Must be module-level (not a lambda or closure) so it can be pickled.

    Uses the analytical path (genetic_data=None) to stay fast on serverless
    environments — stochastic=True still randomises population dynamics.

    args: (model_num, Ne, generations, lambda_val, env_sigma,
           catastrophe_prob, catastrophe_magnitude, max_novel_alleles)
    max_novel_alleles is pre-computed from CSV data by the caller so that
    the generic analytical Na formula uses the correct per-model ceiling
    rather than the hardcoded default.
    """
    model_num, Ne, generations, lambda_val, env_sigma, catastrophe_prob, catastrophe_magnitude, max_novel_alleles = args
    return run_model(model_num, Ne, generations, lambda_val, stochastic=True, genetic_data=None,
                     env_sigma=env_sigma, catastrophe_prob=catastrophe_prob,
                     catastrophe_magnitude=catastrophe_magnitude,
                     max_novel_alleles=max_novel_alleles)


def _pad(arr, n_points):
    arr = np.array(arr, dtype=float)[:n_points]
    if len(arr) < n_points:
        arr = np.pad(arr, (0, n_points - len(arr)), 'edge')
    return arr


def _compute_summary(all_runs, first_result, completed, total, extinction_threshold):
    """Compute aggregate stats from the accumulated replicate lists."""
    pct_levels = [5, 25, 50, 75, 95]
    summary_metrics = {}
    for m in _METRICS:
        data = np.array(all_runs[m])          # (completed, n_points)
        pcts = np.percentile(data, pct_levels, axis=0)
        summary_metrics[m] = {
            'mean': np.mean(data, axis=0).tolist(),
            'std':  np.std(data,  axis=0).tolist(),
            'p5':   pcts[0].tolist(),
            'p25':  pcts[1].tolist(),
            'p50':  pcts[2].tolist(),
            'p75':  pcts[3].tolist(),
            'p95':  pcts[4].tolist(),
        }

    pop_data = np.array(all_runs['population'])
    extinction_prob = (pop_data < extinction_threshold).mean(axis=0).tolist()

    return {
        'type': 'complete' if completed == total else 'progress',
        'completed': completed,
        'total': total,
        'model_number': first_result.get('model_number'),
        'model_name':   first_result.get('model_name', ''),
        'n_replicates': completed,
        'generations':  first_result.get('generations', []),
        'years':        first_result.get('years', []),
        'extinction_threshold':    extinction_threshold,
        'metrics':                 summary_metrics,
        'extinction_probability':  extinction_prob,
        'extinction_probability_final': extinction_prob[-1],
        'initial':    first_result.get('initial', {}),
        'parameters': first_result.get('parameters', {}),
    }


def iter_monte_carlo(model_num, n_replicates, Ne, generations, lambda_val,
                     extinction_threshold=100, batch_size=5,
                     env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40,
                     max_novel_alleles=None):
    """
    Generator that yields partial Monte Carlo results as replicates complete.

    Tries ProcessPoolExecutor(max_workers=2) for speed; falls back to sequential
    if the environment doesn't allow subprocess spawning (e.g. Vercel sandbox).

    Each yielded dict has the same shape as the final run_monte_carlo() result,
    with an added 'completed' / 'total' / 'type' key for progress tracking.

    max_novel_alleles: per-model novel allele ceiling pre-computed from CSV data by
        the API layer (api/index.py). Passed through to _run_generic() so that Na
        uses the correct empirical cap instead of the hardcoded 2.7 default.
    """
    n_points = generations + 1
    all_runs = {m: [] for m in _METRICS}
    first_result = None
    args = (model_num, Ne, generations, lambda_val, env_sigma, catastrophe_prob,
            catastrophe_magnitude, max_novel_alleles)

    def _accumulate(result):
        nonlocal first_result
        if first_result is None:
            first_result = result
        for m in _METRICS:
            all_runs[m].append(_pad(result.get(m, []), n_points))

    def _maybe_yield(completed):
        if completed % batch_size == 0 or completed == n_replicates:
            yield _compute_summary(all_runs, first_result, completed, n_replicates, extinction_threshold)

    # --- attempt parallel execution ---
    try:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        with ProcessPoolExecutor(max_workers=2) as executor:
            futures = [executor.submit(_run_single_replicate, args) for _ in range(n_replicates)]
            completed = 0
            for future in as_completed(futures):
                _accumulate(future.result())
                completed += 1
                yield from _maybe_yield(completed)
        return  # done
    except Exception:
        # Reset and fall through to sequential
        all_runs = {m: [] for m in _METRICS}
        first_result = None

    # --- sequential fallback ---
    for i in range(n_replicates):
        _accumulate(_run_single_replicate(args))
        yield from _maybe_yield(i + 1)


def run_monte_carlo(model_num, n_replicates, Ne, generations, lambda_val,
                    genetic_data=None, extinction_threshold=100,
                    env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40,
                    max_novel_alleles=None):
    """
    Blocking Monte Carlo run. Returns full summary dict.

    NOTE: genetic_data is accepted for API compatibility but ignored —
    the streaming path always uses the analytical model for speed.
    For a fully empirical MC run use iter_monte_carlo() via the SSE endpoint.

    max_novel_alleles: per-model novel allele ceiling pre-computed from CSV data.
    """
    n_points = generations + 1
    all_runs = {m: np.zeros((n_replicates, n_points)) for m in _METRICS}
    first_result = None

    for i in range(n_replicates):
        result = _run_single_replicate((model_num, Ne, generations, lambda_val,
                                        env_sigma, catastrophe_prob, catastrophe_magnitude,
                                        max_novel_alleles))
        if first_result is None:
            first_result = result
        for m in _METRICS:
            all_runs[m][i] = _pad(result.get(m, []), n_points)

    pct_levels = [5, 25, 50, 75, 95]
    summary_metrics = {}
    for m in _METRICS:
        data = all_runs[m]
        pcts = np.percentile(data, pct_levels, axis=0)
        summary_metrics[m] = {
            'mean': np.mean(data, axis=0).tolist(),
            'std':  np.std(data,  axis=0).tolist(),
            'p5':   pcts[0].tolist(),
            'p25':  pcts[1].tolist(),
            'p50':  pcts[2].tolist(),
            'p75':  pcts[3].tolist(),
            'p95':  pcts[4].tolist(),
        }

    pop_data = all_runs['population']
    extinction_prob = (pop_data < extinction_threshold).mean(axis=0).tolist()

    return {
        'model_number': model_num,
        'model_name':   first_result.get('model_name', f'Model {model_num}'),
        'n_replicates': n_replicates,
        'generations':  first_result.get('generations', list(range(n_points))),
        'years':        first_result.get('years', [g * 26 for g in range(n_points)]),
        'extinction_threshold':    extinction_threshold,
        'metrics':                 summary_metrics,
        'extinction_probability':  extinction_prob,
        'extinction_probability_final': extinction_prob[-1],
        'initial':    first_result.get('initial', {}),
        'parameters': first_result.get('parameters', {}),
    }
