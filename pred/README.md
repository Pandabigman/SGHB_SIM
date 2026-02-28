# Southern Ground Hornbill Population Viability Simulator

A population genetics and viability analysis (PVA) tool for Southern Ground Hornbills, using real microsatellite data from 199 wild and 145 captive individuals across 14 loci. Compares 6 conservation supplementation strategies over 260–2600 years (10–100 generations).

---

## Project Structure

```
pred/
├── api/index.py          # Flask app — endpoints, data loading, model dispatch
├── models/               # base.py + model1-6.py + monte_carlo.py + __init__.py
├── genetics/             # metrics.py, population.py, supplementation.py,
│                         #   simulation.py (Mendelian breeding), parsing.py
├── data/
│   ├── blue.csv          # 199 wild birds (EC, Kruger, KZN, Limpopo) — 14 loci × 2 alleles
│   └── red.csv           # 145 captive birds (PAAZA, AZA, EAZA) — same format
├── templates/index.html  # Single-page UI
└── static/js/app.js      # Plotly charts + SSE streaming client
```

---

## Quick Start

```bash
cd pred
.venv/bin/python api/index.py
# → http://localhost:5001
```

---

## The 6 Conservation Scenarios

| # | Name | Supplementation | Birds/gen | Source pool |
|---|------|-----------------|-----------|-------------|
| 1 | Baseline | None | 0 | All 4 wild sites |
| 2 | Population Loss | None | 0 | Kruger + Limpopo only (EC + KZN extinct) |
| 3 | Low Supplementation | SA captive | 4 | PAAZA (~90 birds) |
| 4 | High Supplementation | SA captive | 10 | PAAZA (~90 birds) |
| 5 | International Mix | Global captive | 4 | PAAZA + AZA + EAZA (~145 birds) |
| 6 | Mixed SA Source | Multi-source | 10 (6+2+2) | 6 PAAZA captive + 2 KZN wild + 2 EC wild |

Model 2 loses the alleles unique to EC and KZN populations — these are identified from the real CSV data and permanently removed. Model 6 includes independent Mendelian breeding populations for each wild source (KZN and EC), not just CSV sampling.

---

## Genetic Calculation Approach

### The problem this branch solves

Earlier code ran Wright-Fisher simulation at Ne ≈ 199 (the CSV sample size, not the true population Ne) for all genetic metrics in supplementation models. Drift at Ne=199 is far stronger than drift at the true Ne=500, so conservation scenarios appeared *genetically worse* than the no-intervention baseline — a simulation artefact, not biology.

### The fix: analytical formulas at corrected Ne

**Deterministic mode** — Models 3–6 use analytical formulas throughout, evaluated at Ne corrected for gene flow:

| Metric | Formula | Notes |
|--------|---------|-------|
| Ho, He | `Ht = H0 × (1 − 1/2Ne)^t` | Ne boosted by gene flow: `Ne_eff = Ne + 0.5 × cumulative_birds` |
| F | `Ft = 1 − (1 − 1/2Ne)^t` | Same boosted Ne |
| Na | drift loss + cumulative gene flow gain | See below |

**Stochastic mode** — when real CSV data is loaded, Ho/He/F are extracted from a Wright-Fisher breeding simulation run at the corrected Ne (not the sample size). Na remains analytical in all modes.

**Without CSV data** — all metrics analytical regardless of stochastic/deterministic setting.

### Na with gene flow

```
At = A0 × exp(−t / 4Ne)                   # drift loss of original alleles
   + Σ(g=1..t) [novel/gen × exp(−(t−g) / 4Ne)]  # retained novel alleles, capped
```

Novel allele introduction rate = `birds_per_gen × max_novel / captive_pool_size`, where `max_novel` is computed from actual CSV data via `compute_max_novel_alleles()` in [models/base.py](models/base.py):

- PAAZA only (~90 birds): **1.29 novel alleles/locus**
- Full international pool (~145 birds): **2.64 novel alleles/locus**

The cumulative sum is capped at `max_novel` — the finite captive allele pool prevents unbounded Na growth.

### Key demographic formulas

```
# Inbreeding depression on population growth
N(t) = N(t−1) × λ × exp(−3.14 × F(t))

# Released captive bird survival (per generation since release)
survival = exp(−3.14 × F × 0.15) × 0.95
# — 15% hybrid vigor (heterosis), 95% captive survival rate vs wild
```

B = 3.14 lethal equivalents is the bird-specific value (O'Grady et al. 2006).

---

## Simulation Modes

### Deterministic (default)
Pure analytical computation, random seed = 42. Fast, reproducible, no variance.

### Stochastic
Adds randomness each generation:

- **Environmental stochasticity**: `λ_adj += Normal(0, σ_env)`, default σ = 0.06 (tuned to SGHB cooperative breeding biology)
- **Demographic stochasticity**: `N(t) ~ Normal(μ=expected_N, σ=sqrt(26 × expected_N))` — the sqrt(26) factor accounts for cumulative annual variance within a 26-year generation

### Catastrophe model (optional, within stochastic mode)
Each generation: `if rand() < catastrophe_prob → N × (1 − catastrophe_magnitude)`
Models disease outbreaks, drought, etc. Off by default.

### Monte Carlo
Runs N stochastic replicates and aggregates results. Uses the analytical path (no breeding simulation per replicate) for speed. Returns per-generation mean, std, and percentiles (p5/p25/p75/p95) for all metrics, plus **quasi-extinction probability** P(N < threshold) per generation.

Parallelised with `ProcessPoolExecutor(max_workers=2)`, falling back to sequential when subprocesses are unavailable (e.g., Vercel serverless).

---

## API Reference

### `POST /api/simulate`

```json
// Request
{
  "model": 1,                    // 1–6
  "Ne": 500,                     // 375–625
  "generations": 50,             // 10–100
  "lambda": 1.0,                 // 0.5–1.5
  "stochastic": false,
  "env_sigma": 0.06,             // Environmental noise σ (0.01–0.30)
  "catastrophe_prob": 0.0,       // Per-generation catastrophe probability
  "catastrophe_magnitude": 0.40  // Fraction of population killed if catastrophe occurs
}

// Response (per-generation arrays)
{ "generations", "years", "Ho", "He", "F", "Na", "population",
  "initial": {Ho, He, Na, N}, "parameters": {...} }
```

### `POST /api/monte-carlo`

```json
// Request (same base params as /simulate, plus):
{
  "n_replicates": 100,           // 10–500
  "extinction_threshold": 100    // N below which = quasi-extinct
}

// Response
{
  "metrics": {
    "Ho": { "mean", "std", "p5", "p25", "p50", "p75", "p95" },
    // He, F, Na, population — same structure
  },
  "extinction_probability": [...],       // per generation
  "extinction_probability_final": 0.25
}
```

### `GET /api/monte-carlo/stream`
Same query params as the POST endpoint. Returns an SSE stream with `{"type":"partial", "completed":N, "total":M, "metrics":{...}}` events during computation, then a final `{"type":"complete", ...}` event. Used by the frontend for the real-time progress bar.

### `GET /api/data/info`
Returns loaded genetic data statistics: sample sizes, Ho/He/Na/FIS per population, lost allele count (Model 2), novel allele counts per captive source.

---

## Visualisations

All charts are Plotly.js 2.27.0 interactive line plots.

**Standard mode (4 charts):** Heterozygosity (Ho) · Inbreeding coefficient (F) · Allelic richness (Na) · Population size (N)

**Monte Carlo mode:** Same 4 charts with shaded confidence bands (50% and 90% CIs around the median) + a 5th **quasi-extinction probability** chart (P(N < threshold) per generation, fills red).

---

## Data Sources

- **Wild** (`blue.csv`): 199 individuals — Eastern Cape (48), Kruger NP (52), KwaZulu-Natal (51), Limpopo (48)
- **Captive** (`red.csv`): 145 individuals — PAAZA/SA zoos (~90), AZA/USA-Canada (~28), EAZA/Europe (~27)
- **Loci (14):** Buco2/4/9/11/16/18/19/21/24/25, GHB19/20/21/26

---

## References

- O'Grady et al. (2006). *Biological Conservation* 133:42–51 — lethal equivalents B = 3.14 for birds
- Wright (1931). *Genetics* 16:97–159 — heterozygosity loss under drift
- Mills & Allendorf (1996). *Conservation Biology* 10:1509–1518 — one-migrant-per-generation rule
