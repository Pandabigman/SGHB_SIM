"""
Base model class and shared genetic calculations for SGHB simulation.
"""

import numpy as np
from abc import ABC, abstractmethod

from genetics import (
    LOCI,
    CaptiveBreedingParams,
    simulate_supplementation_effect_breeding,
    simulate_mixed_source_supplementation,
)
from utils.cache import get_cache_key, get_cached_simulation


# Biological constants
LETHAL_EQUIVALENTS_BIRDS = 3.14  # O'Grady et al. 2006
GENERATION_TIME_YEARS = 26


def calculate_heterozygosity_loss(H0, Ne, t):
    """Wright's equation: Ht = H0 * (1 - 1/(2*Ne))^t"""
    return H0 * np.power(1 - 1 / (2 * Ne), t)


def calculate_inbreeding(Ne, t):
    """Inbreeding coefficient: F = 1 - (1 - 1/(2*Ne))^t"""
    return 1 - np.power(1 - 1 / (2 * Ne), t)


def calculate_allelic_diversity(A0, Ne, t):
    """
    Allelic diversity loss via genetic drift only
    At = A0 * exp(-t/(4*Ne))

    Note: This formula ONLY models allele loss and does not account for:
    - Gene flow (immigration adding novel alleles)
    - Mutation (negligible over short timescales)

    For supplementation scenarios, use calculate_allelic_diversity_with_geneflow()
    """
    At = A0 * np.exp(-t / (4 * Ne))
    return np.maximum(At, 2.0)


def calculate_allelic_diversity_with_geneflow(A0, Ne, t, novel_alleles_per_gen=0,
                                              max_novel_alleles=None):
    """
    Allelic diversity accounting for both drift loss AND gene flow gain

    This model combines:
    1. Drift-driven loss of initial alleles: A0 * exp(-t/(4*Ne))
    2. Gain from gene flow: Sum of novel alleles added each generation,
       adjusted for subsequent drift loss, capped at a finite pool

    Args:
        A0: Initial mean alleles per locus
        Ne: Effective population size (can vary with gene flow)
        t: Time in generations (scalar or array)
        novel_alleles_per_gen: Mean novel alleles added per locus per generation
        max_novel_alleles: Hard ceiling on total novel alleles that can ever be
            introduced (the finite captive allele pool). Defaults to 2.7 — the
            empirical surplus of the full PAAZA+AZA+EAZA pool over wild:
            combined Na 9.07 − wild Na 6.43 = 2.64, rounded up to 2.7.
            Pass a higher value for international mixes with more distinct alleles.

    Returns:
        Expected allele count (scalar or array matching t)
    """
    if max_novel_alleles is None:
        max_novel_alleles = 2.7  # SGHB global captive pool surplus over wild

    # Handle both scalar and array inputs
    is_scalar = np.isscalar(t)
    t_array = np.atleast_1d(t)

    results = np.zeros_like(t_array, dtype=float)

    for idx, time_point in enumerate(t_array):
        # Component 1: Drift loss of original alleles
        At_drift = A0 * np.exp(-time_point / (4 * Ne))

        # Component 2: Gene flow gain (cumulative across generations)
        allele_gain = 0.0
        if novel_alleles_per_gen > 0 and time_point > 0:
            # For each past generation, calculate alleles added and retained
            for gen in range(1, int(time_point) + 1):
                # Alleles added at generation 'gen' experience drift for (time_point - gen) generations
                time_since_addition = time_point - gen
                # Retention probability follows drift model
                alleles_retained = novel_alleles_per_gen * np.exp(-time_since_addition / (4 * Ne))
                allele_gain += alleles_retained

        # Cap at the finite pool of novel alleles available in the captive source
        allele_gain = min(allele_gain, max_novel_alleles)

        # Total alleles = original (after drift) + novel (after drift)
        At_total = At_drift + allele_gain

        # Biological minimum: at least 2 alleles per locus (one per chromosome)
        results[idx] = max(At_total, 2.0)

    return results[0] if is_scalar else results


def calculate_population_size(N0, lambda_val, t):
    """Population size: Nt = N0 * lambda^t"""
    return N0 * np.power(lambda_val, t)


def calculate_population_size_with_inbreeding(
    N0, lambda_val, F_array, lethal_equivalents=LETHAL_EQUIVALENTS_BIRDS, stochastic=False,
    env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40
):
    """
    Population size with inbreeding depression

    Survival is reduced by: s = exp(-B * F)
    where B = lethal equivalents (typical range: 3-12 for mammals, 2-6 for birds)

    For Southern Ground Hornbill, using B = 3.14 (moderate inbreeding depression for birds)
    Based on O'Grady et al. 2006 review showing birds average ~3.14 lethal equivalents

    Args:
        N0: Initial population size
        lambda_val: Base growth rate (without inbreeding)
        F_array: Array of inbreeding coefficients over time
        lethal_equivalents: Number of lethal equivalents (default 3.14)
        stochastic: If True, add demographic and environmental stochasticity
        env_sigma: Environmental stochasticity std dev per generation (default 0.06).
            SGHB-tuned: cooperative breeding + high adult survival (~93-96%/yr) buffers
            generation-level variance (Morris & Doak 2002; Reed et al. 2002 report
            annual sigma_e ~0.05-0.10 for birds; generation-level ~0.05-0.08 for SGHB).
        catastrophe_prob: Per-generation probability of a catastrophic event (default 0.0).
            Models disease outbreaks, severe drought, etc. Set >0 to enable.
        catastrophe_magnitude: Fraction of population killed in a catastrophe (default 0.40).
            e.g. 0.40 = 40% of individuals die in the event.

    Returns: Array of population sizes incorporating inbreeding depression
    """
    N = np.zeros(len(F_array))
    N[0] = N0

    for t in range(1, len(F_array)):
        # Calculate survival reduction due to inbreeding from the previous generation
        inbreeding_survival = np.exp(-lethal_equivalents * F_array[t - 1])

        # Base adjusted growth rate
        lambda_adjusted = lambda_val * inbreeding_survival

        if stochastic:
            # ENVIRONMENTAL STOCHASTICITY: Random variation in growth rate
            # sigma=0.06 reflects SGHB's buffered life history; high-end would be ~0.10
            lambda_stochastic = lambda_adjusted + np.random.normal(0, env_sigma)
            # Ensure lambda stays positive
            lambda_stochastic = max(lambda_stochastic, 0.5)

            # DEMOGRAPHIC STOCHASTICITY: variance scaled by generation length (26 yrs)
            # A single per-generation draw underestimates cumulative annual variance by ~sqrt(26).
            # Normal(mu, sqrt(26*mu)) approximates a Galton-Watson process over 26 annual steps.
            expected_N = N[t - 1] * lambda_stochastic
            demo_sigma = np.sqrt(26 * max(expected_N, 0))
            N[t] = max(0, int(round(np.random.normal(expected_N, demo_sigma))))

            # CATASTROPHE: stochastic mass-mortality event (disease, drought, poisoning)
            if catastrophe_prob > 0 and np.random.random() < catastrophe_prob:
                N[t] = N[t] * (1 - catastrophe_magnitude)
        else:
            # Deterministic calculation
            N[t] = N[t - 1] * lambda_adjusted

        # Allow true extinction (N=0); do not rescue with artificial floor
        N[t] = max(N[t], 0)

    return N


def calculate_released_birds_optimized(birds_per_gen, F_array_rescued, captive_survival_multiplier=0.95, lethal_equivalents=LETHAL_EQUIVALENTS_BIRDS):
    """
    OPTIMIZED: Calculate released bird populations across all generations

    This vectorized approach replaces the nested loop O(n²) with O(n) complexity

    Args:
        birds_per_gen: Number of birds released per generation
        F_array_rescued: Array of inbreeding coefficients (with genetic rescue)
        captive_survival_multiplier: Base survival rate for captive birds
        lethal_equivalents: Lethal equivalents for inbreeding depression

    Returns: Array of released bird counts per generation
    """
    n_gens = len(F_array_rescued)
    N_released = np.zeros(n_gens)

    # For each release cohort (starts at 1 — no birds released before gen 1)
    for release_gen in range(1, n_gens):
        # Calculate how long each cohort survives (vectorized)
        future_gens = np.arange(release_gen, n_gens)
        time_since_release = future_gens - release_gen

        # Calculate average F for each time point
        F_cumsum = np.cumsum(F_array_rescued)
        avg_F_values = np.zeros(len(future_gens))

        for idx, gen in enumerate(future_gens):
            if gen == release_gen:
                avg_F_values[idx] = F_array_rescued[release_gen]
            else:
                avg_F_values[idx] = (F_cumsum[gen] - F_cumsum[release_gen - 1]) / (gen - release_gen + 1)

        # Vectorized survival calculation
        # Released birds experience 15% of normal inbreeding depression (hybrid vigor)
        survival_rate = np.exp(-lethal_equivalents * avg_F_values * 0.15) * captive_survival_multiplier

        # Calculate survivors at each time point
        cohort_survivors = birds_per_gen * (survival_rate ** time_since_release)

        # Add to total released population
        N_released[future_gens] += cohort_survivors

    return N_released


def calculate_genetic_rescue_effect(F_array, birds_per_gen, N0):
    """
    Apply genetic rescue effect to inbreeding coefficients.

    Gene flow from supplementation reduces effective inbreeding.

    Args:
        F_array: Original inbreeding coefficient array
        birds_per_gen: Number of birds added per generation
        N0: Base population size

    Returns: Rescued F array with reduced inbreeding
    """
    F_array_rescued = F_array.copy()
    for gen in range(len(F_array)):
        if gen > 0:
            cumulative_birds = birds_per_gen * gen
            total_pop_estimate = N0 + cumulative_birds
            gene_flow_proportion = cumulative_birds / total_pop_estimate
            F_array_rescued[gen] = F_array[gen] * (1 - gene_flow_proportion * 0.5)
    return F_array_rescued


class BaseModel(ABC):
    """Abstract base class for all population models."""

    model_number: int = 0
    model_name: str = "Base Model"

    @abstractmethod
    def run(self, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None,
            env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40):
        """
        Run the simulation model.

        Args:
            Ne: Effective population size
            generations: Number of generations to simulate
            lambda_val: Population growth rate
            stochastic: Whether to include stochasticity
            genetic_data: Loaded genetic data (optional)
            env_sigma: Environmental stochasticity std dev per generation
            catastrophe_prob: Per-generation probability of catastrophic event
            catastrophe_magnitude: Fraction of population killed in catastrophe

        Returns: Dict with simulation results
        """
        pass

    def _build_result(self, t, Ho_vals, He_vals, F_vals, Na_vals, N_vals, Ne,
                      initial_values, parameters, **extra_data):
        """Build standardized result dictionary."""
        result = {
            "model_number": self.model_number,
            "model_name": self.model_name,
            "generations": t.tolist() if hasattr(t, 'tolist') else list(t),
            "years": (t * GENERATION_TIME_YEARS).tolist() if hasattr(t, 'tolist') else [x * GENERATION_TIME_YEARS for x in t],
            "Ho": Ho_vals.tolist() if hasattr(Ho_vals, 'tolist') else list(Ho_vals),
            "He": He_vals.tolist() if hasattr(He_vals, 'tolist') else list(He_vals),
            "F": F_vals.tolist() if hasattr(F_vals, 'tolist') else list(F_vals),
            "Na": Na_vals.tolist() if hasattr(Na_vals, 'tolist') else list(Na_vals),
            "population": N_vals.tolist() if hasattr(N_vals, 'tolist') else list(N_vals),
            "Ne": Ne,
            "initial": initial_values,
            "parameters": parameters,
        }
        # Add any extra data (like effective_Ne, captive_F_mean, etc.)
        result.update(extra_data)
        return result
