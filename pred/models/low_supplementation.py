"""
Model 3: Low Supplementation - Add 4 SA Captive birds per generation.
"""

import numpy as np

from genetics import LOCI, CaptiveBreedingParams, simulate_supplementation_effect_breeding
from utils.cache import get_cache_key, get_cached_simulation

from .base import (
    BaseModel,
    calculate_population_size_with_inbreeding,
    calculate_released_birds_optimized,
    calculate_genetic_rescue_effect,
    calculate_heterozygosity_loss,
    calculate_inbreeding,
    calculate_allelic_diversity_with_geneflow,
    compute_max_novel_alleles,
)


class LowSupplementationModel(BaseModel):
    """Model 3: Low supplementation with 4 SA captive birds per generation."""

    model_number = 3
    model_name = "Low Supplementation (+4 SA Captive/gen)"

    def run(self, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None,
            env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40,
            max_novel_alleles=None):
        """
        Run low supplementation model with dynamic captive breeding.

        Uses Mendelian inheritance to create new genetic combinations each generation.
        """
        birds_per_gen = 4

        if genetic_data:
            return self._run_with_csv_data(Ne, generations, lambda_val, stochastic, genetic_data, birds_per_gen,
                                           env_sigma, catastrophe_prob, catastrophe_magnitude)
        else:
            return self._run_generic(Ne, generations, lambda_val, stochastic, birds_per_gen,
                                     env_sigma, catastrophe_prob, catastrophe_magnitude,
                                     max_novel_alleles=max_novel_alleles)

    def _run_with_csv_data(self, Ne, generations, lambda_val, stochastic, genetic_data, birds_per_gen,
                           env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40):
        """Run model using real CSV data with dynamic breeding simulation."""
        wild_df = genetic_data["dataframes"]["wild_all"]
        paaza_df = genetic_data["dataframes"]["paaza"]

        # Configure captive breeding parameters
        captive_params = CaptiveBreedingParams(
            target_population_size=len(paaza_df),
            captive_Ne=35.0,
            breeding_success_rate=0.7,
            offspring_per_pair=1.5
        )

        seed = None if stochastic else 42

        # Use breeding simulation with caching (provides Ho/He/F in stochastic mode)
        cache_key = get_cache_key(3, birds_per_gen, generations, Ne, stochastic) + "_breeding"
        results = get_cached_simulation(
            cache_key,
            simulate_supplementation_effect_breeding,
            stochastic,
            wild_df, paaza_df,
            birds_per_gen=birds_per_gen,
            generations=generations,
            loci_list=LOCI,
            base_Ne=Ne,
            captive_params=captive_params,
            seed=seed
        )

        effective_Ne_vals = [r["effective_Ne"] for r in results]
        captive_F_mean = [r.get("captive_F_mean", 0) for r in results]

        # Population dynamics
        N0 = 2500
        captive_survival_multiplier = 0.95
        t = np.arange(0, generations + 1)

        # --- Genetic metrics: analytical to correct Ne mismatch ---
        # The breeding simulation uses wild_pop_size=len(wild_df)=199 as the carrying
        # capacity, so Wright-Fisher drift operates at effective Ne≈199 rather than
        # the user-supplied Ne (default 500). This would make conservation models
        # appear genetically *worse* than the baseline, which is biologically wrong.
        # Fix: compute all genetic metrics analytically at the intended Ne, calibrated
        # from the actual CSV initial values and empirical novel-allele counts.
        A0  = genetic_data["summaries"]["wild_all"]["Na"]
        H0  = genetic_data["summaries"]["wild_all"]["Ho"]
        He0 = genetic_data["summaries"]["wild_all"]["He"]

        max_novel = compute_max_novel_alleles(genetic_data["novel_alleles"], ["paaza"])
        # Rate: expected novel alleles introduced per locus per generation
        # = probability a random released bird carries a novel allele × birds_per_gen
        paaza_size = len(paaza_df)
        novel_per_gen = birds_per_gen * max_novel / paaza_size

        # Na — always analytical (prevents Ne≈199 sample-size bias in all modes)
        Na_vals = calculate_allelic_diversity_with_geneflow(A0, Ne, t, novel_per_gen, max_novel)

        if not stochastic:
            # Deterministic: analytical Ho/He/F with gene-flow-boosted effective Ne
            Ho_vals, He_vals, F_vals = [], [], []
            for gen in t:
                eff_Ne = Ne + birds_per_gen * gen * 0.5
                Ho_vals.append(calculate_heterozygosity_loss(H0, eff_Ne, gen))
                He_vals.append(calculate_heterozygosity_loss(He0, eff_Ne, gen))
                F_vals.append(calculate_inbreeding(eff_Ne, gen))
            F_array = np.array(F_vals)
            fis_baseline = F_array[0]
            F_array = np.maximum(F_array - fis_baseline, 0.0)
            F_array_rescued = F_array
        else:
            # Stochastic single-run: simulation Ho/He/F (preserves genetic variance)
            Ho_vals = [r["Ho"] for r in results]
            He_vals = [r["He"] for r in results]
            F_vals  = [r["FIS"] for r in results]
            fis_baseline = F_vals[0]
            F_array = np.maximum(np.array(F_vals) - fis_baseline, 0.0)
            F_array_rescued = F_array

        # Calculate population sizes
        N_wild = calculate_population_size_with_inbreeding(
            N0, lambda_val, F_array_rescued, stochastic=stochastic,
            env_sigma=env_sigma, catastrophe_prob=catastrophe_prob, catastrophe_magnitude=catastrophe_magnitude
        )
        N_released = calculate_released_birds_optimized(birds_per_gen, F_array_rescued, captive_survival_multiplier)
        N = N_wild + N_released

        return self._build_result(
            t=t,
            Ho_vals=Ho_vals,
            He_vals=He_vals,
            F_vals=F_vals,
            Na_vals=Na_vals,
            N_vals=N,
            Ne=Ne,
            initial_values={"Ho": Ho_vals[0], "He": He_vals[0], "Na": float(Na_vals[0]), "N": N[0]},
            parameters={
                "Ne_base": Ne,
                "Ne_final": effective_Ne_vals[-1],
                "N0": int(N[0]),
                "lambda": lambda_val,
                "data_source": "CSV_breeding_simulation",
                "supplementation": f"{birds_per_gen} South African captive birds per generation",
                "supplementation_source": "SA Captive (South African Zoos and Aquaria)",
                "breeding_model": "dynamic (Mendelian inheritance, new combinations each gen)",
                "novel_alleles_added": len(genetic_data["novel_alleles"]["paaza"]),
                "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)",
            },
            effective_Ne=effective_Ne_vals,
            captive_F_mean=captive_F_mean,
        )

    def _run_generic(self, Ne, generations, lambda_val, stochastic, birds_per_gen,
                     env_sigma=0.06, catastrophe_prob=0.0, catastrophe_magnitude=0.40,
                     max_novel_alleles=None):
        """Run generic model when CSV data not available."""
        H0 = 0.502
        He0 = 0.568
        A0 = 6.429
        N0 = 2500

        t = np.arange(0, generations + 1)

        # Per-model novel allele ceiling and captive pool size.
        # When max_novel_alleles is provided (e.g. from Monte Carlo caller), use it;
        # otherwise fall back to a conservative PAAZA-specific default (1.5 vs the
        # combined-pool value of 2.7 that previously inflated Na to ~9).
        if max_novel_alleles is None:
            max_novel_alleles = 1.5  # PAAZA-only surplus over wild (conservative)
        captive_pool_size = 90       # Approximate PAAZA population size
        novel_alleles_per_gen = birds_per_gen * max_novel_alleles / captive_pool_size

        Ho_vals = []
        He_vals = []
        F_vals = []
        Na_vals = []
        effective_Ne_vals = []

        for gen in t:
            birds_added = birds_per_gen * gen
            effective_Ne = Ne + (birds_added * 0.5)
            effective_Ne_vals.append(effective_Ne)

            Ho = calculate_heterozygosity_loss(H0, effective_Ne, gen)
            He = calculate_heterozygosity_loss(He0, effective_Ne, gen)
            F = 1 - np.power(1 - 1 / (2 * effective_Ne), gen)
            Na = calculate_allelic_diversity_with_geneflow(
                A0, effective_Ne, gen, novel_alleles_per_gen, max_novel_alleles
            )

            Ho_vals.append(Ho)
            He_vals.append(He)
            F_vals.append(F)
            Na_vals.append(Na)

        F_array = np.array(F_vals)
        F_array_rescued = calculate_genetic_rescue_effect(F_array, birds_per_gen, N0)

        captive_survival_multiplier = 0.95
        N_wild = calculate_population_size_with_inbreeding(
            N0, lambda_val, F_array_rescued, stochastic=stochastic,
            env_sigma=env_sigma, catastrophe_prob=catastrophe_prob, catastrophe_magnitude=catastrophe_magnitude
        )
        N_released = calculate_released_birds_optimized(birds_per_gen, F_array_rescued, captive_survival_multiplier)
        N_vals = N_wild + N_released

        return self._build_result(
            t=t,
            Ho_vals=Ho_vals,
            He_vals=He_vals,
            F_vals=F_vals,
            Na_vals=Na_vals,
            N_vals=N_vals,
            Ne=Ne,
            initial_values={"Ho": H0, "He": He0, "Na": A0, "N": N0},
            parameters={
                "Ne_base": Ne,
                "Ne_final": effective_Ne_vals[-1],
                "N0": N0,
                "lambda": lambda_val,
                "data_source": "generic",
                "supplementation": f"{birds_per_gen} SA Captive birds per generation",
                "novel_alleles_per_gen": novel_alleles_per_gen,
                "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)",
            },
            effective_Ne=effective_Ne_vals,
        )
