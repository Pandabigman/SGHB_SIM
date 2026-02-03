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
)


class LowSupplementationModel(BaseModel):
    """Model 3: Low supplementation with 4 SA captive birds per generation."""

    model_number = 3
    model_name = "Low Supplementation (+4 SA Captive/gen)"

    def run(self, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None):
        """
        Run low supplementation model with dynamic captive breeding.

        Uses Mendelian inheritance to create new genetic combinations each generation.
        """
        birds_per_gen = 4

        if genetic_data:
            return self._run_with_csv_data(Ne, generations, lambda_val, stochastic, genetic_data, birds_per_gen)
        else:
            return self._run_generic(Ne, generations, lambda_val, stochastic, birds_per_gen)

    def _run_with_csv_data(self, Ne, generations, lambda_val, stochastic, genetic_data, birds_per_gen):
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

        # Use breeding simulation with caching
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

        # Extract data
        Ho = [r["Ho"] for r in results]
        He = [r["He"] for r in results]
        Na = [r["Na"] for r in results]
        F = [r["FIS"] for r in results]
        effective_Ne_vals = [r["effective_Ne"] for r in results]
        captive_F_mean = [r.get("captive_F_mean", 0) for r in results]

        # Population dynamics
        N0 = 2500
        captive_survival_multiplier = 0.95

        t = np.arange(0, generations + 1)
        F_array = np.zeros(len(t))
        for gen in range(len(t)):
            Ne_current = effective_Ne_vals[gen]
            F_array[gen] = 1 - np.power(1 - 1 / (2 * Ne_current), gen)

        # Apply genetic rescue effect
        F_array_rescued = calculate_genetic_rescue_effect(F_array, birds_per_gen, N0)

        # Calculate population sizes
        N_wild = calculate_population_size_with_inbreeding(N0, lambda_val, F_array_rescued, stochastic=stochastic)
        N_released = calculate_released_birds_optimized(birds_per_gen, F_array_rescued, captive_survival_multiplier)
        N = N_wild + N_released

        return self._build_result(
            t=t,
            Ho_vals=Ho,
            He_vals=He,
            F_vals=F,
            Na_vals=Na,
            N_vals=N,
            Ne=Ne,
            initial_values={"Ho": Ho[0], "He": He[0], "Na": Na[0], "N": N[0]},
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

    def _run_generic(self, Ne, generations, lambda_val, stochastic, birds_per_gen):
        """Run generic model when CSV data not available."""
        H0 = 0.502
        He0 = 0.568
        A0 = 6.429
        N0 = 2500

        t = np.arange(0, generations + 1)

        # Estimate novel alleles
        novel_alleles_per_gen = birds_per_gen * 0.15

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
            Na = calculate_allelic_diversity_with_geneflow(A0, effective_Ne, gen, novel_alleles_per_gen)

            Ho_vals.append(Ho)
            He_vals.append(He)
            F_vals.append(F)
            Na_vals.append(Na)

        F_array = np.array(F_vals)
        F_array_rescued = calculate_genetic_rescue_effect(F_array, birds_per_gen, N0)

        captive_survival_multiplier = 0.95
        N_wild = calculate_population_size_with_inbreeding(N0, lambda_val, F_array_rescued, stochastic=stochastic)
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
