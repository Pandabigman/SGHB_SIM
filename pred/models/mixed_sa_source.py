"""
Model 6: Mixed South African Sourcing with wild population depletion.
"""

import numpy as np

from genetics import LOCI, CaptiveBreedingParams, simulate_mixed_source_supplementation
from utils.cache import get_cache_key, get_cached_simulation

from .base import (
    BaseModel,
    calculate_population_size_with_inbreeding,
    calculate_released_birds_optimized,
    calculate_genetic_rescue_effect,
    calculate_heterozygosity_loss,
    calculate_allelic_diversity_with_geneflow,
)


class MixedSASourceModel(BaseModel):
    """Model 6: Mixed SA sourcing with captive + wild translocation."""

    model_number = 6
    model_name = "Mixed SA Source (+6 Captive +2 KZN +2 EC/gen)"

    def run(self, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None):
        """
        Run mixed SA source model with wild population depletion.

        Sources per generation:
        - 6 birds from SA Captive (breeding PAAZA population)
        - 2 birds from KwaZulu-Natal wild (translocated, depletes source)
        - 2 birds from Eastern Cape wild (translocated, depletes source)
        """
        if genetic_data:
            return self._run_with_csv_data(Ne, generations, lambda_val, stochastic, genetic_data)
        else:
            return self._run_generic(Ne, generations, lambda_val, stochastic)

    def _run_with_csv_data(self, Ne, generations, lambda_val, stochastic, genetic_data):
        """Run model using real CSV data with wild depletion."""
        wild_df = genetic_data["dataframes"]["wild_all"]
        paaza_df = genetic_data["dataframes"]["paaza"]
        kzn_df = genetic_data["dataframes"]["wild_kzn"]
        ec_df = genetic_data["dataframes"]["wild_ec"]

        # Configure captive breeding parameters
        captive_params = CaptiveBreedingParams(
            target_population_size=len(paaza_df),
            captive_Ne=35.0,
            breeding_success_rate=0.7,
            offspring_per_pair=1.5
        )

        seed = None if stochastic else 42
        birds_per_gen = 10  # 6 captive + 2 KZN + 2 EC

        # Use mixed source simulation with caching
        cache_key = get_cache_key(6, birds_per_gen, generations, Ne, stochastic) + "_mixed_sa"
        results = get_cached_simulation(
            cache_key,
            simulate_mixed_source_supplementation,
            stochastic,
            wild_df, paaza_df, kzn_df, ec_df,
            captive_birds_per_gen=6,
            kzn_birds_per_gen=2,
            ec_birds_per_gen=2,
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
        kzn_remaining = [r.get("kzn_source_remaining", 0) for r in results]
        ec_remaining = [r.get("ec_source_remaining", 0) for r in results]

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
                "data_source": "CSV_mixed_source_simulation",
                "supplementation": "10 birds/gen (6 SA captive + 2 KZN wild + 2 EC wild)",
                "supplementation_sources": {
                    "captive": "6 from PAAZA (breeding population)",
                    "kzn_wild": "2 from KwaZulu-Natal (translocated, depletes source)",
                    "ec_wild": "2 from Eastern Cape (translocated, depletes source)"
                },
                "breeding_model": "dynamic captive + wild translocation with depletion",
                "source_depletion": "enabled (wild sources decrease over time)",
                "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)",
            },
            effective_Ne=effective_Ne_vals,
            captive_F_mean=captive_F_mean,
            kzn_source_remaining=kzn_remaining,
            ec_source_remaining=ec_remaining,
        )

    def _run_generic(self, Ne, generations, lambda_val, stochastic):
        """Run generic model when CSV data not available."""
        H0 = 0.502
        He0 = 0.568
        A0 = 6.429
        N0 = 2500
        birds_per_gen = 10

        t = np.arange(0, generations + 1)

        # Higher novel alleles for mixed SA sources
        novel_alleles_per_gen = birds_per_gen * 0.20

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
                "supplementation": f"{birds_per_gen} Mixed SA birds per generation",
                "novel_alleles_per_gen": novel_alleles_per_gen,
                "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)",
            },
            effective_Ne=effective_Ne_vals,
        )
