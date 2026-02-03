"""
Model 2: Population Loss - Only Kruger + Limpopo remain.
"""

import numpy as np

from .base import (
    BaseModel,
    calculate_heterozygosity_loss,
    calculate_inbreeding,
    calculate_allelic_diversity,
    calculate_population_size_with_inbreeding,
)


class PopulationLossModel(BaseModel):
    """Model 2: Scenario where Eastern Cape and KwaZulu-Natal populations are lost."""

    model_number = 2
    model_name = "Population Loss (Kruger + Limpopo Only)"

    def run(self, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None):
        """
        Run population loss model showing impact of losing EC and KZN populations.

        Uses real CSV data to show actual allele loss.
        """
        if genetic_data:
            summary = genetic_data["summaries"]["wild_no_ec_kzn"]
            H0 = summary["Ho"]
            He0 = summary["He"]
            A0 = summary["Na"]
            N0 = summary["sample_size"] * (2500 / 199)

            # Get actual lost alleles
            lost_alleles = genetic_data["lost_alleles_ec_kzn"]
            lost_count = sum(len(alleles) for alleles in lost_alleles.values())
            data_source = "CSV"
        else:
            H0 = 0.498
            He0 = 0.565
            A0 = 6.1
            N0 = 2000
            lost_count = 0
            data_source = "default"

        # Scale Ne proportionally
        Ne_scaled = int(Ne * (N0 / 2500))
        t = np.arange(0, generations + 1)

        # Calculate genetic metrics
        Ho_vals = calculate_heterozygosity_loss(H0, Ne_scaled, t)
        He_vals = calculate_heterozygosity_loss(He0, Ne_scaled, t)
        F_vals = calculate_inbreeding(Ne_scaled, t)
        Na_vals = calculate_allelic_diversity(A0, Ne_scaled, t)

        # Calculate population size WITH inbreeding depression
        N_vals = calculate_population_size_with_inbreeding(N0, lambda_val, F_vals, stochastic=stochastic)

        return self._build_result(
            t=t,
            Ho_vals=Ho_vals,
            He_vals=He_vals,
            F_vals=F_vals,
            Na_vals=Na_vals,
            N_vals=N_vals,
            Ne=Ne_scaled,
            initial_values={"Ho": H0, "He": He0, "Na": A0, "N": N0},
            parameters={
                "Ne": Ne_scaled,
                "N0": int(N0),
                "lambda": lambda_val,
                "data_source": data_source,
                "populations": ["Kruger NP", "Limpopo"],
                "lost_populations": ["Eastern Cape", "KwaZulu-Natal"],
                "alleles_lost": lost_count if genetic_data else "unknown",
                "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)",
            }
        )
