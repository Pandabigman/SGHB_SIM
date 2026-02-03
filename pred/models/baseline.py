"""
Model 1: Baseline - All wild populations (no supplementation).
"""

import numpy as np

from .base import (
    BaseModel,
    calculate_heterozygosity_loss,
    calculate_inbreeding,
    calculate_allelic_diversity,
    calculate_population_size_with_inbreeding,
)


class BaselineModel(BaseModel):
    """Model 1: Baseline scenario with all wild populations."""

    model_number = 1
    model_name = "Baseline (All Wild Populations)"

    def run(self, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None):
        """
        Run baseline model using real CSV data if available.

        This model represents the current wild population without any intervention.
        """
        if genetic_data:
            # Use real data
            summary = genetic_data["summaries"]["wild_all"]
            H0 = summary["Ho"]
            He0 = summary["He"]
            A0 = summary["Na"]
            N0 = summary["sample_size"] * (2500 / 199)  # Scale to census
            data_source = "CSV"
        else:
            # Fallback to defaults
            H0 = 0.502
            He0 = 0.568
            A0 = 6.429
            N0 = 2500
            data_source = "default"

        t = np.arange(0, generations + 1)

        # Calculate genetic metrics
        Ho_vals = calculate_heterozygosity_loss(H0, Ne, t)
        He_vals = calculate_heterozygosity_loss(He0, Ne, t)
        F_vals = calculate_inbreeding(Ne, t)
        Na_vals = calculate_allelic_diversity(A0, Ne, t)

        # Calculate population size WITH inbreeding depression
        N_vals = calculate_population_size_with_inbreeding(N0, lambda_val, F_vals, stochastic=stochastic)

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
                "Ne": Ne,
                "N0": int(N0),
                "lambda": lambda_val,
                "data_source": data_source,
                "populations": ["Eastern Cape", "Kruger NP", "KwaZulu-Natal", "Limpopo"],
                "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)",
            }
        )
