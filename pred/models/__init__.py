"""
Models package for SGHB population simulation.

Each model represents a different conservation scenario:
- Model 1: Baseline (all wild populations, no intervention)
- Model 2: Population Loss (Kruger + Limpopo only)
- Model 3: Low Supplementation (4 SA captive birds/gen)
- Model 4: High Supplementation (10 SA captive birds/gen)
- Model 5: International Mix (4 pooled SA/AZA/EAZA birds/gen)
- Model 6: Mixed SA Source (6 captive + 2 KZN + 2 EC birds/gen)
"""

from .base import (
    BaseModel,
    calculate_heterozygosity_loss,
    calculate_inbreeding,
    calculate_allelic_diversity,
    calculate_allelic_diversity_with_geneflow,
    calculate_population_size,
    calculate_population_size_with_inbreeding,
    calculate_released_birds_optimized,
    calculate_genetic_rescue_effect,
    LETHAL_EQUIVALENTS_BIRDS,
    GENERATION_TIME_YEARS,
)

from .baseline import BaselineModel
from .population_loss import PopulationLossModel
from .low_supplementation import LowSupplementationModel
from .high_supplementation import HighSupplementationModel
from .international_mix import InternationalMixModel
from .mixed_sa_source import MixedSASourceModel


# Model registry for easy lookup
MODEL_REGISTRY = {
    1: BaselineModel(),
    2: PopulationLossModel(),
    3: LowSupplementationModel(),
    4: HighSupplementationModel(),
    5: InternationalMixModel(),
    6: MixedSASourceModel(),
}


def run_model(model_num, Ne, generations, lambda_val=1.0, stochastic=False, genetic_data=None):
    """
    Run the specified model.

    Args:
        model_num: Model number (1-6)
        Ne: Effective population size
        generations: Number of generations to simulate
        lambda_val: Population growth rate
        stochastic: Whether to include stochasticity
        genetic_data: Loaded genetic data (optional)

    Returns: Dict with simulation results
    """
    if model_num not in MODEL_REGISTRY:
        raise ValueError(f"Model {model_num} not found. Available models: {list(MODEL_REGISTRY.keys())}")
    return MODEL_REGISTRY[model_num].run(Ne, generations, lambda_val, stochastic, genetic_data)


def get_model(model_num):
    """Get a model instance by number."""
    if model_num not in MODEL_REGISTRY:
        raise ValueError(f"Model {model_num} not found.")
    return MODEL_REGISTRY[model_num]


__all__ = [
    # Base
    'BaseModel',
    'calculate_heterozygosity_loss',
    'calculate_inbreeding',
    'calculate_allelic_diversity',
    'calculate_allelic_diversity_with_geneflow',
    'calculate_population_size',
    'calculate_population_size_with_inbreeding',
    'calculate_released_birds_optimized',
    'calculate_genetic_rescue_effect',
    'LETHAL_EQUIVALENTS_BIRDS',
    'GENERATION_TIME_YEARS',
    # Models
    'BaselineModel',
    'PopulationLossModel',
    'LowSupplementationModel',
    'HighSupplementationModel',
    'InternationalMixModel',
    'MixedSASourceModel',
    # Registry
    'MODEL_REGISTRY',
    'run_model',
    'get_model',
]
