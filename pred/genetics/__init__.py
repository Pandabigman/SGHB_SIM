"""
Genetics package for SGHB population simulation.

This package provides functions for:
- Parsing microsatellite CSV data
- Calculating genetic metrics (Ho, He, Na, Fis)
- Population structure analysis
- Captive breeding simulation with Mendelian inheritance
- Supplementation effect modeling
"""

# Constants
from .constants import (
    LOCI,
    DEFAULT_H0,
    DEFAULT_HE0,
    DEFAULT_A0,
    DEFAULT_N0,
    LETHAL_EQUIVALENTS_BIRDS,
    GENERATION_TIME_YEARS
)

# CSV parsing
from .parsing import (
    parse_microsatellite_csv,
    get_population_subset
)

# Genetic metrics
from .metrics import (
    get_genotypes_for_locus,
    get_all_alleles_for_locus,
    calculate_observed_heterozygosity,
    calculate_allele_frequencies,
    calculate_expected_heterozygosity,
    calculate_expected_heterozygosity_population,
    calculate_allelic_richness,
    calculate_fis
)

# Population structure
from .population import (
    get_population_allele_pool,
    get_unique_alleles_per_locus,
    identify_lost_alleles,
    identify_novel_alleles,
    calculate_population_summary,
    load_and_process_genetic_data
)

# Data classes
from .classes import (
    Bird,
    CaptiveBreedingParams,
    CaptivePopulation
)

# Breeding simulation
from .simulation import (
    initialize_captive_population,
    create_offspring,
    breed_captive_population,
    sample_supplementation_birds,
    birds_to_dataframe
)

# Supplementation simulation
from .supplementation import (
    add_individuals_to_population,
    simulate_supplementation_effect,
    simulate_supplementation_effect_breeding,
    simulate_mixed_source_supplementation
)

__all__ = [
    # Constants
    'LOCI',
    'DEFAULT_H0',
    'DEFAULT_HE0',
    'DEFAULT_A0',
    'DEFAULT_N0',
    'LETHAL_EQUIVALENTS_BIRDS',
    'GENERATION_TIME_YEARS',
    # Parsing
    'parse_microsatellite_csv',
    'get_population_subset',
    # Metrics
    'get_genotypes_for_locus',
    'get_all_alleles_for_locus',
    'calculate_observed_heterozygosity',
    'calculate_allele_frequencies',
    'calculate_expected_heterozygosity',
    'calculate_expected_heterozygosity_population',
    'calculate_allelic_richness',
    'calculate_fis',
    # Population
    'get_population_allele_pool',
    'get_unique_alleles_per_locus',
    'identify_lost_alleles',
    'identify_novel_alleles',
    'calculate_population_summary',
    'load_and_process_genetic_data',
    # Classes
    'Bird',
    'CaptiveBreedingParams',
    'CaptivePopulation',
    # Simulation
    'initialize_captive_population',
    'create_offspring',
    'breed_captive_population',
    'sample_supplementation_birds',
    'birds_to_dataframe',
    # Supplementation
    'add_individuals_to_population',
    'simulate_supplementation_effect',
    'simulate_supplementation_effect_breeding',
    'simulate_mixed_source_supplementation',
]
