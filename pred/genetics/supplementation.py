"""
Supplementation simulation functions.
"""

import numpy as np
import pandas as pd
from typing import Dict, List

from .constants import LOCI
from .classes import CaptiveBreedingParams, CaptivePopulation
from .metrics import (
    calculate_observed_heterozygosity,
    calculate_expected_heterozygosity_population,
    calculate_allelic_richness,
    calculate_fis
)
from .simulation import (
    initialize_captive_population,
    breed_captive_population,
    sample_supplementation_birds,
    birds_to_dataframe
)


def add_individuals_to_population(recipient_df, donor_df, n_individuals):
    """
    OPTIMIZED: Add n random individuals from donor to recipient population
    Uses replace=True to allow resampling and avoid copy overhead

    Returns: Combined dataframe
    """
    # OPTIMIZED: Always sample with replacement to avoid exhausting small donor pools
    # and avoid expensive copy operations
    donor_sample = donor_df.sample(n=n_individuals, replace=True)

    # OPTIMIZED: Use pd.concat with copy=False where possible
    combined = pd.concat([recipient_df, donor_sample], ignore_index=True, copy=False)
    return combined


def simulate_supplementation_effect(wild_df, captive_df, birds_per_gen, generations, loci_list=None, base_Ne=500):
    """
    Simulate adding captive birds to wild population over generations
    Track actual genetic changes INCLUDING effective population size (Ne)

    Args:
        wild_df: Wild population dataframe
        captive_df: Captive population dataframe
        birds_per_gen: Number of birds added per generation
        generations: Number of generations to simulate
        loci_list: List of loci to analyze
        base_Ne: Base effective population size (default 500)

    Returns: List of genetic metrics per generation (includes effective_Ne)
    """
    if loci_list is None:
        loci_list = LOCI

    results = []
    current_population = wild_df.copy()

    for gen in range(generations + 1):
        # Calculate current metrics
        Ho = calculate_observed_heterozygosity(current_population, loci_list)
        He = calculate_expected_heterozygosity_population(current_population, loci_list)
        Na = calculate_allelic_richness(current_population, loci_list)
        FIS = calculate_fis(Ho, He)
        pop_size = len(current_population)

        # Calculate effective Ne based on gene flow
        # Gene flow increases Ne: each immigrant contributes ~0.5 to effective size
        cumulative_immigrants = birds_per_gen * gen
        effective_Ne = base_Ne + (cumulative_immigrants * 0.5)

        results.append({
            'generation': gen,
            'Ho': Ho,
            'He': He,
            'Na': Na,
            'FIS': FIS,
            'population_size': pop_size,
            'effective_Ne': effective_Ne
        })

        # Add birds for next generation (except at last generation)
        if gen < generations:
            current_population = add_individuals_to_population(
                current_population, captive_df, birds_per_gen
            )

    return results


def simulate_supplementation_effect_breeding(
    wild_df: pd.DataFrame,
    initial_captive_df: pd.DataFrame,
    birds_per_gen: int,
    generations: int,
    loci_list: List[str] = None,
    base_Ne: float = 500,
    captive_params: CaptiveBreedingParams = None,
    seed: int = None
) -> List[Dict]:
    """
    Simulate supplementation with dynamic captive breeding population.

    Key difference from original:
    - Captive population breeds each generation, creating new genetic combinations
    - Supplementation birds are sampled from CURRENT captive population
    - Tracks captive population genetics separately
    """
    if loci_list is None:
        loci_list = LOCI

    rng = np.random.default_rng(seed)

    if captive_params is None:
        captive_params = CaptiveBreedingParams()

    # Initialize captive population from CSV (founders)
    captive_pop = initialize_captive_population(initial_captive_df, loci_list)

    results = []
    current_wild_population = wild_df.copy()

    for gen in range(generations + 1):
        # Calculate current wild population metrics
        Ho = calculate_observed_heterozygosity(current_wild_population, loci_list)
        He = calculate_expected_heterozygosity_population(current_wild_population, loci_list)
        Na = calculate_allelic_richness(current_wild_population, loci_list)
        FIS = calculate_fis(Ho, He)
        pop_size = len(current_wild_population)

        # Calculate effective Ne with gene flow
        cumulative_immigrants = birds_per_gen * gen
        effective_Ne = base_Ne + (cumulative_immigrants * 0.5)

        # Track captive population metrics
        captive_birds_list = list(captive_pop.birds.values())
        captive_F_mean = np.mean([b.inbreeding_coefficient for b in captive_birds_list]) if captive_birds_list else 0.0

        results.append({
            'generation': gen,
            'Ho': Ho,
            'He': He,
            'Na': Na,
            'FIS': FIS,
            'population_size': pop_size,
            'effective_Ne': effective_Ne,
            'captive_F_mean': captive_F_mean,
            'captive_pop_size': len(captive_pop.birds)
        })

        # Advance simulation (except at last generation)
        if gen < generations:
            # 1. Breed captive population for next generation
            captive_pop = breed_captive_population(captive_pop, captive_params, rng, loci_list)

            # 2. Sample supplementation birds from CURRENT captive population
            supplementation_birds = sample_supplementation_birds(captive_pop, birds_per_gen, rng)

            # 3. Convert to DataFrame and add to wild population
            supplementation_df = birds_to_dataframe(supplementation_birds, loci_list)
            if len(supplementation_df) > 0:
                current_wild_population = pd.concat(
                    [current_wild_population, supplementation_df],
                    ignore_index=True
                )

    return results


def simulate_mixed_source_supplementation(
    wild_df: pd.DataFrame,
    captive_df: pd.DataFrame,
    kzn_df: pd.DataFrame,
    ec_df: pd.DataFrame,
    captive_birds_per_gen: int,
    kzn_birds_per_gen: int,
    ec_birds_per_gen: int,
    generations: int,
    loci_list: List[str] = None,
    base_Ne: float = 500,
    captive_params: CaptiveBreedingParams = None,
    kzn_params: CaptiveBreedingParams = None,
    ec_params: CaptiveBreedingParams = None,
    seed: int = None
) -> List[Dict]:
    """
    Model 6: Mixed source supplementation with breeding wild populations.

    Sources per generation:
    - captive_birds_per_gen from breeding PAAZA population
    - kzn_birds_per_gen from breeding KZN wild population
    - ec_birds_per_gen from breeding EC wild population

    All source populations breed each generation via Mendelian inheritance,
    producing new recombinant genotypes over time.
    """
    if loci_list is None:
        loci_list = LOCI

    rng = np.random.default_rng(seed)

    if captive_params is None:
        captive_params = CaptiveBreedingParams()

    # Initialize captive breeding population
    captive_pop = initialize_captive_population(captive_df, loci_list)

    # Initialize KZN and EC as breeding populations
    if kzn_params is None:
        kzn_params = CaptiveBreedingParams(
            target_population_size=max(len(kzn_df), 50),
            captive_Ne=15.0,
            breeding_success_rate=0.5,
            offspring_per_pair=1.2
        )
    if ec_params is None:
        ec_params = CaptiveBreedingParams(
            target_population_size=max(len(ec_df), 30),
            captive_Ne=10.0,
            breeding_success_rate=0.5,
            offspring_per_pair=1.2
        )

    kzn_pop = initialize_captive_population(
        kzn_df, loci_list, id_prefix='FK',
        effective_Ne=kzn_params.captive_Ne,
        target_population_size=kzn_params.target_population_size
    )
    ec_pop = initialize_captive_population(
        ec_df, loci_list, id_prefix='FE',
        effective_Ne=ec_params.captive_Ne,
        target_population_size=ec_params.target_population_size
    )

    results = []
    current_wild_population = wild_df.copy()

    total_birds_per_gen = captive_birds_per_gen + kzn_birds_per_gen + ec_birds_per_gen

    for gen in range(generations + 1):
        # Calculate current wild population metrics
        Ho = calculate_observed_heterozygosity(current_wild_population, loci_list)
        He = calculate_expected_heterozygosity_population(current_wild_population, loci_list)
        Na = calculate_allelic_richness(current_wild_population, loci_list)
        FIS = calculate_fis(Ho, He)
        pop_size = len(current_wild_population)

        # Calculate effective Ne with gene flow
        cumulative_immigrants = total_birds_per_gen * gen
        effective_Ne = base_Ne + (cumulative_immigrants * 0.5)

        # Track captive and source population metrics
        captive_birds_list = list(captive_pop.birds.values())
        captive_F_mean = np.mean([b.inbreeding_coefficient for b in captive_birds_list]) if captive_birds_list else 0.0

        kzn_birds_list = list(kzn_pop.birds.values())
        kzn_F_mean = np.mean([b.inbreeding_coefficient for b in kzn_birds_list]) if kzn_birds_list else 0.0

        ec_birds_list = list(ec_pop.birds.values())
        ec_F_mean = np.mean([b.inbreeding_coefficient for b in ec_birds_list]) if ec_birds_list else 0.0

        results.append({
            'generation': gen,
            'Ho': Ho,
            'He': He,
            'Na': Na,
            'FIS': FIS,
            'population_size': pop_size,
            'effective_Ne': effective_Ne,
            'captive_F_mean': captive_F_mean,
            'captive_pop_size': len(captive_pop.birds),
            'kzn_F_mean': kzn_F_mean,
            'kzn_pop_size': len(kzn_pop.birds),
            'ec_F_mean': ec_F_mean,
            'ec_pop_size': len(ec_pop.birds),
        })

        # Advance simulation (except at last generation)
        if gen < generations:
            supplementation_dfs = []

            # 1. Breed captive population and sample
            captive_pop = breed_captive_population(captive_pop, captive_params, rng, loci_list)
            captive_birds = sample_supplementation_birds(captive_pop, captive_birds_per_gen, rng)
            captive_sample_df = birds_to_dataframe(captive_birds, loci_list)
            if len(captive_sample_df) > 0:
                supplementation_dfs.append(captive_sample_df)

            # 2. Breed KZN wild population and sample
            kzn_pop = breed_captive_population(
                kzn_pop, kzn_params, rng, loci_list,
                origin='KwaZulu-Natal', id_prefix='KZN'
            )
            if kzn_birds_per_gen > 0 and len(kzn_pop.birds) > 0:
                kzn_birds = sample_supplementation_birds(kzn_pop, kzn_birds_per_gen, rng)
                kzn_sample_df = birds_to_dataframe(kzn_birds, loci_list)
                if len(kzn_sample_df) > 0:
                    supplementation_dfs.append(kzn_sample_df)

            # 3. Breed EC wild population and sample
            ec_pop = breed_captive_population(
                ec_pop, ec_params, rng, loci_list,
                origin='Eastern Cape', id_prefix='EC'
            )
            if ec_birds_per_gen > 0 and len(ec_pop.birds) > 0:
                ec_birds = sample_supplementation_birds(ec_pop, ec_birds_per_gen, rng)
                ec_sample_df = birds_to_dataframe(ec_birds, loci_list)
                if len(ec_sample_df) > 0:
                    supplementation_dfs.append(ec_sample_df)

            # Combine all supplementation sources
            if supplementation_dfs:
                combined_supplementation = pd.concat(supplementation_dfs, ignore_index=True)
                current_wild_population = pd.concat(
                    [current_wild_population, combined_supplementation],
                    ignore_index=True
                )

    return results
