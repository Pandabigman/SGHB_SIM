"""
Genetic Data Parser and Calculator
Processes microsatellite CSV data to calculate real genetic metrics
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter



# CSV PARSING


# Microsatellite loci (14 total)
LOCI = ['Buco4', 'Buco11', 'Buco2', 'Buco9', 'GHB21', 'GHB19', 'GHB26', 
        'GHB20', 'Buco16', 'Buco18', 'Buco19', 'Buco21', 'Buco24', 'Buco25']


def parse_microsatellite_csv(csv_path):
    """
    Parse microsatellite CSV with format:
    Code, Site, Status, Locus1_allele1, Locus1_allele2, Locus2_allele1, ...
    
    Returns: DataFrame with processed genetic data
    """
    df = pd.read_csv(csv_path)
    
    # Clean column names (remove extra spaces)
    df.columns = df.columns.str.strip()
    
    # Standardize site names
    site_mapping = {
        'Eastern Cape province': 'Eastern Cape',
        'Kruger National Park': 'Kruger',
        'KwaZulu-Natal province': 'KwaZulu-Natal',
        'Limpopo province': 'Limpopo'
    }
    df['Site'] = df['Site'].replace(site_mapping)
    
    return df


def get_population_subset(df, populations):
    """
    Extract specific populations from dataframe
    
    Args:
        df: Full dataframe
        populations: List of population names (e.g., ['Eastern Cape', 'Kruger'])
    
    Returns: Subset dataframe
    """
    return df[df['Site'].isin(populations)].copy()



# GENOTYPE EXTRACTION


def get_genotypes_for_locus(df, locus, loci_list=LOCI):
    """
    OPTIMIZED: Extract both alleles for a specific locus from all individuals
    Uses vectorized operations instead of iterrows()

    Returns: List of tuples [(allele1, allele2), ...]
    """
    # Find column indices for this locus
    # CSV format: each locus has 2 consecutive columns
    locus_idx = loci_list.index(locus)
    col_start = 3 + (locus_idx * 2)  # Skip Code, Site, Status columns

    # OPTIMIZED: Use vectorized column access instead of row iteration
    allele1_col = df.iloc[:, col_start].replace(0, np.nan)
    allele2_col = df.iloc[:, col_start + 1].replace(0, np.nan)

    # Convert to numeric, coercing errors to NaN
    allele1_col = pd.to_numeric(allele1_col, errors='coerce')
    allele2_col = pd.to_numeric(allele2_col, errors='coerce')

    # Filter out rows with missing data
    valid_mask = allele1_col.notna() & allele2_col.notna()

    # Create genotype tuples using zip (much faster than loop)
    genotypes = list(zip(allele1_col[valid_mask].values, allele2_col[valid_mask].values))

    return genotypes


def get_all_alleles_for_locus(df, locus, loci_list=LOCI):
    """
    Get list of all alleles (not genotypes) for a locus
    
    Returns: List of alleles [162, 167, 180, 162, ...]
    """
    genotypes = get_genotypes_for_locus(df, locus, loci_list)
    alleles = []
    for a1, a2 in genotypes:
        alleles.extend([a1, a2])
    return alleles



# GENETIC METRIC CALCULATIONS


def calculate_observed_heterozygosity(df, loci_list=LOCI):
    """
    Calculate observed heterozygosity from actual genotypes
    Ho = proportion of heterozygous individuals
    
    Returns: Mean Ho across all loci
    """
    ho_per_locus = []
    
    for locus in loci_list:
        genotypes = get_genotypes_for_locus(df, locus, loci_list)
        
        if len(genotypes) == 0:
            continue
        
        # Count heterozygotes (allele1 != allele2)
        heterozygotes = sum(1 for a1, a2 in genotypes if a1 != a2)
        ho = heterozygotes / len(genotypes)
        ho_per_locus.append(ho)
    
    return np.mean(ho_per_locus) if ho_per_locus else 0.0


def calculate_allele_frequencies(df, locus, loci_list=LOCI):
    """
    Calculate allele frequencies for a locus
    
    Returns: Dict {allele: frequency}
    """
    alleles = get_all_alleles_for_locus(df, locus, loci_list)
    
    if len(alleles) == 0:
        return {}
    
    allele_counts = Counter(alleles)
    total = len(alleles)
    
    return {allele: count/total for allele, count in allele_counts.items()}


def calculate_expected_heterozygosity(allele_frequencies):
    """
    Calculate expected heterozygosity from allele frequencies
    He = 1 - Σ(pi²)
    
    Args:
        allele_frequencies: Dict {allele: frequency}
    
    Returns: Expected heterozygosity
    """
    if not allele_frequencies:
        return 0.0
    
    sum_squared = sum(freq**2 for freq in allele_frequencies.values())
    return 1 - sum_squared


def calculate_expected_heterozygosity_population(df, loci_list=LOCI):
    """
    Calculate mean expected heterozygosity across all loci
    
    Returns: Mean He
    """
    he_per_locus = []
    
    for locus in loci_list:
        freqs = calculate_allele_frequencies(df, locus, loci_list)
        if freqs:
            he = calculate_expected_heterozygosity(freqs)
            he_per_locus.append(he)
    
    return np.mean(he_per_locus) if he_per_locus else 0.0


def calculate_allelic_richness(df, loci_list=LOCI):
    """
    Calculate mean number of alleles per locus
    
    Returns: Mean Na
    """
    na_per_locus = []
    
    for locus in loci_list:
        alleles = get_all_alleles_for_locus(df, locus, loci_list)
        unique_alleles = set(alleles)
        na_per_locus.append(len(unique_alleles))
    
    return np.mean(na_per_locus) if na_per_locus else 0.0


def calculate_fis(Ho, He):
    """
    Calculate inbreeding coefficient
    FIS = (He - Ho) / He
    """
    if He == 0:
        return 0.0
    return (He - Ho) / He



# POPULATION GENETIC STRUCTURE


def get_population_allele_pool(df, loci_list=LOCI):
    """
    Get complete allele pool for a population
    
    Returns: Dict {locus: [list of all alleles]}
    """
    allele_pool = {}
    
    for locus in loci_list:
        alleles = get_all_alleles_for_locus(df, locus, loci_list)
        allele_pool[locus] = alleles
    
    return allele_pool


def get_unique_alleles_per_locus(df, loci_list=LOCI):
    """
    Get set of unique alleles for each locus
    
    Returns: Dict {locus: set(unique_alleles)}
    """
    unique_alleles = {}
    
    for locus in loci_list:
        alleles = get_all_alleles_for_locus(df, locus, loci_list)
        unique_alleles[locus] = set(alleles)
    
    return unique_alleles


def identify_lost_alleles(full_population_df, reduced_population_df, loci_list=LOCI):
    """
    Identify alleles present in full population but lost in reduced population
    
    Returns: Dict {locus: set(lost_alleles)}
    """
    full_alleles = get_unique_alleles_per_locus(full_population_df, loci_list)
    reduced_alleles = get_unique_alleles_per_locus(reduced_population_df, loci_list)
    
    lost_alleles = {}
    for locus in loci_list:
        lost = full_alleles[locus] - reduced_alleles[locus]
        if lost:
            lost_alleles[locus] = lost
    
    return lost_alleles


def identify_novel_alleles(recipient_df, donor_df, loci_list=LOCI):
    """
    Identify alleles in donor population not present in recipient
    These are "novel" alleles that supplementation would add
    
    Returns: Dict {locus: set(novel_alleles)}
    """
    recipient_alleles = get_unique_alleles_per_locus(recipient_df, loci_list)
    donor_alleles = get_unique_alleles_per_locus(donor_df, loci_list)
    
    novel_alleles = {}
    for locus in loci_list:
        novel = donor_alleles[locus] - recipient_alleles[locus]
        if novel:
            novel_alleles[locus] = novel
    
    return novel_alleles



# SUPPLEMENTATION SIMULATION


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


def simulate_supplementation_effect(wild_df, captive_df, birds_per_gen, generations, loci_list=LOCI, base_Ne=500):
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
    results = []
    current_population = wild_df.copy()
    initial_wild_size = len(wild_df)

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



# SUMMARY STATISTICS


def calculate_population_summary(df, population_name, loci_list=LOCI):
    """
    Calculate complete genetic summary for a population
    
    Returns: Dict with all metrics
    """
    Ho = calculate_observed_heterozygosity(df, loci_list)
    He = calculate_expected_heterozygosity_population(df, loci_list)
    Na = calculate_allelic_richness(df, loci_list)
    FIS = calculate_fis(Ho, He)
    
    # Get allele frequencies for all loci
    allele_freqs = {}
    for locus in loci_list:
        freqs = calculate_allele_frequencies(df, locus, loci_list)
        allele_freqs[locus] = freqs
    
    return {
        'name': population_name,
        'sample_size': len(df),
        'Ho': float(Ho),
        'He': float(He),
        'Na': float(Na),
        'FIS': float(FIS),
        'allele_frequencies': allele_freqs,
        'unique_alleles_per_locus': {
            locus: len(freqs) for locus, freqs in allele_freqs.items()
        }
    }



# DATA LOADER


def load_and_process_genetic_data(wild_csv_path, captive_csv_path):
    """
    Load both CSVs and calculate all genetic metrics
    
    Returns: Dict with all population data
    """
    # Load CSVs
    wild_df = parse_microsatellite_csv(wild_csv_path)
    captive_df = parse_microsatellite_csv(captive_csv_path)
    
    # Separate populations
    populations = {
        'wild_all': wild_df[wild_df['Site'].isin(['Eastern Cape', 'Kruger', 'KwaZulu-Natal', 'Limpopo'])],
        'wild_ec': wild_df[wild_df['Site'] == 'Eastern Cape'],
        'wild_kruger': wild_df[wild_df['Site'] == 'Kruger'],
        'wild_kzn': wild_df[wild_df['Site'] == 'KwaZulu-Natal'],
        'wild_limpopo': wild_df[wild_df['Site'] == 'Limpopo'],
        'wild_no_ec_kzn': wild_df[wild_df['Site'].isin(['Kruger', 'Limpopo'])],
        'paaza': captive_df[captive_df['Site'] == 'PAAZA'],
        'aza': captive_df[captive_df['Site'] == 'AZA'],
        'eaza': captive_df[captive_df['Site'] == 'EAZA']
    }
    
    # Calculate summaries for each population
    summaries = {}
    for pop_name, pop_df in populations.items():
        if len(pop_df) > 0:
            summaries[pop_name] = calculate_population_summary(pop_df, pop_name)
    
    # Identify lost alleles (Model 2)
    lost_alleles = identify_lost_alleles(
        populations['wild_all'], 
        populations['wild_no_ec_kzn']
    )
    
    # Identify novel alleles from captive populations (Models 3-5)
    novel_from_paaza = identify_novel_alleles(populations['wild_all'], populations['paaza'])
    novel_from_aza = identify_novel_alleles(populations['wild_all'], populations['aza'])
    novel_from_eaza = identify_novel_alleles(populations['wild_all'], populations['eaza'])
    
    return {
        'dataframes': populations,
        'summaries': summaries,
        'lost_alleles_ec_kzn': lost_alleles,
        'novel_alleles': {
            'paaza': novel_from_paaza,
            'aza': novel_from_aza,
            'eaza': novel_from_eaza
        },
        'loci': LOCI
    }