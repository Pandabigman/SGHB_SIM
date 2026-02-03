"""
Genetic metric calculations (Ho, He, Na, Fis).
"""

import numpy as np
import pandas as pd
from collections import Counter

from .constants import LOCI


def get_genotypes_for_locus(df, locus, loci_list=None):
    """
    OPTIMIZED: Extract both alleles for a specific locus from all individuals
    Uses vectorized operations instead of iterrows()

    Returns: List of tuples [(allele1, allele2), ...]
    """
    if loci_list is None:
        loci_list = LOCI

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


def get_all_alleles_for_locus(df, locus, loci_list=None):
    """
    Get list of all alleles (not genotypes) for a locus

    Returns: List of alleles [162, 167, 180, 162, ...]
    """
    if loci_list is None:
        loci_list = LOCI

    genotypes = get_genotypes_for_locus(df, locus, loci_list)
    alleles = []
    for a1, a2 in genotypes:
        alleles.extend([a1, a2])
    return alleles


def calculate_observed_heterozygosity(df, loci_list=None):
    """
    Calculate observed heterozygosity from actual genotypes
    Ho = proportion of heterozygous individuals

    Returns: Mean Ho across all loci
    """
    if loci_list is None:
        loci_list = LOCI

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


def calculate_allele_frequencies(df, locus, loci_list=None):
    """
    Calculate allele frequencies for a locus

    Returns: Dict {allele: frequency}
    """
    if loci_list is None:
        loci_list = LOCI

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


def calculate_expected_heterozygosity_population(df, loci_list=None):
    """
    Calculate mean expected heterozygosity across all loci

    Returns: Mean He
    """
    if loci_list is None:
        loci_list = LOCI

    he_per_locus = []

    for locus in loci_list:
        freqs = calculate_allele_frequencies(df, locus, loci_list)
        if freqs:
            he = calculate_expected_heterozygosity(freqs)
            he_per_locus.append(he)

    return np.mean(he_per_locus) if he_per_locus else 0.0


def calculate_allelic_richness(df, loci_list=None):
    """
    Calculate mean number of alleles per locus

    Returns: Mean Na
    """
    if loci_list is None:
        loci_list = LOCI

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
