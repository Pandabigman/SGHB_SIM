"""
Population genetic structure functions and data loading.
"""

from .constants import LOCI
from .parsing import parse_microsatellite_csv
from .metrics import (
    get_all_alleles_for_locus,
    calculate_observed_heterozygosity,
    calculate_expected_heterozygosity_population,
    calculate_allelic_richness,
    calculate_fis,
    calculate_allele_frequencies
)


def get_population_allele_pool(df, loci_list=None):
    """
    Get complete allele pool for a population

    Returns: Dict {locus: [list of all alleles]}
    """
    if loci_list is None:
        loci_list = LOCI

    allele_pool = {}

    for locus in loci_list:
        alleles = get_all_alleles_for_locus(df, locus, loci_list)
        allele_pool[locus] = alleles

    return allele_pool


def get_unique_alleles_per_locus(df, loci_list=None):
    """
    Get set of unique alleles for each locus

    Returns: Dict {locus: set(unique_alleles)}
    """
    if loci_list is None:
        loci_list = LOCI

    unique_alleles = {}

    for locus in loci_list:
        alleles = get_all_alleles_for_locus(df, locus, loci_list)
        unique_alleles[locus] = set(alleles)

    return unique_alleles


def identify_lost_alleles(full_population_df, reduced_population_df, loci_list=None):
    """
    Identify alleles present in full population but lost in reduced population

    Returns: Dict {locus: set(lost_alleles)}
    """
    if loci_list is None:
        loci_list = LOCI

    full_alleles = get_unique_alleles_per_locus(full_population_df, loci_list)
    reduced_alleles = get_unique_alleles_per_locus(reduced_population_df, loci_list)

    lost_alleles = {}
    for locus in loci_list:
        lost = full_alleles[locus] - reduced_alleles[locus]
        if lost:
            lost_alleles[locus] = lost

    return lost_alleles


def identify_novel_alleles(recipient_df, donor_df, loci_list=None):
    """
    Identify alleles in donor population not present in recipient
    These are "novel" alleles that supplementation would add

    Returns: Dict {locus: set(novel_alleles)}
    """
    if loci_list is None:
        loci_list = LOCI

    recipient_alleles = get_unique_alleles_per_locus(recipient_df, loci_list)
    donor_alleles = get_unique_alleles_per_locus(donor_df, loci_list)

    novel_alleles = {}
    for locus in loci_list:
        novel = donor_alleles[locus] - recipient_alleles[locus]
        if novel:
            novel_alleles[locus] = novel

    return novel_alleles


def calculate_population_summary(df, population_name, loci_list=None):
    """
    Calculate complete genetic summary for a population

    Returns: Dict with all metrics
    """
    if loci_list is None:
        loci_list = LOCI

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
