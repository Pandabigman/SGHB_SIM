"""
CSV parsing functions for microsatellite genetic data.
"""

import pandas as pd


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
    df['Site'] = df['Site'].str.strip().replace(site_mapping)

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
