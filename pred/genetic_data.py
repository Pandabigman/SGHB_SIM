"""
Genetic Data Parser and Calculator
Processes microsatellite CSV data to calculate real genetic metrics
"""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional
from uuid import uuid4



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



# DYNAMIC CAPTIVE BREEDING SIMULATION


@dataclass
class Bird:
    """Represents a single diploid bird with genotype data"""
    id: str
    genotype: Dict[str, Optional[Tuple[float, float]]]  # {locus: (allele1, allele2) or None}
    sex: str  # 'M' or 'F'
    generation: int  # 0 = founder from CSV
    origin: str  # 'PAAZA', 'AZA', 'EAZA', or 'captive_bred'
    parents: Optional[Tuple[str, str]]  # (sire_id, dam_id) or None for founders
    inbreeding_coefficient: float  # Individual F value


@dataclass
class CaptiveBreedingParams:
    """Parameters controlling captive breeding simulation"""
    target_population_size: int = 140
    captive_Ne: float = 50.0
    breeding_success_rate: float = 0.7
    offspring_per_pair: float = 1.5
    max_breeding_generations: int = 5  # How long a bird can breed


@dataclass
class CaptivePopulation:
    """Manages a breeding captive population"""
    birds: Dict[str, Bird]
    pedigree: Dict[str, Tuple[str, str]]  # {bird_id: (sire_id, dam_id)}
    current_generation: int
    target_population_size: int
    effective_Ne: float
    mean_inbreeding_history: List[float] = field(default_factory=list)


def initialize_captive_population(captive_df: pd.DataFrame, loci_list: List[str] = None) -> CaptivePopulation:
    """
    Initialize captive population from CSV data.
    CSV birds become founders (generation 0).
    """
    if loci_list is None:
        loci_list = LOCI

    birds = {}

    for idx, row in captive_df.iterrows():
        bird_id = f"F_{row['Code']}"  # F = founder
        origin = row['Site'] if 'Site' in row else 'Unknown'

        # Extract genotype for all loci
        genotype = {}
        for locus in loci_list:
            locus_idx = loci_list.index(locus)
            col_start = 3 + (locus_idx * 2)

            allele1 = row.iloc[col_start]
            allele2 = row.iloc[col_start + 1]

            # Handle missing data (0 values)
            if allele1 == 0 or allele2 == 0 or pd.isna(allele1) or pd.isna(allele2):
                genotype[locus] = None
            else:
                genotype[locus] = (float(allele1), float(allele2))

        # Assign sex (alternating for founders)
        sex = 'M' if idx % 2 == 0 else 'F'

        bird = Bird(
            id=bird_id,
            genotype=genotype,
            sex=sex,
            generation=0,
            origin=origin,
            parents=None,
            inbreeding_coefficient=0.0  # Founders assumed non-inbred
        )
        birds[bird_id] = bird

    return CaptivePopulation(
        birds=birds,
        pedigree={},
        current_generation=0,
        target_population_size=len(birds),
        effective_Ne=50.0,
        mean_inbreeding_history=[0.0]
    )


def create_offspring(
    sire: Bird,
    dam: Bird,
    generation: int,
    pedigree: Dict[str, Tuple[str, str]],
    rng: np.random.Generator,
    loci_list: List[str] = None
) -> Bird:
    """
    Create offspring via Mendelian inheritance.
    For each locus, randomly select one allele from each parent.
    """
    if loci_list is None:
        loci_list = LOCI

    offspring_genotype = {}

    for locus in loci_list:
        sire_geno = sire.genotype.get(locus)
        dam_geno = dam.genotype.get(locus)

        # Handle missing data - if either parent missing, offspring missing
        if sire_geno is None or dam_geno is None:
            offspring_genotype[locus] = None
        else:
            # Mendel's Law: random allele from each parent
            sire_allele = sire_geno[rng.integers(0, 2)]
            dam_allele = dam_geno[rng.integers(0, 2)]
            offspring_genotype[locus] = (sire_allele, dam_allele)

    # Assign sex randomly (1:1 ratio)
    sex = 'M' if rng.random() < 0.5 else 'F'

    # Calculate inbreeding coefficient
    # F_offspring = kinship(sire, dam) - approximate using parental F
    parent_F_avg = (sire.inbreeding_coefficient + dam.inbreeding_coefficient) / 2
    # Simplified: offspring F increases slightly each generation
    offspring_F = min(parent_F_avg + 0.01, 1.0)

    return Bird(
        id=f"CB_G{generation}_{uuid4().hex[:8]}",
        genotype=offspring_genotype,
        sex=sex,
        generation=generation,
        origin='captive_bred',
        parents=(sire.id, dam.id),
        inbreeding_coefficient=offspring_F
    )


def breed_captive_population(
    captive_pop: CaptivePopulation,
    params: CaptiveBreedingParams,
    rng: np.random.Generator,
    loci_list: List[str] = None
) -> CaptivePopulation:
    """
    Simulate one generation of captive breeding.

    Algorithm:
    1. Select breeding pairs (random pairing)
    2. For each pair, produce offspring via Mendelian inheritance
    3. Apply mortality/culling to maintain target population size
    4. Update pedigree
    """
    if loci_list is None:
        loci_list = LOCI

    # Get breeding-age birds
    males = [b for b in captive_pop.birds.values()
             if b.sex == 'M' and (captive_pop.current_generation - b.generation) < params.max_breeding_generations]
    females = [b for b in captive_pop.birds.values()
               if b.sex == 'F' and (captive_pop.current_generation - b.generation) < params.max_breeding_generations]

    if len(males) == 0 or len(females) == 0:
        # No breeding possible - return as is with incremented generation
        return CaptivePopulation(
            birds=captive_pop.birds.copy(),
            pedigree=captive_pop.pedigree.copy(),
            current_generation=captive_pop.current_generation + 1,
            target_population_size=captive_pop.target_population_size,
            effective_Ne=captive_pop.effective_Ne,
            mean_inbreeding_history=captive_pop.mean_inbreeding_history.copy()
        )

    # Random pairing
    rng.shuffle(males)
    rng.shuffle(females)
    n_pairs = min(len(males), len(females))
    pairs = list(zip(males[:n_pairs], females[:n_pairs]))

    # Produce offspring
    offspring = []
    new_pedigree = captive_pop.pedigree.copy()

    for sire, dam in pairs:
        if rng.random() < params.breeding_success_rate:
            n_offspring = rng.poisson(params.offspring_per_pair)
            for _ in range(n_offspring):
                chick = create_offspring(
                    sire, dam,
                    captive_pop.current_generation + 1,
                    new_pedigree,
                    rng, loci_list
                )
                offspring.append(chick)
                new_pedigree[chick.id] = (sire.id, dam.id)

    # Combine existing birds and offspring
    all_birds = list(captive_pop.birds.values()) + offspring

    # Cull to maintain target population size (remove oldest first)
    if len(all_birds) > params.target_population_size:
        # Sort by generation (oldest first), then randomly within generation
        all_birds.sort(key=lambda b: (b.generation, rng.random()))
        # Keep youngest birds up to target size
        all_birds = all_birds[-params.target_population_size:]

    # Calculate mean inbreeding
    mean_F = np.mean([b.inbreeding_coefficient for b in all_birds]) if all_birds else 0.0
    new_history = captive_pop.mean_inbreeding_history + [mean_F]

    return CaptivePopulation(
        birds={b.id: b for b in all_birds},
        pedigree=new_pedigree,
        current_generation=captive_pop.current_generation + 1,
        target_population_size=params.target_population_size,
        effective_Ne=captive_pop.effective_Ne,
        mean_inbreeding_history=new_history
    )


def sample_supplementation_birds(
    captive_pop: CaptivePopulation,
    n_birds: int,
    rng: np.random.Generator
) -> List[Bird]:
    """
    Sample n birds from current captive population for supplementation.
    Samples WITHOUT replacement, preferring birds with lower inbreeding.
    """
    available_birds = list(captive_pop.birds.values())

    if n_birds >= len(available_birds):
        return available_birds.copy()

    # Weight selection by inverse inbreeding (less inbred = more likely selected)
    weights = np.array([1.0 / (1.0 + b.inbreeding_coefficient) for b in available_birds])
    weights /= weights.sum()

    selected_indices = rng.choice(
        len(available_birds),
        size=n_birds,
        replace=False,
        p=weights
    )

    return [available_birds[i] for i in selected_indices]


def birds_to_dataframe(birds: List[Bird], loci_list: List[str] = None) -> pd.DataFrame:
    """
    Convert list of Bird objects back to DataFrame format for compatibility
    with existing genetic calculation functions.
    """
    if loci_list is None:
        loci_list = LOCI

    if not birds:
        return pd.DataFrame()

    rows = []

    for bird in birds:
        row = {
            'Code': bird.id,
            'Site': bird.origin,
            'Status': 'Captive'
        }
        rows.append(row)

    df = pd.DataFrame(rows)

    # Add genotype columns in the correct positions
    for locus_idx, locus in enumerate(loci_list):
        col1_values = []
        col2_values = []

        for bird in birds:
            geno = bird.genotype.get(locus)
            if geno is None:
                col1_values.append(0)
                col2_values.append(0)
            else:
                col1_values.append(geno[0])
                col2_values.append(geno[1])

        # Insert at correct column positions (after Code, Site, Status)
        df.insert(3 + locus_idx * 2, f'Locus{locus_idx}_1', col1_values)
        df.insert(3 + locus_idx * 2 + 1, f'Locus{locus_idx}_2', col2_values)

    return df


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
    seed: int = None
) -> List[Dict]:
    """
    Model 6: Mixed source supplementation with wild population depletion.

    Sources per generation:
    - captive_birds_per_gen from breeding PAAZA population
    - kzn_birds_per_gen from KZN wild (depletes source)
    - ec_birds_per_gen from E Cape wild (depletes source)

    Wild source populations are sampled WITHOUT replacement, realistically
    modeling removal of birds from already-threatened source populations.
    """
    if loci_list is None:
        loci_list = LOCI

    rng = np.random.default_rng(seed)

    if captive_params is None:
        captive_params = CaptiveBreedingParams()

    # Initialize captive breeding population
    captive_pop = initialize_captive_population(captive_df, loci_list)

    # Create mutable copies of wild source populations (for depletion)
    kzn_source = kzn_df.copy()
    ec_source = ec_df.copy()

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
            'kzn_source_remaining': len(kzn_source),
            'ec_source_remaining': len(ec_source)
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

            # 2. Sample from KZN wild (with depletion)
            n_kzn = min(kzn_birds_per_gen, len(kzn_source))
            if n_kzn > 0:
                kzn_sample = kzn_source.sample(n=n_kzn, replace=False)
                supplementation_dfs.append(kzn_sample)
                # Remove sampled birds from source (depletion)
                kzn_source = kzn_source.drop(kzn_sample.index)

            # 3. Sample from E Cape wild (with depletion)
            n_ec = min(ec_birds_per_gen, len(ec_source))
            if n_ec > 0:
                ec_sample = ec_source.sample(n=n_ec, replace=False)
                supplementation_dfs.append(ec_sample)
                # Remove sampled birds from source (depletion)
                ec_source = ec_source.drop(ec_sample.index)

            # Combine all supplementation sources
            if supplementation_dfs:
                combined_supplementation = pd.concat(supplementation_dfs, ignore_index=True)
                current_wild_population = pd.concat(
                    [current_wild_population, combined_supplementation],
                    ignore_index=True
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