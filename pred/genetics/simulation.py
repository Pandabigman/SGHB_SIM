"""
Captive breeding simulation functions with Mendelian inheritance.
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, List
from uuid import uuid4

from .constants import LOCI
from .classes import Bird, CaptiveBreedingParams, CaptivePopulation


def initialize_captive_population(
    captive_df: pd.DataFrame,
    loci_list: List[str] = None,
    id_prefix: str = 'F',
    effective_Ne: float = 50.0,
    target_population_size: int = None
) -> CaptivePopulation:
    """
    Initialize captive population from CSV data.
    CSV birds become founders (generation 0).
    """
    if loci_list is None:
        loci_list = LOCI

    birds = {}

    for idx, row in captive_df.iterrows():
        bird_id = f"{id_prefix}_{row['Code']}"
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
        target_population_size=target_population_size if target_population_size is not None else len(birds),
        effective_Ne=effective_Ne,
        mean_inbreeding_history=[0.0]
    )


def create_offspring(
    sire: Bird,
    dam: Bird,
    generation: int,
    pedigree: Dict[str, Tuple[str, str]],
    rng: np.random.Generator,
    loci_list: List[str] = None,
    origin: str = 'captive_bred',
    id_prefix: str = 'CB'
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
        id=f"{id_prefix}_G{generation}_{uuid4().hex[:8]}",
        genotype=offspring_genotype,
        sex=sex,
        generation=generation,
        origin=origin,
        parents=(sire.id, dam.id),
        inbreeding_coefficient=offspring_F
    )


def breed_captive_population(
    captive_pop: CaptivePopulation,
    params: CaptiveBreedingParams,
    rng: np.random.Generator,
    loci_list: List[str] = None,
    origin: str = 'captive_bred',
    id_prefix: str = 'CB'
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
                    rng, loci_list,
                    origin=origin, id_prefix=id_prefix
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
