"""
Data classes for captive breeding simulation.
"""

from dataclasses import dataclass, field
from typing import Dict, Tuple, List, Optional


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
