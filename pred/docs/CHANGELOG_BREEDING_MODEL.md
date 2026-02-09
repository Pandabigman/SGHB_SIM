# Changelog: Dynamic Captive Breeding & Model 6

**Date**: February 2026

---

## Summary

Two major changes were implemented to address a fundamental flaw in the genetic simulation and to add a new conservation scenario requested by collaborators.

---

## Change 1: Dynamic Captive Breeding Simulation

### The Problem

The original simulation sampled supplementation birds from a **fixed pool of 140 captive birds** (from `red.csv`) using `replace=True`. This meant:

1. **Same birds repeatedly sampled**: Over 50 generations adding 10 birds/gen = 500 samples from only 140 birds
2. **No new genetic combinations**: After initial alleles were introduced, re-adding the same birds didn't increase diversity
3. **Unrealistic genetics**: Real captive populations breed each generation, creating new allele combinations through recombination
4. **Potential inbreeding underestimate**: Adding the same lineages repeatedly could increase relatedness in the wild population (not modeled)

### The Solution

Implemented a **Mendelian inheritance model** where the captive population breeds each generation:

```
Generation 0: Initialize from CSV (140 founders)
Generation 1: Breed → new offspring with recombined genetics
Generation 2: Breed → more new combinations
...
Generation N: Sample from CURRENT population (not original CSV)
```

### Key Components Added

| Component | Purpose |
|-----------|---------|
| `Bird` dataclass | Represents individual birds with genotype, sex, parents, inbreeding coefficient |
| `CaptivePopulation` dataclass | Manages breeding population with pedigree tracking |
| `CaptiveBreedingParams` dataclass | Configuration (breeding success rate, offspring per pair, etc.) |
| `breed_captive_population()` | Simulates one generation of random mating |
| `create_offspring()` | Mendelian inheritance - random allele from each parent per locus |
| `simulate_supplementation_effect_breeding()` | New simulation function using breeding model |

### Breeding Algorithm

```python
for each generation:
    1. Select breeding pairs (random pairing of males and females)
    2. For each pair with breeding success (70% chance):
       - Produce offspring (Poisson distributed, mean 1.5)
       - Each offspring inherits one random allele per locus from each parent
    3. Cull population to target size (remove oldest first)
    4. Sample supplementation birds from current population
    5. Add to wild population
```

### Why This Matters

- **New genetic combinations each generation**: Recombination creates novel genotypes even from the same founding alleles
- **Realistic captive population dynamics**: Population size maintained, inbreeding tracked
- **More conservative diversity estimates**: Can't infinitely re-use the same genetic material

---

## Change 2: Model 6 - Mixed South African Sourcing



### Model 6 Specification

| Source | Birds/Generation | Notes |
|--------|------------------|-------|
| PAAZA Captive | 6 | From breeding captive population |
| KwaZulu-Natal Wild | 2 | Translocated, depletes source |
| Eastern Cape Wild | 2 | Translocated, depletes source |
| **Total** | **10** | Same as Model 4 |

### Key Features

1. **Wild source depletion**: KZN and EC populations are sampled WITHOUT replacement
   - Realistically models removal from already-threatened populations
   - Source populations shrink over time
   - Model tracks `kzn_source_remaining` and `ec_source_remaining`

2. **Mixed genetic contribution**:
   - Captive birds: Managed breeding, lower inbreeding
   - Wild birds: Higher genetic diversity, but removes individuals from wild

3. **Practical conservation scenario**: Tests intra-South African gene flow without international logistics (unlike Model 5)

### Trade-offs This Model Reveals

| Advantage | Disadvantage |
|-----------|--------------|
| Higher genetic diversity from wild sources | Depletes already small wild populations |
| No international logistics | Removes breeding individuals from source |
| Wild birds may be better adapted | Limited source pool (~50 birds each) |


## Parameters

### Captive Breeding Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `breeding_success_rate` | 0.7 | 70% of pairs successfully breed |
| `offspring_per_pair` | 1.5 | Mean offspring per successful pair (Poisson) |
| `max_breeding_generations` | 5 | Birds can breed for ~5 generations |
| `captive_Ne` | 35-70 | Effective population size for captive pop |

### Model 6 Source Allocation

| Source | Initial Size | Birds/Gen | Exhausted After |
|--------|--------------|-----------|-----------------|
| PAAZA Captive | ~70 | 6 | Never (breeds) |
| KZN Wild | ~51 | 2 | ~25 generations |
| EC Wild | ~48 | 2 | ~24 generations |

---

## Verification

To test the changes:

1. **Run all 6 models** via the API with default parameters
2. **Compare Model 3/4 old vs new**: New breeding model should show slower diversity gain after initial generations (can't re-use same genetic material)
3. **Check Model 6 depletion**: `kzn_source_remaining` and `ec_source_remaining` should decrease each generation
4. **Verify Mendelian inheritance**: Offspring alleles should only contain alleles present in parents

---

## Future Considerations

1. **Minimum kinship pairing**: Could implement managed breeding to minimize inbreeding in captive population
2. **Wild population breeding**: Currently wild sources are static; could model their breeding too
3. **Sensitivity analysis**: Test different source allocations (e.g., 4+3+3 vs 6+2+2)
4. **Lambda adjustment**: Literature suggests λ ≈ 0.88 for unmanaged wild populations (currently defaults to 1.0)
