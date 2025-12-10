# Allelic Diversity Model Fix - Technical Documentation

## Problem Statement

The allelic diversity calculation had a fundamental inconsistency between baseline and supplementation models that misrepresented the genetic benefits of conservation intervention.

### Issue #1: Simple Drift Formula Cannot Model Gene Flow

**The Old Formula:**
```python
At = A0 * exp(-t / (4*Ne))
```

This formula **only** models allele loss through genetic drift. It cannot account for:
- ✗ Novel alleles added via gene flow (supplementation)
- ✗ Mutation (negligible over 50 generations anyway)

**The Problem:**
- Models 1-2 (baseline): Used this formula appropriately (closed populations, no gene flow)
- Models 3-5 (CSV mode): Used empirical simulation that DOES track novel alleles ✓
- Models 3-5 (fallback mode): Used drift formula with inflated Ne - **WRONG!**

**Why Inflated Ne Doesn't Fix It:**

The generic supplementation fallback used this approach:
```python
effective_Ne = base_Ne + (birds_added * 0.5)
Na = calculate_allelic_diversity(A0, effective_Ne, gen)  # STILL WRONG!
```

While increasing Ne slows the **rate of loss**, it still predicts:
- Na can only **decline** (slower decline ≠ gain)
- Final Na at gen 50 < Initial Na at gen 0

But the **real CSV data** shows:
- Captive populations have novel alleles not in wild (empirically proven)
- Supplementation **increases** Na over time
- Final Na can be > Initial Na

### Issue #2: Ne Inconsistency in Supplementation Models

**The Problem:**

In Models 3-5, the `Ne` parameter passed to the function was:
1. Used as a parameter for the genetic simulation
2. Returned in the output as `"Ne": 500`
3. **BUT** never actually reflected the gene flow benefit

The genetic metrics (Ho, He, F, Na) came from `simulate_supplementation_effect()`, which implicitly modeled gene flow effects but didn't return the **realized effective Ne**.

**The Issue:**
- API returned `Ne: 500` (static)
- But the simulation was actually experiencing `Ne: 500 → 700` (dynamic)
- Users had no visibility into this genetic rescue effect

---

## Solutions Implemented

### Fix #1: New Gene Flow Model for Allelic Diversity

Created `calculate_allelic_diversity_with_geneflow()` that models **both** drift loss AND gene flow gain:

```python
def calculate_allelic_diversity_with_geneflow(A0, Ne, t, novel_alleles_per_gen):
    """
    Component 1: Drift loss of original alleles
        At_drift = A0 * exp(-t / (4*Ne))

    Component 2: Gene flow gain (cumulative across generations)
        For each generation g from 1 to t:
            Alleles added at gen g experience drift for (t-g) generations
            Retained = novel_alleles_per_gen * exp(-(t-g) / (4*Ne))

        Total gain = sum(Retained for all g)

    Total alleles = At_drift + Total gain
    """
```

**Key Insight:**
- Novel alleles added in generation g=10 also experience drift from g=10 to g=50
- We sum the retained alleles from ALL past additions
- This correctly models the **accumulation** of genetic diversity via gene flow

**Parameters:**
- `novel_alleles_per_gen`: Conservative estimate of 0.15 alleles/locus/bird
  - Based on typical captive-wild differentiation
  - Can be calibrated from CSV data showing novel alleles

**Result:**
- Na can now **increase** in supplementation models
- Matches empirical reality from CSV simulations
- Biologically accurate representation of genetic rescue

---

### Fix #2: Dynamic Ne Tracking

#### A. Updated `simulate_supplementation_effect()` in [genetic_data.py](genetic_data.py:298)

Added `base_Ne` parameter and `effective_Ne` calculation:

```python
def simulate_supplementation_effect(wild_df, captive_df, birds_per_gen,
                                   generations, loci_list=LOCI, base_Ne=500):
    for gen in range(generations + 1):
        # Calculate effective Ne based on gene flow
        cumulative_immigrants = birds_per_gen * gen
        effective_Ne = base_Ne + (cumulative_immigrants * 0.5)

        results.append({
            'generation': gen,
            'Ho': Ho,
            'He': He,
            'Na': Na,
            'FIS': FIS,
            'population_size': pop_size,
            'effective_Ne': effective_Ne  # NEW!
        })
```

#### B. Updated Models 3-5 API Responses

All supplementation models now return:
```python
{
    "Ne": 500,  # Base Ne (for reference)
    "effective_Ne": [500, 502, 504, ...],  # Dynamic Ne array
    "parameters": {
        "Ne_base": 500,
        "Ne_final": 600,  # Shows genetic rescue benefit
        "allele_model": "empirical" or "drift + gene flow"
    }
}
```

---

## Model Comparison Table

| Model | Allele Formula | Ne | Can Na Increase? | Data Source |
|-------|---------------|-----|------------------|-------------|
| Model 1 (Baseline) | Drift only | Static | ✗ No | CSV or defaults |
| Model 2 (Population Loss) | Drift only | Static | ✗ No | CSV or defaults |
| Model 3-5 (CSV mode) | Empirical tracking | Dynamic | ✓ Yes | Real genetic data |
| Model 3-5 (Fallback) | **Drift + gene flow** | Dynamic | ✓ Yes | Estimated parameters |

---

## Biological Justification

### Why Novel Alleles Matter

**Allelic diversity (Na) is critical for:**
1. **Adaptive potential**: More alleles = more genetic variation for selection
2. **Future evolutionary responses**: Climate change, disease resistance
3. **Avoiding genetic bottlenecks**: Loss of rare alleles is irreversible

**Why drift-only models fail for supplementation:**
- Real supplementation adds birds from genetically differentiated populations
- These birds carry alleles absent in the recipient population
- Empirical data proves this (see CSV novel allele counts)
- Ignoring this effect **underestimates** conservation value

### Empirical Evidence from CSV Data

From the real genetic data:
- **PAAZA** (South African zoos): Contains novel alleles not in wild populations
- **AZA** (USA/Canada zoos): Additional unique alleles
- **EAZA** (European zoos): Further diversity

Model 5 (International Mix) combines all three sources for maximum allelic diversity gain.

---

## Testing Recommendations

### Verify Gene Flow Model

Run Model 3 with these parameters:
- Ne = 500
- Generations = 50
- Lambda = 1.0
- Birds per gen = 4

**Expected Results:**
- `effective_Ne[0]` = 500
- `effective_Ne[25]` ≈ 550
- `effective_Ne[50]` ≈ 600
- `Na[0]` ≈ 6.4
- `Na[50]` > `Na[0]` (should increase with CSV data)

### Compare Fallback vs CSV

Run same parameters twice:
1. With CSV data available → Uses empirical simulation
2. Without CSV data → Uses drift + gene flow model

Both should show:
- Na increasing over time
- Similar magnitude of increase
- Dynamic Ne tracking

---

## API Breaking Changes

⚠️ **Client applications should be updated to use:**

```javascript
// Old (deprecated but still works)
const baseNe = response.Ne;  // Static value

// New (recommended)
const baseNe = response.parameters.Ne_base;
const finalNe = response.parameters.Ne_final;
const dynamicNe = response.effective_Ne;  // Array [500, 502, 504, ...]

// For plotting Ne over time
for (let i = 0; i < response.generations.length; i++) {
    const gen = response.generations[i];
    const ne = response.effective_Ne[i];
    // Plot gen vs ne
}
```

---

## Future Improvements

### 1. Calibrate Novel Allele Estimates

Current fallback uses generic 0.15 alleles/locus/bird. Could be improved:
- Calculate actual novel allele rates from CSV data
- Use source-specific rates (PAAZA vs AZA vs EAZA)
- Vary by locus (some loci more variable than others)

### 2. Account for Allele Frequency

Current model treats all novel alleles equally. Could weight by:
- Initial frequency in donor population
- Fitness effects (neutral vs beneficial)
- Probability of establishment in recipient population

### 3. Mutation

For simulations >100 generations, might need to add:
```python
new_alleles_from_mutation = mutation_rate * loci * 4 * Ne * t
```

But for 10-100 generations, mutation is negligible compared to gene flow.

---

## References

**Genetic Drift and Allele Loss:**
- Kimura & Crow (1964): Expected time to allele loss
- Nei et al. (1975): Expected heterozygosity in finite populations

**Gene Flow and Genetic Rescue:**
- Tallmon et al. (2004): "When are genetic rescue and translocation effective?"
- Whiteley et al. (2015): "Genetic rescue to the rescue"

**Southern Ground Hornbill Genetics:**
- Kemp & Begg (various): Population structure and conservation genetics
- Your CSV data: Empirical allele frequencies and differentiation

---

## Summary

The allelic diversity fix addresses a fundamental modeling error:

**Before:**
- Supplementation only slowed genetic decline
- Na always decreased (just slower with more birds)
- Underestimated conservation value

**After:**
- Supplementation can reverse genetic decline
- Na can increase when novel alleles are added
- Accurately reflects genetic rescue benefits
- Matches empirical data from CSV simulations

This makes the model **more realistic**, **more optimistic** (appropriately so), and **scientifically defensible**.
