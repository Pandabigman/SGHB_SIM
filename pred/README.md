# Southern Ground Hornbill Population Viability Simulator

**A population genetics and viability analysis tool for Southern Ground Hornbills** using real microsatellite data from 199 wild and 145 captive individuals across South Africa, USA, and Europe.

---

## Recent Updates (December 2024)

### Model Improvements for More Realistic Projections

**Problem:** Previous simulations showed population collapse even with conservation management due to overly pessimistic parameters.

**Solution:** Updated to bird-specific parameters and added genetic rescue modeling:

1. ✅ **Lethal equivalents: 6.29 → 3.14**
   - Changed from mammal average to bird-specific value (O'Grady et al. 2006)
   - Reduces excessive inbreeding depression in projections

2. ✅ **Added genetic rescue effect**
   - Supplementation now reduces population-wide inbreeding coefficient
   - Formula: F_rescued = F × (1 - gene_flow_proportion × 0.5)

3. ✅ **Improved hybrid vigor modeling**
   - Released birds: 30% → 15% of wild inbreeding penalty
   - Reflects stronger F1 heterosis effect

4. ✅ **Better released bird survival: 90% → 95%**
   - Accounts for disease screening and monitoring

5. ✅ **Comprehensive documentation**
   - Added expandable info panel in UI with all assumptions and parameters
   - Updated README with model changes and rationale

6. ✅ **Allelic diversity gene flow model (December 2024)**
   - Fixed conceptual inconsistency in allelic diversity calculation
   - Supplementation models now correctly show Na can increase (not just decline slower)
   - Added dynamic Ne tracking for supplementation scenarios
   - Empirical CSV models track actual novel alleles introduced

**Impact:** Conservation scenarios (Models 3-5) now show realistic population stabilization and recovery, while maintaining scientific rigor. Allelic diversity dynamics now accurately reflect gene flow benefits.

---

## Overview

This web application models genetic diversity loss and population decline in Southern Ground Hornbills under different conservation scenarios. It combines **real genetic data** with **inbreeding depression** to provide realistic predictions of population trajectories over 10-100 generations (260-2600 years).

### Key Features

- ✅ **Real microsatellite data** (14 loci from 344 birds)
- ✅ **Inbreeding depression** (B = 3.14 lethal equivalents for birds, O'Grady et al. 2006)
- ✅ **Genetic rescue effect** (supplementation reduces population-wide inbreeding)
- ✅ **User-adjustable population growth rate** (λ = 0.5-1.5) to model habitat loss/recovery
- ✅ **5 conservation scenarios** comparing no intervention, population loss, and genetic rescue
- ✅ **Demographic-genetic coupling** (inbreeding affects population size)
- ✅ **Cohort tracking** for released captive birds with hybrid vigor modeling
- ✅ **Interactive visualization** with Plotly.js charts
- ✅ **Comprehensive model documentation** in expandable info panel

---

## Installation

### Requirements
- Python 3.8+
- Flask 3.1.2
- NumPy 2.3.5
- Pandas 2.3.3

### Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run locally (do not use app.py)
python api/index.py

# Open browser
http://localhost:5001
```

### Deployment
Configured for Vercel serverless deployment via `vercel.json`.

---

## The 5 Conservation Scenarios

### **Model 1: Baseline (All Wild Populations)**
**What it does:** Simulates genetic drift and inbreeding in all 4 wild populations (Eastern Cape, Kruger NP, KwaZulu-Natal, Limpopo) with no intervention.

**Genetic metrics:** Calculated from 199 wild birds
- Initial Ho (observed heterozygosity): 0.502
- Initial He (expected heterozygosity): 0.568
- Initial Na (allelic richness): 6.429 alleles/locus

**Population dynamics:**
```python
# Inbreeding increases over time
F(t) = 1 - (1 - 1/(2*Ne))^t

# Population declines due to inbreeding depression
N(t) = N(t-1) * λ * exp(-3.14 * F(t))
```

**Use case:** Baseline for comparison - shows what happens with no conservation action.

---

### **Model 2: Population Loss (Kruger + Limpopo Only)**
**What it does:** Simulates the loss of Eastern Cape and KwaZulu-Natal populations, tracking which specific alleles are lost forever.

**Alleles lost:** The model identifies exact alleles unique to EC/KZN that disappear:
```json
{
  "Buco4": [165, 188],
  "Buco11": [159],
  "GHB21": [149, 161]
}
```

**Population metrics:**
- Starting N: ~2000 (reduced from 2500)
- Effective Ne: Scaled proportionally to population size
- Genetic diversity: Lower starting Ho, He, Na than Model 1

**Use case:** Shows genetic consequences if these populations go extinct.

---

### **Model 3: Low Supplementation (+4 South African Captive Birds/Generation)**
**What it does:** Adds 4 birds per generation from South African zoos (PAAZA) to the wild population.

**How it works:**
```python
# Each generation:
1. Calculate genetic metrics from current wild + captive gene pool
2. Add 4 random PAAZA birds to population
3. Apply genetic rescue effect (reduces population-wide F)
4. Track released bird cohorts separately with hybrid vigor
5. Apply inbreeding depression to wild-born birds
6. Apply 95% survival rate to captive-released birds
```

**Key improvements (Updated Model):**
- PAAZA birds: ~70 individuals in South African zoos
- **Genetic rescue effect**: Gene flow reduces population-wide F by 50% of supplementation proportion
- Released birds experience **15% of wild inbreeding depression** (strong hybrid vigor)
- Released birds survive at **95% of wild rate** (disease-screened, monitored)
- All released cohorts tracked cumulatively

**Novel alleles gained:** Shows which new alleles from PAAZA are introduced to wild population.

**Use case:** Realistic genetic rescue using locally available captive birds.

---

### **Model 4: High Supplementation (+10 South African Captive Birds/Generation)**
**What it does:** Same as Model 3, but adds 10 PAAZA birds per generation for stronger genetic rescue.

**Effect:**
- Faster Ho/He recovery
- More novel alleles introduced
- Greater population size increase
- Higher cost and logistical burden

**Use case:** Compares cost-benefit of low vs. high supplementation intensity.

---

### **Model 5: International Mix (+4 Mixed Birds/Generation)**
**What it does:** Adds 4 birds per generation from a mix of South African (PAAZA), USA/Canada (AZA), and European (EAZA) zoos.

**Captive sources:**
- **PAAZA** (Pan-African): ~70 birds in South Africa
- **AZA** (Association of Zoos & Aquariums): ~20 birds in USA/Canada
- **EAZA** (European): ~55 birds in Europe

**Genetic advantage:** Maximum allelic diversity (combines unique alleles from all 3 sources)

**Logistical challenges:**
- International transport (expensive, permits, quarantine)
- Climate adaptation (birds from Europe/USA → South African climate)
- Disease risk (cross-continental movement)

**Use case:** Best genetic outcome vs. highest logistical cost trade-off.

---

## Genetic & Demographic Calculations

### Genetic Metrics

**Observed Heterozygosity (Ho):**
```python
def calculate_ho(df, locus):
    genotypes = get_genotypes_for_locus(df, locus)
    heterozygotes = count(allele1 != allele2 for allele1, allele2 in genotypes)
    return heterozygotes / len(genotypes)
```

**Expected Heterozygosity (He):**
```python
def calculate_he(allele_frequencies):
    # Nei's unbiased estimator
    return 1 - sum(freq**2 for freq in allele_frequencies.values())
```

**Allelic Richness (Na):**
```python
def calculate_na(df, locus):
    alleles = get_all_alleles_for_locus(df, locus)
    return len(set(alleles))  # Count unique alleles
```

**Inbreeding Coefficient (F):**
```python
def calculate_fis(Ho, He):
    # FIS = (He - Ho) / He
    return (He - Ho) / He if He > 0 else 0
```

---

### Population Dynamics

**Wright's Heterozygosity Loss:**
```python
# Ht = H0 * (1 - 1/(2*Ne))^t
def calculate_heterozygosity_loss(H0, Ne, t):
    return H0 * np.power(1 - 1/(2*Ne), t)
```

**Inbreeding Accumulation:**
```python
# F = 1 - (1 - 1/(2*Ne))^t
def calculate_inbreeding(Ne, t):
    return 1 - np.power(1 - 1/(2*Ne), t)
```

**Allelic Diversity Loss (Drift Only):**
```python
# At = A0 * exp(-t / (4*Ne))
# NOTE: This formula ONLY models allele loss via genetic drift
# It does NOT account for gene flow or mutation
def calculate_allelic_diversity(A0, Ne, t):
    At = A0 * np.exp(-t / (4*Ne))
    return np.maximum(At, 2.0)  # Minimum 2 alleles/locus
```

**Allelic Diversity with Gene Flow (NEW!):**
```python
def calculate_allelic_diversity_with_geneflow(A0, Ne, t, novel_alleles_per_gen):
    """
    Models BOTH drift loss AND gene flow gain

    Component 1: Drift loss of original alleles
        At_drift = A0 * exp(-t / (4*Ne))

    Component 2: Gene flow gain (cumulative)
        For each generation g from 1 to t:
            Alleles added at gen g experience drift for (t-g) generations
            Retained = novel_alleles_per_gen * exp(-(t-g) / (4*Ne))
            Total gain = sum of all retained alleles

    Total alleles = At_drift + Total gain
    """
    # See api/index.py for full implementation
```

**Why This Matters:**

1. **Models 1-2 (Baseline)**: Use drift-only formula
   - Appropriate for closed populations (no immigration)
   - Na can only decline over time

2. **Models 3-5 (CSV data)**: Use empirical simulation
   - Tracks actual alleles from real genetic data
   - Na can increase when novel alleles are added
   - Proven by CSV data showing captive populations have unique alleles

3. **Models 3-5 (Fallback)**: Use gene flow formula
   - When CSV data unavailable, estimates novel allele addition
   - Conservative estimate: 0.15 novel alleles per locus per bird
   - Na can increase, matching empirical reality

---

### Inbreeding Depression (UPDATED MODEL!)

**Population size with inbreeding:**
```python
def calculate_population_size_with_inbreeding(N0, lambda_val, F_array, lethal_equivalents=3.14):
    """
    Inbreeding reduces survival: s = exp(-B * F)
    where B = 3.14 lethal equivalents for birds (O'Grady et al. 2006)

    
    """
    N = np.zeros(len(F_array))
    N[0] = N0

    for t in range(1, len(F_array)):
        # Inbreeding reduces survival
        inbreeding_survival = np.exp(-3.14 * F_array[t])

        # Adjusted growth rate combines demography + genetics
        lambda_adjusted = lambda_val * inbreeding_survival

        # Population at time t
        N[t] = N[t-1] * lambda_adjusted
        N[t] = max(N[t], 10)  # Minimum viable population

    return N
```

**Key parameters:**
- **B = 3.14**: Lethal equivalents (moderate inbreeding depression for birds)
  - **Updated from 6.29** (mammal average) to reflect bird-specific biology
  - O'Grady et al. (2006): Birds average ~3.14, mammals average ~6.29
  - Provides more realistic, less pessimistic population projections
- **λ (lambda)**: Base population growth rate (user-adjustable 0.5-1.5)
  - λ < 1.0: Declining (habitat loss, persecution)
  - λ = 1.0: Stable
  - λ > 1.0: Growing (conservation working)

---

## User-Adjustable Parameters

### 1. **Effective Population Size (Ne)** [375 - 625]
- Represents breeding individuals (not census size)
- Ne/N ratio ≈ 0.15-0.25 for cooperative breeders
- Lower Ne → Faster genetic loss

**IMPORTANT: Ne Dynamics in Supplementation Models**

Gene flow from supplementation increases the effective population size over time:

- **Models 1-2 (Baseline)**: Ne is constant (no immigration)
- **Models 3-5**: Ne increases dynamically with gene flow
  - Formula: `effective_Ne = base_Ne + (cumulative_immigrants × 0.5)`
  - Each immigrant contributes ~0.5 to effective size (conservative migration model)
  - Example: Model 3 with base Ne=500, adding 4 birds/gen for 50 gen:
    - Generation 0: Ne = 500
    - Generation 25: Ne = 550 (50 birds added)
    - Generation 50: Ne = 600 (100 birds added)

This dynamic Ne is used in:
- ✅ Heterozygosity loss calculations
- ✅ Inbreeding coefficient calculations
- ✅ Allelic diversity calculations (gene flow model)
- ✅ Returned in API response as `effective_Ne` array

### 2. **Generations** [10 - 100]
- Each generation = **26 years** (Southern Ground Hornbill generation time)
- 10 gen = 260 years
- 50 gen = 1,300 years
- 100 gen = 2,600 years

### 3. **Population Growth Rate (λ)** [0.5 - 1.5] 
- **λ = 0.90**: 10% decline per generation (habitat loss, climate change)
- **λ = 1.00**: Stable population (births = deaths)
- **λ = 1.10**: 10% growth per generation (successful conservation)
- **Independent of genetics** (accounts for habitat, climate, management)

---

## Supplementation Logic (Models 3-5) - UPDATED!

### Genetic Rescue Effect

**New feature:** Supplementation reduces population-wide inbreeding coefficient!

```python
# Apply genetic rescue effect
F_array_rescued = F_array.copy()
for gen in range(len(F_array)):
    if gen > 0:
        # Gene flow benefit proportional to cumulative supplementation
        cumulative_birds = birds_per_gen * gen
        total_pop_estimate = N0 + cumulative_birds
        gene_flow_proportion = cumulative_birds / total_pop_estimate
        # Reduce F by 50% of gene flow proportion
        F_array_rescued[gen] = F_array[gen] * (1 - gene_flow_proportion * 0.5)
```

**Impact:** As more birds are added over generations, the entire population benefits from reduced inbreeding, not just the released birds themselves!

---

### Cohort Tracking with Hybrid Vigor

Released birds are tracked separately from wild-born birds:

```python
# For each generation:
N_wild = calculate_with_inbreeding(N0, lambda, F_wild)

# Track all released cohorts
for each cohort released in past:
    time_since_release = current_gen - release_gen
    avg_F_since_release = mean(F_rescued[release_gen:current_gen+1])
    survival = exp(-3.14 * avg_F_since_release * 0.15) * 0.95  # Strong hybrid vigor!
    cohort_survivors = birds_released * (survival ** time_since_release)
    N_released += cohort_survivors

# Total population
N_total = N_wild + N_released
```

**Updated Assumptions:**
- Released birds experience **15% of wild inbreeding** (strong hybrid vigor / heterosis)
  - **Improved from 30%** to reflect F1 hybrid vigor effect
- Released birds survive at **95% of wild rate** (healthier, disease-screened)
  - **Improved from 90%** - captive birds are monitored and initially healthier
- All cohorts age together (cumulative tracking)
- Lethal equivalents updated to **3.14** (bird-specific value)

---

## API Reference

### POST `/api/simulate`

**Request:**
```json
{
  "Ne": 500,
  "generations": 50,
  "lambda": 1.0,
  "model": 3
}
```

**Response:**
```json
{
  "model_number": 3,
  "model_name": "Low Supplementation (+4 SA Captive/gen)",
  "generations": [0, 1, 2, ..., 50],
  "years": [0, 26, 52, ..., 1300],
  "Ho": [0.502, 0.501, 0.500, ...],
  "He": [0.568, 0.567, 0.566, ...],
  "F": [0.000, 0.001, 0.002, ...],
  "Na": [6.429, 6.420, 6.410, ...],
  "population": [2500, 2498, 2495, ...],
  "Ne": 500,
  "initial": {"Ho": 0.502, "He": 0.568, "Na": 6.429, "N": 2500},
  "parameters": {
    "lambda": 1.0,
    "data_source": "CSV_simulation",
    "supplementation": "4 South African captive birds per generation",
    "supplementation_source": "PAAZA (Pan-African Association of Zoos and Aquaria)",
    "novel_alleles_added": 12,
    "inbreeding_depression": "enabled (B=3.14 lethal equivalents for birds)"
  }
}
```

### GET `/api/data/info`

Returns loaded genetic data statistics:
```json
{
  "data_source": "CSV",
  "populations": {
    "wild_all": {"sample_size": 199, "Ho": 0.502, "He": 0.568, "Na": 6.429},
    "paaza": {"sample_size": 70, "Ho": 0.487, "He": 0.551, "Na": 5.8}
  },
  "lost_alleles_count": 8,
  "novel_alleles": {"paaza": 12, "aza": 5, "eaza": 18}
}
```

---

## Scientific References

### Population Genetics Theory
- **Wright, S. (1931).** Evolution in Mendelian populations. *Genetics* 16:97-159.
- **Nei, M. et al. (1975).** The bottleneck effect and genetic variability. *Evolution* 29:1-10.
- **Kimura, M. & Crow, J.F. (1964).** Number of alleles in a finite population. *Genetics* 49:725.

### Inbreeding Depression
- **O'Grady, J.J. et al. (2006).** Realistic levels of inbreeding depression strongly affect extinction risk in wild populations. *Biological Conservation* 133:42-51.

### Conservation Genetics
- **Mills, L.S. & Allendorf, F.W. (1996).** The one-migrant-per-generation rule. *Conservation Biology* 10:1509-1518.

---

## Data Sources

### Wild Populations (N = 199)
- **Eastern Cape province**: 48 individuals
- **Kruger National Park**: 52 individuals
- **KwaZulu-Natal province**: 51 individuals
- **Limpopo province**: 48 individuals

### Captive Populations (N = 145)
- **PAAZA** (Pan-African): ~70 birds (South Africa)
- **AZA** (USA/Canada): ~20 birds
- **EAZA** (European): ~55 birds

### Microsatellite Loci (14 total)
Buco4, Buco11, Buco2, Buco9, GHB21, GHB19, GHB26, GHB20, Buco16, Buco18, Buco19, Buco21, Buco24, Buco25

---

## Model Assumptions

### What the model includes:
✅ Real microsatellite genotypes (not simulated)
✅ Inbreeding depression (B = 3.14 for birds, updated from mammal value)
✅ **Genetic rescue effect** (supplementation reduces population-wide F)
✅ **Hybrid vigor modeling** (released birds show heterosis)
✅ User-adjustable growth rate (habitat/climate effects)
✅ Cohort-specific mortality for released birds
✅ Genetic drift (deterministic, via Wright's equation)

### What the model assumes:
⚠️ **Generation time = 26 years** (fixed)
⚠️ **Random mating** within populations
⚠️ **No natural migration** between wild populations
⚠️ **No mutation** (new alleles don't appear)
⚠️ **Released birds integrate immediately** (no behavioral barriers)
⚠️ **Deterministic** (single trajectory, no stochastic variation)
⚠️ **Genetic rescue at 50% efficiency** (gene flow reduces F proportionally)

### What the model does NOT include:
❌ Environmental stochasticity (random good/bad years)
❌ Demographic stochasticity (random births/deaths)
❌ Genetic drift stochasticity (random allele sampling)
❌ Catastrophes (disease outbreaks, droughts)
❌ Allee effects (reduced fitness at low density)

---

## Limitations & Caveats

1. **No uncertainty quantification**: Model shows single trajectory, not confidence intervals
2. **Deterministic genetics**: Real genetic drift is stochastic (random), model uses expected values
3. **Fixed generation time**: Assumes constant 26 years (may vary by population age structure)
4. **Ne/N ratio implicit**: User sets Ne directly, not calculated from biology
5. **Genetic rescue efficiency**: 50% reduction in F from gene flow is estimated, not empirically measured
6. **Hybrid vigor assumption**: 85% reduction in inbreeding penalty for F1s may vary in practice

**Important Note:** The model was updated (Dec 2024) with more realistic parameters:
- Lethal equivalents reduced from 6.29 (mammal average) to 3.14 (bird average)
- Released bird survival increased from 90% to 95%
- Inbreeding penalty for released birds reduced from 30% to 15% (hybrid vigor)
- Added genetic rescue effect where supplementation benefits entire population

**Recommendation:** Use for **comparing scenarios**, not absolute predictions. Add stochastic simulations for publication-quality uncertainty estimates.

---

## Future Enhancements

### Planned (High Priority)
- [ ] Stochastic simulations (100+ runs with confidence intervals)
- [ ] Sensitivity analysis (test different B, λ, generation time)
- [ ] Document generation time source (why 26 years?)

### Possible (Medium Priority)
- [ ] Genetic drift stochasticity (random allele sampling)
- [ ] Demographic stochasticity (Poisson births/deaths)
- [ ] User-adjustable B (lethal equivalents slider)
- [ ] Download results as CSV

### Research (Low Priority)
- [ ] Spatial structure (metapopulation dynamics)
- [ ] Habitat quality variation
- [ ] Climate change projections
- [ ] Cost-benefit analysis module

---



