# GHB_SIM — Formula Reference Table

All mathematical relationships used in the simulation, where they are applied, and what they drive.

---

## Genetic Diversity Metrics
*(Computed from real microsatellite data each generation — `genetics/metrics.py`)*

| Formula Name | Equation | Used For | Output |
|---|---|---|---|
| **Observed Heterozygosity** | `Ho = (1/L) * Σ_l [ proportion of heterozygous genotypes at locus l ]` | Measures actual genetic variation in sampled population | Ho chart line |
| **Expected Heterozygosity** (Nei) | `He = (1/L) * Σ_l [ 1 - Σ_i p_i² ]` where `p_i` = allele frequency | What Ho would be under random mating; measures allele diversity | He chart line |
| **Allelic Richness** | `Na = (1/L) * Σ_l [ number of unique alleles at locus l ]` | Number of distinct alleles per locus; irreversible loss metric | Na chart line |
| **FIS (Fixation Index)** | `FIS = (He - Ho) / He` | Departure from Hardy-Weinberg; proxy for recent inbreeding | Inbreeding chart line (Models 3–6 CSV) |

---

## Analytical Genetic Projections
*(Applied without real data — `models/base.py`, used by all `_run_generic` fallbacks and Models 1–2)*

| Formula Name | Equation | Used For | Source |
|---|---|---|---|
| **Wright's Heterozygosity Loss** | `H(t) = H₀ × (1 - 1/2Ne)^t` | Projects how Ho/He decline from drift alone over time | Wright 1931 |
| **Wright's Inbreeding Accumulation** | `F(t) = 1 - (1 - 1/2Ne)^t` | Projects inbreeding coefficient accumulation from drift | Wright 1931 |
| **Allelic Diversity — Drift Only** | `A(t) = A₀ × exp(-t / 4Ne)`, min 2 | Projects allele loss from drift; used in Models 1–2 and generic fallbacks | Maruyama & Fuerst 1985 |
| **Allelic Diversity — With Gene Flow** | `A(t) = A₀ × exp(-t/4Ne) + Σ_{g=1}^{t} [ novel_per_gen × exp(-(t-g)/4Ne) ]` | Drift loss of original alleles plus cumulative retention of novel alleles added each generation | Models 3–6 generic fallback |

---

## Population Dynamics
*(Applied to census-scale population — `models/base.py`, all models)*

| Formula Name | Equation | Used For | Notes |
|---|---|---|---|
| **Inbreeding Depression Survival** | `s(t) = exp(-B × F(t))` where `B = 3.14` lethal equivalents | Survival reduction factor applied to each generation's growth rate | O'Grady et al. 2006; `B = 3.14` is bird-specific average |
| **Population Growth with Inbreeding** | `N(t) = N(t-1) × λ × exp(-B × F(t-1))` | Projects census population size each generation accounting for lambda and inbreeding | Core demographic model; `N₀ = 2500` |
| **Stochastic Environmental Variation** | `λ_eff = λ_adjusted + Normal(0, σ=0.10)`, clipped at min 0.5 | Adds year-to-year variation in growth rate when stochastic mode is on | Only active when stochastic=True |
| **Stochastic Demographic Variation** | `N(t) ~ Poisson(N(t-1) × λ_eff)` for N < 1000, else Normal | Adds integer-level randomness for small populations | Only active when stochastic=True |

---

## Effective Population Size Adjustments
*(Applied during genetic simulation — `genetics/supplementation.py`)*

| Formula Name | Equation | Used For | Notes |
|---|---|---|---|
| **Ne with Gene Flow** | `Ne_eff(t) = Ne_base + (birds_per_gen × t × 0.5)` | Increases effective Ne as immigrants accumulate (each contributes ~0.5 to Ne) | Used to compute effective_Ne output; in old code also fed F_array |

---

## Supplementation & Released Bird Survival
*(Applied to released bird cohort — `models/base.py`)*

| Formula Name | Equation | Used For | Notes |
|---|---|---|---|
| **Released Bird Cohort Survival** | `survivors(t) = birds_per_gen × (survival_rate(t))^(t - release_gen)` | Tracks how many birds from each release cohort survive through time | Runs for every cohort released from gen 1 onwards |
| **Released Bird Survival Rate** | `survival_rate = exp(-B × F_avg × 0.15) × 0.95` | Survival per generation for captive-origin birds | 0.15 = hybrid vigour (only 15% of inbreeding penalty); 0.95 = captive health screening discount |
| **Total Released Population** | `N_released(t) = Σ_{g=1}^{t} survivors_from_cohort_g(t)` | All surviving released cohorts summed at each generation | Added to N_wild to get total population on chart |
| **Total Population** | `N(t) = N_wild(t) + N_released(t)` | Final census count shown on population chart | Released birds tracked separately due to different survival rules |

---

## Inbreeding F Sources by Model
*(Which F value drives the population survival formula `exp(-B × F)`)*

| Model | F Source | Method |
|---|---|---|
| **1 — Baseline** | Wright's formula: `F(t) = 1 - (1 - 1/2Ne)^t` | Analytical, Ne fixed |
| **2 — Population Loss** | Wright's formula with scaled Ne: `Ne_scaled = Ne × (N0/2500)` | Analytical, smaller Ne for reduced population |
| **3 — Low Supplementation (CSV)** | `delta-FIS(t) = max(FIS_sim(t) - FIS_sim(0), 0)` | From real genetic simulation |
| **4 — High Supplementation (CSV)** | `delta-FIS(t) = max(FIS_sim(t) - FIS_sim(0), 0)` | From real genetic simulation |
| **5 — International Mix (CSV)** | `delta-FIS(t) = max(FIS_sim(t) - FIS_sim(0), 0)` | From real genetic simulation |
| **6 — Mixed SA Source (CSV)** | `delta-FIS(t) = max(FIS_sim(t) - FIS_sim(0), 0)` | From real genetic simulation |
| **3–6 generic fallback** | Wright's formula + genetic rescue discount | Analytical approximation, no CSV |

---

## Captive Breeding (Mendelian Engine)
*(`genetics/simulation.py` — runs each generation for Models 3–6)*

| Formula Name | Equation / Rule | Used For |
|---|---|---|
| **Mendelian Inheritance** | For each locus: `offspring_allele = random choice from {parent_allele_1, parent_allele_2}` | Creates new recombinant genotypes each generation |
| **Offspring Inbreeding** | `F_offspring = min((F_sire + F_dam)/2 + 0.01, 1.0)` | Simplified inbreeding accumulation per generation in captive pedigree |
| **Breeding Success** | `Bernoulli(0.70)` per pair | Whether a given pair produces offspring in a given generation |
| **Clutch Size** | `n_offspring ~ Poisson(1.5)` per successful pair | Number of chicks produced per breeding pair per generation |
| **Breeding Eligibility** | `current_gen - bird.generation < max_breeding_generations (5)` | Birds older than 5 generations excluded from pairing |
| **Captive Cull** | Remove oldest birds if `len(all_birds) > target_population_size` | Maintains captive population at stable size |
| **Supplementation Selection** | `weight(bird) = 1 / (1 + F_bird)`, sample without replacement | Prefers least-inbred individuals for release |

---

## Allele Tracking
*(`genetics/population.py`)*

| Calculation | Definition | Used For |
|---|---|---|
| **Lost Alleles** | Alleles present in all-wild CSV but absent in Kruger+Limpopo subset | Quantifies genetic cost of EC+KZN extinction in Model 2 |
| **Novel Alleles** | Alleles present in captive CSV but absent in wild CSV | Quantifies potential genetic gain from supplementation; used in Na gene flow formula |

---

## Constants

| Symbol | Value | Meaning | Source |
|---|---|---|---|
| `B` | 3.14 | Lethal equivalents for birds | O'Grady et al. 2006 |
| `Ne` | 500 (default) | Effective population size | ~20% of census N = 2500 |
| `N₀` | 2500 | Census population size | Wild SGH estimate |
| `λ` | 1.0 (default) | Annual population growth rate | User-adjustable 0.95–1.05 |
| `T` | 26 years | Generation time | SGH biology |
| `σ_env` | 0.10 | Environmental stochasticity SD | Stochastic mode only |
| Hybrid vigour | 0.15 | Fraction of inbreeding penalty applied to F1 captive releases | Conservation genetics standard |
| Captive survival | 0.95 | Survival multiplier for captive-origin birds | Health screening advantage |
| Breeding success | 0.70 | Probability a pair breeds successfully per generation | PAAZA program estimate |
| Offspring per pair | Poisson(1.5) | Mean clutch size per successful pair | SGH reproductive biology |
