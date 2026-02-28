# Southern Ground Hornbill Population Viability Simulator

A population genetics and viability analysis (PVA) tool for Southern Ground Hornbills, using real microsatellite data from 199 wild and 145 captive individuals across 14 loci.

---

## The 6 Conservation Scenarios

| # | Name | Description |
|---|------|-------------|
| 1 | Baseline | All 4 wild populations, no intervention |
| 2 | Population Loss | EC + KZN go extinct; Kruger + Limpopo only |
| 3 | Low Supplementation | +4 SA captive birds/generation |
| 4 | High Supplementation | +10 SA captive birds/generation |
| 5 | International Mix | +4 birds/gen from SA + AZA + EAZA sources |
| 6 | Mixed SA Source | +6 captive + 2 KZN + 2 EC wild birds/gen (wild source depletes) |

Each generation = 26 years. Default Ne = 500, N₀ = 2500.

---

## Genetic Calculation Approach

### The problem this branch solves

Earlier versions ran Wright-Fisher simulation at Ne ≈ 199 (the CSV sample size, not the true population Ne) for Ho, He, and F in supplementation models. This made conservation scenarios appear *worse* than the no-intervention baseline — an artefact of simulating genetic drift at the wrong scale, not a biological signal.

### The fix

**Deterministic mode** — Models 3–6 now use analytical formulas throughout:

- **Ho / He**: Wright's heterozygosity loss at corrected Ne = 500 from the empirical initial values in the CSV.
- **F**: Drift-based accumulation `F(t) = 1 − (1 − 1/2Ne)^t`, with Ne adjusted upward by gene flow (`Ne_eff = base_Ne + cumulative_immigrants × 0.5`).
- **Na**: `calculate_allelic_diversity_with_geneflow()` — models drift loss of existing alleles *plus* cumulative retention of novel alleles introduced each generation. Novel allele rate is derived from actual CSV data (`compute_max_novel_alleles()`), not a fixed constant.

**Stochastic mode** — Models 3–6 keep Na analytical (same formula); Ho, He, and F come from a Wright-Fisher simulation at the *corrected* Ne, providing realistic variance without the scale error.

**Models 1–2** are unaffected: Ho/He/Na are always analytical (closed populations, no gene flow term).

### Novel allele calibration

`compute_max_novel_alleles()` in [models/base.py](models/base.py) computes mean novel alleles per locus from the real CSV genotypes for each captive source:

- PAAZA (SA captive, ~70 birds): **1.29 novel alleles/locus**
- International pool (~140 birds): **2.64 novel alleles/locus**

The introduction rate per generation scales as `birds_per_gen × max_novel / captive_pool_size` and is pre-computed in [api/index.py](api/index.py) before being threaded through to each model run.

---

## Monte Carlo Mode

Running N stochastic replicates reveals quasi-extinction risk that a single deterministic trajectory cannot.

- **Endpoint**: `POST /api/monte-carlo`
- **Engine**: `run_monte_carlo()` in [models/monte_carlo.py](models/monte_carlo.py)
- Each replicate uses `stochastic=True`:
  - Environmental stochasticity: λ += Normal(0, 0.10) per generation
  - Demographic stochasticity: N ~ Poisson(N·λ) for small N, Normal for large N
- Returns per-generation mean, std, and percentiles (p5/p25/p75/p95)
- Returns **quasi-extinction probability** P(N < threshold) per generation
- Frontend renders confidence bands via Plotly `fill='tonexty'` and an extinction probability chart panel

---

## Key Parameters

| Parameter | Range | Notes |
|-----------|-------|-------|
| Ne | 375–625 | Effective population size; increases dynamically with gene flow in models 3–6 |
| Generations | 10–100 | Each = 26 years |
| λ (lambda) | 0.5–1.5 | Base growth rate, independent of genetics |
| B | 3.14 (fixed) | Lethal equivalents for birds (O'Grady et al. 2006) |
| Hybrid vigor | 15% | Released birds experience 15% of wild inbreeding depression |
| Captive survival | 95% | Released birds survive at 95% of wild rate |

---

## Quick Start

```bash
cd pred
.venv/bin/python api/index.py
# → http://localhost:5001
```

---

## References

- O'Grady et al. (2006). *Biological Conservation* 133:42–51. — lethal equivalents B = 3.14 for birds
- Wright (1931). *Genetics* 16:97–159. — heterozygosity loss, drift
- Mills & Allendorf (1996). *Conservation Biology* 10:1509–1518. — one-migrant-per-generation rule
