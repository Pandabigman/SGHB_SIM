# Southern Ground Hornbill Population Genetics Simulator

A React-based application for modeling genetic diversity and population viability of Southern Ground Hornbills across 5 conservation scenarios.





## 🔧 Installation

```bash
# Install dependencies
npm install react recharts lucide-react

# Or with yarn
yarn add react recharts lucide-react
```



## 📊 Models

1. **Model 1**: Baseline (all wild populations)
2. **Model 2**: Population loss (Kruger + Limpopo only)
3. **Model 3**: Low supplementation (+4 PAAZA birds/generation)
4. **Model 4**: High supplementation (+10 PAAZA birds/generation)
5. **Model 5**: International mix (+4 PAAZA/AZA/EAZA birds/generation)

## 🧬 Genetic Equations

### Heterozygosity Loss (Wright, 1931)
```
Ht = H0 × (1 - 1/(2Ne))^t
```

### Inbreeding Accumulation
```
F = 1 - (1 - 1/(2Ne))^t
```

### Allelic Diversity (Nei et al., 1975)
```
At = A0 × exp(-t/(4Ne))
```

## 🎲 Deterministic vs Stochastic Simulations

### Deterministic Mode (Default)
- Uses **pure mathematical equations**
- Same inputs → **always same outputs**
- No randomness
- Good for: comparing scenarios, understanding trends

### Stochastic Mode (Optional Checkbox)
- Adds **random variation** to simulate real-world uncertainty
- Each run produces **different results**
- Includes:
  - Demographic stochasticity (±5% heterozygosity variation)
  - Environmental variation (±15% population fluctuation)
  - Random allele loss (±10% variation)
  - Inbreeding variation (±3%)


 By default, the model is deterministic (like a calculator). Enable stochastic mode to see variation between runs.

- **Ne Range**: 20-200 (default: 100)
- **Generations**: 10-100 (default: 50)
- **Generation Time**: 15 years
- **Wild Population**: ~450 birds
  
📐 Leslie Matrix Demographic Model
Age-Structured Population Growth
The model uses a proper Leslie matrix to calculate population growth rate (λ):

When To Use:

Default (λ = 1.0): Assumes stable population, focuses on genetic drift
Leslie Matrix: Calculates realistic λ from age structure (~0.97-1.02)

Result: More realistic population projections that account for age structure and delayed breeding.

Ne Range: 20-200 (default: 100)
Generations: 10-100 (default: 50)
Generation Time: 15 years
Wild Population: ~450 birds

## 📚 Key References

1. Wright, S. (1931). Evolution in Mendelian populations. *Genetics* 16(97)
2. Nei, M. et al. (1975). The bottleneck effect and genetic variability. *Evolution* 29(1)
3. Mills, L.S. & Allendorf, F.W. (1996). The one-migrant-per-generation rule. *Conservation Biology* 10(6)
4. Kemp, A.C. & Webster, P.J. (2008). Family Bucerotidae (hornbills). *Handbook of Birds of the World*
5. Taylor, M.R. (2015). The Eskom Red Data Book of Birds of South Africa

##  Customization


### Adjust Demographics
Edit `demographics` in `constants.js`:
```javascript
export const demographics = {
  breedingAge: 10,
  generationTime: 15,
  // ...
};
```

### Add New Models
1. Add model description in `getModelName()` in `constants.js`
2. Add switch case in `runSimulation()` in main component
3. Add color to `modelColors`

