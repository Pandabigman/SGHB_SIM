# Southern Ground Hornbill Population Genetics Simulator

Python Flask application with **microsatellite CSV data processing** for accurate genetic modeling.



### ✅ **All 5 Models Complete!**
1. ✅ Model 1: Baseline (all populations)
2. ✅ Model 2: Population loss (shows actual alleles lost)
3. ✅ Model 3: +4 PAAZA birds/gen (real genetic rescue)
4. ✅ Model 4: +10 PAAZA birds/gen (stronger rescue)
5. ✅ Model 5: +4 mixed birds/gen (maximum diversity)

### ✅ **CSV Data Processing**
- Parses actual microsatellite genotypes
- Calculates Ho, He, Na from raw alleles
- Tracks specific alleles gained/lost
- Models realistic supplementation effects

---



## Installation and Start

### 1. Install & Run
```bash
# Install dependencies
pip install -r requirements.txt

# Run application
python app.py
```

### 2. Open Browser
Navigate to: `http://localhost:5000`

---


## 🎯 What Each Model Does

### **Model 1: Baseline** ✅

 Calculates real Ho from 199 genotypes

```python
# Calculates from actual data
Ho = count_heterozygotes / total_individuals
He = 1 - sum(allele_freq²)
Na = mean_unique_alleles_per_locus
```

### **Model 2: Population Loss** ✅
**Without CSV:** Generic 1% Ho reduction
**With CSV:** Shows actual alleles lost

```python
# Removes EC + KZN individuals
wild_reduced = wild_df[~Site.isin(['Eastern Cape', 'KwaZulu-Natal'])]

# Identifies lost alleles
lost = {
  'Buco4': {165, 188},  # These alleles gone forever!
  'Buco11': {159},
  ...
}
```

### **Model 3: +4 PAAZA Birds** ✅
**Without CSV:** Generic Ne increase
**With CSV:** Adds 4 random PAAZA individuals each generation

```python
# Generation 0: Wild only
population = wild_df.copy()

# Generation 1: Add 4 PAAZA
birds = paaza_df.sample(4)
population = concat([population, birds])
calculate_metrics(population)  # Real Ho, He, Na

# Generation 2: Add 4 more
birds = paaza_df.sample(4)
population = concat([population, birds])
calculate_metrics(population)

# ... continues for all generations
```

**Result:** Shows REAL genetic rescue effect from actual alleles!

### **Model 4: +10 PAAZA Birds** ✅
Same as Model 3 but 10 birds per generation → Stronger effect

### **Model 5: +4 Mixed Birds** ✅
 Samples from PAAZA + AZA + EAZA

```python
# Mix all captive sources
mixed = concat([paaza_df, aza_df, eaza_df])

# Add 4 random from mixed pool
birds = mixed.sample(4)  # Maximum diversity!
```

**Result:** Best genetic outcome (most novel alleles)



## 🧮  Genetic Calculations

### **Real Observed Heterozygosity:**
```python
def calculate_ho(df, locus):
    genotypes = extract_genotypes(df, locus)
    heterozygotes = count(allele1 != allele2)
    return heterozygotes / total
```

### **Real Expected Heterozygosity:**
```python
def calculate_he(allele_frequencies):
    return 1 - sum(freq² for freq in frequencies)
```

### **Allelic Richness:**
```python
def calculate_na(df, locus):
    alleles = extract_all_alleles(df, locus)
    return len(set(alleles))
```

### **Supplementation Simulation:**
```python
def simulate_supplementation(wild, captive, n_birds, generations):
    population = wild.copy()
    results = []
    
    for gen in range(generations):
        # Calculate current metrics
        Ho = calculate_ho(population)
        He = calculate_he(population)
        Na = calculate_na(population)
        results.append({Ho, He, Na})
        
        # Add birds for next generation
        new_birds = captive.sample(n_birds)
        population = concat([population, new_birds])
    
    return results
```

---

## 🔬 What can be tracked

### 1. **Specific Alleles Lost** (Model 2)
```json
{
  "lost_alleles": {
    "Buco4": [165, 188],
    "Buco11": [159],
    "GHB21": [149, 161]
  },
  "total_lost": 5
}
```

### 2. **Novel Alleles Gained** (Models 3-5)
```json
{
  "novel_from_paaza": {
    "Buco2": [206],
    "GHB20": [183, 185]
  },
  "total_novel": 12
}
```

### 3. **Per-Locus Diversity**
```json
{
  "Buco4": {
    "Ho": 0.523,
    "He": 0.642,
    "Na": 8,
    "alleles": [162, 167, 180, 182, ...]
  }
}
```

---

## 📡 API Endpoints

### **POST `/api/simulate`**
Run any model (1-5)

```bash
curl -X POST http://localhost:5000/api/simulate \
  -H "Content-Type: application/json" \
  -d '{"Ne": 100, "generations": 50, "model": 3}'
```

**Response includes:**
- Genetic metrics per generation
- Data source (CSV vs default)
- Novel alleles count (for Models 3-5)
- Lost alleles count (for Model 2)

### **GET `/api/data/info`**
Check if CSV data loaded

```bash
curl http://localhost:5000/api/data/info
```
---


## ✅ Complete Feature List

### Backend
- ✅ All 5 models implemented
- ✅ Real CSV data parsing
- ✅ Per-locus genetic calculations
- ✅ Allele frequency tracking
- ✅ Supplementation simulation
- ✅ Allele loss identification
- ✅ Novel allele detection
- ✅ RESTful API

### Frontend
- ✅ Model selector tabs
- ✅ Interactive parameter sliders
- ✅ 4 Plotly charts (Ho, F, Na, N)
- ✅ Summary statistics
- ✅ Responsive design
- ✅ Error handling
- ✅ Loading states

### Data Processing
- ✅ 14 microsatellite loci
- ✅ Handle missing data (0 values)
- ✅ Calculate real Ho, He, Na
- ✅ Track specific alleles
- ✅ Realistic supplementation
- ✅ Generation-by-generation simulation

---

## 🧪 Testing

### Test All Models:
```bash
# Model 1
curl -X POST http://localhost:5000/api/simulate \
  -d '{"Ne": 100, "generations": 50, "model": 1}' \
  -H "Content-Type: application/json"

# Model 2
curl -X POST http://localhost:5000/api/simulate \
  -d '{"Ne": 100, "generations": 50, "model": 2}' \
  -H "Content-Type: application/json"

# Model 3
curl -X POST http://localhost:5000/api/simulate \
  -d '{"Ne": 100, "generations": 50, "model": 3}' \
  -H "Content-Type: application/json"

# Model 4
curl -X POST http://localhost:5000/api/simulate \
  -d '{"Ne": 100, "generations": 50, "model": 4}' \
  -H "Content-Type: application/json"

# Model 5
curl -X POST http://localhost:5000/api/simulate \
  -d '{"Ne": 100, "generations": 50, "model": 5}' \
  -H "Content-Type: application/json"
```

### Check Data Source:
```bash
curl http://localhost:5000/api/data/info | python -m json.tool
```

---

## 🎓 Scientific Basis

### Formulas Used:
- **Heterozygosity:** Wright (1931)
- **Inbreeding:** Nei et al. (1975)
- **Allelic diversity:** Kimura & Crow (1964)
- **Gene flow:** Mills & Allendorf (1996)

### CSV Data Provides:
- Real genotypes → Actual Ho calculation
- Allele frequencies → Accurate He calculation
- Specific alleles → Realistic drift simulation
- Individual birds → True supplementation effects

---
