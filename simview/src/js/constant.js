// Genetic data from CSV files - initial heterozygosity values
export const initialGenetics = {
  wild: {
    Ho: 0.502,
    He: 0.568,
    Na: 6.429,
    Ne_effective: 3.143,
    FIS: 0.121,
    sampleSize: 199, // SAE(8) + SAG(154) + SAK(14) + SAL(23)
    censusSize: 450 // Estimated total wild population
  },
  wildNoECKZN: {
    // Model 2: Only Kruger (154) + Limpopo (23) = 177 samples
    // Losing EC + KZN means losing genetic diversity
    Ho: 0.498, // Slightly lower due to loss of population diversity
    He: 0.565, // Expected heterozygosity also reduced
    Na: 6.1, // Fewer alleles (lost some unique alleles from EC/KZN)
    sampleSize: 177,
    censusSize: 398 // ~52 fewer birds (proportional to samples lost: 22/199 * 450)
  },
  paaza: {
    Ho: 0.483,
    He: 0.565,
    Na: 7.071,
    Ne_effective: 3.201,
    FIS: 0.139,
    sampleSize: 70
  },
  aza: {
    Ho: 0.515,
    He: 0.534,
    Na: 5.429,
    Ne_effective: 2.666,
    FIS: 0.051,
    sampleSize: 29
  },
  eaza: {
    Ho: 0.564,
    He: 0.582,
    Na: 5.571,
    Ne_effective: 3.100,
    FIS: 0.041,
    sampleSize: 46
  }
};

// Demographic parameters based on provided breeding information
export const demographics = {
  breedingAge: 10, // years
  generationTime: 15, // years
  breedingProportion: 0.20, // 20% of population breeds
  chicksPerBreeding: 1 / 3, // 1 chick every 3 years
  groupSize: 5.5, // Average 3-9, median ~5.5
  survivalAdult: 0.92, // Assumed from literature (high adult survival)
  survivalJuvenile: 0.65, // Assumed (lower juvenile survival)
  survivalSubadult: 0.85, // Assumed intermediate
  
  // Age class structure for Leslie matrix
  ageClasses: [
    { name: 'Juvenile', minAge: 0, maxAge: 4, survival: 0.65, fecundity: 0 },
    { name: 'Subadult', minAge: 5, maxAge: 9, survival: 0.85, fecundity: 0 },
    { name: 'Adult', minAge: 10, maxAge: 40, survival: 0.92, fecundity: 0.0667 }
    // Adult fecundity: 20% breed × (1 chick/3 years) = 0.20 × 0.333 = 0.0667
  ]
};

// Model colors for visualization
export const modelColors = {
  model1: '#2563eb', // Blue
  model2: '#dc2626', // Red
  model3: '#16a34a', // Green
  model4: '#9333ea', // Purple
  model5: '#ea580c'  // Orange
};

// Model names and descriptions
export const getModelName = (modelNum) => {
  const names = {
    1: "Baseline (All Wild Populations)",
    2: "Population Loss (Kruger + Limpopo Only)",
    3: "Low Supplementation (+4 PAAZA/gen)",
    4: "High Supplementation (+10 PAAZA/gen)",
    5: "International Supplementation (+4 Mixed/gen)"
  };
  return names[modelNum];
};