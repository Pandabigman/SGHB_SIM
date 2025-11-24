import { demographics } from '../js/constant';

/**
 * Calculate population growth rate (lambda)
 * Based on demographic parameters
 * @returns {number} Lambda (population growth rate)
 */
export const calculateLambda = () => {
  // Simplified Leslie matrix approach
  // With 20% breeding, 1 chick/3 years = 0.333 * 0.20 = 0.0667 offspring per capita per year
  // Annual fecundity considering juvenile survival
  const annualFecundity = demographics.breedingProportion * 
                          demographics.chicksPerBreeding * 
                          demographics.survivalJuvenile;
  
  // Generation-based growth (simplified)
  // Assuming stable population in wild (lambda ≈ 1.0)
  return 1.0; // Stable population assumption
};

/**
 * Calculate heterozygosity at generation t using Wright's equation
 * WITH STOCHASTIC VARIATION (optional)
 * Ht = H0 * (1 - 1/(2*Ne))^t
 * @param {number} H0 - Initial heterozygosity
 * @param {number} Ne - Effective population size
 * @param {number} t - Number of generations
 * @param {boolean} stochastic - Add random demographic variation
 * @returns {number} Heterozygosity at generation t
 */
export const calculateGeneticDiversity = (H0, Ne, t, stochastic = false) => {
  let Ht = H0 * Math.pow(1 - 1 / (2 * Ne), t);
  
  if (stochastic && t > 0) {
    // Add demographic stochasticity (±5% random variation per generation)
    const variance = 0.05;
    const randomFactor = 1 + (Math.random() - 0.5) * 2 * variance;
    Ht = Ht * randomFactor;
    Ht = Math.max(0, Math.min(1, Ht)); // Clamp between 0 and 1
  }
  
  return Ht;
};

/**
 * Calculate inbreeding coefficient at generation t
 * WITH STOCHASTIC VARIATION (optional)
 * F = 1 - (1 - 1/(2*Ne))^t
 * @param {number} Ne - Effective population size
 * @param {number} t - Number of generations
 * @param {boolean} stochastic - Add random variation
 * @returns {number} Inbreeding coefficient at generation t
 */
export const calculateInbreeding = (Ne, t, stochastic = false) => {
  let F = 1 - Math.pow(1 - 1 / (2 * Ne), t);
  
  if (stochastic && t > 0) {
    // Add stochasticity (±3% random variation)
    const variance = 0.03;
    const randomFactor = 1 + (Math.random() - 0.5) * 2 * variance;
    F = F * randomFactor;
    F = Math.max(0, Math.min(1, F)); // Clamp between 0 and 1
  }
  
  return F;
};

/**
 * Calculate allelic richness at generation t (simplified exponential decay)
 * WITH STOCHASTIC VARIATION (optional)
 * At = A0 * exp(-t / (4*Ne))
 * @param {number} A0 - Initial number of alleles
 * @param {number} Ne - Effective population size
 * @param {number} t - Number of generations
 * @param {boolean} stochastic - Add random allele loss
 * @returns {number} Number of alleles at generation t
 */
export const calculateAlleles = (A0, Ne, t, stochastic = false) => {
  let At = A0 * Math.exp(-t / (4 * Ne));
  
  if (stochastic && t > 0) {
    // Alleles can be randomly lost (especially rare ones)
    // Add ±10% variation (allele loss is more stochastic)
    const variance = 0.10;
    const randomFactor = 1 + (Math.random() - 0.5) * 2 * variance;
    At = At * randomFactor;
  }
  
  return Math.max(At, 2); // Minimum 2 alleles
};

/**
 * Project population size at generation t
 * WITH STOCHASTIC VARIATION (optional)
 * @param {number} N0 - Initial population size
 * @param {number} lambda - Population growth rate
 * @param {number} t - Number of generations
 * @param {number|null} carryingCapacity - Optional carrying capacity for logistic growth
 * @param {boolean} stochastic - Add demographic stochasticity
 * @returns {number} Population size at generation t
 */
export const projectPopulationSize = (N0, lambda, t, carryingCapacity = null, stochastic = false) => {
  // Exponential growth: Nt = N0 * lambda^t
  let Nt = N0 * Math.pow(lambda, t);
  
  // If carrying capacity is set, use logistic growth
  if (carryingCapacity) {
    // Logistic: Nt = K / (1 + ((K - N0) / N0) * exp(-r*t))
    const r = Math.log(lambda);
    Nt = carryingCapacity / (1 + ((carryingCapacity - N0) / N0) * Math.exp(-r * t));
  }
  
  if (stochastic && t > 0) {
    // Environmental stochasticity: good years vs bad years
    // ±15% variation in population size per generation
    const variance = 0.15;
    const randomFactor = 1 + (Math.random() - 0.5) * 2 * variance;
    Nt = Nt * randomFactor;
    Nt = Math.max(10, Nt); // Minimum 10 birds (quasi-extinction threshold)
  }
  
  return Math.round(Nt);
};

/**
 * Estimate initial population size from sample sizes and literature
 * Literature suggests ~400-500 wild birds in South Africa
 * @param {number} sampleSize - Sample size for specific population
 * @param {number} totalSamples - Total samples across all populations
 * @param {number} totalPopulation - Estimated total population from literature
 * @returns {number} Estimated population size
 */
export const estimatePopulationSize = (sampleSize, totalSamples = 199, totalPopulation = 450) => {
  // Proportional estimation
  return Math.round((sampleSize / totalSamples) * totalPopulation);
};

/**
 * Calculate effective Ne based on census population size
 * Ne/N ratio typically 0.1-0.3 for most species
 * Using 0.2 as reasonable estimate for cooperative breeders
 * @param {number} censusSize - Census population size
 * @param {number} ratio - Ne/N ratio (default 0.2)
 * @returns {number} Effective population size
 */
export const calculateEffectiveNe = (censusSize, ratio = 0.2) => {
  return Math.round(censusSize * ratio);
};