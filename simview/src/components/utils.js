import { demographics } from '../js/constant';

/**
 * Calculate dominant eigenvalue of a square matrix
 * Uses power iteration method for computational efficiency
 * @param {number[][]} matrix - Square matrix
 * @param {number} maxIterations - Maximum iterations (default 1000)
 * @param {number} tolerance - Convergence tolerance (default 1e-6)
 * @returns {number} Dominant eigenvalue
 */
const dominantEigenvalue = (matrix, maxIterations = 1000, tolerance = 1e-6) => {
  const n = matrix.length;
  
  // Start with random vector
  let vector = Array(n).fill(0).map(() => Math.random());
  
  // Normalize
  const magnitude = Math.sqrt(vector.reduce((sum, v) => sum + v * v, 0));
  vector = vector.map(v => v / magnitude);
  
  let eigenvalue = 0;
  let prevEigenvalue = 0;
  
  for (let iter = 0; iter < maxIterations; iter++) {
    // Multiply matrix by vector
    const newVector = matrix.map(row => 
      row.reduce((sum, val, i) => sum + val * vector[i], 0)
    );
    
    // Calculate eigenvalue (Rayleigh quotient)
    eigenvalue = newVector.reduce((sum, v, i) => sum + v * vector[i], 0) /
                 vector.reduce((sum, v) => sum + v * v, 0);
    
    // Normalize new vector
    const newMagnitude = Math.sqrt(newVector.reduce((sum, v) => sum + v * v, 0));
    vector = newVector.map(v => v / newMagnitude);
    
    // Check convergence
    if (Math.abs(eigenvalue - prevEigenvalue) < tolerance) {
      break;
    }
    
    prevEigenvalue = eigenvalue;
  }
  
  return eigenvalue;
};

/**
 * Build Leslie matrix from age-structured demographic parameters
 * Leslie matrix structure:
 *   [f1  f2  f3]   First row: fecundities
 *   [s1  0   0 ]   Subdiagonal: survival rates
 *   [0   s2  s3]   All other entries: 0
 * @returns {number[][]} Leslie matrix
 */
const buildLeslieMatrix = () => {
  const ages = demographics.ageClasses;
  const n = ages.length;
  
  // Initialize matrix with zeros
  const matrix = Array(n).fill(0).map(() => Array(n).fill(0));
  
  // First row: fecundities (births per individual in each age class)
  ages.forEach((ageClass, i) => {
    matrix[0][i] = ageClass.fecundity;
  });
  
  // Subdiagonal: survival rates (probability of surviving to next age class)
  for (let i = 0; i < n - 1; i++) {
    matrix[i + 1][i] = ages[i].survival;
  }
  
  // Last diagonal: adult survival (staying in adult class)
  matrix[n - 1][n - 1] = ages[n - 1].survival;
  
  return matrix;
};

/**
 * Calculate population growth rate (lambda) using Leslie matrix eigenvalue analysis
 * This is the PROPER demographic method for age-structured populations
 * 
 * λ > 1: Population growing
 * λ = 1: Population stable
 * λ < 1: Population declining
 * 
 * @param {boolean} useLeslieMatrix - If true, calculate from Leslie matrix. If false, assume λ=1.0
 * @returns {number} Lambda (population growth rate per generation)
 */
export const calculateLambda = (useLeslieMatrix = false) => {
  if (!useLeslieMatrix) {
    // Simple assumption: stable population (lambda = 1.0)
    // Rationale: Simplifies genetic drift analysis by removing demographic confounds
    return 1.0;
  }
  
  // Build age-structured Leslie matrix
  const leslieMatrix = buildLeslieMatrix();
  
  // The dominant eigenvalue of the Leslie matrix is the population growth rate
  // This is the mathematically correct way to calculate λ for age-structured populations
  const annualLambda = dominantEigenvalue(leslieMatrix);
  
  // Convert annual λ to generation-based λ
  // λ_generation = (λ_annual)^generation_time
  const generationLambda = Math.pow(annualLambda, demographics.generationTime);
  
  // Sanity check: clamp to biologically reasonable range
  // Long-lived birds rarely have |λ - 1| > 0.2 per generation
  const clampedLambda = Math.max(0.80, Math.min(1.20, generationLambda));
  
  return clampedLambda;
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