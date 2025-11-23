

// --- Constants ---

// Genetic data from CSV files - initial heterozygosity values
const initialGenetics = {
    wild: { Ho: 0.502, He: 0.568, Na: 6.429, Ne_effective: 3.143, FIS: 0.121, sampleSize: 199 },
    wildNoECKZN: { Ho: 0.502, He: 0.568, Na: 6.2, sampleSize: 177 },
    paaza: { Ho: 0.483, He: 0.565, Na: 7.071, Ne_effective: 3.201, FIS: 0.139, sampleSize: 70 },
    aza: { Ho: 0.515, He: 0.534, Na: 5.429, Ne_effective: 2.666, FIS: 0.051, sampleSize: 29 },
    eaza: { Ho: 0.564, He: 0.582, Na: 5.571, Ne_effective: 3.100, FIS: 0.041, sampleSize: 46 }
};

// Demographic parameters based on provided breeding information
const demographics = {
    breedingAge: 10, // years
    generationTime: 15, // years
    breedingProportion: 0.20, // 20% of population breeds
    chicksPerBreeding: 1 / 3, // 1 chick every 3 years
    groupSize: 5.5, // Average 3-9, median ~5.5
    survivalAdult: 0.92, // Assumed from literature (high adult survival)
    survivalJuvenile: 0.65, // Assumed (lower juvenile survival)
    survivalSubadult: 0.85 // Assumed intermediate
};

// --- Calculation Functions ---

const calculateLambda = () => 1.0; // Stable population assumption

const calculateGeneticDiversity = (H0, Ne, t) => H0 * Math.pow(1 - 1 / (2 * Ne), t);

const calculateInbreeding = (Ne, t) => 1 - Math.pow(1 - 1 / (2 * Ne), t);

const calculateAlleles = (A0, Ne, t) => {
    const At = A0 * Math.exp(-t / (4 * Ne));
    return Math.max(At, 2); // Minimum 2 alleles
};

const projectPopulationSize = (N0, lambda, t, carryingCapacity = null) => {
    let Nt = N0 * Math.pow(lambda, t);
    if (carryingCapacity) {
        const r = Math.log(lambda);
        Nt = carryingCapacity / (1 + ((carryingCapacity - N0) / N0) * Math.exp(-r * t));
    }
    return Math.round(Nt);
};

const estimatePopulationSize = (sampleSize, totalSamples = 199, totalPopulation = 450) => {
    return Math.round((sampleSize / totalSamples) * totalPopulation);
};

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

// --- Main Simulation Runner ---

export const runSimulationModels = ({ ne, generations }) => {
    // Initial population estimates
    const N_wild = estimatePopulationSize(199, 199, 450);
    const N_wildNoECKZN = estimatePopulationSize(177, 199, 450);
    const lambda = calculateLambda();

    const models = [];

    for (let model = 1; model <= 5; model++) {
        const modelData = [];

        for (let t = 0; t <= generations; t++) {
            let H0, A0, N0, suppAdded = 0;

            switch (model) {
                case 1: // Baseline wild
                    H0 = initialGenetics.wild.Ho;
                    A0 = initialGenetics.wild.Na;
                    N0 = N_wild;
                    break;
                case 2: // Wild without EC and KZN
                    H0 = initialGenetics.wildNoECKZN.Ho;
                    A0 = initialGenetics.wildNoECKZN.Na;
                    N0 = N_wildNoECKZN;
                    break;
                case 3: // Wild + 4 PAAZA birds/generation
                    H0 = initialGenetics.wild.Ho;
                    A0 = initialGenetics.wild.Na;
                    N0 = N_wild;
                    suppAdded = 4 * t;
                    break;
                case 4: // Wild + 10 PAAZA birds/generation
                    H0 = initialGenetics.wild.Ho;
                    A0 = initialGenetics.wild.Na;
                    N0 = N_wild;
                    suppAdded = 10 * t;
                    break;
                case 5: // Wild + 4 mixed captive birds/generation
                    H0 = initialGenetics.wild.Ho;
                    A0 = initialGenetics.wild.Na;
                    N0 = N_wild;
                    suppAdded = 4 * t;
                    break;
            }

            let effectiveNe = ne;
            if (model >= 3) {
                const immigrantContribution = suppAdded * 0.5;
                effectiveNe = ne + immigrantContribution;
            }

            const Ho = calculateGeneticDiversity(H0, effectiveNe, t);
            const He = calculateGeneticDiversity(initialGenetics.wild.He, effectiveNe, t);
            const F = calculateInbreeding(effectiveNe, t);
            const Na = calculateAlleles(A0, effectiveNe, t);
            const popSize = projectPopulationSize(N0, lambda, t) + suppAdded;

            modelData.push({
                generation: t,
                year: t * demographics.generationTime,
                Ho: Ho.toFixed(4),
                He: He.toFixed(4),
                F: F.toFixed(4),
                Na: Na.toFixed(2),
                popSize: popSize,
                model: model
            });
        }

        models.push({ modelNumber: model, name: getModelName(model), data: modelData });
    }

    return models;
};