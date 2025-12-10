from flask import Flask, render_template, request, jsonify
import numpy as np
import os
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add parent directory to path to import genetic_data
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genetic_data import (
    load_and_process_genetic_data,
    simulate_supplementation_effect,
    calculate_observed_heterozygosity,
    calculate_expected_heterozygosity_population,
    calculate_allelic_richness,
    LOCI
)

# Get base directory for path resolution
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

app = Flask(__name__,
            template_folder=os.path.join(BASE_DIR, 'templates'),
            static_folder=os.path.join(BASE_DIR, 'static'))


# LOAD GENETIC DATA FROM CSVs


# Paths to CSV files
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WILD_CSV = os.path.join(BASE_DIR, 'data', 'blue.csv')
CAPTIVE_CSV = os.path.join(BASE_DIR, 'data', 'red.csv')

# Global variable to store processed genetic data
GENETIC_DATA = None

def initialize_genetic_data():
    """Load and process CSV data on startup"""
    global GENETIC_DATA

    if os.path.exists(WILD_CSV) and os.path.exists(CAPTIVE_CSV):
        logger.info("Loading genetic data from CSVs...")
        GENETIC_DATA = load_and_process_genetic_data(WILD_CSV, CAPTIVE_CSV)
        logger.info(f"Loaded {len(GENETIC_DATA['summaries'])} populations")
        logger.info(f"Wild population Ho: {GENETIC_DATA['summaries']['wild_all']['Ho']:.4f}")
    else:
        logger.warning("CSV files not found. Using default values.")
        GENETIC_DATA = None



# GENETIC CALCULATION FUNCTIONS


def calculate_heterozygosity_loss(H0, Ne, t):
    """Wright's equation: Ht = H0 * (1 - 1/(2*Ne))^t"""
    return H0 * np.power(1 - 1/(2*Ne), t)


def calculate_inbreeding(Ne, t):
    """Inbreeding coefficient: F = 1 - (1 - 1/(2*Ne))^t"""
    return 1 - np.power(1 - 1/(2*Ne), t)


def calculate_allelic_diversity(A0, Ne, t):
    """Allelic diversity: At = A0 * exp(-t/(4*Ne))"""
    At = A0 * np.exp(-t / (4*Ne))
    return np.maximum(At, 2.0)


def calculate_population_size(N0, lambda_val, t):
    """Population size: Nt = N0 * lambda^t"""
    return N0 * np.power(lambda_val, t)


def calculate_population_size_with_inbreeding(N0, lambda_val, F_array, lethal_equivalents=3.14):
    """
    Population size with inbreeding depression

    Survival is reduced by: s = exp(-B * F)
    where B = lethal equivalents (typical range: 3-12 for mammals, 2-6 for birds)

    For Southern Ground Hornbill, using B = 3.14 (moderate inbreeding depression for birds)
    Based on O'Grady et al. 2006 review showing birds average ~3.14 lethal equivalents
    Previous value of 6.29 was too severe and more appropriate for mammals

    Args:
        N0: Initial population size
        lambda_val: Base growth rate (without inbreeding)
        F_array: Array of inbreeding coefficients over time
        lethal_equivalents: Number of lethal equivalents (default 3.14)

    Returns: Array of population sizes incorporating inbreeding depression
    """
    N = np.zeros(len(F_array))
    N[0] = N0

    for t in range(1, len(F_array)):
        # Calculate survival reduction due to inbreeding
        inbreeding_survival = np.exp(-lethal_equivalents * F_array[t])

        # Adjusted growth rate
        lambda_adjusted = lambda_val * inbreeding_survival

        # Population size at time t
        N[t] = N[t-1] * lambda_adjusted

        # Prevent population from going extinct (minimum viable)
        N[t] = max(N[t], 10)

    return N



# MODEL IMPLEMENTATIONS


def run_model_1(Ne, generations, lambda_val=1.0):
    """
    Model 1: Baseline - All wild populations
    Uses real CSV data if available
    """
    if GENETIC_DATA:
        # Use real data
        summary = GENETIC_DATA['summaries']['wild_all']
        H0 = summary['Ho']
        He0 = summary['He']
        A0 = summary['Na']
        N0 = summary['sample_size'] * (2500 / 199)  # Scale to census
    else:
        # Fallback to defaults
        H0 = 0.502
        He0 = 0.568
        A0 = 6.429
        N0 = 2500
    t = np.arange(0, generations + 1)

    # Calculate genetic metrics
    Ho_vals = calculate_heterozygosity_loss(H0, Ne, t)
    He_vals = calculate_heterozygosity_loss(He0, Ne, t)
    F_vals = calculate_inbreeding(Ne, t)
    Na_vals = calculate_allelic_diversity(A0, Ne, t)

    # Calculate population size WITH inbreeding depression
    N_vals = calculate_population_size_with_inbreeding(N0, lambda_val, F_vals)

    return {
        'model_number': 1,
        'model_name': 'Baseline (All Wild Populations)',
        'generations': t.tolist(),
        'years': (t * 26).tolist(),
        'Ho': Ho_vals.tolist(),
        'He': He_vals.tolist(),
        'F': F_vals.tolist(),
        'Na': Na_vals.tolist(),
        'population': N_vals.tolist(),
        'Ne': Ne,
        'initial': {'Ho': H0, 'He': He0, 'Na': A0, 'N': N0},
        'parameters': {
            'Ne': Ne,
            'N0': int(N0),
            'lambda': lambda_val,
            'data_source': 'CSV' if GENETIC_DATA else 'default',
            'populations': ['Eastern Cape', 'Kruger NP', 'KwaZulu-Natal', 'Limpopo'],
            'inbreeding_depression': 'enabled (B=3.14 lethal equivalents for birds)'
        }
    }


def run_model_2(Ne, generations, lambda_val=1.0):
    """
    Model 2: Population Loss - Only Kruger + Limpopo
    Uses real CSV data to show actual allele loss
    """
    if GENETIC_DATA:
        summary = GENETIC_DATA['summaries']['wild_no_ec_kzn']
        H0 = summary['Ho']
        He0 = summary['He']
        A0 = summary['Na']
        N0 = summary['sample_size'] * (2500 / 199)

        # Get actual lost alleles
        lost_alleles = GENETIC_DATA['lost_alleles_ec_kzn']
        lost_count = sum(len(alleles) for alleles in lost_alleles.values())
    else:
        H0 = 0.498
        He0 = 0.565
        A0 = 6.1
        N0 = 2000
        lost_count = 0

    # Scale Ne proportionally
    Ne_scaled = int(Ne * (N0 / 2500))
    t = np.arange(0, generations + 1)

    # Calculate genetic metrics
    Ho_vals = calculate_heterozygosity_loss(H0, Ne_scaled, t)
    He_vals = calculate_heterozygosity_loss(He0, Ne_scaled, t)
    F_vals = calculate_inbreeding(Ne_scaled, t)
    Na_vals = calculate_allelic_diversity(A0, Ne_scaled, t)

    # Calculate population size WITH inbreeding depression
    N_vals = calculate_population_size_with_inbreeding(N0, lambda_val, F_vals)

    return {
        'model_number': 2,
        'model_name': 'Population Loss (Kruger + Limpopo Only)',
        'generations': t.tolist(),
        'years': (t * 26).tolist(),
        'Ho': Ho_vals.tolist(),
        'He': He_vals.tolist(),
        'F': F_vals.tolist(),
        'Na': Na_vals.tolist(),
        'population': N_vals.tolist(),
        'Ne': Ne_scaled,
        'initial': {'Ho': H0, 'He': He0, 'Na': A0, 'N': N0},
        'parameters': {
            'Ne': Ne_scaled,
            'N0': int(N0),
            'lambda': lambda_val,
            'data_source': 'CSV' if GENETIC_DATA else 'default',
            'populations': ['Kruger NP', 'Limpopo'],
            'lost_populations': ['Eastern Cape', 'KwaZulu-Natal'],
            'alleles_lost': lost_count if GENETIC_DATA else 'unknown',
            'inbreeding_depression': 'enabled (B=3.14 lethal equivalents for birds)'
        }
    }


def run_model_3(Ne, generations, lambda_val=1.0):
    """
    Model 3: Low Supplementation - Add 4 PAAZA birds per generation
    Uses real CSV data to track actual genetic contribution
    """
    if GENETIC_DATA:
        # Use real CSV simulation
        wild_df = GENETIC_DATA['dataframes']['wild_all']
        paaza_df = GENETIC_DATA['dataframes']['paaza']

        results = simulate_supplementation_effect(
            wild_df, paaza_df,
            birds_per_gen=4,
            generations=generations,
            loci_list=LOCI
        )

        # Extract data
        Ho = [r['Ho'] for r in results]
        He = [r['He'] for r in results]
        Na = [r['Na'] for r in results]
        F = [r['FIS'] for r in results]
        sample_sizes = [r['population_size'] for r in results]

        # CRITICAL FIX: Keep initial wild population size separate from genetic samples
        # The genetic simulation adds samples for tracking alleles, not census population
        N0 = 2500  # Fixed census size for wild population

        # Extract F array from genetic simulation
        F_array = np.array(F)

        # Track supplemented birds parameters
        birds_per_gen = 4
        captive_survival_multiplier = 0.95

        # CRITICAL FIX: Apply genetic rescue effect BEFORE calculating wild population dynamics
        # Gene flow from supplementation reduces population-wide inbreeding
        F_array_rescued = F_array.copy()
        for gen in range(len(F_array)):
            if gen > 0:
                # Gene flow benefit proportional to cumulative supplementation
                cumulative_birds = birds_per_gen * gen
                total_pop_estimate = N0 + cumulative_birds
                gene_flow_proportion = cumulative_birds / total_pop_estimate
                # Reduce F by 50% of gene flow proportion (genetic rescue)
                F_array_rescued[gen] = F_array[gen] * (1 - gene_flow_proportion * 0.5)

        # Apply inbreeding depression to wild-born population using RESCUED F values
        N_wild = calculate_population_size_with_inbreeding(N0, lambda_val, F_array_rescued)

        # CRITICAL FIX: Don't double-count supplemented birds
        # Birds are already in genetic simulation; demographic tracking is for display only
        # Track them separately to show their contribution but they don't add to total
        N_released = np.zeros(len(N_wild))

        for gen in range(len(N_wild)):
            # Track all cohorts of released birds for visualization
            for release_gen in range(gen + 1):
                time_since_release = gen - release_gen
                # Released birds experience much lower inbreeding depression (15% vs wild)
                # This represents hybrid vigor / heterosis benefit
                avg_F_since_release = np.mean(F_array_rescued[release_gen:gen+1])
                survival = np.exp(-3.14 * avg_F_since_release * 0.15) * captive_survival_multiplier
                cohort_survivors = birds_per_gen * (survival ** time_since_release)
                N_released[gen] += cohort_survivors

        # Total population = wild-born + released survivors
        # NOTE: In reality, released birds are integrated into wild population
        # This separation is for tracking genetic rescue contribution
        N = N_wild + N_released

        t = np.arange(0, generations + 1)

        return {
            'model_number': 3,
            'model_name': 'Low Supplementation (+4 SA Captive/gen)',
            'generations': t.tolist(),
            'years': (t * 26).tolist(),
            'Ho': Ho,
            'He': He,
            'F': F,
            'Na': Na,
            'population': N.tolist(),
            'Ne': Ne,
            'initial': {'Ho': Ho[0], 'He': He[0], 'Na': Na[0], 'N': N[0]},
            'parameters': {
                'Ne': Ne,
                'N0': int(N[0]),
                'lambda': lambda_val,
                'data_source': 'CSV_simulation',
                'supplementation': '4 South African captive birds per generation',
                'supplementation_source': 'PAAZA (Pan-African Association of Zoos and Aquaria)',
                'novel_alleles_added': len(GENETIC_DATA['novel_alleles']['paaza']),
                'inbreeding_depression': 'enabled (B=3.14 lethal equivalents for birds)'
            }
        }
    else:
        # Fallback: Generic supplementation model
        return run_generic_supplementation_model(Ne, generations, lambda_val,
                                                  birds_per_gen=4,
                                                  model_num=3,
                                                  source='SA Captive')


def run_model_4(Ne, generations, lambda_val=1.0):
    """
    Model 4: High Supplementation - Add 10 PAAZA birds per generation
    """
    if GENETIC_DATA:
        wild_df = GENETIC_DATA['dataframes']['wild_all']
        paaza_df = GENETIC_DATA['dataframes']['paaza']

        results = simulate_supplementation_effect(
            wild_df, paaza_df,
            birds_per_gen=10,
            generations=generations,
            loci_list=LOCI
        )

        Ho = [r['Ho'] for r in results]
        He = [r['He'] for r in results]
        Na = [r['Na'] for r in results]
        F = [r['FIS'] for r in results]
        sample_sizes = [r['population_size'] for r in results]

        # CRITICAL FIX: Keep initial wild population size separate from genetic samples
        N0 = 2500  # Fixed census size for wild population

        # Extract F array from genetic simulation
        F_array = np.array(F)

        # Track supplemented birds parameters
        birds_per_gen = 10
        captive_survival_multiplier = 0.95

        # CRITICAL FIX: Apply genetic rescue effect BEFORE calculating wild population dynamics
        F_array_rescued = F_array.copy()
        for gen in range(len(F_array)):
            if gen > 0:
                cumulative_birds = birds_per_gen * gen
                total_pop_estimate = N0 + cumulative_birds
                gene_flow_proportion = cumulative_birds / total_pop_estimate
                F_array_rescued[gen] = F_array[gen] * (1 - gene_flow_proportion * 0.5)

        # Apply inbreeding depression to wild-born population using RESCUED F values
        N_wild = calculate_population_size_with_inbreeding(N0, lambda_val, F_array_rescued)

        # CRITICAL FIX: Don't double-count supplemented birds
        N_released = np.zeros(len(N_wild))

        for gen in range(len(N_wild)):
            # Track all cohorts of released birds for visualization
            for release_gen in range(gen + 1):
                time_since_release = gen - release_gen
                avg_F_since_release = np.mean(F_array_rescued[release_gen:gen+1])
                survival = np.exp(-3.14 * avg_F_since_release * 0.15) * captive_survival_multiplier
                cohort_survivors = birds_per_gen * (survival ** time_since_release)
                N_released[gen] += cohort_survivors

        # Total population = wild-born + released survivors
        N = N_wild + N_released

        t = np.arange(0, generations + 1)

        return {
            'model_number': 4,
            'model_name': 'High Supplementation (+10 SA Captive/gen)',
            'generations': t.tolist(),
            'years': (t * 26).tolist(),
            'Ho': Ho,
            'He': He,
            'F': F,
            'Na': Na,
            'population': N.tolist(),
            'Ne': Ne,
            'initial': {'Ho': Ho[0], 'He': He[0], 'Na': Na[0], 'N': N[0]},
            'parameters': {
                'Ne': Ne,
                'N0': int(N[0]),
                'lambda': lambda_val,
                'data_source': 'CSV_simulation',
                'supplementation': '10 South African captive birds per generation',
                'supplementation_source': 'PAAZA (Pan-African Association of Zoos and Aquaria)',
                'inbreeding_depression': 'enabled (B=3.14 lethal equivalents for birds)'
            }
        }
    else:
        return run_generic_supplementation_model(Ne, generations, lambda_val,
                                                  birds_per_gen=10,
                                                  model_num=4,
                                                  source='SA Captive')


def run_model_5(Ne, generations, lambda_val=1.0):
    """
    Model 5: International Mix - Add 4 mixed birds (PAAZA/AZA/EAZA) per generation
    """
    if GENETIC_DATA:
        wild_df = GENETIC_DATA['dataframes']['wild_all']

        # Mix captive populations
        import pandas as pd
        paaza_df = GENETIC_DATA['dataframes']['paaza']
        aza_df = GENETIC_DATA['dataframes']['aza']
        eaza_df = GENETIC_DATA['dataframes']['eaza']
        mixed_captive = pd.concat([paaza_df, aza_df, eaza_df], ignore_index=True)

        results = simulate_supplementation_effect(
            wild_df, mixed_captive,
            birds_per_gen=4,
            generations=generations,
            loci_list=LOCI
        )

        Ho = [r['Ho'] for r in results]
        He = [r['He'] for r in results]
        Na = [r['Na'] for r in results]
        F = [r['FIS'] for r in results]

        # CRITICAL FIX: Keep initial wild population size separate from genetic samples
        N0 = 2500  # Fixed census size for wild population

        # Extract F array from genetic simulation
        F_array = np.array(F)

        # Track supplemented birds parameters
        birds_per_gen = 4
        captive_survival_multiplier = 0.95

        # CRITICAL FIX: Apply genetic rescue effect BEFORE calculating wild population dynamics
        F_array_rescued = F_array.copy()
        for gen in range(len(F_array)):
            if gen > 0:
                cumulative_birds = birds_per_gen * gen
                total_pop_estimate = N0 + cumulative_birds
                gene_flow_proportion = cumulative_birds / total_pop_estimate
                F_array_rescued[gen] = F_array[gen] * (1 - gene_flow_proportion * 0.5)

        # Apply inbreeding depression to wild-born population using RESCUED F values
        N_wild = calculate_population_size_with_inbreeding(N0, lambda_val, F_array_rescued)

        # CRITICAL FIX: Don't double-count supplemented birds
        N_released = np.zeros(len(N_wild))

        for gen in range(len(N_wild)):
            # Track all cohorts of released birds for visualization
            for release_gen in range(gen + 1):
                time_since_release = gen - release_gen
                avg_F_since_release = np.mean(F_array_rescued[release_gen:gen+1])
                survival = np.exp(-3.14 * avg_F_since_release * 0.15) * captive_survival_multiplier
                cohort_survivors = birds_per_gen * (survival ** time_since_release)
                N_released[gen] += cohort_survivors

        # Total population = wild-born + released survivors
        N = N_wild + N_released

        t = np.arange(0, generations + 1)

        # Count total novel alleles
        novel_count = sum(
            len(GENETIC_DATA['novel_alleles'][src])
            for src in ['paaza', 'aza', 'eaza']
        )

        return {
            'model_number': 5,
            'model_name': 'International Mix (+4 Mixed/gen)',
            'generations': t.tolist(),
            'years': (t * 26).tolist(),
            'Ho': Ho,
            'He': He,
            'F': F,
            'Na': Na,
            'population': N.tolist(),
            'Ne': Ne,
            'initial': {'Ho': Ho[0], 'He': He[0], 'Na': Na[0], 'N': N[0]},
            'parameters': {
                'Ne': Ne,
                'N0': int(N[0]),
                'lambda': lambda_val,
                'data_source': 'CSV_simulation',
                'supplementation': '4 mixed birds (South African + USA + European zoos) per generation',
                'supplementation_sources': 'PAAZA (South Africa) + AZA (USA/Canada) + EAZA (Europe)',
                'novel_alleles_total': novel_count,
                'inbreeding_depression': 'enabled (B=3.14 lethal equivalents for birds)'
            }
        }
    else:
        return run_generic_supplementation_model(Ne, generations, lambda_val,
                                                  birds_per_gen=4,
                                                  model_num=5,
                                                  source='Mixed')


def run_generic_supplementation_model(Ne, generations, lambda_val, birds_per_gen, model_num, source):
    """
    Generic supplementation model when CSV data not available
    """
    H0 = 0.502
    He0 = 0.568
    A0 = 6.429
    N0 = 2500

    t = np.arange(0, generations + 1)

    # Simple model: Ne increases with supplementation
    Ho_vals = []
    He_vals = []
    F_vals = []
    Na_vals = []

    for gen in t:
        # Cumulative birds added
        birds_added = birds_per_gen * gen

        # Effective Ne increases with gene flow
        effective_Ne = Ne + (birds_added * 0.5)

        # Calculate metrics
        Ho = calculate_heterozygosity_loss(H0, effective_Ne, gen)
        He = calculate_heterozygosity_loss(He0, effective_Ne, gen)
        F = calculate_inbreeding(effective_Ne, gen)
        Na = calculate_allelic_diversity(A0, effective_Ne, gen)

        Ho_vals.append(Ho)
        He_vals.append(He)
        F_vals.append(F)
        Na_vals.append(Na)

    # Extract F array
    F_array = np.array(F_vals)

    # Track supplemented birds parameters
    captive_survival_multiplier = 0.95

    # CRITICAL FIX: Apply genetic rescue effect BEFORE calculating wild population dynamics
    F_array_rescued = F_array.copy()
    for gen in range(len(F_array)):
        if gen > 0:
            cumulative_birds = birds_per_gen * gen
            total_pop_estimate = N0 + cumulative_birds
            gene_flow_proportion = cumulative_birds / total_pop_estimate
            F_array_rescued[gen] = F_array[gen] * (1 - gene_flow_proportion * 0.5)

    # Calculate population size WITH inbreeding depression using RESCUED F values
    N_wild = calculate_population_size_with_inbreeding(N0, lambda_val, F_array_rescued)

    # CRITICAL FIX: Don't double-count supplemented birds
    N_released = np.zeros(len(N_wild))

    for gen in range(len(t)):
        # Track all cohorts of released birds
        for release_gen in range(gen + 1):
            time_since_release = gen - release_gen
            # Released birds experience much lower inbreeding depression (15% vs wild, improved from 30%)
            # This represents hybrid vigor / heterosis benefit
            avg_F_since_release = np.mean(F_array_rescued[release_gen:gen+1])
            survival = np.exp(-3.14 * avg_F_since_release * 0.15) * captive_survival_multiplier
            cohort_survivors = birds_per_gen * (survival ** time_since_release)
            N_released[gen] += cohort_survivors

    # Total population = wild-born + released survivors
    N_vals = N_wild + N_released

    return {
        'model_number': model_num,
        'model_name': f'Supplementation (+{birds_per_gen} {source}/gen)',
        'generations': t.tolist(),
        'years': (t * 26).tolist(),
        'Ho': Ho_vals,
        'He': He_vals,
        'F': F_vals,
        'Na': Na_vals,
        'population': N_vals.tolist(),
        'Ne': Ne,
        'initial': {'Ho': H0, 'He': He0, 'Na': A0, 'N': N0},
        'parameters': {
            'Ne': Ne,
            'N0': N0,
            'lambda': lambda_val,
            'data_source': 'generic',
            'supplementation': f'{birds_per_gen} {source} birds per generation',
            'inbreeding_depression': 'enabled (B=3.14 lethal equivalents for birds)'
        }
    }



# API ROUTES


@app.route('/')
def index():
    """Serve the main page"""
    return render_template('index.html')


@app.route('/api/simulate', methods=['POST'])
def simulate():
    """Run simulation for specified model"""
    try:
        data = request.get_json()

        Ne = int(data.get('Ne', 500))
        generations = int(data.get('generations', 50))
        lambda_val = float(data.get('lambda', 1.0))
        model = int(data.get('model', 1))

        # Validate inputs
        if not 375 <= Ne <= 625:
            return jsonify({'error': 'Ne must be between 375 and 625'}), 400
        if not 10 <= generations <= 100:
            return jsonify({'error': 'Generations must be between 10 and 100'}), 400
        if not 0.5 <= lambda_val <= 1.5:
            return jsonify({'error': 'Lambda must be between 0.5 and 1.5'}), 400
        if not 1 <= model <= 5:
            return jsonify({'error': 'Model must be between 1 and 5'}), 400

        # Run appropriate model
        if model == 1:
            results = run_model_1(Ne, generations, lambda_val)
        elif model == 2:
            results = run_model_2(Ne, generations, lambda_val)
        elif model == 3:
            results = run_model_3(Ne, generations, lambda_val)
        elif model == 4:
            results = run_model_4(Ne, generations, lambda_val)
        elif model == 5:
            results = run_model_5(Ne, generations, lambda_val)

        return jsonify(results)

    except Exception as e:
        import traceback
        traceback.print_exc()
        logger.error(f"Simulation error: {str(e)}")
        return jsonify({'error': str(e)}), 500


@app.route('/api/data/info', methods=['GET'])
def get_data_info():
    """Get information about genetic data"""
    if GENETIC_DATA:
        return jsonify({
            'data_source': 'CSV',
            'populations': {
                name: {
                    'sample_size': summary['sample_size'],
                    'Ho': summary['Ho'],
                    'He': summary['He'],
                    'Na': summary['Na'],
                    'FIS': summary['FIS']
                }
                for name, summary in GENETIC_DATA['summaries'].items()
            },
            'lost_alleles_count': sum(
                len(alleles) for alleles in GENETIC_DATA['lost_alleles_ec_kzn'].values()
            ),
            'novel_alleles': {
                'paaza': sum(len(a) for a in GENETIC_DATA['novel_alleles']['paaza'].values()),
                'aza': sum(len(a) for a in GENETIC_DATA['novel_alleles']['aza'].values()),
                'eaza': sum(len(a) for a in GENETIC_DATA['novel_alleles']['eaza'].values())
            }
        })
    else:
        return jsonify({
            'data_source': 'default',
            'message': 'Using default genetic parameters. Upload CSVs for real data.'
        })


# Initialize data on import
initialize_genetic_data()

# For Vercel serverless functions
app = app

if __name__ == '__main__':
    
    app.run(debug=True, port=5001) # Using a different port to avoid conflict with app.py
