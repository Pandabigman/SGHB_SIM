from flask import Flask, render_template, request, jsonify
import numpy as np
import os
from genetic_data import (
    load_and_process_genetic_data,
    simulate_supplementation_effect,
    calculate_observed_heterozygosity,
    calculate_expected_heterozygosity_population,
    calculate_allelic_richness,
    LOCI
)

app = Flask(__name__)


# LOAD GENETIC DATA FROM CSVs


# Paths to CSV files
WILD_CSV = 'data/blue.csv'
CAPTIVE_CSV = 'data/red.csv'

# Global variable to store processed genetic data
GENETIC_DATA = None

def initialize_genetic_data():
    """Load and process CSV data on startup"""
    global GENETIC_DATA
    
    if os.path.exists(WILD_CSV) and os.path.exists(CAPTIVE_CSV):
        print("Loading genetic data from CSVs...")
        GENETIC_DATA = load_and_process_genetic_data(WILD_CSV, CAPTIVE_CSV)
        print(f"Loaded {len(GENETIC_DATA['summaries'])} populations")
        print(f"Wild population Ho: {GENETIC_DATA['summaries']['wild_all']['Ho']:.4f}")
    else:
        print("WARNING: CSV files not found. Using default values.")
        GENETIC_DATA = None



# GENETIC CALCULATION FUNCTIONS (for models without CSV)


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



# MODEL IMPLEMENTATIONS


def run_model_1(Ne, generations):
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
        N0 = summary['sample_size'] * (450 / 199)  # Scale to census
    else:
        # Fallback to defaults
        H0 = 0.502
        He0 = 0.568
        A0 = 6.429
        N0 = 450
    
    lambda_val = 1.0
    t = np.arange(0, generations + 1)
    
    return {
        'model_number': 1,
        'model_name': 'Baseline (All Wild Populations)',
        'generations': t.tolist(),
        'years': (t * 15).tolist(),
        'Ho': calculate_heterozygosity_loss(H0, Ne, t).tolist(),
        'He': calculate_heterozygosity_loss(He0, Ne, t).tolist(),
        'F': calculate_inbreeding(Ne, t).tolist(),
        'Na': calculate_allelic_diversity(A0, Ne, t).tolist(),
        'population': calculate_population_size(N0, lambda_val, t).tolist(),
        'Ne': Ne,
        'initial': {'Ho': H0, 'He': He0, 'Na': A0, 'N': N0},
        'parameters': {
            'Ne': Ne,
            'N0': int(N0),
            'lambda': lambda_val,
            'data_source': 'CSV' if GENETIC_DATA else 'default',
            'populations': ['Eastern Cape', 'Kruger NP', 'KwaZulu-Natal', 'Limpopo']
        }
    }


def run_model_2(Ne, generations):
    """
    Model 2: Population Loss - Only Kruger + Limpopo
    Uses real CSV data to show actual allele loss
    """
    if GENETIC_DATA:
        summary = GENETIC_DATA['summaries']['wild_no_ec_kzn']
        H0 = summary['Ho']
        He0 = summary['He']
        A0 = summary['Na']
        N0 = summary['sample_size'] * (450 / 199)
        
        # Get actual lost alleles
        lost_alleles = GENETIC_DATA['lost_alleles_ec_kzn']
        lost_count = sum(len(alleles) for alleles in lost_alleles.values())
    else:
        H0 = 0.498
        He0 = 0.565
        A0 = 6.1
        N0 = 398
        lost_count = 0
    
    # Scale Ne proportionally
    Ne_scaled = int(Ne * (N0 / 450))
    lambda_val = 1.0
    t = np.arange(0, generations + 1)
    
    return {
        'model_number': 2,
        'model_name': 'Population Loss (Kruger + Limpopo Only)',
        'generations': t.tolist(),
        'years': (t * 15).tolist(),
        'Ho': calculate_heterozygosity_loss(H0, Ne_scaled, t).tolist(),
        'He': calculate_heterozygosity_loss(He0, Ne_scaled, t).tolist(),
        'F': calculate_inbreeding(Ne_scaled, t).tolist(),
        'Na': calculate_allelic_diversity(A0, Ne_scaled, t).tolist(),
        'population': calculate_population_size(N0, lambda_val, t).tolist(),
        'Ne': Ne_scaled,
        'initial': {'Ho': H0, 'He': He0, 'Na': A0, 'N': N0},
        'parameters': {
            'Ne': Ne_scaled,
            'N0': int(N0),
            'lambda': lambda_val,
            'data_source': 'CSV' if GENETIC_DATA else 'default',
            'populations': ['Kruger NP', 'Limpopo'],
            'lost_populations': ['Eastern Cape', 'KwaZulu-Natal'],
            'alleles_lost': lost_count if GENETIC_DATA else 'unknown'
        }
    }


def run_model_3(Ne, generations):
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
        N = [r['population_size'] for r in results]
        
        # Scale population to census size
        N = [n * (450 / 199) for n in N]
        
        t = np.arange(0, generations + 1)
        
        return {
            'model_number': 3,
            'model_name': 'Low Supplementation (+4 PAAZA/gen)',
            'generations': t.tolist(),
            'years': (t * 15).tolist(),
            'Ho': Ho,
            'He': He,
            'F': F,
            'Na': Na,
            'population': N,
            'Ne': Ne,
            'initial': {'Ho': Ho[0], 'He': He[0], 'Na': Na[0], 'N': N[0]},
            'parameters': {
                'Ne': Ne,
                'N0': int(N[0]),
                'lambda': 1.0,
                'data_source': 'CSV_simulation',
                'supplementation': '4 PAAZA birds per generation',
                'novel_alleles_added': len(GENETIC_DATA['novel_alleles']['paaza'])
            }
        }
    else:
        # Fallback: Generic supplementation model
        return run_generic_supplementation_model(Ne, generations, 
                                                  birds_per_gen=4, 
                                                  model_num=3,
                                                  source='PAAZA')


def run_model_4(Ne, generations):
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
        N = [r['population_size'] * (450 / 199) for r in results]
        
        t = np.arange(0, generations + 1)
        
        return {
            'model_number': 4,
            'model_name': 'High Supplementation (+10 PAAZA/gen)',
            'generations': t.tolist(),
            'years': (t * 15).tolist(),
            'Ho': Ho,
            'He': He,
            'F': F,
            'Na': Na,
            'population': N,
            'Ne': Ne,
            'initial': {'Ho': Ho[0], 'He': He[0], 'Na': Na[0], 'N': N[0]},
            'parameters': {
                'Ne': Ne,
                'N0': int(N[0]),
                'lambda': 1.0,
                'data_source': 'CSV_simulation',
                'supplementation': '10 PAAZA birds per generation'
            }
        }
    else:
        return run_generic_supplementation_model(Ne, generations,
                                                  birds_per_gen=10,
                                                  model_num=4,
                                                  source='PAAZA')


def run_model_5(Ne, generations):
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
        N = [r['population_size'] * (450 / 199) for r in results]
        
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
            'years': (t * 15).tolist(),
            'Ho': Ho,
            'He': He,
            'F': F,
            'Na': Na,
            'population': N,
            'Ne': Ne,
            'initial': {'Ho': Ho[0], 'He': He[0], 'Na': Na[0], 'N': N[0]},
            'parameters': {
                'Ne': Ne,
                'N0': int(N[0]),
                'lambda': 1.0,
                'data_source': 'CSV_simulation',
                'supplementation': '4 mixed birds (PAAZA/AZA/EAZA) per generation',
                'novel_alleles_total': novel_count
            }
        }
    else:
        return run_generic_supplementation_model(Ne, generations,
                                                  birds_per_gen=4,
                                                  model_num=5,
                                                  source='Mixed')


def run_generic_supplementation_model(Ne, generations, birds_per_gen, model_num, source):
    """
    Generic supplementation model when CSV data not available
    """
    H0 = 0.502
    He0 = 0.568
    A0 = 6.429
    N0 = 450
    lambda_val = 1.0
    
    t = np.arange(0, generations + 1)
    
    # Simple model: Ne increases with supplementation
    Ho_vals = []
    He_vals = []
    F_vals = []
    Na_vals = []
    N_vals = []
    
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
        N = N0 + birds_added
        
        Ho_vals.append(Ho)
        He_vals.append(He)
        F_vals.append(F)
        Na_vals.append(Na)
        N_vals.append(N)
    
    return {
        'model_number': model_num,
        'model_name': f'Supplementation (+{birds_per_gen} {source}/gen)',
        'generations': t.tolist(),
        'years': (t * 15).tolist(),
        'Ho': Ho_vals,
        'He': He_vals,
        'F': F_vals,
        'Na': Na_vals,
        'population': N_vals,
        'Ne': Ne,
        'initial': {'Ho': H0, 'He': He0, 'Na': A0, 'N': N0},
        'parameters': {
            'Ne': Ne,
            'N0': N0,
            'lambda': lambda_val,
            'data_source': 'generic',
            'supplementation': f'{birds_per_gen} {source} birds per generation'
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
        
        Ne = int(data.get('Ne', 100))
        generations = int(data.get('generations', 50))
        model = int(data.get('model', 1))
        
        # Validate inputs
        if not 20 <= Ne <= 200:
            return jsonify({'error': 'Ne must be between 20 and 200'}), 400
        if not 10 <= generations <= 100:
            return jsonify({'error': 'Generations must be between 10 and 100'}), 400
        if not 1 <= model <= 5:
            return jsonify({'error': 'Model must be between 1 and 5'}), 400
        
        # Run appropriate model
        if model == 1:
            results = run_model_1(Ne, generations)
        elif model == 2:
            results = run_model_2(Ne, generations)
        elif model == 3:
            results = run_model_3(Ne, generations)
        elif model == 4:
            results = run_model_4(Ne, generations)
        elif model == 5:
            results = run_model_5(Ne, generations)
        
        return jsonify(results)
    
    except Exception as e:
        import traceback
        traceback.print_exc()
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


if __name__ == '__main__':
    initialize_genetic_data()
    app.run(debug=True, port=5000)