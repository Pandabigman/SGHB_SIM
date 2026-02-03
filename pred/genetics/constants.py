"""
Genetic constants and default parameters for SGHB simulation.
"""

# Microsatellite loci (14 total)
LOCI = ['Buco4', 'Buco11', 'Buco2', 'Buco9', 'GHB21', 'GHB19', 'GHB26',
        'GHB20', 'Buco16', 'Buco18', 'Buco19', 'Buco21', 'Buco24', 'Buco25']

# Default genetic parameters (fallback when CSV not available)
DEFAULT_H0 = 0.502
DEFAULT_HE0 = 0.568
DEFAULT_A0 = 6.429
DEFAULT_N0 = 2500

# Biological constants
LETHAL_EQUIVALENTS_BIRDS = 3.14  # O'Grady et al. 2006
GENERATION_TIME_YEARS = 26
