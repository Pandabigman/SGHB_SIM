"""
Caching utilities for expensive genetic simulations.
"""

import logging

logger = logging.getLogger(__name__)

# Global cache for simulation results
SIMULATION_CACHE = {}


def get_cache_key(model_num, birds_per_gen, generations, base_Ne, stochastic=False):
    """
    Generate cache key for genetic simulation results.

    Args:
        model_num: Model number
        birds_per_gen: Birds added per generation
        generations: Number of generations
        base_Ne: Base effective population size
        stochastic: Whether simulation is stochastic

    Returns: String cache key
    """
    stoch_suffix = "_stoch" if stochastic else ""
    return f"model{model_num}_birds{birds_per_gen}_gen{generations}_Ne{base_Ne}{stoch_suffix}"


def get_cached_simulation(cache_key, simulation_func, stochastic=False, *args, **kwargs):
    """
    Cache wrapper for expensive genetic simulations.

    Args:
        cache_key: Unique identifier for this simulation
        simulation_func: The function to call if cache miss
        stochastic: If True, skip caching (each run should be unique)
        *args, **kwargs: Arguments to pass to simulation_func

    Returns: Simulation results (from cache or freshly computed)
    """
    # IMPORTANT: Don't cache stochastic simulations - they should vary each run
    if stochastic:
        logger.info(f"Stochastic mode - skipping cache")
        return simulation_func(*args, **kwargs)

    if cache_key in SIMULATION_CACHE:
        logger.info(f"Cache HIT for {cache_key}")
        return SIMULATION_CACHE[cache_key]

    logger.info(f"Cache MISS for {cache_key}, computing...")
    result = simulation_func(*args, **kwargs)
    SIMULATION_CACHE[cache_key] = result
    return result


def clear_cache():
    """Clear the simulation cache."""
    global SIMULATION_CACHE
    SIMULATION_CACHE = {}
    logger.info("Simulation cache cleared")
