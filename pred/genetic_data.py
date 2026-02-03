"""
Backwards compatibility module.

All functionality has been moved to the genetics/ package.
This file re-exports everything for existing imports.

Usage:
    # Old import (still works):
    from genetic_data import LOCI, load_and_process_genetic_data

    # New import (preferred):
    from genetics import LOCI, load_and_process_genetic_data
"""

from genetics import *
