"""
Analytical constants shared across the workflow, with their units.

These are the values the published analysis used. Changing any of them changes a published
number, so they are declared once here rather than repeated in each script.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""

# Permutation convention. A complete unmodified run of a script reproduces its published
# probabilities exactly; see README.md for the limitation this wording is careful about.
NPERM = 100000          # relabellings per permutation test
SEED = 20260827         # generator seed

# Catchment containment cut-offs, percent of upstream catchment area inside the domain.
# India is densely sampled with small catchments; Australia is sparsely sampled with large
# ones, so the Australian cut-off is relaxed.
IN_THR_PCT = 50.0
AUS_THR_PCT = 25.0

# CI chondrite reference concentrations, ppm, for the europium anomaly and the normalised
# lanthanum/ytterbium ratio.
CI_CHONDRITE_PPM = dict(La=.237, Yb=.170, Sm=.148, Eu=.0580, Gd=.199)

# The six weathering-robust provenance ratios used to match catchments.
MATCH_RATIOS = ["Th/Sc", "La/Sc", "Th/Co", "EuEu", "La/Yb_n", "Nb/Y"]

# The eight pathfinder elements whose transfer between matched catchments is tested.
PATHFINDER_ELEMENTS = ["Zr", "Hf", "Ti", "Ce", "Nd", "Pr", "Dy", "P"]

# Oxide to element conversion. The survey reports Ti as TiO2 and P as P2O5 in weight per cent;
# the analysis uses element concentrations in mg/kg. The factors are the element's mass fraction
# of the oxide (Ti 47.867/79.866, P 2x30.974/141.944, Al 2x26.982/101.961) and one weight per
# cent is ten thousand mg/kg.
WT_PCT_TO_MG_KG = 1e4
TI_MASS_FRACTION_OF_TIO2 = 0.5995
P_MASS_FRACTION_OF_P2O5 = 0.4364
AL_MASS_FRACTION_OF_AL2O3 = 0.5293
