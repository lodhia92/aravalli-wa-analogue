"""
Shared statistics, composition and survey-reading code for the Aravalli-Western Australia
analogue study.

Every script in scripts/ imports its rank statistics, its false-discovery-rate correction, its
provenance ratios and its survey column lookup from here, so that the matching, the rank
statistics and the correction cannot drift between analyses. Before this package existed the
same helpers were defined independently in up to ten files.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
