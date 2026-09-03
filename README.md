# Critical mineral pathfinder targeting in the Aravalli basement

Analysis code for the study that matches north-west Indian and Western Australian stream-sediment
catchments on weathering-robust provenance ratios, tests whether critical mineral pathfinder
element signatures transfer between the matched catchments, and scores the Aravalli catchments
with the signatures that transfer.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au

## What this repository contains

Every analysis behind the published tables and the reported statistics, from the drainage
selection through to the sensitivity tests. It contains no figure code: the manuscript figures
are published with the paper and supplied to the journal separately. It does not redistribute
any of the source datasets. It does carry the published tables themselves, and four small
derived files that the pipeline reads and no script here writes; `data/README.md` lists those
and says where each came from. The run order below gives each step, the script that performs it
and what that script produces.

## Before you start

Install the requirements and put the input data in place:

```
pip install -r requirements.txt
python scripts/verify_reproduction.py
```

`data/README.md` lists every dataset, its provider and the path the code expects. The
verification script confirms that your copy reproduces the published result: twenty analogue
pairs, a light rare-earth transfer correlation of 0.696 and a heavy rare-earth correlation of
0.506. Run it first. If it passes, everything downstream is running on the same data as the paper.

## Run order

Stage 1, study areas and drainage

| Script | Produces |
|---|---|
| `prepare_ngcm.py` | the Indian geochemistry table from the survey packages |
| `extract_sandmata_polygon.py` | the Sandmata Complex domain polygon, from a source figure that is not redistributed here; the polygon it writes is supplied instead |
| `ngcm_coverage.py` | Indian sample coverage, per state and per domain |
| `drainage_containment.py` | for every sample, the fraction of its upstream catchment inside a domain |
| `pair_catchments.py` | upstream catchment polygons for the analogue samples |

Stage 2, matching

| Script | Produces |
|---|---|
| `provenance_ratios.py` | the four drainage-defined sample sets and their provenance ratios |
| `analogue_pairing.py` | **the twenty analogue pairs** |
| `pair_validation.py` | per-element and per-fingerprint correlations across the pairs |

Stage 3, transfer and scoring

| Script | Produces |
|---|---|
| `fingerprint_transfer.py` | **Tables 2 and 7**: transfer correlations, corrected probabilities, resampling |
| `catchment_scores.py` | per-catchment pathfinder scores and top-decile flags |

Stage 4, catchment geology

| Script | Produces |
|---|---|
| `catchment_geology.py` | **Table 3**: area fraction of each catchment underlain by each mapped unit, in `results/catchment_geology.csv` |
| `terrane_attribution.py` | tectonic unit of every Australian sample and pair member |
| `ree_occurrences.py` | documented rare-earth occurrences inside the analogue catchments |

Stage 5, sensitivity and validation

| Script | Tests |
|---|---|
| `matching_sensitivity.py` | whether the matching algorithm manufactures the transfer |
| `matcher_independence.py` | matchers built only from elements absent from the target minerals |
| `standardisation_invariance.py` | dependence of the result on the within-survey standardisation |
| `grain_size_morphometry.py` | grain size, catchment area and transport distance |
| `containment_sensitivity.py` | the two catchment containment thresholds |
| `lithological_mixing.py` | mixed lithologies within the catchments |
| `score_threshold_sensitivity.py` | the top-decile rule used to flag high-scoring catchments |
| `hree_transfer.py` | broader heavy rare-earth fingerprints |
| `broadened_fingerprints.py` | broader fingerprints for the remaining host minerals |
| `mineralogical_validation.py` | fingerprints against measured heavy-mineral grain counts |
| `alternative_fingerprints.py` | alternative element sets against the same grain counts |
| `data_quality.py` | analytical methods, detection limits, censored and missing values |

Stage 6, table content

| Script | Produces |
|---|---|
| `build_table04.py` | Table 4, sample identifiers and coordinates |
| `build_table05.py` | Table 5, Indian raw geochemistry |
| `build_table06.py` | Table 6, Australian raw geochemistry |
| `build_table08.py` | Table 8, matching sensitivity |
| `build_table09.py` | Table 9, mineralogical validation |

`matching_sensitivity.py` is also imported as a library by the four scripts that follow it, so
that the matching, the rank statistics and the false-discovery-rate correction cannot drift
between analyses.

## Conventions

Permutation probabilities use 100 000 relabellings with a fixed seed, so a rerun reproduces the
published probabilities exactly. Multiple comparisons are corrected by the Benjamini-Hochberg
procedure across the five fingerprints, which is the family the paper reports. Values below
detection and values not reported are treated as missing and are never replaced by a substitute
concentration; the matching and the scoring use only samples carrying the complete element suite.

## Licence

Code: MIT, see LICENSE. The input datasets are not covered by that licence and remain under the
terms of their providers, listed in `data/README.md`.
