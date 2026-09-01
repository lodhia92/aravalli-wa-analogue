"""
Appendix figure: the analysis workflow, keyed to the scripts in the public repository.

Boxes carry a short code and a two or three word descriptor. The key on the right gives the
code, the script filename and what the script produces, so the boxes stay small enough for the
text to remain legible at print size.

Author: Bhavik Harish Lodhia, Curtin University
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ----------------------------------------------------------------------------------
# content
# ----------------------------------------------------------------------------------
# stage: (title, ncols, [(code, box label, script, produces)])
# a script of "" marks a dataset box, which is drawn in the flow but not listed in the key
STAGES = [
    (
        "0  Datasets",
        5,
        [
            ("NGCM", "Geological Survey\nof India", "", ""),
            ("NGSA", "Geoscience\nAustralia", "", ""),
            ("HMMA", "Heavy Mineral Map\nof Australia", "", ""),
            ("HYDRO", "HydroSHEDS\nHydroBASINS L12", "", ""),
            ("GSWA", "Geological Survey of\nWestern Australia", "", ""),
        ],
    ),
    (
        "1  Inputs and study areas",
        5,
        [
            (
                "D1",
                "Indian\ngeochemistry",
                "prepare_ngcm.py",
                "one table per sample from the survey packages",
            ),
            ("D2", "Domain\npolygon", "extract_sandmata_polygon.py", "Sandmata Complex boundary"),
            ("D3", "Sample\ncoverage", "ngcm_coverage.py", "Indian coverage by state and domain"),
            (
                "D4",
                "Drainage\ncontainment",
                "drainage_containment.py",
                "catchment fraction inside each domain",
            ),
            (
                "D5",
                "Catchment\npolygons",
                "pair_catchments.py",
                "upstream catchments of the analogue samples",
            ),
        ],
    ),
    (
        "2  Analogue matching",
        3,
        [
            (
                "M1",
                "Provenance\nratios",
                "provenance_ratios.py",
                "the four drainage-defined sample sets",
            ),
            ("M2", "Analogue\npairing", "analogue_pairing.py", "the twenty analogue pairs"),
            (
                "M3",
                "Pair\nvalidation",
                "pair_validation.py",
                "per-element correlations across the pairs",
            ),
        ],
    ),
    (
        "3  Transfer and scoring",
        2,
        [
            (
                "T1",
                "Fingerprint\ntransfer",
                "fingerprint_transfer.py",
                "Tables 2 and 7; correlations and resampling",
            ),
            (
                "T2",
                "Catchment\nscores",
                "catchment_scores.py",
                "per-catchment scores and top-decile flags",
            ),
        ],
    ),
    (
        "4  Catchment geology",
        3,
        [
            (
                "G1",
                "Catchment\ngeology",
                "catchment_geology.py",
                "Table 3; area fraction per mapped unit",
            ),
            (
                "G2",
                "Terrane\nattribution",
                "terrane_attribution.py",
                "tectonic unit of every Australian sample",
            ),
            (
                "G3",
                "Rare-earth\noccurrences",
                "ree_occurrences.py",
                "occurrences inside the analogue catchments",
            ),
        ],
    ),
    (
        "5  Sensitivity and validation",
        4,
        [
            (
                "S1",
                "Matching\nsensitivity",
                "matching_sensitivity.py",
                "matching variants and random control",
            ),
            (
                "S2",
                "Matcher\nindependence",
                "matcher_independence.py",
                "matchers free of the target minerals",
            ),
            (
                "S3",
                "Standardisation",
                "standardisation_invariance.py",
                "dependence on within-survey scaling",
            ),
            (
                "S4",
                "Grain size and\nmorphometry",
                "grain_size_morphometry.py",
                "fraction, area and transport distance",
            ),
            (
                "S5",
                "Containment\nthresholds",
                "containment_sensitivity.py",
                "the two catchment thresholds",
            ),
            (
                "S6",
                "Lithological\nmixing",
                "lithological_mixing.py",
                "mixed lithologies within catchments",
            ),
            ("S7", "Score\nthreshold", "score_threshold_sensitivity.py", "the top-decile rule"),
            (
                "S8",
                "Heavy rare-earth\ntransfer",
                "hree_transfer.py",
                "broader heavy rare-earth fingerprints",
            ),
            (
                "S9",
                "Broadened\nfingerprints",
                "broadened_fingerprints.py",
                "broader sets for the other hosts",
            ),
            (
                "S10",
                "Mineralogical\nvalidation",
                "mineralogical_validation.py",
                "fingerprints against measured grains",
            ),
            (
                "S11",
                "Alternative\nfingerprints",
                "alternative_fingerprints.py",
                "alternative sets against the same grains",
            ),
            ("S12", "Data\nquality", "data_quality.py", "methods, detection limits, missing data"),
        ],
    ),
    (
        "6  Table content",
        5,
        [
            ("B1", "Table 4", "build_table04.py", "sample identifiers and coordinates"),
            ("B2", "Table 5", "build_table05.py", "Indian raw geochemistry"),
            ("B3", "Table 6", "build_table06.py", "Australian raw geochemistry"),
            ("B4", "Table 8", "build_table08.py", "matching sensitivity"),
            ("B5", "Table 9", "build_table09.py", "mineralogical validation"),
        ],
    ),
]

HILITE = set()

# ----------------------------------------------------------------------------------
# palette
# ----------------------------------------------------------------------------------
INK = "#1a1a1a"
RULE = "#9aa4ad"
SRC_FILL = "#eceff1"
SRC_EDGE = "#7d8a95"
BOX_FILL = "#ffffff"
BOX_EDGE = "#4a5560"
HI_FILL = "#dfe8ef"
HI_EDGE = "#1f3a4d"
BAND = "#f5f7f8"
KEYRULE = "#c7ced4"

FS_TITLE = 11.5
FS_STAGE = 10.5
FS_CODE = 11
FS_LABEL = 8.2
FS_SRC = 8.0
FS_KEY = 7.6
FS_KEYH = 8.6

# ----------------------------------------------------------------------------------
# layout
# ----------------------------------------------------------------------------------
FIGW, FIGH = 16.0, 10.4
fig = plt.figure(figsize=(FIGW, FIGH), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

LEFT, RIGHT = 0.035, 0.560  # flow panel
KEY_L, KEY_R = 0.600, 0.985  # key panel
TOP = 0.955
BOT = 0.035

# vertical budget: source row + 6 stages (stage 5 takes 3 rows of boxes)
GRID = 5
rows_per_stage = [1, 1, 1, 1, 1, 3, 1]
SRC_H = 0.072
GAP = 0.020
STAGE_HEAD = 0.028
BOX_H = 0.072
BOX_VGAP = 0.012


def stage_height(nr):
    return STAGE_HEAD + nr * BOX_H + (nr - 1) * BOX_VGAP + 0.012


total = sum(stage_height(n) + GAP for n in rows_per_stage) - GAP
scale = (TOP - BOT) / total
SRC_H *= scale
GAP *= scale
STAGE_HEAD *= scale
BOX_H *= scale
BOX_VGAP *= scale


def rbox(x, y, w, h, fc, ec, lw=1.0, r=0.008):
    ax.add_patch(
        FancyBboxPatch(
            (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=3
        )
    )


def arrow(x, y0, y1):
    ax.add_patch(
        FancyArrowPatch(
            (x, y0),
            (x, y1),
            arrowstyle="-|>",
            mutation_scale=11,
            lw=1.1,
            color=RULE,
            zorder=2,
            shrinkA=0,
            shrinkB=0,
        )
    )


cursor = TOP

for si, (title, ncol, items) in enumerate(STAGES):
    nr = rows_per_stage[si]
    h = STAGE_HEAD + nr * BOX_H + (nr - 1) * BOX_VGAP + 0.012 * scale
    top_y = cursor
    bot_y = cursor - h
    ax.add_patch(
        FancyBboxPatch(
            (LEFT - 0.008, bot_y),
            (RIGHT - LEFT) + 0.016,
            h,
            boxstyle="round,pad=0,rounding_size=0.006",
            fc=BAND,
            ec="none",
            zorder=1,
        )
    )
    ax.text(
        LEFT,
        top_y - STAGE_HEAD * 0.55,
        title,
        fontsize=FS_STAGE,
        fontweight="bold",
        color="#2c3a46",
        va="center",
        zorder=4,
    )

    gx = 0.010
    bw = (RIGHT - LEFT - gx * (GRID - 1)) / GRID
    for k, (code, label, _, _) in enumerate(items):
        r, c = divmod(k, GRID)
        bx = LEFT + c * (bw + gx)
        by = top_y - STAGE_HEAD - (r + 1) * BOX_H - r * BOX_VGAP
        data = items[k][2] == ""
        rbox(
            bx,
            by,
            bw,
            BOX_H,
            SRC_FILL if data else BOX_FILL,
            SRC_EDGE if data else BOX_EDGE,
            0.9 if data else 1.0,
        )
        ax.text(
            bx + bw / 2,
            by + BOX_H * 0.71,
            code,
            ha="center",
            va="center",
            fontsize=FS_CODE - (0.5 if data else 0),
            fontweight="bold",
            color=INK,
            zorder=4,
        )
        ax.text(
            bx + bw / 2,
            by + BOX_H * 0.27,
            label,
            ha="center",
            va="center",
            fontsize=FS_SRC if data else FS_LABEL,
            color="#41505c" if data else "#3b4750",
            linespacing=1.22,
            zorder=4,
        )
    cursor = bot_y - GAP

# ---- key
ky = TOP + 0.012

CODE_X = KEY_L
NAME_X = KEY_L + 0.038
DESC_X = KEY_L + 0.192
LINE_H = 0.0207
HEAD_H = 0.026

for title, _, items in STAGES:
    if all(i[2] == "" for i in items):
        continue
    ky -= HEAD_H
    ax.text(CODE_X, ky, title, fontsize=FS_KEYH, fontweight="bold", color="#2c3a46", va="center")
    ky -= 0.004
    ax.plot([KEY_L, KEY_R], [ky, ky], color=KEYRULE, lw=0.5)
    for code, _, script, produces in items:
        ky -= LINE_H
        hot = code in HILITE
        ax.text(
            CODE_X,
            ky,
            code,
            fontsize=FS_KEY,
            fontweight="bold",
            color=HI_EDGE if hot else INK,
            va="center",
        )
        ax.text(
            NAME_X,
            ky,
            script,
            fontsize=FS_KEY,
            family="DejaVu Sans Mono",
            color="#1f3a4d" if hot else "#243039",
            va="center",
        )
        ax.text(DESC_X, ky, produces, fontsize=FS_KEY, color="#4b5763", va="center")
    ky -= 0.006

ky -= 0.010
ky -= 0.018
ax.text(
    CODE_X,
    ky,
    "Utilities, not part of the flow",
    fontsize=FS_KEYH,
    fontweight="bold",
    color="#2c3a46",
    va="center",
)
ky -= 0.004
ax.plot([KEY_L, KEY_R], [ky, ky], color=KEYRULE, lw=0.5)
ky -= 0.018
ax.text(CODE_X, ky, "V", fontsize=FS_KEY, fontweight="bold", color=INK, va="center")
ax.text(
    NAME_X,
    ky,
    "verify_reproduction.py",
    fontsize=FS_KEY,
    family="DejaVu Sans Mono",
    color="#243039",
    va="center",
)
ax.text(
    DESC_X,
    ky,
    "self-test: 20 pairs, rho = 0.696 and 0.506",
    fontsize=FS_KEY,
    color="#4b5763",
    va="center",
)
ky -= LINE_H
ax.text(
    NAME_X, ky, "paths.py", fontsize=FS_KEY, family="DejaVu Sans Mono", color="#243039", va="center"
)
ax.text(
    DESC_X, ky, "resolves every external dataset", fontsize=FS_KEY, color="#4b5763", va="center"
)

fig.savefig("/tmp/fig/workflow.png", dpi=200, facecolor="white")
print("wrote /tmp/fig/workflow.png")
