"""
Build Table 5: raw stream-sediment geochemistry of the twenty north-west Indian
(National Geochemical Mapping) analogue samples, one column per analogue pair.

Values are reported as supplied by the survey. A value the survey did not report, whether
because the element was not determined or because the result was below detection, is written
as an en dash and is never replaced by a substitute concentration.

Inputs : results/analogue_pairs.csv, the NGCM table resolved by paths.NGCM_TABLE
Output : results/table05_cells.tsv

Author: Bhavik Harish Lodhia, Curtin University
"""

import csv
import os

import pandas as pd

import paths

OUT = os.path.join(paths.RESULTS, "table05_cells.tsv")
DASH = "–"
OXIDES = ["SiO2", "TiO2", "Al2O3", "Fe2O3", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5"]
TRACE = [
    "Sc",
    "V",
    "Cr",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "Pb",
    "Th",
    "U",
]


def label(el):
    return f"{el} (wt%)" if el in OXIDES else f"{el} (ppm)"


def fmt(v):
    if pd.isna(v):
        return DASH
    return f"{v:.2f}" if abs(v) < 100 else f"{v:.1f}"


def main():
    paths.require("NGCM_TABLE")
    pairs = pd.read_csv(os.path.join(paths.RESULTS, "analogue_pairs.csv"))
    ar = pd.read_csv(paths.NGCM_TABLE, low_memory=False)
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)
    cols = []
    for r in pairs.itertuples():
        lat, lon = r.india_sid.split("_")[-2:]
        k = f"{round(float(lat), 4)}_{round(float(lon), 4)}"
        m = ar[ar.k == k]
        assert len(m) == 1, f"expected one NGCM row for {k}, found {len(m)}"
        cols.append(m.iloc[0])
    rows = [["Variable"] + [str(i) for i in range(1, len(cols) + 1)]]
    for el in OXIDES + TRACE:
        rows.append([label(el)] + [fmt(pd.to_numeric(c.get(el), errors="coerce")) for c in cols])
    with open(OUT, "w", newline="", encoding="utf8") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)
    print(f"wrote {OUT} ({len(rows) - 1} variables x {len(cols)} samples)")


if __name__ == "__main__":
    main()
