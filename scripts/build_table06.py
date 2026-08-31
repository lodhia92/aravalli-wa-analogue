"""
Build Table 6: raw stream-sediment geochemistry of the twenty Western Australian
(National Geochemical Survey of Australia) analogue samples, matched to Table 5.

The survey analyses three grain-size fractions at two sampling depths. The value reported here for
each site is the median across those analyses, which is the same value the transfer analysis uses
and is what Section 2.1 of the manuscript describes. Note that the survey's Bulk rows carry no
geochemistry, so a table built from the Bulk sample alone would be empty.

Total iron is reported by the survey under the header FeT rather than Fe. It is resolved here
explicitly, because a header match on Fe alone returns no column and leaves the row empty.

Inputs : results/analogue_pairs.csv, the NGSA table resolved by paths.NGSA
Output : results/table06_cells.tsv

Author: Bhavik Harish Lodhia, Curtin University
"""
import csv
import os

import pandas as pd

import paths

csv.field_size_limit(10 ** 7)
OUT = os.path.join(paths.RESULTS, "table06_cells.tsv")
DASH = "–"
ELS = ["Si", "Ti", "Al", "Fe", "Mn", "Mg", "Ca", "Na", "K", "P", "Sc", "V", "Cr", "Co", "Ni",
       "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
       "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Yb", "Lu", "Hf", "Ta", "Pb", "Th", "U"]
STEM = {"Fe": "FeT"}          # the survey reports total iron under FeT


def find(header, el):
    stem = STEM.get(el, el)
    for meth in ("XRF", "ICP-MS"):
        for i, c in enumerate(header):
            if c.strip().startswith(f"{stem} {meth}"):
                return i
    return None


def fmt(v):
    if pd.isna(v):
        return DASH
    return f"{v:.2f}" if abs(v) < 100 else f"{v:.1f}"


def main():
    paths.require("NGSA")
    pairs = pd.read_csv(os.path.join(paths.RESULTS, "analogue_pairs.csv"))
    sites = [str(r.aus_sid).split("_")[-1].split(".")[0] for r in pairs.itertuples()]
    with open(paths.NGSA, encoding="latin-1") as fh:
        header = [c.strip() for c in list(csv.reader(fh))[11]]
    idx = {el: find(header, el) for el in ELS}
    missing = [el for el, i in idx.items() if i is None]
    assert not missing, f"no total-suite column for {missing}"
    use = sorted(set(idx.values()) | {0, 7})
    ng = pd.read_csv(paths.NGSA, header=None, skiprows=12, usecols=use, encoding="latin-1",
                     low_memory=False)
    ng.columns = [header[i] for i in use]
    ng["SITEID"] = pd.to_numeric(ng[header[0]], errors="coerce").astype("Int64").astype(str)
    ng["SAMPLEID"] = ng[header[7]].astype(str).str.strip()
    num = ng.copy()
    for el in ELS:
        num[header[idx[el]]] = pd.to_numeric(num[header[idx[el]]], errors="coerce")
    med = num.groupby("SITEID")[[header[idx[el]] for el in ELS]].median()
    cols = []
    for s in sites:
        assert s in med.index, f"site {s} not found in the survey table"
        cols.append(med.loc[s])
    rows = [["Variable"] + [str(i) for i in range(1, len(cols) + 1)]]
    for el in ELS:
        name = header[idx[el]]
        rows.append([f"{el} (ppm)"] +
                    [fmt(c[name]) for c in cols])
    with open(OUT, "w", newline="", encoding="utf8") as fh:
        csv.writer(fh, delimiter="\t").writerows(rows)
    print(f"wrote {OUT} ({len(rows) - 1} variables x {len(cols)} samples)")


if __name__ == "__main__":
    main()
