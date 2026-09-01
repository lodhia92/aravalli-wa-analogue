"""
Build Table 4: identifiers and coordinates of the forty analogue samples.

India samples are identified by latitude and longitude, because the National Geochemical Mapping
release carries no stable public sample identifier for these composites. Australia samples are
identified by their National Geochemical Survey of Australia site and sample identifiers.

Input : results/analogue_pairs.csv
Output: results/table04_cells.tsv

Author: Bhavik Harish Lodhia, Curtin University
"""

import csv
import os

import pandas as pd

import paths

OUT = os.path.join(paths.RESULTS, "table04_cells.tsv")
HEAD = ["Pair", "Domain", "India lat", "India lon", "NGSA site", "NGSA sample"]


def main():
    p = pd.read_csv(os.path.join(paths.RESULTS, "analogue_pairs.csv"))
    rows = []
    for n, r in enumerate(p.itertuples(), start=1):
        lat, lon = r.india_sid.split("_")[-2:]
        site = str(r.aus_sid).split("_")[-1].split(".")[0]
        rows.append(
            [str(n), r.domain, f"{float(lat):.3f}", f"{float(lon):.3f}", site, f"{site}001 Bulk"]
        )
    with open(OUT, "w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEAD)
        w.writerows(rows)
    assert len(rows) == 20, f"expected 20 pairs, built {len(rows)}"
    assert all(len(r) == len(HEAD) for r in rows)
    print(f"wrote {OUT} ({len(rows)} pairs, {len(rows) * 2} samples)")


if __name__ == "__main__":
    main()
