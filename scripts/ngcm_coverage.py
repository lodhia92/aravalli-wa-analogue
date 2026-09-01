"""
Compile National Geochemical Mapping sample coverage from the survey packages.

The raw data lives in the directory resolved by paths.NGCM_RAW, organised by state, each
state containing toposheet folders with a *samples.metadata*.xlsx file. Sample coordinates
are in the "Samples" sheet, columns 6 (Latitude WGS84) and 7 (Longitude WGS84), below three
header rows.

Reading ~50+ Excel files at once is slow and can time out, so this runs per state (cache to
results/_cov_<state>.csv) and then combines. Usage:

    # one call per state (each is quick):
    python ngcm_coverage.py --state "2. Rajasthan"
    ... (repeat for every state folder) ...
    # then combine, recount domains, and regenerate the coverage figure:
    python ngcm_coverage.py --combine

--combine writes results/ngcm_all_coverage.csv, prints per-state counts and the counts inside
the Sandmata and Mangalwar domain polygons (so changes vs the last refresh are visible), and
regenerates figures/ngcm_coverage_updated.png.

Author: Bhavik Harish Lodhia, Curtin University
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import warnings
import zipfile

import pandas as pd

import paths

warnings.filterwarnings("ignore")

SKIPPED = []   # workbooks that could not be read, reported at the end

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(HERE)
RAW = paths.NGCM_RAW
RES = os.path.join(HERE, "results")


def load_state(state_dir):
    rows = []
    for f in glob.glob(os.path.join(RAW, state_dir, "**", "*samples.metadata*.xlsx"), recursive=True):
        try:
            df = pd.read_excel(f, sheet_name="Samples", header=None, skiprows=3,
                               usecols=[6, 7], names=["lat", "lon"], engine="openpyxl")
        except (ValueError, KeyError, OSError, zipfile.BadZipFile) as ex:
            # A workbook openpyxl cannot parse, or one with no "Samples" sheet. Named and
            # counted so that a survey file is never dropped from the coverage silently.
            SKIPPED.append(os.path.basename(f))
            print(f"  SKIPPED {os.path.basename(f)}: {type(ex).__name__}: {ex}")
            continue
        df["lat"] = pd.to_numeric(df.lat, errors="coerce")
        df["lon"] = pd.to_numeric(df.lon, errors="coerce")
        df = df.dropna()
        df = df[(df.lat > 5) & (df.lat < 40) & (df.lon > 65) & (df.lon < 100)]
        m = re.search(r"_([0-9]{2}[A-Z])", f)
        df["toposheet"] = m.group(1) if m else "?"
        rows.append(df)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["lat", "lon", "toposheet"])


def combine():
    from shapely.geometry import Point, shape
    from shapely.prepared import prep
    frames = []
    for f in glob.glob(os.path.join(RES, "_cov_*.csv")):
        st = os.path.basename(f).replace("_cov_", "").replace(".csv", "")
        d = pd.read_csv(f); d["state"] = st; frames.append(d)
    if not frames:
        sys.exit("No per-state caches found. Run --state for each state first.")
    A = pd.concat(frames, ignore_index=True)
    A.to_csv(os.path.join(RES, "ngcm_all_coverage.csv"), index=False)
    sand = prep(shape(json.load(open(os.path.join(RES, "sandmata_complex.geojson")))["features"][0]["geometry"]).buffer(0))
    mang = prep(shape(json.load(open(os.path.join(RES, "archaean_mangalwar_domain.geojson")))["features"][0]["geometry"]).buffer(0))

    def cnt(df, pp):
        s = df[(df.lon >= 73.8) & (df.lon <= 76.1) & (df.lat >= 23.9) & (df.lat <= 27.0)]
        return sum(pp.contains(Point(x, y)) for x, y in zip(s.lon, s.lat))

    print("TOTAL located samples:", len(A))
    for st, g in A.groupby("state"):
        print(f"  {st:22s} n={len(g):6d}  Sandmata={cnt(g, sand):5d}  Mangalwar={cnt(g, mang):5d}")
    print("TOTAL in Sandmata:", cnt(A, sand), " in Mangalwar:", cnt(A, mang))
    subprocess.run([sys.executable, os.path.join(HERE, "scripts", "plot_ngcm_coverage.py")])
    print("Wrote results/ngcm_all_coverage.csv and figures/ngcm_coverage_updated.png")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", help='state folder name, e.g. "2. Rajasthan"')
    ap.add_argument("--combine", action="store_true")
    args = ap.parse_args()
    if args.state:
        a = load_state(args.state)
        stn = re.sub(r"^[0-9]+\.\s*", "", args.state).replace(" ", "_")
        a.to_csv(os.path.join(RES, f"_cov_{stn}.csv"), index=False)
        print(f"{stn}: n={len(a)}"
              + (f" lon {a.lon.min():.1f}-{a.lon.max():.1f} lat {a.lat.min():.1f}-{a.lat.max():.1f}" if len(a) else ""))
    elif args.combine:
        combine()
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
