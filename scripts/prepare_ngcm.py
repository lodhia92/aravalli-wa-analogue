"""
Build the Indian geochemistry table from the National Geochemical Mapping survey packages.

The survey supplies each toposheet as a set of workbooks split by analytical method: an X-ray
fluorescence package carrying the major oxides and the higher-abundance trace elements, an
inductively coupled plasma mass spectrometry package carrying the rare earths and the remaining
trace elements, a package covering the elements determined by other methods, and a sample
metadata workbook carrying the coordinates.

This script merges those packages into one row per sample, with LAT, LON and one column per
element or oxide, which is the table every downstream script reads. Elements not determined for
a sample are left empty and are never filled with a substitute value. No spatial selection is
performed here: the drainage containment step selects the samples belonging to each basement
domain.

    python scripts/prepare_ngcm.py --state "2. Rajasthan"     # one state at a time, then
    python scripts/prepare_ngcm.py --combine

Output: the table at paths.NGCM_TABLE

Author: Bhavik Harish Lodhia, Curtin University
"""

import argparse
import glob
import os
import re
import warnings
import zipfile

import pandas as pd

import paths

warnings.filterwarnings("ignore")

SKIPPED = []  # workbooks that could not be read, reported at the end
CACHE = os.path.join(paths.RESULTS, "_ngcm_{}.csv")
PACKAGES = ("*package A*XRF*.xlsx", "*package*ICPMS*.xlsx", "*package B*Other*.xlsx")
ELEMENT = re.compile(r"^[A-Z][a-z]?[0-9]?[A-Za-z0-9]*$")
SKIP = {"Domain", "Comment", "Total", "Datapoint Key", "Aliquot ID", "Spot ID"}


def sheet_elements(path):
    """One row per sample key, one column per element reported in this package."""
    dp = pd.read_excel(path, sheet_name="GC Datapoints", header=2, engine="openpyxl")
    cc = pd.read_excel(path, sheet_name="GC Concentrations", header=2, engine="openpyxl")
    keep = [
        c
        for c in cc.columns
        if str(c).strip() not in SKIP
        and ELEMENT.match(str(c).strip())
        and "Uncertainty" not in str(c)
        and "Measured Mass" not in str(c)
    ]
    out = cc[keep].copy()
    out.columns = [str(c).strip() for c in keep]
    out["sample"] = dp["Sample"].values[: len(out)]
    return out


def load_state(state):
    frames, meta = [], []
    for sheet_dir in sorted(glob.glob(os.path.join(paths.NGCM_RAW, state, "*"))):
        if not os.path.isdir(sheet_dir):
            continue
        toposheet = os.path.basename(sheet_dir)
        parts = []
        for pattern in PACKAGES:
            for f in glob.glob(os.path.join(sheet_dir, pattern)):
                try:
                    parts.append(sheet_elements(f))
                except (ValueError, KeyError, OSError, zipfile.BadZipFile) as ex:
                    # A package workbook openpyxl cannot parse, or one missing a required
                    # sheet. Named and counted so that no survey file is dropped silently.
                    SKIPPED.append(os.path.basename(f))
                    print(f"  SKIPPED {os.path.basename(f)}: {type(ex).__name__}: {ex}")
        if not parts:
            continue
        merged = parts[0]
        for p in parts[1:]:
            merged = merged.merge(p, on="sample", how="outer", suffixes=("", "_dup"))
        merged = merged[[c for c in merged.columns if not c.endswith("_dup")]]
        merged["toposheet"] = toposheet
        frames.append(merged)
        for f in glob.glob(os.path.join(sheet_dir, "*samples.metadata*.xlsx")):
            m = pd.read_excel(
                f,
                sheet_name="Samples",
                header=None,
                skiprows=3,
                usecols=[0, 6, 7],
                names=["sample", "LAT", "LON"],
                engine="openpyxl",
            )
            meta.append(m)
    if not frames:
        raise SystemExit(f"no survey packages found under {os.path.join(paths.NGCM_RAW, state)}")
    d = pd.concat(frames, ignore_index=True)
    if meta:
        mm = pd.concat(meta, ignore_index=True).dropna(subset=["sample"])
        d = d.merge(mm, on="sample", how="left")
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state")
    ap.add_argument("--combine", action="store_true")
    a = ap.parse_args()
    os.makedirs(paths.RESULTS, exist_ok=True)
    if a.state:
        d = load_state(a.state)
        out = CACHE.format(re.sub(r"\W+", "_", a.state))
        d.to_csv(out, index=False)
        print(f"{a.state}: {len(d)} samples -> {out}")
        return
    if a.combine:
        files = sorted(glob.glob(CACHE.format("*")))
        if not files:
            raise SystemExit("no per-state caches; run --state for each state first")
        d = pd.concat([pd.read_csv(f, low_memory=False) for f in files], ignore_index=True)
        for c in d.columns:
            if c not in ("sample", "toposheet"):
                d[c] = pd.to_numeric(d[c], errors="coerce")
        d = d.dropna(subset=["LAT", "LON"])
        os.makedirs(os.path.dirname(paths.NGCM_TABLE), exist_ok=True)
        d.to_csv(paths.NGCM_TABLE, index=False)
        print(f"combined {len(files)} states, {len(d)} samples -> {paths.NGCM_TABLE}")
        return
    ap.error("give --state <name> or --combine")


if __name__ == "__main__":
    main()
