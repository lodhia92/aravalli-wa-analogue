"""
Assemble the four drainage-defined sample sets and test the analogue on
weathering-robust provenance ratios (Th/Sc, Nb/Y, (La/Yb)n).

Pairs:
  Sandmata (India NGCM)  <->  Halls Creek (Australia NGSA)         [Palaeoproterozoic]
  Mangalwar (India NGCM) <->  Youanmi+Narryer (Australia NGSA)     [Archaean]

India geochem: the NGCM table (single-element columns), joined to the drainage selection on
lon/lat. Australia geochem: NGSA_data.csv ICP-MS columns, joined on SITEID.

(La/Yb)n uses CI chondrite (McDonough & Sun 1995): La 0.237, Yb 0.170 ppm.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import os, csv

import numpy as np

import pandas as pd

import paths

def ratios(df):
    out = pd.DataFrame(index=df.index)
    out["Th_Sc"] = df["Th"] / df["Sc"].replace(0, np.nan)
    out["Nb_Y"] = df["Nb"] / df["Y"].replace(0, np.nan)
    out["LaYb_n"] = (df["La"] / df["Yb"].replace(0, np.nan)) / (0.237 / 0.170)
    return out

def india_set(name):
    d = pd.read_csv(os.path.join(DR, f"{name}_contained.csv"))
    d = d[d.pct_in_domain >= THRESH]
    d["k"] = d.lat.round(4).astype(str) + "_" + d.lon.round(4).astype(str)
    m = ar[ar.k.isin(set(d.k))].copy()
    return m

def icp(el):
    for i, c in enumerate(header):
        if c.strip().startswith(f"{el} ICP-MS"):
            return i
    return None

def aus_set(path, pctcol):
    d = pd.read_csv(path); d = d[d[pctcol] >= THRESH]
    ids = set(pd.to_numeric(d["id"] if "id" in d.columns else d["SITEID"], errors="coerce").dropna())
    return ng[ng.SITEID.isin(ids)].copy()

def main():
    global DR, HERE, PROJ, RES, THRESH, _, allr, ar, cols, df, el, f, hc, header, k, mang, ng, nm, ordered, r, rows, sand, sets, summary, usecols, v, yil
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ = os.path.dirname(HERE); RES = os.path.join(HERE, "results"); DR = os.path.join(RES, "drainage")
    PROJ = os.path.dirname(HERE); RES = os.path.join(HERE, "results"); DR = os.path.join(RES, "drainage")
    PROJ = os.path.dirname(HERE); RES = os.path.join(HERE, "results"); DR = os.path.join(RES, "drainage")
    THRESH = 50.0  # primary catchment-containment threshold (%)
    ar = pd.read_csv(paths.NGCM_TABLE)
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)
    sand = india_set("sandmata"); mang = india_set("mangalwar")
    sand = india_set("sandmata"); mang = india_set("mangalwar")
    with open(paths.NGSA, encoding="latin-1") as f:
        header = list(csv.reader(f))[11]  # row 12 = column names
    cols = {el: icp(el) for el in ["Th", "Sc", "Nb", "Y", "La", "Yb"]}
    ordered = sorted(cols.items(), key=lambda kv: kv[1])   # pandas returns usecols in ascending file order
    usecols = [0] + [v for _, v in ordered]
    ng = pd.read_csv(paths.NGSA, header=None,
                     skiprows=12, usecols=usecols, encoding="latin-1", low_memory=False)
    ng.columns = ["SITEID"] + [k for k, _ in ordered]
    for el in cols: ng[el] = pd.to_numeric(ng[el], errors="coerce")
    ng["SITEID"] = pd.to_numeric(ng["SITEID"], errors="coerce")
    hc = aus_set(os.path.join(RES, "ngsa_draining_hallscreek.csv"), "pct_in_hallscreek")
    yil = aus_set(os.path.join(DR, "yilgarn_contained.csv"), "pct_in_domain")
    sets = {"Sandmata (India)": sand, "Halls Creek (Aus)": hc,
            "Mangalwar (India)": mang, "Youanmi+Narryer (Aus)": yil}
    rows = []
    for nm, df in sets.items():
        r = ratios(df)
        rows.append({"set": nm, "n": len(df),
                     "Th/Sc med": round(r.Th_Sc.median(), 3), "Nb/Y med": round(r.Nb_Y.median(), 3),
                     "(La/Yb)n med": round(r.LaYb_n.median(), 2),
                     "Th/Sc IQR": f"{r.Th_Sc.quantile(.25):.2f}-{r.Th_Sc.quantile(.75):.2f}",
                     "Nb/Y IQR": f"{r.Nb_Y.quantile(.25):.2f}-{r.Nb_Y.quantile(.75):.2f}",
                     "(La/Yb)n IQR": f"{r.LaYb_n.quantile(.25):.1f}-{r.LaYb_n.quantile(.75):.1f}"})
    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(RES, "provenance_summary.csv"), index=False)
    print(f"threshold: catchment >= {THRESH}% in domain")
    print(summary.to_string(index=False))
    allr = []
    for nm, df in sets.items():
        r = ratios(df); r["set"] = nm; allr.append(r)
    pd.concat(allr).to_csv(os.path.join(RES, "provenance_ratio_values.csv"), index=False)
    print("\nwrote results/provenance_summary.csv and provenance_ratio_values.csv")

if __name__ == "__main__":
    main()
