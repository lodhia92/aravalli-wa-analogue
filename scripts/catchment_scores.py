"""Per-sample pathfinder scores for every drainage-selected catchment.

For every drainage-selected sample (India: Sandmata + Mangalwar at >=50% catchment
containment; Australia: WA Palaeoproterozoic + Youanmi/Narryer at >=25%) compute one score
per pathfinder fingerprint:
    monazite_LREE  = mean standardised log (Ce, Nd, Pr)
    xenotime_HREE  = standardised log Dy
    zircon_ZrHf    = mean standardised log (Zr, Hf)
    tioxide_Ti     = standardised log Ti
Standardisation is within-survey (z-score of ln(element) over all drainage-selected
catchments of that survey), identical to pair_validation.py, so the scores and the pair
validation use one method. Transfer status (robust 20 pairs, FDR): monazite and xenotime
significant (q<0.05); zircon and Ti-oxide not significant.

Outputs: results/catchment_scores.csv with per-domain top-decile flags and
robust-pair numbers. Verification: reproduces the pair_validation.py per-mineral rho on the
robust-20 pairs from the exported scores.

Author: Bhavik Harish Lodhia, Curtin University
"""

import csv
import json
import os

import numpy as np
import pandas as pd

import paths
from aravalli_wa.constants import AUS_THR_PCT, TI_MASS_FRACTION_OF_TIO2, WT_PCT_TO_MG_KG
from aravalli_wa.survey import indian_pool


def india(name):
    m = indian_pool(ar, pd.read_csv(os.path.join(DR, f"{name}_contained.csv")), name)
    m = m.rename(columns={"LAT": "lat", "LON": "lon"})
    m["survey"] = "NGCM"
    m["domain"] = {"sandmata": "Sandmata", "mangalwar": "Mangalwar"}[name]
    return m[["sid", "survey", "domain", "lat", "lon"] + ELEMS]


def findcol(el):
    for meth in ("ICP-MS", "XRF"):
        for i, c in enumerate(H):
            if c.strip().startswith(f"{el} {meth}"):
                return i


def aus(name, label, domain):
    d = pd.read_csv(f"{DR}/{name}_contained.csv")
    d = d[d.pct_in_domain >= AUS_THR_PCT]
    d["id"] = pd.to_numeric(d["id"], errors="coerce")
    m = ng[ng.index.isin(set(d.id.dropna()))].copy().reset_index()
    m = m.merge(d[["id", "lat", "lon"]].drop_duplicates("id"), left_on="SITEID", right_on="id")
    m["sid"] = [f"{label}_{i}" for i in m.SITEID]
    m["survey"] = "NGSA"
    m["domain"] = domain
    return m[["sid", "survey", "domain", "lat", "lon"] + ELEMS]


def zscore(pool):
    x = np.log(pool[ELEMS].where(pool[ELEMS] > 0))
    return x.mean(), x.std(ddof=0)


def score(df, mu, sd):
    z = (np.log(df[ELEMS].where(df[ELEMS] > 0)) - mu) / sd
    out = df[["sid", "survey", "domain", "lat", "lon"]].copy()
    for name, grp in FP.items():
        out[name] = z[grp].mean(axis=1)
    return out


def main():
    global \
        DR, \
        ELEMS, \
        FP, \
        H, \
        HERE, \
        NG, \
        PROJ, \
        RES, \
        _, \
        allsc, \
        ar, \
        c, \
        cidx, \
        d, \
        dom, \
        el, \
        expected, \
        f, \
        g, \
        k, \
        mang, \
        match, \
        mu_au, \
        mu_in, \
        name, \
        ng, \
        ngcm_all, \
        ngsa_all, \
        ok, \
        ordered, \
        pno, \
        r, \
        rho, \
        robust, \
        rp, \
        s, \
        sand, \
        sd_au, \
        sd_in, \
        thr, \
        use, \
        v, \
        wapp, \
        x, \
        y, \
        yiln
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ = os.path.dirname(HERE)
    RES = f"{HERE}/results"
    DR = f"{RES}/drainage"
    ELEMS = ["Ce", "Nd", "Pr", "Dy", "Zr", "Hf", "Ti"]
    FP = {
        "monazite_LREE": ["Ce", "Nd", "Pr"],
        "xenotime_HREE": ["Dy"],
        "zircon_ZrHf": ["Zr", "Hf"],
        "tioxide_Ti": ["Ti"],
    }
    ar = pd.read_csv(paths.NGCM_TABLE)
    for c in ["TiO2", "Zr", "Hf", "Ce", "Nd", "Pr", "Dy", "LAT", "LON"]:
        ar[c] = pd.to_numeric(ar[c], errors="coerce")
    ar["Ti"] = ar["TiO2"] * WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)
    sand = india("sandmata")
    mang = india("mangalwar")
    NG = paths.NGSA
    with open(NG, encoding="latin-1") as f:
        H = list(csv.reader(f))[11]
    cidx = {el: findcol(el) for el in ELEMS}
    cidx = {k: v for k, v in cidx.items() if v is not None}
    ordered = sorted(cidx.items(), key=lambda kv: kv[1])
    use = [0] + [v for _, v in ordered]
    ng = pd.read_csv(
        NG, header=None, skiprows=12, usecols=use, encoding="latin-1", low_memory=False
    )
    ng.columns = ["SITEID"] + [k for k, _ in ordered]
    for c in ng.columns:
        ng[c] = pd.to_numeric(ng[c], errors="coerce")
    ng = ng.groupby("SITEID").median(numeric_only=True)
    wapp = aus("wa_palaeoprot", "WA_PP", "WA_Palaeoproterozoic")
    yiln = aus("yilgarn", "Y+N", "Youanmi+Narryer")
    ngcm_all = pd.concat([sand, mang], ignore_index=True)
    ngsa_all = pd.concat([wapp, yiln], ignore_index=True)
    mu_in, sd_in = zscore(ngcm_all)
    mu_au, sd_au = zscore(ngsa_all)
    allsc = pd.concat(
        [score(ngcm_all, mu_in, sd_in), score(ngsa_all, mu_au, sd_au)], ignore_index=True
    )
    for name in FP:
        allsc[f"top10_{name}"] = False
        for dom, g in allsc.groupby("domain"):
            thr = g[name].quantile(0.90)
            allsc.loc[g.index[g[name] >= thr], f"top10_{name}"] = True
    rp = pd.read_csv(f"{RES}/analogue_pairs.csv")
    pno = {}
    for f in json.load(open(f"{RES}/pair_catchments_india.geojson"))["features"]:
        pno[f["properties"]["sid"]] = f["properties"]["pair_no"]
    for f in json.load(open(f"{RES}/pair_catchments_australia.geojson"))["features"]:
        pno[f["properties"]["sid"]] = f["properties"]["pair_no"]
    allsc["pair_no"] = allsc.sid.map(pno).astype("Int64")
    allsc.to_csv(f"{RES}/catchment_scores.csv", index=False)
    print("wrote catchment_scores.csv:", len(allsc), "samples")
    print(allsc.groupby(["survey", "domain"]).size())
    robust = pd.concat(
        [
            rp[rp.domain == d].sort_values(["mnn", "dist"], ascending=[False, True]).head(10)
            for d in ["Palaeoproterozoic", "Archaean"]
        ]
    )
    s = allsc.set_index("sid")
    expected = {
        "monazite_LREE": 0.696,
        "xenotime_HREE": 0.506,
        "zircon_ZrHf": 0.313,
        "tioxide_Ti": -0.240,
    }
    print("\nverification vs pair_per_element_validation.csv (robust-20 Spearman rho):")
    ok = True
    for name in FP:
        x = [s.loc[r.india_sid, name] for _, r in robust.iterrows()]
        y = [s.loc[r.aus_sid, name] for _, r in robust.iterrows()]
        rho = pd.Series(x).rank().corr(pd.Series(y).rank())
        match = abs(rho - expected[name]) < 0.005
        ok &= match
        print(
            f"  {name:15} rho={rho:+.3f}  expected {expected[name]:+.3f}  {'OK' if match else 'MISMATCH'}"
        )
    print("VERIFICATION", "PASSED" if ok else "FAILED")


if __name__ == "__main__":
    main()
