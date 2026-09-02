"""
Analytical methods, detection limits, censored values and missing data, per element and per pool.

Counts only. No permutations, no matching, no scoring. Paths, containment thresholds and
column selection are copied from matching_sensitivity.py so that the pools counted here are
exactly the pools the published analysis uses.

Outputs
  results/data_quality.csv   one row per element per pool
  results/analytical_methods.csv        analytical method, unit and detection limit per element

Author: Bhavik Harish Lodhia, Curtin University
"""

import csv
import os

import pandas as pd

import paths
from aravalli_wa.constants import (
    AUS_THR_PCT,
    P_MASS_FRACTION_OF_P2O5,
    TI_MASS_FRACTION_OF_TIO2,
    WT_PCT_TO_MG_KG,
)
from aravalli_wa.survey import indian_pool


def india(name):
    contained = pd.read_csv(os.path.join(DR, f"{name}_contained.csv"))
    return indian_pool(ar, contained, name).drop(columns="sid")


def findcol(el):
    for meth in ("ICP-MS", "XRF"):
        for i, c in enumerate(H):
            if c.strip().startswith(f"{el} {meth}"):
                return i, c.strip()
    return None, None


def aus(path):
    d = pd.read_csv(path)
    d = d[d.pct_in_domain >= AUS_THR_PCT]
    ids = set(pd.to_numeric(d["id"], errors="coerce").dropna())
    return ngn[ngn.index.isin(ids)].copy(), cens_site[cens_site.index.isin(ids)].copy()


def main():
    global \
        DR, \
        ELS, \
        H, \
        HERE, \
        MATCH_EL, \
        NG, \
        PATH, \
        PROJ, \
        RES, \
        _, \
        a1, \
        a2, \
        ar, \
        c, \
        c1, \
        c2, \
        cens_flag, \
        cens_rows, \
        cens_site, \
        cidx, \
        cols, \
        el, \
        f, \
        i, \
        k, \
        m, \
        meth, \
        n, \
        name_of, \
        ngcm, \
        ngn, \
        ngsa, \
        ngsa_c, \
        out, \
        raw, \
        rows, \
        v
    csv.field_size_limit(10**7)
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ = os.path.dirname(os.path.dirname(HERE))
    RES = os.path.join(HERE, "results")
    DR = os.path.join(RES, "drainage")
    MATCH_EL = ["Th", "Sc", "Co", "La", "Eu", "Sm", "Gd", "Yb", "Nb", "Y"]
    PATH = ["Zr", "Hf", "Ti", "Ce", "Nd", "Pr", "Dy", "P"]
    ELS = MATCH_EL + PATH
    ar = pd.read_csv(paths.NGCM_TABLE, low_memory=False)
    for c in ELS + ["TiO2", "P2O5", "LAT", "LON"]:
        if c in ar.columns:
            ar[c] = pd.to_numeric(ar[c], errors="coerce")
    ar["Ti"] = ar["TiO2"] * WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2
    ar["P"] = ar["P2O5"] * WT_PCT_TO_MG_KG * P_MASS_FRACTION_OF_P2O5
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)
    ngcm = pd.concat([india("sandmata"), india("mangalwar")])
    NG = paths.NGSA
    with open(NG, encoding="latin-1") as f:
        H = list(csv.reader(f))[11]
    cidx = {el: findcol(el) for el in ELS}
    cols = sorted(i for i, _ in cidx.values() if i is not None)
    raw = pd.read_csv(
        NG,
        header=None,
        skiprows=12,
        usecols=[0] + cols,
        encoding="latin-1",
        low_memory=False,
        dtype=str,
    )
    raw.columns = ["SITEID"] + [H[i].strip() for i in cols]
    name_of = {el: cidx[el][1] for el in ELS}
    cens_rows = {
        el: int(raw[name_of[el]].astype(str).str.strip().str.startswith("<").sum()) for el in ELS
    }
    ngn = raw.copy()
    for c in ngn.columns:
        ngn[c] = pd.to_numeric(ngn[c], errors="coerce")
    ngn = ngn.groupby("SITEID").median(numeric_only=True)
    ngn = ngn.rename(columns={v: k for k, v in name_of.items()})
    cens_flag = raw.copy()
    for el in ELS:
        cens_flag[el + "_c"] = cens_flag[name_of[el]].astype(str).str.strip().str.startswith("<")
    cens_site = cens_flag.groupby("SITEID")[[el + "_c" for el in ELS]].any()
    a1, c1 = aus(os.path.join(DR, "wa_palaeoprot_contained.csv"))
    a2, c2 = aus(os.path.join(DR, "yilgarn_contained.csv"))
    ngsa, ngsa_c = pd.concat([a1, a2]), pd.concat([c1, c2])
    rows = []
    for el in ELS:
        n = len(ngcm)
        m = int(ngcm[el].isna().sum())
        rows.append(
            dict(
                pool="India NGCM drainage-selected",
                element=el,
                n=n,
                missing=m,
                pct_missing=round(100 * m / n, 1),
                censored="not distinguishable",
                pct_censored="",
            )
        )
    for el in ELS:
        n = len(ngsa)
        m = int(ngsa[el].isna().sum())
        c = int(ngsa_c[el + "_c"].sum())
        rows.append(
            dict(
                pool="Australia NGSA drainage-selected",
                element=el,
                n=n,
                missing=m,
                pct_missing=round(100 * m / n, 1),
                censored=c,
                pct_censored=round(100 * c / n, 1),
            )
        )
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "data_quality.csv"), index=False)
    meth = pd.DataFrame(
        [
            dict(element=el, ngsa_column=name_of[el], censored_analyses_all_ngsa=cens_rows[el])
            for el in ELS
        ]
    )
    meth.to_csv(os.path.join(RES, "analytical_methods.csv"), index=False)
    print("India pool n =", len(ngcm), " Australia pool n =", len(ngsa))
    print()
    print(out.to_string(index=False))
    print()
    print(meth.to_string(index=False))


if __name__ == "__main__":
    main()
