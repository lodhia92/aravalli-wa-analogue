"""
Independent mineralogical validation of the transferred fingerprints against measured
heavy-mineral grain counts.

Question: do the per-mineral fingerprint scores used in this paper, computed from NGSA
stream-sediment chemistry, track the heavy minerals actually counted in the SAME samples by
the Heavy Mineral Map of Australia?

HMMA mineralogy is measured on the NGSA samples themselves, so the join is by Site_ID and is
exact. Scores are the within-survey standardised log-abundances of Equation 5, using the same
NGSA standardisation pool as the matching (all drainage-selected Australian catchments).

Tested at three scales: the 20 analogue pair sites, the 85 drainage-selected Australian sites,
and all NGSA sites carrying HMMA mineralogy.

Output: results/mineralogical_validation.csv, results/mineralogical_site_scores.csv
Run:    python scripts/mineralogical_validation.py

Author: Bhavik Harish Lodhia, Curtin University
"""

import csv
import os

import numpy as np
import pandas as pd

import paths
from aravalli_wa.constants import AUS_THR_PCT
from aravalli_wa.stats import perm_p

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
DR = os.path.join(RES, "drainage")
NG = paths.NGSA
HM = paths.HMMA
CACHE = os.path.join(RES, "hmma_grain_counts.csv")
NPERM = 100000

# fingerprint -> (elements used in the paper, HMMA mineral column stem)
FP = {
    "Monazite (LREE)": (["Ce", "Nd", "Pr"], "Monazite"),
    "Xenotime (HREE)": (["Dy"], "Xenotime-Y"),
    "Zircon": (["Zr", "Hf"], "Zircon"),
    "Ti oxides": (["Ti"], "Ilmenite"),
    "Apatite": (["P"], "Apatite"),
}
EXTRA_HM = ["Allanite", "Rutile", "Florencite", "Huttonite"]


def spearman(x, y):
    xr = pd.Series(x).rank().values
    yr = pd.Series(y).rank().values
    xr = xr - xr.mean()
    yr = yr - yr.mean()
    d = np.sqrt((xr**2).sum() * (yr**2).sum())
    return float((xr * yr).sum() / d) if d > 0 else np.nan


def load_ngsa():
    with open(NG, encoding="latin-1") as f:
        H = list(csv.reader(f))[11]
    ELEM = dict(
        Ce="ICP-MS",
        Nd="ICP-MS",
        Pr="ICP-MS",
        Dy="ICP-MS",
        Zr="ICP-MS",
        Hf="ICP-MS",
        Ti="XRF",
        P="XRF",
    )
    cidx = {}
    for el, meth in ELEM.items():
        for i, c in enumerate(H):
            if c.strip().startswith("%s %s" % (el, meth)):
                cidx[el] = i
                break
    ordered = sorted(cidx.items(), key=lambda kv: kv[1])
    ng = pd.read_csv(
        NG,
        header=None,
        skiprows=12,
        usecols=[0] + [v for _, v in ordered],
        encoding="latin-1",
        low_memory=False,
    )
    ng.columns = ["SITEID"] + [k for k, _ in ordered]
    for c in ng.columns:
        ng[c] = pd.to_numeric(ng[c], errors="coerce")
    return ng.groupby("SITEID").median(numeric_only=True).reset_index()


def load_hmma():
    if os.path.exists(CACHE):
        return pd.read_csv(CACHE)
    want = ["Site_ID"]
    stems = [v[1] for v in FP.values()] + EXTRA_HM
    for s in stems:
        want += ["%s (obs)" % s, "%s (pmo)" % s]
    a = pd.read_excel(HM, sheet_name="HMMA Dataset v.1.0", skiprows=11, header=0)
    keep = [c for c in want if c in a.columns]
    a = a[keep].copy()
    a = a.rename(columns={"Site_ID": "SITEID"})
    a.to_csv(CACHE, index=False)
    return a


def main():
    ng = load_ngsa()
    hm = load_hmma()
    print("NGSA sites %d; HMMA sites %d" % (len(ng), len(hm)))

    # drainage-selected Australian pool = the standardisation pool used for matching
    ids = set()
    for f in ["wa_palaeoprot_contained.csv", "yilgarn_contained.csv"]:
        d = pd.read_csv(os.path.join(DR, f))
        ids |= set(pd.to_numeric(d[d.pct_in_domain >= AUS_THR_PCT]["id"], errors="coerce").dropna())
    pool = ng[ng.SITEID.isin(ids)].copy()
    print("drainage-selected Australian pool: %d sites" % len(pool))

    # within-survey standardised log-abundance over the pool (Equation 5)
    z = pd.DataFrame({"SITEID": pool.SITEID.values})
    for el in ["Ce", "Nd", "Pr", "Dy", "Zr", "Hf", "Ti", "P"]:
        v = pd.to_numeric(pool[el], errors="coerce")
        x = np.log(v.where(v > 0))
        z[el] = ((x - x.mean()) / x.std(ddof=0)).values
    for name, (els, _) in FP.items():
        z[name] = z[els].mean(axis=1)

    pairs = pd.read_csv(os.path.join(RES, "analogue_pairs.csv"))
    pairs = pairs[pairs.mnn.astype(str).str.lower() == "true"]
    pair_ids = set(float(str(s).split("_")[-1]) for s in pairs.aus_sid)

    j = z.merge(hm, on="SITEID", how="inner")
    j["is_pair"] = j.SITEID.isin(pair_ids)
    j.to_csv(os.path.join(RES, "mineralogical_site_scores.csv"), index=False)
    print(
        "joined to HMMA: %d of %d pool sites; pair sites present: %d of 20"
        % (len(j), len(pool), int(j.is_pair.sum()))
    )

    rows = []
    scales = [("analogue pairs", j[j.is_pair]), ("drainage-selected pool", j)]
    for label, d in scales:
        for name, (els, stem) in FP.items():
            pcol, ocol = "%s (pmo)" % stem, "%s (obs)" % stem
            if pcol not in d.columns:
                continue
            m = d[[name, pcol, ocol]].apply(pd.to_numeric, errors="coerce").dropna()
            if len(m) < 6:
                continue
            rho = spearman(m[name].values, m[pcol].values)
            p = perm_p(m[name].values, m[pcol].values, rho)
            rows.append(
                dict(
                    scale=label,
                    fingerprint=name,
                    elements="+".join(els),
                    hmma_mineral=stem,
                    n=len(m),
                    pct_sites_present=round(100 * (m[ocol] > 0).mean(), 1),
                    rho=round(rho, 3),
                    perm_p=round(p, 4),
                )
            )
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "mineralogical_validation.csv"), index=False)
    print()
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
