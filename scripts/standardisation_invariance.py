"""
How much of the transfer result depends on the within-survey standardisation?

Standardising separately within each survey means the analysis compares relative rank rather than
absolute abundance. This script establishes the consequences of that choice by measurement.

The transfer statistic is a Spearman rank correlation, which is invariant under any strictly
increasing transformation applied within a survey. Three of the five fingerprints are single
elements (Dy, Ti, P), so for those the standardisation cannot change the reported correlation at
all. Only the two multi-element fingerprints, monazite (Ce, Nd, Pr) and zircon (Zr, Hf), can move,
and only through the relative weight the transformation gives their elements when they are averaged.
This script demonstrates that rather than asserting it, by recomputing the transfer under five
treatments of the same data on the same twenty pairs:

  1. within-survey z-score of log abundance                (the published method)
  2. raw log abundance, no standardisation at all          (identical to normalising both surveys
                                                            to any common external reference, which
                                                            only shifts each element by a constant)
  3. pooled standardisation over both surveys combined     (survey offsets preserved, not removed)
  4. robust within-survey standardisation, median and MAD
  5. within-survey percentile rank

Pair construction, data loading and the permutation procedure are those of fingerprint_transfer.py and
matching_sensitivity.py: 100 000 permutations, seed 20260827. Treatment 1 must reproduce monazite
rho = 0.696 and xenotime rho = 0.506.

Also reports the median inter-survey offset factor for each pathfinder element over the twenty
pairs, which is the absolute information the standardisation removes.

Outputs
  results/standardisation_invariance.csv   one row per treatment per fingerprint
  results/survey_element_offsets.csv           median Indian-to-Australian abundance ratio per element
Run
  python scripts/standardisation_invariance.py

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import csv
import os

import numpy as np
import pandas as pd
import paths

csv.field_size_limit(10 ** 7)
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
DR = os.path.join(RES, "drainage")
IN_THR, AUS_THR = 50.0, 25.0
CI = dict(La=.237, Yb=.170, Sm=.148, Eu=.0580, Gd=.199)
MATCH = ["Th/Sc", "La/Sc", "Th/Co", "EuEu", "La/Yb_n", "Nb/Y"]
PATH = ["Zr", "Hf", "Ti", "Ce", "Nd", "Pr", "Dy", "P"]
MIN = {"monazite (Ce,Nd,Pr)": ["Ce", "Nd", "Pr"], "xenotime (Dy)": ["Dy"],
       "zircon (Zr,Hf)": ["Zr", "Hf"], "Ti-oxide (Ti)": ["Ti"], "apatite (P)": ["P"]}
LAB = list(MIN)
NPERM = 100000
rng = np.random.default_rng(20260827)
_BANK = {}


def avg_rank(A):
    A = np.atleast_2d(A)
    m, n = A.shape
    order = np.argsort(A, axis=1, kind="mergesort")
    ranks = np.empty((m, n), float)
    ranks[np.arange(m)[:, None], order] = np.arange(1, n + 1)[None, :]
    S = np.take_along_axis(A, order, axis=1)
    for i in range(m):
        j = 0
        while j < n:
            k = j
            while k + 1 < n and S[i, k + 1] == S[i, j]:
                k += 1
            if k > j:
                ranks[i, order[i, j:k + 1]] = (j + k + 2) / 2.0
            j = k + 1
    return ranks


def corr_rows(Z, z):
    Zc = Z - Z.mean(1, keepdims=True)
    zc = z - z.mean()
    den = np.sqrt((Zc ** 2).sum(1) * (zc ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (Zc @ zc) / den, np.nan)


def spearman_perm(a, b):
    ra, rb = avg_rank(a)[0], avg_rank(b)[0]
    rho = float(corr_rows(ra[None, :], rb)[0])
    n = len(rb)
    if n not in _BANK:
        _BANK[n] = rng.random((NPERM, n)).argsort(axis=1)
    null = np.abs(corr_rows(rb[_BANK[n]], ra))
    return rho, (np.sum(null >= abs(rho)) + 1) / (NPERM + 1)


def bh(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    m = len(ps)
    q = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(o[::-1]):
        prev = min(prev, ps[idx] * m / (m - rank))
        q[idx] = prev
    return q


def ratios(df):
    o = pd.DataFrame(index=df.index)
    o["Th/Sc"] = df.Th / df.Sc
    o["La/Sc"] = df.La / df.Sc
    o["Th/Co"] = df.Th / df.Co
    o["EuEu"] = (df.Eu / CI["Eu"]) / np.sqrt((df.Sm / CI["Sm"]) * (df.Gd / CI["Gd"]))
    o["La/Yb_n"] = (df.La / df.Yb) / (CI["La"] / CI["Yb"])
    o["Nb/Y"] = df.Nb / df.Y
    return o


def logr(df):
    r = ratios(df).replace([np.inf, -np.inf], np.nan)
    return np.log(r.where(r > 0)).replace([np.inf, -np.inf], np.nan)


ar = pd.read_csv(paths.NGCM_TABLE)
for c in ["Th", "Sc", "Co", "La", "Eu", "Sm", "Gd", "Yb", "Nb", "Y", "TiO2", "P2O5",
          "Zr", "Hf", "Ce", "Nd", "Pr", "Dy", "LAT", "LON"]:
    ar[c] = pd.to_numeric(ar[c], errors="coerce")
ar["Ti"] = ar["TiO2"] * 1e4 * 0.5995
ar["P"] = ar["P2O5"] * 1e4 * 0.4364
ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)


def india(name):
    d = pd.read_csv(os.path.join(DR, f"{name}_contained.csv"))
    d = d[d.pct_in_domain >= IN_THR]
    d["k"] = d.lat.round(4).astype(str) + "_" + d.lon.round(4).astype(str)
    m = ar[ar.k.isin(set(d.k))].copy()
    m["sid"] = name + "_" + m.k
    return m


sand, mang = india("sandmata"), india("mangalwar")
NG = paths.NGSA
with open(NG, encoding="latin-1") as f:
    H = list(csv.reader(f))[11]


def findcol(el):
    for meth in ("ICP-MS", "XRF"):
        for i, c in enumerate(H):
            if c.strip().startswith(f"{el} {meth}"):
                return i


cidx = {el: findcol(el) for el in ["Th", "Sc", "Nb", "Y", "La", "Yb", "Co", "Eu", "Sm",
                                   "Gd"] + PATH}
ordered = sorted(((k, v) for k, v in cidx.items() if v is not None), key=lambda kv: kv[1])
ng = pd.read_csv(NG, header=None, skiprows=12, usecols=[0] + [v for _, v in ordered],
                 encoding="latin-1", low_memory=False)
ng.columns = ["SITEID"] + [k for k, _ in ordered]
for c in ng.columns:
    ng[c] = pd.to_numeric(ng[c], errors="coerce")
ng = ng.groupby("SITEID").median(numeric_only=True)


def aus(path, label):
    d = pd.read_csv(path)
    d = d[d.pct_in_domain >= AUS_THR]
    ids = set(pd.to_numeric(d["id"], errors="coerce").dropna())
    m = ng[ng.index.isin(ids)].copy()
    m["sid"] = [f"{label}_{i}" for i in m.index]
    return m.reset_index()


wapp, yiln = aus(os.path.join(DR, "wa_palaeoprot_contained.csv"), "WA_PP"), \
    aus(os.path.join(DR, "yilgarn_contained.csv"), "Y+N")
ngcm_all, ngsa_all = pd.concat([sand, mang]), pd.concat([wapp, yiln])
mi, si = logr(ngcm_all).mean(), logr(ngcm_all).std(ddof=0)
ma, sa = logr(ngsa_all).mean(), logr(ngsa_all).std(ddof=0)

# ---------------------------------------------------------------- the twenty published pairs
rows = []
for idf, adf, dom in [(sand, wapp, "Palaeoproterozoic"), (mang, yiln, "Archaean")]:
    zi = ((logr(idf) - mi) / si)[MATCH].dropna()
    za = ((logr(adf) - ma) / sa)[MATCH].dropna()
    Iv, Av = zi.values, za.values
    for a in range(len(Av)):
        di = np.sqrt(((Iv - Av[a]) ** 2).sum(1))
        bi = int(di.argmin())
        da = np.sqrt(((Av - Iv[bi]) ** 2).sum(1))
        rows.append(dict(domain=dom, india_sid=idf.loc[zi.index[bi], "sid"],
                         aus_sid=adf.loc[za.index[a], "sid"], dist=float(di[bi]),
                         mnn=(int(da.argmin()) == a)))
P = pd.DataFrame(rows)
pairs = pd.concat([P[P.domain == d].sort_values(["mnn", "dist"], ascending=[False, True]).head(10)
                   for d in ("Palaeoproterozoic", "Archaean")]).reset_index(drop=True)
print("pairs built: %d" % len(pairs), flush=True)

ic, acn = ngcm_all.set_index("sid"), ngsa_all.set_index("sid")
LI = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0))
LA = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0))
li = np.log(ic.loc[pairs.india_sid, PATH].where(ic.loc[pairs.india_sid, PATH] > 0))
la = np.log(acn.loc[pairs.aus_sid, PATH].where(acn.loc[pairs.aus_sid, PATH] > 0))


def mad(x):
    m = np.nanmedian(x, axis=0)
    return np.nanmedian(np.abs(x - m), axis=0) * 1.4826


def transform(name):
    """Return the per-element transformed values of the paired sites, Indian then Australian."""
    if name == "within-survey z-score of log abundance (published)":
        return ((li - LI.mean()) / LI.std(ddof=0)), ((la - LA.mean()) / LA.std(ddof=0))
    if name == "raw log abundance, no standardisation":
        return li.copy(), la.copy()
    if name == "pooled standardisation over both surveys":
        both = pd.concat([LI, LA])
        return ((li - both.mean()) / both.std(ddof=0)), ((la - both.mean()) / both.std(ddof=0))
    if name == "robust within-survey standardisation (median, MAD)":
        return (((li - np.nanmedian(LI, axis=0)) / mad(LI.values)),
                ((la - np.nanmedian(LA, axis=0)) / mad(LA.values)))
    if name == "within-survey percentile rank":
        out = []
        for pool, sub in ((LI, li), (LA, la)):
            r = pd.DataFrame(index=sub.index, columns=PATH, dtype=float)
            for e in PATH:
                v = pool[e].dropna().values
                r[e] = [np.nan if not np.isfinite(x) else (np.sum(v <= x) / (len(v) + 1.0))
                        for x in sub[e].values]
            out.append(r)
        return out[0], out[1]
    raise ValueError(name)


TREATMENTS = ["within-survey z-score of log abundance (published)",
              "raw log abundance, no standardisation",
              "pooled standardisation over both surveys",
              "robust within-survey standardisation (median, MAD)",
              "within-survey percentile rank"]

out = []
for t in TREATMENTS:
    TI, TA = transform(t)
    res = {}
    for lab in LAB:
        els = MIN[lab]
        with np.errstate(invalid="ignore"):
            a = np.nanmean(TI[els].values, axis=1)
            b = np.nanmean(TA[els].values, axis=1)
        ok = np.isfinite(a) & np.isfinite(b)
        res[lab] = spearman_perm(a[ok], b[ok]) + (int(ok.sum()),)
    qs = dict(zip(LAB, bh([res[l][1] for l in LAB])))
    for lab in LAB:
        r, p, n = res[lab]
        out.append(dict(treatment=t, mineral=lab, n_elements=len(MIN[lab]), n=n,
                        rho=round(r, 3), perm_p=round(p, 4), q=round(float(qs[lab]), 4),
                        transfers="Yes" if qs[lab] < 0.05 else "No"))
    print("  %-52s monazite %.3f  xenotime %.3f  zircon %.3f"
          % (t[:52], res[LAB[0]][0], res[LAB[1]][0], res[LAB[2]][0]), flush=True)

O = pd.DataFrame(out)
O.to_csv(os.path.join(RES, "standardisation_invariance.csv"), index=False)

# ---------------------------------------------------------------- absolute inter-survey offsets
off = []
for e in PATH:
    ratio = np.exp(li[e].values) / np.exp(la[e].values)
    ratio = ratio[np.isfinite(ratio)]
    off.append(dict(element=e, n_pairs=len(ratio),
                    median_india_over_australia=round(float(np.median(ratio)), 2),
                    p25=round(float(np.percentile(ratio, 25)), 2),
                    p75=round(float(np.percentile(ratio, 75)), 2)))
pd.DataFrame(off).to_csv(os.path.join(RES, "survey_element_offsets.csv"), index=False)

pd.set_option("display.width", 250)
print("\n==== transfer under each treatment ====")
print(O.to_string(index=False))
print("\n==== median Indian / Australian abundance over the twenty pairs ====")
print(pd.DataFrame(off).to_string(index=False))
print("\nwrote results/standardisation_invariance.csv and results/survey_element_offsets.csv")
