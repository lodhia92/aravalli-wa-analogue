"""
Does the fingerprint transfer depend on the matching ratios sharing mineral hosts with the
fingerprints?

This tests the transfer using provenance variables that do not reside in the targeted
rare-earth-bearing minerals. The earlier sensitivity test (dropping (La/Yb)n and Nb/Y from the six
ratios) is CONFOUNDED: it removes the host-shared ratios AND degrades the matcher at the same time,
so a fall in transfer cannot be attributed to either cause. This script separates them by building
matchers entirely from elements that sit in NO targeted mineral, and by reporting MATCH QUALITY
alongside the transfer statistic so the two effects are visible separately.

Ratio sets (all matched, standardised and paired by the analogue_pairing.py method, unchanged):
  SIX    Th/Sc, La/Sc, Th/Co, Eu/Eu*, (La/Yb)n, Nb/Y      original, reported as the control
  SET_A  Cr/Ni, V/Cr, Sc/Co, Ga/Al                        strict: low-solubility elements only
  SET_C  Cr/Ni, V/Cr, Sc/Co, Ga/Al, Ba/Sr                 adds feldspar-hosted discrimination

None of Cr, Ni, V, Sc, Co, Ga, Al, Ba or Sr occurs in monazite, xenotime, zircon, rutile,
ilmenite or apatite, so SET_A and SET_C are independent of the fingerprints at the mineral-host
level, not only at the element level.

Caveat carried into the output: the NGSA chromium and nickel values used here are the
inductively coupled plasma mass spectrometry determinations; a partial digestion under-reports
refractory chromite, so Cr/Ni and V/Cr carry a method caveat that Th/Sc and La/Sc do not.

Outputs
  results/matcher_independence.csv    per ratio set: match quality and per-mineral transfer
  results/matcher_independence_pairs.csv      the pair set produced by each ratio set
Run
  python scripts/matcher_independence.py

Author: Bhavik Harish Lodhia, Curtin University
"""
import os, csv

import numpy as np, pandas as pd

import paths

from aravalli_wa.stats import zscore_elem
from aravalli_wa.constants import AL_MASS_FRACTION_OF_AL2O3, P_MASS_FRACTION_OF_P2O5, TI_MASS_FRACTION_OF_TIO2, WT_PCT_TO_MG_KG

def ratios(df):
    o = pd.DataFrame(index=df.index)
    o["Th/Sc"] = df.Th / df.Sc
    o["La/Sc"] = df.La / df.Sc
    o["Th/Co"] = df.Th / df.Co
    o["EuEu"] = (df.Eu / CI["Eu"]) / np.sqrt((df.Sm / CI["Sm"]) * (df.Gd / CI["Gd"]))
    o["La/Yb_n"] = (df.La / df.Yb) / (CI["La"] / CI["Yb"])
    o["Nb/Y"] = df.Nb / df.Y
    o["Cr/Ni"] = df.Cr / df.Ni
    o["V/Cr"] = df.V / df.Cr
    o["Sc/Co"] = df.Sc / df.Co
    o["Ga/Al"] = df.Ga / df.Al
    o["Ba/Sr"] = df.Ba / df.Sr
    return o

def india(name):
    d = pd.read_csv(os.path.join(DR, f"{name}_contained.csv"))
    d = d[d.pct_in_domain >= IN_THR]
    d["k"] = d.lat.round(4).astype(str) + "_" + d.lon.round(4).astype(str)
    m = ar[ar.k.isin(set(d.k))].copy()
    m["sid"] = name + "_" + m.k
    return m

def findcol(el):
    for meth in ("ICP-MS", "XRF"):
        for i, c in enumerate(H):
            if c.strip().startswith(f"{el} {meth}"):
                return i

def aus(path, label):
    d = pd.read_csv(path)
    d = d[d.pct_in_domain >= AUS_THR]
    ids = set(pd.to_numeric(d["id"], errors="coerce").dropna())
    m = ng[ng.index.isin(ids)].copy()
    m["sid"] = [f"{label}_{i}" for i in m.index]
    return m.reset_index()

def logratios(df):
    r = ratios(df).replace([np.inf, -np.inf], np.nan)
    return np.log(r.where(r > 0)).replace([np.inf, -np.inf], np.nan)

def build(MATCH):
    out = []
    for idf, adf, dom in [(sand, wapp, "Palaeoproterozoic"), (mang, yiln, "Archaean")]:
        zi = ((logratios(idf) - mi) / si)[MATCH].dropna()
        za = ((logratios(adf) - ma) / sa)[MATCH].dropna()
        Iv, Av = zi.values, za.values
        if not len(Iv) or not len(Av):
            continue
        for a in range(len(Av)):
            di = np.sqrt(((Iv - Av[a]) ** 2).sum(1)); bi = int(di.argmin())
            da = np.sqrt(((Av - Iv[bi]) ** 2).sum(1)); mnn = (int(da.argmin()) == a)
            out.append(dict(domain=dom, india_sid=idf.loc[zi.index[bi], "sid"],
                            aus_sid=adf.loc[za.index[a], "sid"],
                            dist=float(di[bi]), dist_scaled=float(di[bi]) / np.sqrt(len(MATCH)),
                            mnn=mnn))
    P = pd.DataFrame(out)
    return pd.concat([P[P.domain == d].sort_values(["mnn", "dist"], ascending=[False, True]).head(10)
                      for d in ("Palaeoproterozoic", "Archaean")]), P

def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b); a, b = a[ok], b[ok]
    if len(a) < 4:
        return np.nan, np.nan, len(a)
    ra, rb = pd.Series(a).rank().values, pd.Series(b).rank().values
    rho = float(np.corrcoef(ra, rb)[0, 1])
    n = len(rb)
    B = rb[rng.random((NPERM, n)).argsort(axis=1)]      # NPERM relabellings at once
    Bc = B - B.mean(1, keepdims=True)
    ac = ra - ra.mean()
    den = np.sqrt((Bc ** 2).sum(1) * (ac ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        null = np.abs(np.where(den > 0, (Bc @ ac) / den, np.nan))
    return rho, (np.sum(null >= abs(rho)) + 1) / (NPERM + 1), n

def bh(ps):
    ps = np.array(ps, float); o = np.argsort(ps); m = len(ps); q = np.empty(m); prev = 1.0
    for rank, idx in enumerate(o[::-1]):
        i = m - rank; prev = min(prev, ps[idx] * m / i); q[idx] = prev
    return q

def main():
    global AUS0, AUS_THR, CI, DR, H, HERE, IN_THR, KEY0, MATCH, MIN, NG, NG_ELEMS, NPERM, NUM, P, PATH, PROJ, RES, ROB0, SETS, _, acn, allpairs, ar, c, cidx, da, di, e, el, f, grp, ic, k, key, l, label, ma, mang, mi, missing, mu_au, mu_in, n, name, ng, ngcm_all, ngsa_all, ordered, out, p, q, qs, r, rho, rng, rob, rows, sa, sand, sd_au, sd_in, si, v, wapp, xs, xs_ys, yiln, ys
    csv.field_size_limit(10 ** 7)
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ = os.path.dirname(os.path.dirname(HERE))
    RES = os.path.join(HERE, "results"); DR = os.path.join(RES, "drainage")
    IN_THR, AUS_THR = 50.0, 25.0
    CI = dict(La=.237, Yb=.170, Sm=.148, Eu=.0580, Gd=.199)
    PATH = ["Zr", "Hf", "Ti", "Ce", "Nd", "Pr", "Dy", "P"]
    MIN = {"zircon (Zr,Hf)": ["Zr", "Hf"], "monazite (Ce,Nd,Pr)": ["Ce", "Nd", "Pr"],
           "Ti-oxide (Ti)": ["Ti"], "xenotime (Dy)": ["Dy"], "apatite (P)": ["P"]}
    SETS = {
        "SIX (original)":  ["Th/Sc", "La/Sc", "Th/Co", "EuEu", "La/Yb_n", "Nb/Y"],
        "SET_A (strict)":  ["Cr/Ni", "V/Cr", "Sc/Co", "Ga/Al"],
        "SET_C (5-ratio)": ["Cr/Ni", "V/Cr", "Sc/Co", "Ga/Al", "Ba/Sr"],
    }
    NPERM = 100000
    rng = np.random.default_rng(20260827)
    ar = pd.read_csv(paths.NGCM_TABLE)
    NUM = ["Th", "Sc", "Co", "La", "Eu", "Sm", "Gd", "Yb", "Nb", "Y", "TiO2", "P2O5",
           "Zr", "Hf", "Ce", "Nd", "Pr", "Dy", "Cr", "Ni", "V", "Ga", "Ba", "Sr", "Al2O3", "LAT", "LON"]
    for c in NUM:
        ar[c] = pd.to_numeric(ar[c], errors="coerce")
    ar["Ti"] = ar["TiO2"] * WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2
    ar["P"] = ar["P2O5"] * WT_PCT_TO_MG_KG * P_MASS_FRACTION_OF_P2O5
    ar["Al"] = ar["Al2O3"] * WT_PCT_TO_MG_KG * AL_MASS_FRACTION_OF_AL2O3          # oxide wt% -> element mg/kg
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)
    sand, mang = india("sandmata"), india("mangalwar")
    NG = paths.NGSA
    with open(NG, encoding="latin-1") as f:
        H = list(csv.reader(f))[11]
    NG_ELEMS = ["Th", "Sc", "Nb", "Y", "La", "Yb", "Co", "Eu", "Sm", "Gd",
                "Cr", "Ni", "V", "Ga", "Ba", "Sr", "Al"] + PATH
    cidx = {el: findcol(el) for el in NG_ELEMS}
    missing = [k for k, v in cidx.items() if v is None]
    assert not missing, "NGSA columns not found: %s" % missing
    ordered = sorted(cidx.items(), key=lambda kv: kv[1])
    ng = pd.read_csv(NG, header=None, skiprows=12, usecols=[0] + [v for _, v in ordered],
                     encoding="latin-1", low_memory=False)
    ng.columns = ["SITEID"] + [k for k, _ in ordered]
    for c in ng.columns:
        ng[c] = pd.to_numeric(ng[c], errors="coerce")
    ng = ng.groupby("SITEID").median(numeric_only=True)
    wapp, yiln = aus(os.path.join(DR, "wa_palaeoprot_contained.csv"), "WA_PP"), \
                 aus(os.path.join(DR, "yilgarn_contained.csv"), "Y+N")
    ngcm_all, ngsa_all = pd.concat([sand, mang]), pd.concat([wapp, yiln])
    print("pools: NGCM %d, NGSA %d" % (len(ngcm_all), len(ngsa_all)))
    mi, si = logratios(ngcm_all).mean(), logratios(ngcm_all).std(ddof=0)
    ma, sa = logratios(ngsa_all).mean(), logratios(ngsa_all).std(ddof=0)
    mu_in, sd_in = zscore_elem(ngcm_all, PATH)
    mu_au, sd_au = zscore_elem(ngsa_all, PATH)
    ic, acn = ngcm_all.set_index("sid"), ngsa_all.set_index("sid")
    ROB0, _ = build(SETS["SIX (original)"])
    KEY0 = set(zip(ROB0.india_sid, ROB0.aus_sid))
    AUS0 = set(ROB0.aus_sid)
    rows = []
    allpairs = []
    for name, MATCH in SETS.items():
        rob, P = build(MATCH)
        rob = rob.copy(); rob["ratio_set"] = name
        allpairs.append(rob)
        key = set(zip(rob.india_sid, rob.aus_sid))
        xs_ys = {}
        for label, grp in MIN.items():
            xs, ys = [], []
            for _, r in rob.iterrows():
                if r.india_sid not in ic.index or r.aus_sid not in acn.index:
                    continue
                di, da = ic.loc[r.india_sid], acn.loc[r.aus_sid]
                xs.append(np.nanmean([(np.log(di[e]) - mu_in[e]) / sd_in[e] for e in grp]))
                ys.append(np.nanmean([(np.log(da[e]) - mu_au[e]) / sd_au[e] for e in grp]))
            xs_ys[label] = spearman(xs, ys)
        qs = bh([xs_ys[l][1] for l in MIN])
        for (label, (rho, p, n)), q in zip(xs_ys.items(), qs):
            rows.append(dict(ratio_set=name, n_ratios=len(MATCH), mineral=label,
                             n_pairs=n, rho=round(rho, 3), p=round(p, 4), q=round(q, 4),
                             transfers="Yes" if q < 0.05 else "No",
                             n_mnn=int(P.mnn.sum()), n_candidates=len(P),
                             median_dist_scaled=round(rob.dist_scaled.median(), 3),
                             pairs_shared_with_original=len(key & KEY0),
                             aus_sites_shared_with_original=len(set(rob.aus_sid) & AUS0)))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "matcher_independence.csv"), index=False)
    pd.concat(allpairs).to_csv(os.path.join(RES, "matcher_independence_pairs.csv"), index=False)
    pd.set_option("display.width", 200)
    print("\n==== MATCH QUALITY ====")
    print(out.drop_duplicates("ratio_set")[["ratio_set", "n_ratios", "n_mnn", "n_candidates",
          "median_dist_scaled", "pairs_shared_with_original", "aus_sites_shared_with_original"]]
          .to_string(index=False))
    print("\n==== PER-MINERAL TRANSFER ====")
    print(out[["ratio_set", "mineral", "n_pairs", "rho", "p", "q", "transfers"]].to_string(index=False))
    print("\nwrote results/matcher_independence.csv and matcher_independence_pairs.csv")

if __name__ == "__main__":
    main()
