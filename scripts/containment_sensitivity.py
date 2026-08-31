"""How do the analogue pairs and the transferred fingerprints change with the catchment containment
thresholds, and does terrane crossing on the Australian side carry the transfer?

The Indian and Australian sample sets use different containment thresholds, 50 per cent and 25 per
cent, because the two sets of catchments differ by about an order of magnitude in area. Australian
catchments can therefore cross terrane boundaries. This script recomputes the pairs and the
transfer across a grid of thresholds to establish whether that asymmetry carries the result.

Two parts.

--part grid
    Rebuilds the whole analysis, standardisation and matching included, at twelve threshold
    combinations: the Indian threshold varied with Australia held at 25, the Australian threshold
    varied with India held at 50, and both raised together. Reports pool sizes, pair count, pairs
    shared with the published twenty, and the Spearman correlation with its Benjamini-Hochberg q
    for all five fingerprints. Flags any combination where the Australian pool in a domain has
    fallen to the number of pairs requested, because at that point the top-ten rule is no longer
    selecting anything.

--part confound
    Treats containment as the continuous variable it is rather than as a cut. pct_in_domain is the
    fraction of a sample's catchment lying inside its domain, so it measures terrane crossing
    directly. Tests whether it correlates with the fingerprint scores, then partials it out of the
    transfer over the published pairs, exactly as grain_size_morphometry.py does for catchment area.

The matching, the statistics and the correction are imported from matching_sensitivity.py so they
cannot drift; only the loader is reimplemented here, with the thresholds as arguments and without
the cache. The published combination (50, 25) is asserted to return the published correlations,
monazite 0.696 and xenotime 0.506.

Convention shared with the rest of the study: 100 000 permutations, seed 20260827, inherited from
matching_sensitivity.py.

Inputs : the NGCM and NGSA tables resolved by paths.py
         results/drainage/{sandmata,mangalwar,wa_palaeoprot,yilgarn}_contained.csv
Outputs: results/containment_grid.csv
         results/containment_vs_scores.csv
         results/containment_covariate.csv

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import argparse
import csv
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matching_sensitivity as core  # matching, statistics, BH correction
import paths
from aravalli_wa.constants import P_MASS_FRACTION_OF_P2O5, TI_MASS_FRACTION_OF_TIO2, WT_PCT_TO_MG_KG

PROJ = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
DR = os.path.join(RES, "drainage")
DOMS, LAB, MIN = core.DOMS, core.LAB, core.MIN
MATCH, PATH = core.MATCH, core.PATH
PUB_IN, PUB_AU = 50.0, 25.0
PUB_RHO = {"monazite (Ce,Nd,Pr)": 0.696, "xenotime (Dy)": 0.506}

_RAW = {}
_D = {}


# ---------------------------------------------------------------- loading

def load_raw():
    """Read the two geochemistry tables and the four drainage tables once."""
    if _RAW:
        return _RAW
    ar = pd.read_csv(paths.NGCM_TABLE)
    for c in ["Th", "Sc", "Co", "La", "Eu", "Sm", "Gd", "Yb", "Nb", "Y", "TiO2", "P2O5",
              "Zr", "Hf", "Ce", "Nd", "Pr", "Dy", "LAT", "LON"]:
        ar[c] = pd.to_numeric(ar[c], errors="coerce")
    ar["Ti"] = ar["TiO2"] * WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2
    ar["P"] = ar["P2O5"] * WT_PCT_TO_MG_KG * P_MASS_FRACTION_OF_P2O5
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)

    NG = paths.NGSA
    with open(NG, encoding="latin-1") as f:
        H = list(csv.reader(f))[11]

    def findcol(el):
        for meth in ("ICP-MS", "XRF"):
            for i, c in enumerate(H):
                if c.strip().startswith(f"{el} {meth}"):
                    return i

    cidx = {el: findcol(el)
            for el in ["Th", "Sc", "Nb", "Y", "La", "Yb", "Co", "Eu", "Sm", "Gd"] + PATH}
    ordered = sorted(((k, v) for k, v in cidx.items() if v is not None), key=lambda kv: kv[1])
    ng = pd.read_csv(NG, header=None, skiprows=12, usecols=[0] + [v for _, v in ordered],
                     encoding="latin-1", low_memory=False)
    ng.columns = ["SITEID"] + [k for k, _ in ordered]
    for c in ng.columns:
        ng[c] = pd.to_numeric(ng[c], errors="coerce")
    ng = ng.groupby("SITEID").median(numeric_only=True)

    dr = {}
    for name in ["sandmata", "mangalwar", "wa_palaeoprot", "yilgarn"]:
        d = pd.read_csv(os.path.join(DR, f"{name}_contained.csv"))
        if "id" in d.columns:
            d["id"] = pd.to_numeric(d["id"], errors="coerce")
            d = d.dropna(subset=["id"]).drop_duplicates("id")
        else:
            d["k"] = d.lat.round(4).astype(str) + "_" + d.lon.round(4).astype(str)
            d = d.drop_duplicates("k")
        dr[name] = d
    _RAW.update(ar=ar, ng=ng, dr=dr)
    return _RAW


def build_D(in_thr, aus_thr):
    """matching_sensitivity.prepare() with the thresholds as arguments and no cache.

    The standardisation pools are the threshold-filtered pools, as in the published analysis, so
    changing a threshold changes the standardisation too. That is deliberate: each combination is
    a full rebuild of the analysis, not a re-cut of fixed scores. pct_in_domain is carried through
    per site for the confound part.
    """
    key = (in_thr, aus_thr)
    if key in _D:
        return _D[key]
    R = load_raw()
    ar, ng, dr = R["ar"], R["ng"], R["dr"]

    def india(name):
        d = dr[name]
        d = d[d.pct_in_domain >= in_thr]
        pct = dict(zip(d.k, d.pct_in_domain))
        m = ar[ar.k.isin(set(d.k))].copy()
        m["sid"] = name + "_" + m.k
        m["pct"] = m.k.map(pct)
        return m

    def aus(name, label):
        d = dr[name]
        d = d[d.pct_in_domain >= aus_thr]
        pct = dict(zip(d.id, d.pct_in_domain))
        m = ng[ng.index.isin(set(d.id))].copy()
        m["sid"] = [f"{label}_{i}" for i in m.index]
        m["pct"] = [pct[i] for i in m.index]
        return m.reset_index()

    sand, mang = india("sandmata"), india("mangalwar")
    wapp, yiln = aus("wa_palaeoprot", "WA_PP"), aus("yilgarn", "Y+N")
    ngcm_all, ngsa_all = pd.concat([sand, mang]), pd.concat([wapp, yiln])
    mi, si = core.logr(ngcm_all).mean(), core.logr(ngcm_all).std(ddof=0)
    ma, sa = core.logr(ngsa_all).mean(), core.logr(ngsa_all).std(ddof=0)
    mu_in = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0)).mean()
    sd_in = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0)).std(ddof=0)
    mu_au = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0)).mean()
    sd_au = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0)).std(ddof=0)

    def score(df, mu, sd):
        S = {}
        for lab, els in MIN.items():
            Z = np.column_stack([(np.log(df[e].where(df[e] > 0)) - mu[e]) / sd[e] for e in els])
            with np.errstate(invalid="ignore"):
                S[lab] = np.nanmean(Z, axis=1)
        return S

    D = {}
    for dom, idf, adf in [(DOMS[0], sand, wapp), (DOMS[1], mang, yiln)]:
        zi = ((core.logr(idf) - mi) / si)[MATCH].dropna()
        za = ((core.logr(adf) - ma) / sa)[MATCH].dropna()
        isub, asub = idf.loc[zi.index], adf.loc[za.index]
        D[dom] = dict(Iv=zi.values, Av=za.values,
                      isid=isub["sid"].values, asid=asub["sid"].values,
                      ipct=isub["pct"].values.astype(float),
                      apct=asub["pct"].values.astype(float),
                      SI=score(isub, mu_in, sd_in), SA=score(asub, mu_au, sd_au))
    _D[key] = D
    return D


def transfer(D, sel):
    """Spearman rho, permutation p and BH q for every fingerprint under one pair selection."""
    res = {lab: core.spearman_perm(*core.vectors(D, sel, lab)[:2]) for lab in LAB}
    qs = dict(zip(LAB, core.bh([res[l][1] for l in LAB])))
    return res, qs


# ---------------------------------------------------------------- part 1: the grid

GRID = ([("Indian threshold varied, Australia held at 25", t, 25.0) for t in (40., 50., 60., 75.)]
        + [("Australian threshold varied, India held at 50", 50., t) for t in (25., 40., 50., 60.)]
        + [("both raised together", t, t) for t in (40., 50., 60.)])


def run_grid():
    pubD = build_D(PUB_IN, PUB_AU)
    pubsel = core.select(pubD, "published")
    pubk = core.pairkeys(pubD, pubsel)
    pres, _ = transfer(pubD, pubsel)
    for lab, want in PUB_RHO.items():
        got = round(pres[lab][0], 3)
        assert got == want, f"published anchor failed: {lab} {got} != {want}"
    print("published anchor reproduced: monazite %.3f, xenotime %.3f"
          % (pres[LAB[0]][0], pres[LAB[1]][0]), flush=True)

    rows = []
    for part, it, at in GRID:
        D = build_D(it, at)
        sel = core.select(D, "published")
        npairs = sum(len(sel[d][0]) for d in DOMS)
        shared = len(core.pairkeys(D, sel) & pubk)
        au_pool = {d: len(D[d]["Av"]) for d in DOMS}
        in_pool = {d: len(D[d]["Iv"]) for d in DOMS}
        degenerate = min(au_pool.values()) <= core.KPUB
        res, qs = transfer(D, sel)
        for lab in LAB:
            r, p = res[lab]
            rows.append(dict(
                part=part, india_threshold=it, australia_threshold=at, mineral=lab,
                india_pool_palaeoprot=in_pool[DOMS[0]], india_pool_archaean=in_pool[DOMS[1]],
                aus_pool_palaeoprot=au_pool[DOMS[0]], aus_pool_archaean=au_pool[DOMS[1]],
                n_pairs=npairs, shared_with_published=shared,
                australian_pool_exhausted="Yes" if degenerate else "No",
                rho=round(r, 3), perm_p=round(p, 4), q=round(float(qs[lab]), 4),
                transfers="Yes" if qs[lab] < 0.05 else "No"))
        print("  in %3.0f aus %3.0f  aus pool %3d/%3d  pairs %2d  shared %2d%s"
              "  monazite %.3f q %.4f  xenotime %.3f q %.4f"
              % (it, at, au_pool[DOMS[0]], au_pool[DOMS[1]], npairs, shared,
                 "  POOL EXHAUSTED" if degenerate else "",
                 res[LAB[0]][0], qs[LAB[0]], res[LAB[1]][0], qs[LAB[1]]), flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "containment_grid.csv"), index=False)
    print("\nwrote results/containment_grid.csv")


# ---------------------------------------------------------------- part 2: containment as covariate

def corr_pool(s, c):
    """Spearman of a score against containment over a whole eligible pool.

    Permutation at the study convention is used while it is affordable. The Indian pools run to
    thousands of sites, where a 100 000-relabelling bank does not fit in memory and is not needed:
    above 200 points the two-sided asymptotic probability on the rank correlation is accurate to
    well within the precision reported here. The Australian pools, which are the ones this test
    is asking about, are all small enough to keep the permutation test.

    Returns nan when containment does not vary in the pool, which is the case on the Indian side:
    almost every Indian catchment lies entirely inside its domain, so there is nothing to
    correlate against.
    """
    if len(s) < 3 or np.ptp(c) == 0:
        return np.nan, np.nan, "undefined"
    if len(s) <= 200:
        r, p = core.spearman_perm(s, c)
        return r, p, "permutation, 100 000"
    ra, rb = core.avg_rank(s)[0], core.avg_rank(c)[0]
    r = float(core.corr_rows(ra[None, :], rb)[0])
    if not np.isfinite(r) or abs(r) >= 1:
        return r, np.nan, "undefined"
    n = len(s)
    t = r * math.sqrt((n - 2) / (1 - r * r))
    p = math.erfc(abs(t) / math.sqrt(2))          # normal approximation, exact enough at n > 200
    return r, p, "asymptotic"


def run_confound():
    D = build_D(PUB_IN, PUB_AU)
    sel = core.select(D, "published")

    # (a) does containment correlate with the fingerprint scores, over the whole eligible pool?
    rows = []
    for side, vkey, pkey in (("australia", "SA", "apct"), ("india", "SI", "ipct")):
        for scope in DOMS + ["both domains pooled"]:
            doms = DOMS if scope == "both domains pooled" else [scope]
            ps, keep = [], {}
            for lab in LAB:
                s = np.concatenate([D[d][vkey][lab] for d in doms])
                c = np.concatenate([D[d][pkey] for d in doms])
                ok = np.isfinite(s) & np.isfinite(c)
                r, p, how = corr_pool(s[ok], c[ok])
                keep[lab] = (r, p, int(ok.sum()), how)
                ps.append(1.0 if not np.isfinite(p) else p)
            qs = core.bh(ps)
            for j, lab in enumerate(LAB):
                r, p, n, how = keep[lab]
                rows.append(dict(side=side, scope=scope, mineral=lab, n=n,
                                 pct_spread=round(float(np.nanpercentile(
                                     np.concatenate([D[d][pkey] for d in doms]), 75)
                                     - np.nanpercentile(
                                     np.concatenate([D[d][pkey] for d in doms]), 25)), 1),
                                 rho_score_vs_pct_in_domain=("" if not np.isfinite(r)
                                                             else round(r, 3)),
                                 p=("" if not np.isfinite(p) else round(p, 4)), test=how,
                                 q=("" if not np.isfinite(p) else round(float(qs[j]), 4)),
                                 significant=("undefined" if not np.isfinite(r)
                                              else "Yes" if qs[j] < 0.05 else "No")))
            r0 = keep[LAB[0]][0]
            print("  %-9s %-22s monazite %s  n=%d  (%s)"
                  % (side, scope[:22],
                     "undefined, containment does not vary" if not np.isfinite(r0)
                     else "%.3f (q %.4f)" % (r0, qs[0]),
                     keep[LAB[0]][2], keep[LAB[0]][3]), flush=True)
    P = pd.DataFrame(rows)
    P.to_csv(os.path.join(RES, "containment_vs_scores.csv"), index=False)
    print("wrote results/containment_vs_scores.csv\n")

    # (b) partial out containment from the transfer, over the published pairs
    pct_au, pct_in = [], []
    for dom in DOMS:
        ipos, apos = sel[dom]
        pct_in.append(D[dom]["ipct"][ipos])
        pct_au.append(D[dom]["apct"][apos])
    pct_au, pct_in = np.concatenate(pct_au), np.concatenate(pct_in)

    rows = []
    for cname, cov in (("Australian containment", pct_au),
                       ("Indian containment", pct_in),
                       ("both containments", np.column_stack([pct_au, pct_in]))):
        cov2 = cov if cov.ndim > 1 else cov[:, None]
        C = np.column_stack([core.avg_rank(cov2[:, j])[0] for j in range(cov2.shape[1])])
        ps, keep = [], {}
        for lab in LAB:
            a, b, _ = core.vectors(D, sel, lab)
            # vectors() drops non-finite pairs, so rebuild the covariate on the same mask
            aa, bb = [], []
            mask = []
            for gi, dom in enumerate(DOMS):
                ipos, apos = sel[dom]
                aa.append(D[dom]["SI"][lab][ipos])
                bb.append(D[dom]["SA"][lab][apos])
            aa, bb = np.concatenate(aa), np.concatenate(bb)
            mask = np.isfinite(aa) & np.isfinite(bb)
            ra, rb = core.avg_rank(aa[mask])[0], core.avg_rank(bb[mask])[0]
            X = np.column_stack([np.ones(int(mask.sum())), C[mask]])
            resa = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
            resb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
            r = float(core.corr_rows(resa[None, :], resb)[0])
            n = len(resb)
            null = np.abs(core.corr_rows(resb[core.bank(n)], resa))
            p = (np.sum(null >= abs(r)) + 1) / (core.NPERM + 1)
            keep[lab] = (r, p, n)
            ps.append(p)
        qs = core.bh(ps)
        for j, lab in enumerate(LAB):
            r, p, n = keep[lab]
            rows.append(dict(controlling_for=cname, mineral=lab, n=n, partial_rho=round(r, 3),
                             perm_p=round(p, 4), q=round(float(qs[j]), 4),
                             transfers="Yes" if qs[j] < 0.05 else "No"))
        print("  controlling for %-24s monazite %.3f (q %.4f)  xenotime %.3f"
              % (cname, keep[LAB[0]][0], qs[0], keep[LAB[1]][0]), flush=True)
    C2 = pd.DataFrame(rows)
    C2.to_csv(os.path.join(RES, "containment_covariate.csv"), index=False)
    pd.set_option("display.width", 260)
    print("\n" + C2.to_string(index=False))
    print("\nwrote results/containment_covariate.csv")

    # context: how much of each Australian catchment actually sits outside its domain
    print("\ncontainment of the selected Australian pair members, per cent in domain:")
    print("  median %.1f, min %.1f, max %.1f, below 50 per cent: %d of %d"
          % (np.median(pct_au), pct_au.min(), pct_au.max(), (pct_au < 50).sum(), len(pct_au)))
    print("containment of the selected Indian pair members:")
    print("  median %.1f, min %.1f, max %.1f" % (np.median(pct_in), pct_in.min(), pct_in.max()))



# ---------------------------------------------------------------- part 3: is it size, not containment?

POWER_DRAWS = 1000
POWER_SEED = 20260827


def pair_distance(D, sel):
    """Mean Euclidean distance between paired sites in the six-ratio matching space.

    This is the match quality. If raising a threshold empties the pool, the surviving matches get
    worse, and this number rises.
    """
    d = []
    for dom in DOMS:
        ipos, apos = sel[dom]
        d.append(np.linalg.norm(D[dom]["Iv"][ipos] - D[dom]["Av"][apos], axis=1))
    return float(np.concatenate(d).mean())


def subsample(D, sizes, rng):
    """Draw a random Australian sub-pool of the given size in each domain, India left whole."""
    out = {}
    for dom in DOMS:
        n = sizes[dom]
        idx = rng.choice(len(D[dom]["Av"]), size=n, replace=False)
        idx.sort()
        out[dom] = dict(D[dom])
        out[dom]["Av"] = D[dom]["Av"][idx]
        out[dom]["asid"] = D[dom]["asid"][idx]
        out[dom]["apct"] = D[dom]["apct"][idx]
        out[dom]["SA"] = {lab: D[dom]["SA"][lab][idx] for lab in LAB}
    return out


def run_power():
    """Does the transfer collapse because catchments cross terranes, or because the pool empties?

    Raising the Australian containment threshold does two things at once: it removes the
    catchments that cross terrane boundaries, and it shrinks the pool the matching can draw on.
    This separates them. Containment is held at the published 25 per cent, so every site is as
    terrane-crossing as it ever was, and the Australian pool is instead cut at random to the size
    that each raised threshold produces. If a random pool of the same size loses as much
    correlation as the thresholded pool did, the collapse is sample size.
    """
    rng = np.random.default_rng(POWER_SEED)
    Dfull = build_D(PUB_IN, PUB_AU)
    pubsel = core.select(Dfull, "published")
    pubrho = {lab: core.spearman_perm(*core.vectors(Dfull, pubsel, lab)[:2])[0] for lab in LAB}
    pubdist = pair_distance(Dfull, pubsel)
    print("published: monazite %.3f, xenotime %.3f, mean pair distance %.3f"
          % (pubrho[LAB[0]], pubrho[LAB[1]], pubdist), flush=True)

    rows = []
    for at in (40.0, 50.0, 60.0):
        Dt = build_D(PUB_IN, at)
        selt = core.select(Dt, "published")
        sizes = {dom: len(Dt[dom]["Av"]) for dom in DOMS}
        obs = {lab: core.spearman_perm(*core.vectors(Dt, selt, lab)[:2])[0] for lab in LAB}
        obsdist = pair_distance(Dt, selt)

        draws = {lab: np.empty(POWER_DRAWS) for lab in LAB}
        dists = np.empty(POWER_DRAWS)
        for b in range(POWER_DRAWS):
            Ds = subsample(Dfull, sizes, rng)
            sel = core.select(Ds, "published")
            dists[b] = pair_distance(Ds, sel)
            for lab in LAB:
                a, bb, _ = core.vectors(Ds, sel, lab)
                ra, rb = core.avg_rank(a)[0], core.avg_rank(bb)[0]
                draws[lab][b] = float(core.corr_rows(ra[None, :], rb)[0])

        for lab in LAB:
            dr = draws[lab]
            pct = float((dr <= obs[lab]).mean() * 100)
            rows.append(dict(
                australia_threshold=at, mineral=lab,
                aus_pool_palaeoprot=sizes[DOMS[0]], aus_pool_archaean=sizes[DOMS[1]],
                published_rho=round(pubrho[lab], 3),
                thresholded_rho=round(obs[lab], 3),
                random_subsample_median=round(float(np.median(dr)), 3),
                random_subsample_lo=round(float(np.percentile(dr, 2.5)), 3),
                random_subsample_hi=round(float(np.percentile(dr, 97.5)), 3),
                percentile_of_thresholded=round(pct, 1),
                explained_by_pool_size="Yes" if pct >= 2.5 else "No",
                published_pair_distance=round(pubdist, 3),
                thresholded_pair_distance=round(obsdist, 3),
                random_pair_distance_median=round(float(np.median(dists)), 3),
                random_pair_distance_lo=round(float(np.percentile(dists, 2.5)), 3),
                random_pair_distance_hi=round(float(np.percentile(dists, 97.5)), 3)))
        print("  aus %3.0f  pool %3d/%2d  monazite: thresholded %.3f, random pools of the same "
              "size %.3f [%.3f, %.3f], percentile %.1f"
              % (at, sizes[DOMS[0]], sizes[DOMS[1]], obs[LAB[0]],
                 float(np.median(draws[LAB[0]])), float(np.percentile(draws[LAB[0]], 2.5)),
                 float(np.percentile(draws[LAB[0]], 97.5)),
                 float((draws[LAB[0]] <= obs[LAB[0]]).mean() * 100)), flush=True)
        print("        match quality, mean pair distance: published %.3f, thresholded %.3f, "
              "random pools %.3f [%.3f, %.3f]"
              % (pubdist, obsdist, float(np.median(dists)),
                 float(np.percentile(dists, 2.5)), float(np.percentile(dists, 97.5))), flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "containment_random_control.csv"), index=False)
    pd.set_option("display.width", 300)
    print("\n" + out[out.mineral.isin(LAB[:2])].to_string(index=False))
    print("\nwrote results/containment_random_control.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True, choices=["grid", "confound", "power"])
    A = ap.parse_args()
    if A.part == "grid":
        run_grid()
    elif A.part == "confound":
        run_confound()
    else:
        run_power()
