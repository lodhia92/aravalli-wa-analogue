"""
Is the fingerprint transfer a product of the matching algorithm rather than of geology?

Choosing the closest mutual nearest neighbours could manufacture similarity, in which case the
transfer correlation would partly reflect the selection procedure rather than shared provenance.
Four tests separate the two: alternative matching algorithms, a varied number of retained pairs,
changed matching thresholds, and a control population of randomly selected analogues.

All four are implemented here, together with a diagnostic that separates the part of the
correlation carried by the difference between the two domain pools from the part carried by the
pair-level matching. Data loading, standardisation and scoring are copied unchanged from
scripts/fingerprint_transfer.py so that the published pairing reproduces exactly; the published variant
in this script must return rho = 0.696 for monazite and 0.506 for xenotime.

Permutation probabilities use 100 000 relabellings, as in fingerprint_transfer.py. One bank of
permutations is generated per sample size and reused across fingerprints and variants.

Outputs
  results/matching_variants.csv   one row per matching variant per fingerprint (algorithms, counts,
                            thresholds)
  results/matching_random_control.csv    random-pairing control, one row per fingerprint
  results/matching_domain_contrast.csv     domain-contrast diagnostic, one row per fingerprint
  results/_matching_cache.npz     prepared score and ratio arrays, so later parts start quickly
Run
  python scripts/matching_sensitivity.py --part variants
  python scripts/matching_sensitivity.py --part control
  python scripts/matching_sensitivity.py --part domain
  python scripts/matching_sensitivity.py --part summary

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import argparse
import csv
import os

import numpy as np
import pandas as pd
import paths
from aravalli_wa.composition import logr, ratios
from aravalli_wa.stats import avg_rank, bh, corr_rows
from aravalli_wa.constants import P_MASS_FRACTION_OF_P2O5, TI_MASS_FRACTION_OF_TIO2, WT_PCT_TO_MG_KG

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
DOMS = ["Palaeoproterozoic", "Archaean"]
NPERM, NDRAW, KPUB = 100000, 10000, 10
SEED = 20260827
rng = np.random.default_rng(SEED)
_BANK = {}


def bank(n):
    if n not in _BANK:
        _BANK[n] = rng.random((NPERM, n)).argsort(axis=1)
    return _BANK[n]


def spearman_perm(a, b, blocks=None):
    """rho and a two-sided permutation p. blocks restricts relabelling to within groups."""
    ra, rb = avg_rank(a)[0], avg_rank(b)[0]
    rho = float(corr_rows(ra[None, :], rb)[0])
    n = len(rb)
    P = bank(n)
    if blocks is not None:
        P = np.empty((NPERM, n), int)
        for g in np.unique(blocks):
            pos = np.flatnonzero(blocks == g)
            P[:, pos] = pos[rng.random((NPERM, len(pos))).argsort(axis=1)]
    null = np.abs(corr_rows(rb[P], ra))
    return rho, (np.sum(null >= abs(rho)) + 1) / (NPERM + 1)


def hungarian(cost):
    """Minimum-cost one-to-one assignment, Jonker-Volgenant shortest augmenting path.
    cost is (n, m) with n <= m. Verified against scipy.optimize.linear_sum_assignment on
    200 random rectangular matrices before use."""
    cost = np.asarray(cost, float)
    n, m = cost.shape
    u, v = np.zeros(n + 1), np.zeros(m + 1)
    p, way = np.zeros(m + 1, int), np.zeros(m + 1, int)
    for i in range(1, n + 1):
        p[0], j0 = i, 0
        minv = np.full(m + 1, np.inf)
        used = np.zeros(m + 1, bool)
        while True:
            used[j0] = True
            i0 = p[j0]
            js = np.flatnonzero(~used[1:]) + 1
            cur = cost[i0 - 1, js - 1] - u[i0] - v[js]
            better = cur < minv[js]
            minv[js[better]] = cur[better]
            way[js[better]] = j0
            k = int(np.argmin(minv[js]))
            delta, j1 = minv[js[k]], int(js[k])
            uj = np.flatnonzero(used)
            u[p[uj]] += delta
            v[uj] -= delta
            nuj = np.flatnonzero(~used)
            minv[nuj] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    out = np.full(n, -1, int)
    for j in range(1, m + 1):
        if p[j] > 0:
            out[p[j] - 1] = j - 1
    return out


def prepare():
    """Load and standardise exactly as fingerprint_transfer.py, then return per-domain arrays."""
    cache = os.path.join(RES, "_matching_cache.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {k: z[k].item() for k in ["D"]}["D"]
    ar = pd.read_csv(paths.NGCM_TABLE)
    for c in ["Th", "Sc", "Co", "La", "Eu", "Sm", "Gd", "Yb", "Nb", "Y", "TiO2", "P2O5",
              "Zr", "Hf", "Ce", "Nd", "Pr", "Dy", "LAT", "LON"]:
        ar[c] = pd.to_numeric(ar[c], errors="coerce")
    ar["Ti"] = ar["TiO2"] * WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2
    ar["P"] = ar["P2O5"] * WT_PCT_TO_MG_KG * P_MASS_FRACTION_OF_P2O5
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
        zi = ((logr(idf) - mi) / si)[MATCH].dropna()
        za = ((logr(adf) - ma) / sa)[MATCH].dropna()
        isub, asub = idf.loc[zi.index], adf.loc[za.index]
        D[dom] = dict(Iv=zi.values, Av=za.values,
                      isid=isub["sid"].values, asid=asub["sid"].values,
                      SI=score(isub, mu_in, sd_in), SA=score(asub, mu_au, sd_au))
    np.savez(os.path.join(RES, "_matching_cache.npz"), D=np.array(D, dtype=object))
    return D


def candidates(Iv, Av, W=None):
    """For each Australian site: nearest Indian site, its distance, and whether mutual."""
    if W is None:
        dI = np.sqrt(((Iv[None, :, :] - Av[:, None, :]) ** 2).sum(2))
    else:
        Dif = Iv[None, :, :] - Av[:, None, :]
        dI = np.sqrt(np.einsum("aij,jk,aik->ai", Dif, W, Dif))
    bi = dI.argmin(1)
    d = dI[np.arange(len(Av)), bi]
    back = dI.argmin(0)
    mnn = back[bi] == np.arange(len(Av))
    return bi, d, mnn, dI


def select(D, how, k=KPUB, pct=None):
    """Return {domain: (ipos, apos)} for one matching variant."""
    out = {}
    for dom in DOMS:
        Iv, Av = D[dom]["Iv"], D[dom]["Av"]
        if how == "mahalanobis":
            X = np.vstack([Iv, Av])
            W = np.linalg.pinv(np.cov(X.T))
            bi, d, mnn, _ = candidates(Iv, Av, W)
        else:
            bi, d, mnn, dI = candidates(Iv, Av)
        if how in ("published", "mahalanobis"):
            order = np.lexsort((d, ~mnn))[:k]
            out[dom] = (bi[order], order)
        elif how == "oneway":
            order = np.argsort(d)[:k]
            out[dom] = (bi[order], order)
        elif how == "strict_mnn":
            idx = np.flatnonzero(mnn)
            order = idx[np.argsort(d[idx])][:k]
            out[dom] = (bi[order], order)
        elif how == "greedy_unique":
            used_i, keep = set(), []
            for a in np.argsort(d):
                if bi[a] in used_i:
                    continue
                used_i.add(bi[a])
                keep.append(a)
                if len(keep) == k:
                    break
            keep = np.array(keep, int)
            out[dom] = (bi[keep], keep)
        elif how == "optimal":
            C = dI if dI.shape[0] <= dI.shape[1] else dI.T
            asg = hungarian(C)
            if dI.shape[0] <= dI.shape[1]:
                apos, ipos = np.arange(len(asg)), asg
            else:
                apos, ipos = asg, np.arange(len(asg))
            dd = dI[apos, ipos]
            order = np.argsort(dd)[:k]
            out[dom] = (ipos[order], apos[order])
        elif how == "caliper":
            thr = np.percentile(d, pct)
            idx = np.flatnonzero((d <= thr) & mnn)
            idx = idx[np.argsort(d[idx])]
            out[dom] = (bi[idx], idx)
        else:
            raise ValueError(how)
    return out


def vectors(D, sel, lab):
    a, b, blocks = [], [], []
    for gi, dom in enumerate(DOMS):
        ipos, apos = sel[dom]
        a.append(D[dom]["SI"][lab][ipos])
        b.append(D[dom]["SA"][lab][apos])
        blocks.append(np.full(len(ipos), gi))
    a, b, blocks = np.concatenate(a), np.concatenate(b), np.concatenate(blocks)
    ok = np.isfinite(a) & np.isfinite(b)
    return a[ok], b[ok], blocks[ok]


def pairkeys(D, sel):
    return {(D[d]["isid"][i], D[d]["asid"][j]) for d in DOMS for i, j in zip(*sel[d])}


def run_variants(D):
    pub = select(D, "published")
    pubk = pairkeys(D, pub)
    VAR = [("algorithm", "mutual nearest neighbour, top 10 per domain (published)", "published",
            KPUB, None),
           ("algorithm", "one-directional nearest neighbour, mutuality not required", "oneway",
            KPUB, None),
           ("algorithm", "greedy unique, each Indian site used once", "greedy_unique", KPUB, None),
           ("algorithm", "globally optimal one-to-one assignment", "optimal", KPUB, None),
           ("algorithm", "Mahalanobis distance in place of Euclidean", "mahalanobis", KPUB, None),
           ("algorithm", "mutual pairs only, no distance ranking cap", "strict_mnn", 10 ** 6, None)]
    for k in (5, 8, 12, 15, 20, 30):
        VAR.append(("pair count", "top %d pairs per domain" % k, "published", k, None))
    for pct in (25, 50, 75, 100):
        VAR.append(("threshold", "mutual pairs within the %dth distance percentile" % pct,
                    "caliper", KPUB, pct))
    rows = []
    for part, name, how, k, pct in VAR:
        sel = select(D, how, k=k, pct=pct)
        shared = len(pairkeys(D, sel) & pubk)
        npairs = sum(len(sel[d][0]) for d in DOMS)
        res = {lab: spearman_perm(*vectors(D, sel, lab)[:2]) for lab in LAB}
        qs = dict(zip(LAB, bh([res[l][1] for l in LAB])))
        for lab in LAB:
            r, p = res[lab]
            rows.append(dict(part=part, variant=name, mineral=lab,
                             n_selected=npairs, n=len(vectors(D, sel, lab)[0]),
                             shared_with_published=shared, rho=round(r, 3),
                             perm_p=round(p, 4), q=round(float(qs[lab]), 4),
                             transfers="Yes" if qs[lab] < 0.05 else "No"))
        print("  %-58s n=%-3d shared=%-3d monazite rho=%.3f q=%.4f"
              % (name[:58], npairs, shared, res[LAB[0]][0], qs[LAB[0]]), flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "matching_variants.csv"), index=False)
    print("wrote results/matching_variants.csv")


def run_control(D):
    """Control population: analogue pairs drawn at random from the same eligible pools."""
    pub = select(D, "published")
    rows = []
    for lab in LAB:
        a_obs, b_obs, _ = vectors(D, pub, lab)
        rho_obs = float(corr_rows(avg_rank(a_obs), avg_rank(b_obs)[0])[0])
        nulls = []
        blocks = []
        for dom in DOMS:
            npd = len(pub[dom][0])
            ip = np.flatnonzero(np.isfinite(D[dom]["SI"][lab]))
            ap = np.flatnonzero(np.isfinite(D[dom]["SA"][lab]))
            blocks.append((dom, npd, ip, ap))
        done = 0
        while done < NDRAW:
            m = min(2000, NDRAW - done)
            A, B = [], []
            for dom, npd, ip, ap in blocks:
                si = ip[rng.random((m, len(ip))).argpartition(npd, axis=1)[:, :npd]]
                sa = ap[rng.random((m, len(ap))).argpartition(npd, axis=1)[:, :npd]]
                A.append(D[dom]["SI"][lab][si])
                B.append(D[dom]["SA"][lab][sa])
            A, B = np.hstack(A), np.hstack(B)
            RA, RB = avg_rank(A), avg_rank(B)
            Ac, Bc = RA - RA.mean(1, keepdims=True), RB - RB.mean(1, keepdims=True)
            den = np.sqrt((Ac ** 2).sum(1) * (Bc ** 2).sum(1))
            with np.errstate(invalid="ignore", divide="ignore"):
                nulls.append(np.where(den > 0, (Ac * Bc).sum(1) / den, np.nan))
            done += m
            print("    %s %d/%d" % (lab[:20], done, NDRAW), flush=True)
        nl = np.concatenate(nulls)
        nl = nl[np.isfinite(nl)]
        rows.append(dict(mineral=lab, n_pairs=len(a_obs), rho_matched=round(rho_obs, 3),
                         n_draws=len(nl), rho_random_mean=round(float(nl.mean()), 3),
                         rho_random_sd=round(float(nl.std(ddof=1)), 3),
                         rho_random_p50=round(float(np.percentile(nl, 50)), 3),
                         rho_random_p95=round(float(np.percentile(nl, 95)), 3),
                         rho_random_p99=round(float(np.percentile(nl, 99)), 3),
                         p_vs_random=round(float((np.sum(nl >= rho_obs) + 1) / (len(nl) + 1)), 4),
                         frac_random_above_matched=round(float(np.mean(nl >= rho_obs)), 4)))
        print("  %-22s matched %.3f  random mean %.3f  p95 %.3f  p=%.4f"
              % (lab[:22], rho_obs, nl.mean(), np.percentile(nl, 95), rows[-1]["p_vs_random"]),
              flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(RES, "matching_random_control.csv"), index=False)
    print("wrote results/matching_random_control.csv")


def run_domain(D):
    """How much of the correlation is the difference between the two domain pools?"""
    pub = select(D, "published")
    rows = []
    for lab in LAB:
        a, b, g = vectors(D, pub, lab)
        r_all, p_all = spearman_perm(a, b)
        r_within, p_within = spearman_perm(a, b, blocks=g)
        rec = dict(mineral=lab, n=len(a), rho_all=round(r_all, 3), p_all=round(p_all, 4),
                   p_within_domain_permutation=round(p_within, 4))
        for gi, dom in enumerate(DOMS):
            m = g == gi
            if m.sum() >= 4:
                r_d, p_d = spearman_perm(a[m], b[m])
                rec["rho_" + dom[:5]] = round(r_d, 3)
                rec["p_" + dom[:5]] = round(p_d, 4)
                rec["n_" + dom[:5]] = int(m.sum())
            si = D[dom]["SI"][lab]
            sa = D[dom]["SA"][lab]
            rec["poolmean_india_" + dom[:5]] = round(float(np.nanmean(si)), 3)
            rec["poolmean_aus_" + dom[:5]] = round(float(np.nanmean(sa)), 3)
        ac, bc = a.copy(), b.copy()
        for gi in np.unique(g):
            m = g == gi
            ac[m] -= ac[m].mean()
            bc[m] -= bc[m].mean()
        r_c, p_c = spearman_perm(ac, bc, blocks=g)
        rec["rho_domain_centred"] = round(r_c, 3)
        rec["p_domain_centred"] = round(p_c, 4)
        rows.append(rec)
        print("  %-22s all %.3f  centred %.3f (p %.4f)" % (lab[:22], r_all, r_c, p_c), flush=True)
    df = pd.DataFrame(rows)
    df["q_domain_centred"] = np.round(bh(df["p_domain_centred"].values), 4)
    df.to_csv(os.path.join(RES, "matching_domain_contrast.csv"), index=False)
    print("wrote results/matching_domain_contrast.csv")


def summary():
    pd.set_option("display.width", 250)
    for f in ("matching_variants.csv", "matching_random_control.csv", "matching_domain_contrast.csv"):
        p = os.path.join(RES, f)
        if os.path.exists(p):
            print("\n==== %s ====" % f)
            print(pd.read_csv(p).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True,
                    choices=["variants", "control", "domain", "summary", "prepare"])
    A = ap.parse_args()
    if A.part == "summary":
        summary()
    else:
        D = prepare()
        print("prepared: " + ", ".join("%s india %d aus %d"
                                       % (d, len(D[d]["Iv"]), len(D[d]["Av"])) for d in DOMS),
              flush=True)
        if A.part == "variants":
            run_variants(D)
        elif A.part == "control":
            run_control(D)
        elif A.part == "domain":
            run_domain(D)
