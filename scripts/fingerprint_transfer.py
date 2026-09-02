"""
Stability of the fingerprint transfer correlations, and confidence intervals on them.

This script is the single source for every transfer statistic reported in the paper. Tables 2 and 7
are both generated from its output, so the headline probabilities and the resampling results come
from one run and cannot disagree.

Leave-one-out and bootstrap resampling give both the stability of each correlation and its
confidence interval. Two corrections to the earlier analysis are applied here:

1. The false-discovery-rate correction is applied across the FIVE fingerprints, which is the family
   of tests the paper reports and the family described in Section 2.4. The earlier Table 2 took
   its q-values from a correction across the eight individual pathfinder elements, a different
   family, which gave the single xenotime element a weaker correction because monazite contributes
   three elements to that list.
2. The permutation test uses 100 000 relabellings rather than 5 000, so the reported probabilities
   are stable to the third decimal and do not depend on the random seed. Permutations are generated
   vectorised, so this is no slower than the original.

The twenty pairs are NOT twenty independent catchments: Palaeoproterozoic pairs 1 and 4 share an
identical catchment (NGSA sites 2007190003 and 2007190213 lie about 150 m apart) and pairs 5, 6 and
7 are nested on the Gascoyne drainage. The whole analysis is therefore repeated with the duplicated
pair removed.

Outputs
  results/fingerprint_transfer.csv          one row per fingerprint per variant; source for Tables 2 and 7
  results/fingerprint_transfer_jackknife.csv    every leave-one-out refit, for inspection
Run
  python scripts/fingerprint_transfer.py

Author: Bhavik Harish Lodhia, Curtin University
"""

import csv
import os

import numpy as np
import pandas as pd

import paths
from aravalli_wa.composition import logr
from aravalli_wa.constants import (
    AUS_THR_PCT,
    IN_THR_PCT,
    P_MASS_FRACTION_OF_P2O5,
    TI_MASS_FRACTION_OF_TIO2,
    WT_PCT_TO_MG_KG,
)
from aravalli_wa.stats import avg_rank, bh, corr_rows

NPERM, NBOOT = 100000, 10000


def spearman_perm(a, b, nperm=NPERM):
    """rho and a two-sided permutation p; permutations generated vectorised."""
    ra, rb = avg_rank(a)[0], avg_rank(b)[0]
    rho = float(corr_rows(ra[None, :], rb)[0])
    n = len(rb)
    perms = rng.random((nperm, n)).argsort(axis=1)  # nperm random permutations at once
    null = np.abs(corr_rows(rb[perms], ra))
    return rho, (np.sum(null >= abs(rho)) + 1) / (nperm + 1)


def india(name):
    d = pd.read_csv(os.path.join(DR, f"{name}_contained.csv"))
    d = d[d.pct_in_domain >= IN_THR_PCT]
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
    d = d[d.pct_in_domain >= AUS_THR_PCT]
    ids = set(pd.to_numeric(d["id"], errors="coerce").dropna())
    m = ng[ng.index.isin(ids)].copy()
    m["sid"] = [f"{label}_{i}" for i in m.index]
    return m.reset_index()


def vectors(pairs, grp):
    xs, ys = [], []
    for _, r in pairs.iterrows():
        if r.india_sid not in ic.index or r.aus_sid not in acn.index:
            continue
        di, da = ic.loc[r.india_sid], acn.loc[r.aus_sid]
        xs.append(np.nanmean([(np.log(di[e]) - mu_in[e]) / sd_in[e] for e in grp]))
        ys.append(np.nanmean([(np.log(da[e]) - mu_au[e]) / sd_au[e] for e in grp]))
    a, b = np.asarray(xs, float), np.asarray(ys, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return a[ok], b[ok]


def main():
    global \
        Ac, \
        Av, \
        Bc, \
        CI, \
        DR, \
        H, \
        HERE, \
        Iv, \
        LAB, \
        MATCH, \
        MIN, \
        NG, \
        O, \
        P, \
        PATH, \
        PROJ, \
        RA, \
        RB, \
        RES, \
        VARIANTS, \
        _, \
        a, \
        acn, \
        adf, \
        ar, \
        b, \
        bi, \
        bs, \
        c, \
        cidx, \
        d, \
        da, \
        den, \
        detail, \
        di, \
        dom, \
        el, \
        f, \
        full, \
        hi, \
        i, \
        ic, \
        idf, \
        idx, \
        j, \
        jk, \
        k, \
        l, \
        lab, \
        lo, \
        loo_pass, \
        loo_rho, \
        m_, \
        ma, \
        mang, \
        mi, \
        mu_au, \
        mu_in, \
        n, \
        ng, \
        ngcm_all, \
        ngsa_all, \
        ordered, \
        out, \
        p, \
        pairs, \
        ps, \
        qi, \
        qs, \
        r, \
        rng, \
        rob, \
        rows, \
        rs, \
        sa, \
        sand, \
        sd_au, \
        sd_in, \
        se, \
        si, \
        v, \
        vecs, \
        vname, \
        wapp, \
        yiln, \
        z, \
        za, \
        zi
    csv.field_size_limit(10**7)
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    PROJ = os.path.dirname(os.path.dirname(HERE))
    RES = os.path.join(HERE, "results")
    DR = os.path.join(RES, "drainage")
    CI = dict(La=0.237, Yb=0.170, Sm=0.148, Eu=0.0580, Gd=0.199)
    MATCH = ["Th/Sc", "La/Sc", "Th/Co", "EuEu", "La/Yb_n", "Nb/Y"]
    PATH = ["Zr", "Hf", "Ti", "Ce", "Nd", "Pr", "Dy", "P"]
    MIN = {
        "monazite (Ce,Nd,Pr)": ["Ce", "Nd", "Pr"],
        "xenotime (Dy)": ["Dy"],
        "zircon (Zr,Hf)": ["Zr", "Hf"],
        "Ti-oxide (Ti)": ["Ti"],
        "apatite (P)": ["P"],
    }
    rng = np.random.default_rng(20260827)
    ar = pd.read_csv(paths.NGCM_TABLE)
    for c in [
        "Th",
        "Sc",
        "Co",
        "La",
        "Eu",
        "Sm",
        "Gd",
        "Yb",
        "Nb",
        "Y",
        "TiO2",
        "P2O5",
        "Zr",
        "Hf",
        "Ce",
        "Nd",
        "Pr",
        "Dy",
        "LAT",
        "LON",
    ]:
        ar[c] = pd.to_numeric(ar[c], errors="coerce")
    ar["Ti"] = ar["TiO2"] * WT_PCT_TO_MG_KG * TI_MASS_FRACTION_OF_TIO2
    ar["P"] = ar["P2O5"] * WT_PCT_TO_MG_KG * P_MASS_FRACTION_OF_P2O5
    ar["k"] = ar.LAT.round(4).astype(str) + "_" + ar.LON.round(4).astype(str)
    sand, mang = india("sandmata"), india("mangalwar")
    NG = paths.NGSA
    with open(NG, encoding="latin-1") as f:
        H = list(csv.reader(f))[11]
    cidx = {
        el: findcol(el) for el in ["Th", "Sc", "Nb", "Y", "La", "Yb", "Co", "Eu", "Sm", "Gd"] + PATH
    }
    ordered = sorted(((k, v) for k, v in cidx.items() if v is not None), key=lambda kv: kv[1])
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
    ng = ng.groupby("SITEID").median(numeric_only=True)
    wapp = aus(os.path.join(DR, "wa_palaeoprot_contained.csv"), "WA_PP")
    yiln = aus(os.path.join(DR, "yilgarn_contained.csv"), "Y+N")
    ngcm_all, ngsa_all = pd.concat([sand, mang]), pd.concat([wapp, yiln])
    mi, si = logr(ngcm_all).mean(), logr(ngcm_all).std(ddof=0)
    ma, sa = logr(ngsa_all).mean(), logr(ngsa_all).std(ddof=0)
    rows = []
    for idf, adf, dom in [(sand, wapp, "Palaeoproterozoic"), (mang, yiln, "Archaean")]:
        zi = ((logr(idf) - mi) / si)[MATCH].dropna()
        za = ((logr(adf) - ma) / sa)[MATCH].dropna()
        Iv, Av = zi.values, za.values
        for a in range(len(Av)):
            di = np.sqrt(((Iv - Av[a]) ** 2).sum(1))
            bi = int(di.argmin())
            da = np.sqrt(((Av - Iv[bi]) ** 2).sum(1))
            rows.append(
                dict(
                    domain=dom,
                    india_sid=idf.loc[zi.index[bi], "sid"],
                    aus_sid=adf.loc[za.index[a], "sid"],
                    dist=float(di[bi]),
                    mnn=(int(da.argmin()) == a),
                )
            )
    P = pd.DataFrame(rows)
    rob = pd.concat(
        [
            P[P.domain == d].sort_values(["mnn", "dist"], ascending=[False, True]).head(10)
            for d in ("Palaeoproterozoic", "Archaean")
        ]
    ).reset_index(drop=True)
    rob.insert(0, "pair_no", list(range(1, 11)) * 2)
    print("pairs built: %d" % len(rob), flush=True)
    mu_in = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0)).mean()
    sd_in = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0)).std(ddof=0)
    mu_au = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0)).mean()
    sd_au = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0)).std(ddof=0)
    ic, acn = ngcm_all.set_index("sid"), ngsa_all.set_index("sid")
    LAB = list(MIN)
    out, detail = [], []
    VARIANTS = [
        ("all 20 pairs", rob),
        (
            "19 pairs, duplicate catchment removed",
            rob.drop(rob[(rob.domain == "Palaeoproterozoic") & (rob.pair_no == 4)].index),
        ),
    ]
    for vname, pairs in VARIANTS:
        print("\n%s" % vname, flush=True)
        vecs = {lab: vectors(pairs, MIN[lab]) for lab in LAB}
        n = len(vecs[LAB[0]][0])
        full = {lab: spearman_perm(*vecs[lab]) for lab in LAB}
        qs = dict(zip(LAB, bh([full[l][1] for l in LAB])))

        loo_rho = {lab: [] for lab in LAB}
        loo_pass = {lab: 0 for lab in LAB}
        for i in range(n):
            m_ = np.ones(n, bool)
            m_[i] = False
            ps, rs = [], {}
            for lab in LAB:
                a, b = vecs[lab]
                r, p = spearman_perm(a[m_], b[m_])
                ps.append(p)
                rs[lab] = r
            qi = bh(ps)
            for j, lab in enumerate(LAB):
                loo_rho[lab].append(rs[lab])
                if qi[j] < 0.05:
                    loo_pass[lab] += 1
                detail.append(
                    dict(
                        variant=vname,
                        mineral=lab,
                        dropped_pair=i + 1,
                        rho=round(rs[lab], 3),
                        q=round(float(qi[j]), 4),
                    )
                )
            print("  leave-one-out %d/%d" % (i + 1, n), flush=True)

        for lab in LAB:
            a, b = vecs[lab]
            r, p = full[lab]
            jk = np.array(loo_rho[lab], float)
            idx = rng.integers(0, n, size=(NBOOT, n))
            RA, RB = avg_rank(a[idx]), avg_rank(b[idx])
            Ac = RA - RA.mean(1, keepdims=True)
            Bc = RB - RB.mean(1, keepdims=True)
            den = np.sqrt((Ac**2).sum(1) * (Bc**2).sum(1))
            with np.errstate(invalid="ignore", divide="ignore"):
                bs = np.where(den > 0, (Ac * Bc).sum(1) / den, np.nan)
            bs = bs[np.isfinite(bs)]
            lo, hi = np.percentile(bs, [2.5, 97.5])
            z = np.arctanh(np.clip(r, -0.999999, 0.999999))
            se = 1 / np.sqrt(n - 3)
            out.append(
                dict(
                    variant=vname,
                    mineral=lab,
                    n=n,
                    rho=round(r, 3),
                    perm_p=round(p, 4),
                    q=round(qs[lab], 4),
                    jk_min=round(jk.min(), 3),
                    jk_max=round(jk.max(), 3),
                    loo_pass_q05="%d/%d" % (loo_pass[lab], n),
                    boot_lo=round(float(lo), 3),
                    boot_hi=round(float(hi), 3),
                    fisher_lo=round(float(np.tanh(z - 1.96 * se)), 3),
                    fisher_hi=round(float(np.tanh(z + 1.96 * se)), 3),
                )
            )
    O = pd.DataFrame(out)
    O.to_csv(os.path.join(RES, "fingerprint_transfer.csv"), index=False)
    pd.DataFrame(detail).to_csv(
        os.path.join(RES, "fingerprint_transfer_jackknife.csv"), index=False
    )
    pd.set_option("display.width", 250)
    for v in O.variant.unique():
        print("\n==== %s ====" % v)
        print(O[O.variant == v].drop(columns="variant").to_string(index=False))
    print("\nNPERM=%d NBOOT=%d seed=20260827" % (NPERM, NBOOT))
    print("wrote results/fingerprint_transfer.csv and results/fingerprint_transfer_jackknife.csv")


if __name__ == "__main__":
    main()
