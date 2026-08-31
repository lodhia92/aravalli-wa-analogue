"""Sensitivity of the selected high-score set to the threshold rule.

The published maps ring the top decile of the light rare-earth score within each domain. This
script establishes what changes if that rule changes, and whether the score distribution supports
a statistically defined anomaly threshold instead.

This script does four things, all on the published per-sample scores.

1. Applies seven threshold rules to the light rare-earth score within each Indian domain: the
   80th, 85th, 90th (published) and 95th percentiles, and three conventional anomaly rules,
   median + 2 MAD, mean + 2 SD and the Tukey upper fence Q3 + 1.5 IQR. Reports the cut value,
   the number of sites selected, and the number of distinct target areas those sites form.

2. Measures whether the target AREAS move with the threshold, as opposed to simply gaining or
   losing sites. Percentile cuts on one score are nested by construction, so the informative
   quantities are how many distinct areas a stricter rule retains and how much new ground a
   looser rule opens.

3. Tests whether the score distribution contains a break that would define an anomalous
   population independently of any percentile choice. Two tests: the largest gap between
   consecutive order statistics in the upper quarter, against a parametric bootstrap under a
   fitted normal; and the BIC of a two-component Gaussian mixture (plain-numpy EM) against a
   single Gaussian.

4. Quantifies the pooled-versus-per-domain choice: how the published decile would divide
   between Sandmata and Mangalwar if the threshold were set on the two domains pooled.

A target area is a group of selected sites joined by single-linkage clustering at 10 km. The
sample spacing is about 2 km, so this groups sites within one drainage catchment without
merging separated ones. The 5 km and 20 km alternatives are reported as a check that the
conclusion does not depend on that distance.

Convention shared with the rest of the study: 100 000 draws, np.random.default_rng(20260827).

Inputs : results/catchment_scores.csv
         results/drainage/{sandmata,mangalwar}_contained.csv  (catchment cross-check only)
Outputs: results/score_thresholds.csv
         results/score_threshold_stability.csv
         results/score_break_tests.csv

Author: Bhavik Harish Lodhia, Curtin University
"""
import os

import numpy as np

import pandas as pd

NDRAW = 100_000

LINK_KM = 10.0

def haversine_matrix(lat, lon):
    la = np.radians(lat)
    lo = np.radians(lon)
    dla = la[:, None] - la[None, :]
    dlo = lo[:, None] - lo[None, :]
    a = np.sin(dla / 2) ** 2 + np.cos(la[:, None]) * np.cos(la[None, :]) * np.sin(dlo / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))

def single_linkage(lat, lon, km):
    """Union-find single-linkage clustering. Returns an integer label per site."""
    n = len(lat)
    if n == 0:
        return np.array([], dtype=int)
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    d = haversine_matrix(lat, lon)
    ii, jj = np.where(np.triu(d <= km, 1))
    for i, j in zip(ii, jj):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    lab = np.array([find(i) for i in range(n)])
    _, out = np.unique(lab, return_inverse=True)
    return out

def n_areas(g, mask, km=LINK_KM):
    sel = g[mask]
    if len(sel) == 0:
        return 0
    return int(single_linkage(sel.lat.values, sel.lon.values, km).max() + 1)

def rules(s):
    """Threshold rules -> cut value. Percentiles use the published convention exactly:
    pandas quantile with linear interpolation, selection by score >= cut."""
    q3, q1 = s.quantile(0.75), s.quantile(0.25)
    med = s.median()
    mad = (s - med).abs().median() * 1.4826
    return [
        ("Top 20 per cent", s.quantile(0.80)),
        ("Top 15 per cent", s.quantile(0.85)),
        ("Top 10 per cent (published)", s.quantile(0.90)),
        ("Top 5 per cent", s.quantile(0.95)),
        ("Median + 2 MAD", med + 2 * mad),
        ("Mean + 2 SD", s.mean() + 2 * s.std()),
        ("Tukey upper fence", q3 + 1.5 * (q3 - q1)),
    ]

def gap_test(s, rng, ndraw=NDRAW):
    """Largest gap between consecutive order statistics in the upper quarter, against a
    parametric bootstrap under a normal fitted to the same data. Also returns where that gap
    sits and how many sites lie above it, without which the gap cannot be interpreted."""
    x = np.sort(s.values)
    assert np.isfinite(x).all(), "non-finite score reached the gap test"
    n = len(x)
    lo = int(np.ceil(0.75 * n))
    d = np.diff(x[lo:])
    j = int(d.argmax())
    obs = float(d[j])
    at = float(x[lo + j])             # score just below the gap
    n_above = int((x > at).sum())
    mu, sd = x.mean(), x.std(ddof=1)
    null = np.empty(ndraw)
    step = 5000
    for a in range(0, ndraw, step):
        b = min(a + step, ndraw)
        sim = rng.normal(mu, sd, size=(b - a, n))
        sim.sort(axis=1)
        null[a:b] = np.diff(sim[:, lo:], axis=1).max(axis=1)
    p = (1 + (null >= obs).sum()) / (1 + ndraw)
    return obs, at, n_above, float(np.median(null)), float(p)

def mixture_crossover(mu, sd, w):
    """Score at which the two mixture components have equal posterior weight, searched on a
    fine grid between the component means. This is the threshold a two-population reading of
    the data would imply."""
    lo, hi = float(min(mu)), float(max(mu))
    span = float(max(sd)) * 6
    g = np.linspace(lo - span, hi + span, 200001)
    dens = w * np.exp(-0.5 * ((g[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
    hi_i = int(np.argmax(mu))
    diff = dens[:, hi_i] - dens[:, 1 - hi_i]
    cross = np.where(np.diff(np.sign(diff)) != 0)[0]
    # the threshold implied by the mixture is the highest score at which the upper component
    # takes over. If there is none, the mixture never separates the data into two sets.
    up = [i for i in cross if diff[i + 1] > 0]
    if not up:
        return float("nan")
    return float(g[up[-1]])

def em2(x, rng, iters=500):
    """Two-component 1-D Gaussian mixture by EM. Returns BIC."""
    n = len(x)
    q = np.quantile(x, [0.25, 0.75])
    mu = q.copy()
    sd = np.full(2, x.std(ddof=1))
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        p = w * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
        tot = p.sum(1, keepdims=True)
        tot[tot == 0] = 1e-300
        r = p / tot
        nk = r.sum(0)
        w = nk / n
        mu = (r * x[:, None]).sum(0) / nk
        sd = np.sqrt((r * (x[:, None] - mu) ** 2).sum(0) / nk)
        sd = np.maximum(sd, 1e-6)
    p = w * np.exp(-0.5 * ((x[:, None] - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))
    ll = np.log(np.maximum(p.sum(1), 1e-300)).sum()
    k = 5
    return 2 * k * np.log(n) - 2 * ll, mu, sd, w

def main():
    global DOMAINS, DR, FP, FP_ALT, HERE, IN_THR, RES, SEED, bic1, bic2, bk, brk, c, cat, cut, cutp, d, dom, g, gap_at, gdom, gg, i, ind, lab, label, lb, ll1, m, med_null, mine, mu2, n_above_gap, n_above_xover, n_before, n_dom, name, nxt, obs, p, pool, pub, pub_areas, pub_catch, pub_mask, rng, rows, s, sc, sd2, sel, sel_catch, st, stab, sub, thr, v, w2, x, xover
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RES = f"{HERE}/results"
    DR = f"{RES}/drainage"
    SEED = 20260827
    FP = "monazite_LREE"          # the fingerprint carried forward; see Section 3.2
    FP_ALT = "xenotime_HREE"      # shown in Figure 5b, not carried forward
    DOMAINS = ["Sandmata", "Mangalwar"]
    IN_THR = 50.0
    rng = np.random.default_rng(SEED)
    sc = pd.read_csv(f"{RES}/catchment_scores.csv")
    ind = sc[sc.survey == "NGCM"].copy()
    ind["k"] = ind.lat.round(4).astype(str) + "_" + ind.lon.round(4).astype(str)
    cat = []
    for name, dom in [("sandmata", "Sandmata"), ("mangalwar", "Mangalwar")]:
        d = pd.read_csv(f"{DR}/{name}_contained.csv")
        d = d[d.pct_in_domain >= IN_THR]
        d["k"] = d.lat.round(4).astype(str) + "_" + d.lon.round(4).astype(str)
        g = ind[ind.domain == dom].merge(
            d[["k", "catchment_km2"]].drop_duplicates("k"), on="k", how="left")
        assert g.catchment_km2.isna().sum() == 0, f"{dom}: unmatched sites"
        lab = np.full(len(g), -1)
        nxt = 0
        for v, gg in g.groupby("catchment_km2"):
            sub = single_linkage(gg.lat.values, gg.lon.values, LINK_KM)
            lab[gg.index.values] = sub + nxt
            nxt += sub.max() + 1
        g["catch_id"] = [f"{dom}_{i}" for i in lab]
        cat.append(g)
    ind = pd.concat(cat, ignore_index=True)
    n_before = len(ind)
    ind = ind[ind[FP].notna()].reset_index(drop=True)
    print("sites dropped for a missing score:", n_before - len(ind))
    print("sites per domain:", dict(ind.domain.value_counts()))
    print("catchment groups:", dict(ind.groupby("domain").catch_id.nunique()))
    rows = []
    stab = []
    for dom in DOMAINS:
        g = ind[ind.domain == dom].reset_index(drop=True)
        s = g[FP]
        pub_mask = s >= s.quantile(0.90)
        pub_areas = set(single_linkage(g.lat.values[pub_mask.values], g.lon.values[pub_mask.values],
                                       LINK_KM)) if pub_mask.any() else set()
        pub_catch = set(g.catch_id[pub_mask])

        for label, cut in rules(s):
            m = s >= cut
            rows.append(dict(
                domain=dom, rule=label, cut=round(float(cut), 3),
                n_sites=int(m.sum()),
                pct_of_domain=round(100 * m.sum() / len(g), 1),
                n_areas_10km=n_areas(g, m.values, 10.0),
                n_areas_5km=n_areas(g, m.values, 5.0),
                n_areas_20km=n_areas(g, m.values, 20.0),
                n_catchments=int(g.catch_id[m].nunique()),
                median_score=round(float(s[m].median()), 3),
                min_score=round(float(s[m].min()), 3),
            ))
            # stability against the published decile, at catchment level
            sel_catch = set(g.catch_id[m])
            stab.append(dict(
                domain=dom, rule=label,
                catchments_shared_with_decile=len(sel_catch & pub_catch),
                catchments_new_vs_decile=len(sel_catch - pub_catch),
                catchments_lost_vs_decile=len(pub_catch - sel_catch),
                pct_sites_in_decile_catchments=round(
                    100 * g.catch_id[m].isin(pub_catch).mean(), 1) if m.any() else np.nan,
            ))
    thr = pd.DataFrame(rows)
    thr.to_csv(f"{RES}/score_thresholds.csv", index=False)
    st = pd.DataFrame(stab)
    st.to_csv(f"{RES}/score_threshold_stability.csv", index=False)
    print("\n--- thresholds ---")
    print(thr.to_string(index=False))
    print("\n--- stability against the published decile ---")
    print(st.to_string(index=False))
    for dom in DOMAINS:
        g = ind[ind.domain == dom]
        pub = int(g[f"top10_{FP}"].sum())
        mine = int(thr[(thr.domain == dom) & (thr.rule.str.startswith("Top 10"))].n_sites.iloc[0])
        assert pub == mine, f"{dom}: published flag {pub} != recomputed {mine}"
    print("\npublished top-decile flags reproduced exactly")
    brk = []
    for dom in DOMAINS:
        gdom = ind[ind.domain == dom]
        s = gdom[FP]
        obs, gap_at, n_above_gap, med_null, p = gap_test(s, rng)
        x = s.values
        bic2, mu2, sd2, w2 = em2(x, rng)
        xover = mixture_crossover(mu2, sd2, w2)
        n_above_xover = int((s >= xover).sum())
        ll1 = (-0.5 * ((x - x.mean()) / x.std(ddof=1)) ** 2
               - np.log(x.std(ddof=1) * np.sqrt(2 * np.pi))).sum()
        bic1 = 2 * 2 * np.log(len(x)) - 2 * ll1
        brk.append(dict(
            domain=dom, n=len(x),
            largest_upper_gap=round(float(obs), 4),
            gap_at_score=round(gap_at, 3),
            sites_above_gap=n_above_gap,
            median_gap_under_normal=round(med_null, 4),
            gap_p=round(p, 5),
            bic_1component=round(float(bic1), 1),
            bic_2component=round(float(bic2), 1),
            delta_bic=round(float(bic1 - bic2), 1),
            mix_means=f"{mu2[0]:.2f}, {mu2[1]:.2f}",
            mix_sds=f"{sd2[0]:.2f}, {sd2[1]:.2f}",
            mix_weights=f"{w2[0]:.2f}, {w2[1]:.2f}",
            mix_crossover=round(xover, 3),
            sites_above_crossover=n_above_xover,
            pct_above_crossover=round(100 * n_above_xover / len(x), 1),
            skew=round(float(s.skew()), 2),
        ))
    bk = pd.DataFrame(brk)
    bk.to_csv(f"{RES}/score_break_tests.csv", index=False)
    print("\n--- break tests ---")
    print(bk.to_string(index=False))
    pool = ind[ind.domain.isin(DOMAINS)]
    cutp = pool[FP].quantile(0.90)
    sel = pool[pool[FP] >= cutp]
    print("\n--- pooled decile, for comparison with the per-domain decile ---")
    print("pooled cut %.3f, n %d" % (cutp, len(sel)))
    print("pooled split:", dict(sel.domain.value_counts()))
    print("per-domain split:", dict(pool[pool[f"top10_{FP}"]].domain.value_counts()))
    for dom in DOMAINS:
        n_dom = (pool.domain == dom).sum()
        print("  %s: pooled selects %d of %d (%.1f per cent of the domain)"
              % (dom, (sel.domain == dom).sum(), n_dom,
                 100 * (sel.domain == dom).sum() / n_dom))
    print("\n--- Figure 5b fingerprint, same rules, for reference ---")
    for dom in DOMAINS:
        g = ind[ind.domain == dom]
        s = g[FP_ALT]
        print(dom, ", ".join("%s n=%d" % (lb, int((s >= c).sum())) for lb, c in rules(s)))

if __name__ == "__main__":
    main()
