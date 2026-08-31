"""
Grain size, hydraulic sorting and catchment morphometry.

Monazite and xenotime are dense heavy minerals, so an explanation in which zircon and the titanium
oxides fail to transfer through hydraulic sorting while the rare-earth phosphates transfer as bulk
rare-earth character has to be demonstrated rather than asserted. This script tests whether
sediment grain size, catchment area and sediment transport distance correlate with the fingerprint
scores.

Two of those are testable from data already held.

GRAIN SIZE. The National Geochemical Survey of Australia file carries SIX rows per site: three
grain-size fractions (Bulk, <2 mm, <75 um) crossed with two sampling depths (top and bottom of
soil). Every earlier script in this study collapsed all six with a median per site, so the
Australian value used in the paper is a median across fractions and depths rather than a defined
fraction. This script recomputes the transfer on each fraction and each depth separately, in two
modes: holding the twenty published pairs fixed, which isolates the effect on the fingerprint
scores, and rebuilding the matching from the same subset, which is the full recomputation. The
Indian samples are uniformly minus 120 mesh, so the Indian side does not vary.

MORPHOMETRY. Catchment area is carried in results/pair_catchments_{australia,india}.geojson.
Sediment transport distance is approximated by the great-circle distance from the sample site to
the centroid of its own catchment. Relief and drainage gradient are NOT tested: no digital
elevation model covering Western Australia or the Aravalli is held in this project, and the
response says so plainly rather than substituting a proxy for them.

Permutations, seed and scoring follow fingerprint_transfer.py. The published variant must return monazite
rho = 0.696 and xenotime rho = 0.506.

Outputs
  results/grain_size_fractions.csv     one row per grain-size or depth subset per mode per fingerprint
  results/catchment_morphometry.csv   Spearman of each fingerprint score against area and site-to-centroid
                               distance, each side separately
Run
  python scripts/grain_size_morphometry.py --part fractions
  python scripts/grain_size_morphometry.py --part morphometry

Author: Bhavik Harish Lodhia, Curtin University
"""
import argparse
import csv
import json
import os

import numpy as np
import pandas as pd
import paths
from aravalli_wa.composition import logr
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
NPERM = 100000
rng = np.random.default_rng(20260827)
_BANK = {}


def spearman_perm(a, b):
    ra, rb = avg_rank(a)[0], avg_rank(b)[0]
    rho = float(corr_rows(ra[None, :], rb)[0])
    n = len(rb)
    if n not in _BANK:
        _BANK[n] = rng.random((NPERM, n)).argsort(axis=1)
    null = np.abs(corr_rows(rb[_BANK[n]], ra))
    return rho, (np.sum(null >= abs(rho)) + 1) / (NPERM + 1)


# The Indian pool does not vary across the grain-size treatments; it is built once.
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
ngcm_all = pd.concat([sand, mang])
mi, si = logr(ngcm_all).mean(), logr(ngcm_all).std(ddof=0)

# The Australian pool is rebuilt per treatment: two size fractions by two sample depths,
# plus the published median across all of them.
NG = paths.NGSA
with open(NG, encoding="latin-1") as f:
    H = list(csv.reader(f))[11]
GS_COL, DEPTH_COL = 8, 9
assert H[GS_COL].strip().upper() == "GRAIN SIZE" and H[DEPTH_COL].strip().upper() == "DEPTH"


def findcol(el):
    for meth in ("ICP-MS", "XRF"):
        for i, c in enumerate(H):
            if c.strip().startswith(f"{el} {meth}"):
                return i


cidx = {el: findcol(el) for el in ["Th", "Sc", "Nb", "Y", "La", "Yb", "Co", "Eu", "Sm",
                                   "Gd"] + PATH}
ordered = sorted(((k, v) for k, v in cidx.items() if v is not None), key=lambda kv: kv[1])
use = [0, GS_COL, DEPTH_COL] + [v for _, v in ordered]
NGRAW = pd.read_csv(NG, header=None, skiprows=12, usecols=use, encoding="latin-1",
                    low_memory=False)
NGRAW.columns = ["SITEID", "GRAIN_SIZE", "DEPTH"] + [k for k, _ in ordered]
for c in NGRAW.columns:
    if c not in ("GRAIN_SIZE", "DEPTH"):
        NGRAW[c] = pd.to_numeric(NGRAW[c], errors="coerce")
NGRAW["GRAIN_SIZE"] = NGRAW["GRAIN_SIZE"].astype(str).str.strip()
NGRAW["DEPTH"] = NGRAW["DEPTH"].astype(str).str.strip()

# The Bulk rows carry no ICP-MS or XRF geochemistry for any element used here (checked: 0 of
# 2 630 rows), so the "median per site" of the published pipeline is in practice a median of the
# two real fractions across the two depths, four values per element.
SUBSETS = [("all rows, median per site (as published)", None, None),
           ("fine fraction only (<75 um)", "<75 µm", None),
           ("coarse fraction only (<2 mm)", "<2 mm", None),
           ("fine fraction, top of soil only", "<75 µm", "TOS"),
           ("top of soil only, both fractions", None, "TOS"),
           ("bottom of soil only, both fractions", None, "BOS")]


def ngsa(gs=None, depth=None):
    d = NGRAW
    if gs is not None:
        d = d[d.GRAIN_SIZE == gs]
    if depth is not None:
        d = d[d.DEPTH == depth]
    g = d.groupby("SITEID").median(numeric_only=True)
    if g.empty or g[PATH].notna().sum().sum() == 0:
        raise SystemExit("subset gs=%r depth=%r has no geochemistry" % (gs, depth))
    return g


def aus_pool(ng, path, label):
    d = pd.read_csv(path)
    d = d[d.pct_in_domain >= AUS_THR]
    ids = set(pd.to_numeric(d["id"], errors="coerce").dropna())
    m = ng[ng.index.isin(ids)].copy()
    m["sid"] = [f"{label}_{i}" for i in m.index]
    return m.reset_index()


def build(ng):
    """Return the twenty pairs and the two pools for a given Australian subset."""
    wapp = aus_pool(ng, os.path.join(DR, "wa_palaeoprot_contained.csv"), "WA_PP")
    yiln = aus_pool(ng, os.path.join(DR, "yilgarn_contained.csv"), "Y+N")
    ngsa_all = pd.concat([wapp, yiln])
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
            rows.append(dict(domain=dom, india_sid=idf.loc[zi.index[bi], "sid"],
                             aus_sid=adf.loc[za.index[a], "sid"], dist=float(di[bi]),
                             mnn=(int(da.argmin()) == a)))
    P = pd.DataFrame(rows)
    pairs = pd.concat([P[P.domain == d].sort_values(["mnn", "dist"], ascending=[False, True]).head(10)
                       for d in ("Palaeoproterozoic", "Archaean")]).reset_index(drop=True)
    return pairs, ngsa_all


def scores(pairs, ngsa_all):
    mu_in = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0)).mean()
    sd_in = np.log(ngcm_all[PATH].where(ngcm_all[PATH] > 0)).std(ddof=0)
    mu_au = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0)).mean()
    sd_au = np.log(ngsa_all[PATH].where(ngsa_all[PATH] > 0)).std(ddof=0)
    ic, acn = ngcm_all.set_index("sid"), ngsa_all.set_index("sid")
    out = {}
    for lab, els in MIN.items():
        xs, ys = [], []
        for _, r in pairs.iterrows():
            if r.india_sid not in ic.index or r.aus_sid not in acn.index:
                xs.append(np.nan)
                ys.append(np.nan)
                continue
            di, da = ic.loc[r.india_sid], acn.loc[r.aus_sid]
            xs.append(np.nanmean([(np.log(di[e]) - mu_in[e]) / sd_in[e] for e in els]))
            ys.append(np.nanmean([(np.log(da[e]) - mu_au[e]) / sd_au[e] for e in els]))
        out[lab] = (np.asarray(xs, float), np.asarray(ys, float))
    return out


def run_fractions():
    base_pairs, _ = build(ngsa())
    base_key = list(zip(base_pairs.india_sid, base_pairs.aus_sid))
    rows = []
    for name, gs, depth in SUBSETS:
        ng = ngsa(gs, depth)
        pairs_new, pool = build(ng)
        shared = len(set(zip(pairs_new.india_sid, pairs_new.aus_sid)) & set(base_key))
        for mode, pairs in (("published pairs held fixed", base_pairs),
                            ("matching rebuilt from this subset", pairs_new)):
            sc = scores(pairs, pool)
            res = {}
            for lab in LAB:
                a, b = sc[lab]
                ok = np.isfinite(a) & np.isfinite(b)
                res[lab] = spearman_perm(a[ok], b[ok]) + (int(ok.sum()),)
            qs = dict(zip(LAB, bh([res[l][1] for l in LAB])))
            for lab in LAB:
                r, p, n = res[lab]
                rows.append(dict(subset=name, mode=mode, mineral=lab, n=n,
                                 pairs_shared_with_published=(shared if "rebuilt" in mode else 20),
                                 rho=round(r, 3), perm_p=round(p, 4), q=round(float(qs[lab]), 4),
                                 transfers="Yes" if qs[lab] < 0.05 else "No"))
        print("  %-44s shared=%2d  fixed-pairs monazite %.3f  rebuilt monazite %.3f"
              % (name[:44], shared,
                 [x for x in rows if x["subset"] == name
                  and x["mode"].startswith("published") and x["mineral"] == LAB[0]][0]["rho"],
                 [x for x in rows if x["subset"] == name
                  and x["mode"].startswith("matching") and x["mineral"] == LAB[0]][0]["rho"]),
              flush=True)
    O = pd.DataFrame(rows)
    O.to_csv(os.path.join(RES, "grain_size_fractions.csv"), index=False)
    pd.set_option("display.width", 260)
    print("\n" + O.to_string(index=False))
    print("\nwrote results/grain_size_fractions.csv")


def centroid_and_area(feat):
    """Area-weighted centroid of a polygon ring, and its stored area."""
    ring = feat["geometry"]["coordinates"][0]
    x = np.array([p[0] for p in ring], float)
    y = np.array([p[1] for p in ring], float)
    lat0 = np.radians(y.mean())
    X, Y = x * np.cos(lat0), y
    a = X[:-1] * Y[1:] - X[1:] * Y[:-1]
    A = a.sum() / 2.0
    if abs(A) < 1e-12:
        return float(x.mean()), float(y.mean())
    cx = ((X[:-1] + X[1:]) * a).sum() / (6 * A) / np.cos(lat0)
    cy = ((Y[:-1] + Y[1:]) * a).sum() / (6 * A)
    return float(cx), float(cy)


def haversine(lo1, la1, lo2, la2):
    R = 6371.0
    p1, p2 = np.radians(la1), np.radians(la2)
    dp, dl = p2 - p1, np.radians(lo2 - lo1)
    h = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(h))


def run_morphometry():
    pairs, pool = build(ngsa())
    sc = scores(pairs, pool)
    geo = {}
    for side, f in (("india", "pair_catchments_india.geojson"),
                    ("australia", "pair_catchments_australia.geojson")):
        g = json.load(open(os.path.join(RES, f)))
        for ft in g["features"]:
            p = ft["properties"]
            cx, cy = centroid_and_area(ft)
            geo[(side, p["sid"])] = dict(area_km2=float(p["area_km2"]),
                                         dist_km=float(haversine(float(p["lon"]), float(p["lat"]),
                                                                 cx, cy)))
    miss = [k for k in zip(["india"] * 20, pairs.india_sid) if k not in geo] + \
           [k for k in zip(["australia"] * 20, pairs.aus_sid) if k not in geo]
    print("catchment records missing for %d of 40 sites" % len(miss), flush=True)
    V = {}
    for side, sids in (("india", pairs.india_sid), ("australia", pairs.aus_sid)):
        for var in ("area_km2", "dist_km"):
            V[(side, var)] = np.array([geo[(side, s)][var] if (side, s) in geo else np.nan
                                       for s in sids], float)
    rows = []
    for side, idx in (("india", 0), ("australia", 1)):
        for var in ("area_km2", "dist_km"):
            v = V[(side, var)]
            ps, rs, ns = [], {}, {}
            for lab in LAB:
                s_ = sc[lab][idx]
                ok = np.isfinite(s_) & np.isfinite(v)
                r, p = spearman_perm(s_[ok], v[ok])
                ps.append(p)
                rs[lab], ns[lab] = r, int(ok.sum())
            qs = bh(ps)
            for j, lab in enumerate(LAB):
                rows.append(dict(side=side, variable=var, mineral=lab, n=ns[lab],
                                 rho=round(rs[lab], 3), perm_p=round(ps[j], 4),
                                 q=round(float(qs[j]), 4),
                                 significant="Yes" if qs[j] < 0.05 else "No"))
            print("  %-10s %-9s  %s" % (side, var,
                                        "  ".join("%s %.2f" % (l.split(" ")[0][:4], rs[l])
                                                  for l in LAB)), flush=True)
    M = pd.DataFrame(rows)
    M.to_csv(os.path.join(RES, "catchment_morphometry.csv"), index=False)
    A = pd.DataFrame([dict(side=s, sid=k[1], **geo[k]) for s, k in
                      [("india", ("india", x)) for x in pairs.india_sid] +
                      [("australia", ("australia", x)) for x in pairs.aus_sid] if k in geo])
    print("\ncatchment area km2: india median %.0f (%.0f to %.0f), australia median %.0f (%.0f to %.0f)"
          % (A[A.side == "india"].area_km2.median(), A[A.side == "india"].area_km2.min(),
             A[A.side == "india"].area_km2.max(), A[A.side == "australia"].area_km2.median(),
             A[A.side == "australia"].area_km2.min(), A[A.side == "australia"].area_km2.max()))
    print("site to centroid km: india median %.0f, australia median %.0f"
          % (A[A.side == "india"].dist_km.median(), A[A.side == "australia"].dist_km.median()))
    pd.set_option("display.width", 260)
    print("\n" + M.to_string(index=False))
    print("\nwrote results/catchment_morphometry.csv")


def run_confound():
    """Does the transfer survive once Australian catchment area is partialled out?

    Australian catchment area correlates strongly with the Australian fingerprint scores, so the
    India-Australia correlation could in principle be carried by catchment size rather than by
    provenance. Both score vectors are rank-transformed, the rank of Australian catchment area is
    regressed out of each by least squares, and the residuals are correlated. The permutation test
    shuffles one residual vector, which is the correct null for a partial correlation.
    """
    pairs, pool = build(ngsa())
    sc = scores(pairs, pool)
    geo = {}
    for side, f in (("india", "pair_catchments_india.geojson"),
                    ("australia", "pair_catchments_australia.geojson")):
        for ft in json.load(open(os.path.join(RES, f)))["features"]:
            geo[(side, ft["properties"]["sid"])] = float(ft["properties"]["area_km2"])
    area_au = np.array([geo[("australia", s_)] for s_ in pairs.aus_sid], float)
    area_in = np.array([geo[("india", s_)] for s_ in pairs.india_sid], float)
    rows = []
    for cname, cov in (("Australian catchment area", area_au),
                       ("Indian catchment area", area_in),
                       ("both catchment areas", np.column_stack([area_au, area_in]))):
        C = avg_rank(np.atleast_2d(cov.T if cov.ndim > 1 else cov))
        C = C.T if C.shape[0] > 1 else C.reshape(-1, 1)
        X = np.column_stack([np.ones(len(pairs)), C])
        ps, rs = [], {}
        for lab in LAB:
            a, b = sc[lab]
            ok = np.isfinite(a) & np.isfinite(b)
            ra, rb = avg_rank(a[ok])[0], avg_rank(b[ok])[0]
            Xo = X[ok]
            resa = ra - Xo @ np.linalg.lstsq(Xo, ra, rcond=None)[0]
            resb = rb - Xo @ np.linalg.lstsq(Xo, rb, rcond=None)[0]
            r = float(corr_rows(resa[None, :], resb)[0])
            n = len(resb)
            if n not in _BANK:
                _BANK[n] = rng.random((NPERM, n)).argsort(axis=1)
            null = np.abs(corr_rows(resb[_BANK[n]], resa))
            p = (np.sum(null >= abs(r)) + 1) / (NPERM + 1)
            ps.append(p)
            rs[lab] = (r, p, n)
        qs = bh(ps)
        for j, lab in enumerate(LAB):
            r, p, n = rs[lab]
            rows.append(dict(controlling_for=cname, mineral=lab, n=n, partial_rho=round(r, 3),
                             perm_p=round(p, 4), q=round(float(qs[j]), 4),
                             transfers="Yes" if qs[j] < 0.05 else "No"))
        print("  controlling for %-26s monazite %.3f (q %.4f)  xenotime %.3f"
              % (cname, rs[LAB[0]][0], qs[0], rs[LAB[1]][0]), flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(os.path.join(RES, "catchment_area_partial.csv"), index=False)
    pd.set_option("display.width", 260)
    print("\n" + D.to_string(index=False))
    print("\nwrote results/catchment_area_partial.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True, choices=["fractions", "morphometry", "confound"])
    A = ap.parse_args()
    if A.part == "fractions":
        run_fractions()
    elif A.part == "morphometry":
        run_morphometry()
    else:
        run_confound()
