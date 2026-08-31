"""
What geology does each analogue catchment actually drain?

For every one of the twenty analogue pair catchments (Australia and India) this computes the
AREA FRACTION of the catchment underlain by each mapped tectonic unit, rather than relying on
the tectonic unit beneath the sample site or on a single containment percentage.

Area is measured by point counting on a regular lon/lat grid, each point weighted by cos(lat)
so the fractions are true area fractions. Australian units come from the GSWA 1:500 000
tectonic units layer. Grid step is chosen per catchment to give of order 30 000 interior
points.

Output: results/catchment_geology.csv  (one row per catchment per unit)
Run:    python scripts/catchment_geology.py --start 1 --end 20

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import argparse, json, os
import numpy as np
import pandas as pd
import shapefile
import paths
from aravalli_wa.geometry import geom_to_path, shape_to_path

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(os.path.dirname(HERE))
RES  = os.path.join(HERE, "results")
SHP  = paths.TECTONIC
KEEP = ["TECTNAME", "OROGEN", "PROVINCE", "CRATON", "LITHOLOGY", "TECTSETTIN",
        "ERA_FROM", "ERA_TO", "MAX_AGE_MA", "MIN_AGE_MA"]
TARGET_PTS = 30000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--end", type=int, default=20)
    a = ap.parse_args()

    cat = json.load(open(os.path.join(RES, "pair_catchments_australia.geojson")))["features"]
    cat = [f for f in cat if a.start <= f["properties"]["pair_no"] <= a.end]
    print("catchments this run: %d (pairs %d-%d)" % (len(cat), a.start, a.end))

    r = shapefile.Reader(SHP)
    fields = [f[0] for f in r.fields[1:]]
    idx = {k: fields.index(k) for k in KEEP if k in fields}
    area_i = fields.index("SHAPE_STAr")
    n = len(r)
    bb = np.zeros((n, 4))
    for i, s in enumerate(r.iterShapes()):
        bb[i] = s.bbox

    rows = []
    for f in cat:
        p = f["properties"]
        cpath, cbox = geom_to_path(f["geometry"])
        if cpath is None:
            continue
        x0, y0, x1, y1 = cbox
        step = np.sqrt(max((x1 - x0) * (y1 - y0), 1e-9) / TARGET_PTS)
        gx = np.arange(x0 + step / 2, x1, step)
        gy = np.arange(y0 + step / 2, y1, step)
        GX, GY = np.meshgrid(gx, gy)
        pts = np.column_stack([GX.ravel(), GY.ravel()])
        inside = cpath.contains_points(pts)
        pts = pts[inside]
        if len(pts) == 0:
            continue
        w = np.cos(np.radians(pts[:, 1]))
        total = w.sum()

        # smallest containing tectonic unit for each interior point
        best_area = np.full(len(pts), np.inf)
        best_rec = np.full(len(pts), -1, dtype=int)
        cand = np.where((bb[:, 0] <= x1) & (bb[:, 2] >= x0) & (bb[:, 1] <= y1) & (bb[:, 3] >= y0))[0]
        for i in cand:
            sub = np.where((pts[:, 0] >= bb[i, 0]) & (pts[:, 0] <= bb[i, 2]) &
                           (pts[:, 1] >= bb[i, 1]) & (pts[:, 1] <= bb[i, 3]))[0]
            if sub.size == 0:
                continue
            up = shape_to_path(r.shape(i))
            if up is None:
                continue
            ins = up.contains_points(pts[sub])
            if not ins.any():
                continue
            ar = float(r.record(i)[area_i])
            sel = sub[ins]
            better = ar < best_area[sel]
            best_area[sel[better]] = ar
            best_rec[sel[better]] = i

        acc = {}
        for i in np.unique(best_rec):
            m = best_rec == i
            frac = w[m].sum() / total
            if i < 0:
                acc[("unmapped", "", "", "", "", "", "", "", np.nan, np.nan)] = frac
            else:
                rec = r.record(i)
                key = tuple(rec[idx[k]] if k in idx else "" for k in KEEP)
                acc[key] = acc.get(key, 0.0) + frac
        for key, frac in sorted(acc.items(), key=lambda kv: -kv[1]):
            row = dict(zip(KEEP, key))
            row.update(pair_no=p["pair_no"], domain=p["domain"], sid=p["sid"],
                       area_km2=p["area_km2"], n_pts=len(pts), area_frac=round(frac, 4))
            rows.append(row)
        top = sorted(acc.items(), key=lambda kv: -kv[1])[:3]
        print("pair %2d  %-17s %8.0f km2  n=%5d  " % (p["pair_no"], p["domain"], p["area_km2"], len(pts))
              + "; ".join("%s %.0f%%" % (k[0], v * 100) for k, v in top))

    out = os.path.join(RES, "catchment_geology_%02d_%02d.csv" % (a.start, a.end))
    pd.DataFrame(rows)[["pair_no", "domain", "sid", "area_km2", "n_pts", "area_frac"] + KEEP] \
        .to_csv(out, index=False)
    print("wrote", out)


if __name__ == "__main__":
    main()
