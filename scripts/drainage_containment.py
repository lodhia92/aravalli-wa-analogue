"""
Drainage-containment sample selection using HydroBASINS level-12 topology directly.

For each sample: find its L12 basin, trace upstream via NEXT_DOWN, and compute the fraction of
the upstream catchment area (sum of SUB_AREA) that lies within a domain polygon. Samples whose
catchment is dominantly inside the domain are that domain's analogue set.

This uses only HydroBASINS (its own NEXT_DOWN routing) plus shapely/pyshp; it does NOT need
HydroRIVERS or geopandas (geopandas/GDAL do not install in the sandbox). HydroBASINS:
Lehner & Grill (2013), Hydrological Processes 27(15), 2171-2186, doi:10.1002/hyp.9740.

Usage:
  python drainage_containment.py --domain results/sandmata_complex.geojson \
      --samples <csv> --latcol LAT --loncol LON --idcol SITEID \
      --basins data/hydrosheds/hybas_as_lev12/hybas_as_lev12_v1c.shp \
      --region 73.5 24.0 76.0 27.0 --out results/drainage/sandmata_contained.csv

Output CSV has every sample inside the region with its catchment area and the percentage of
that catchment inside the domain (pct_in_domain); filter by threshold afterwards.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import argparse, json, os
import pandas as pd
import shapefile
from shapely.geometry import shape as shp, Point
from shapely.strtree import STRtree
from shapely.prepared import prep


def load_basins(path, region):
    """Read HydroBASINS L12 within region bbox. Returns dict id->(next_down,sub_area), id->geom."""
    r = shapefile.Reader(path)
    fl = [f[0] for f in r.fields[1:]]
    iId, iND, iSub = fl.index("HYBAS_ID"), fl.index("NEXT_DOWN"), fl.index("SUB_AREA")
    info, geoms = {}, {}
    x0, y0, x1, y1 = region
    for sr in r.iterShapeRecords():
        b = sr.shape.bbox
        if b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1:
            continue
        rec = sr.record
        hid = int(rec[iId])
        info[hid] = (int(rec[iND]), float(rec[iSub]))
        geoms[hid] = shp(sr.shape.__geo_interface__).buffer(0)
    return info, geoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--basins", required=True)
    ap.add_argument("--region", nargs=4, type=float, required=True, metavar=("X0", "Y0", "X1", "Y1"))
    ap.add_argument("--latcol", default="LAT")
    ap.add_argument("--loncol", default="LON")
    ap.add_argument("--idcol", default=None)
    ap.add_argument("--encoding", default="utf-8")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    dom = shp(json.load(open(args.domain))["features"][0]["geometry"]).buffer(0)
    pdom = prep(dom)
    info, geoms = load_basins(args.basins, args.region)
    print(f"basins in region: {len(info)}")

    # basins mostly within the domain
    in_dom = set()
    for hid, g in geoms.items():
        if pdom.intersects(g) and g.intersection(dom).area > 0.5 * g.area:
            in_dom.add(hid)
    print(f"basins mostly within domain: {len(in_dom)}")

    # upstream graph
    up = {}
    for hid, (nd, _) in info.items():
        up.setdefault(nd, []).append(hid)

    def upstream(start):
        seen, st = {start}, [start]
        while st:
            cur = st.pop()
            for u in up.get(cur, []):
                if u not in seen:
                    seen.add(u); st.append(u)
        return seen

    # spatial index of basins for point location
    ids = list(geoms.keys())
    tree = STRtree([geoms[i] for i in ids])

    df = pd.read_csv(args.samples, low_memory=False, encoding=args.encoding)
    df["_lat"] = pd.to_numeric(df[args.latcol], errors="coerce")
    df["_lon"] = pd.to_numeric(df[args.loncol], errors="coerce")
    df = df.dropna(subset=["_lat", "_lon"])
    x0, y0, x1, y1 = args.region
    df = df[(df._lon >= x0) & (df._lon <= x1) & (df._lat >= y0) & (df._lat <= y1)]
    print(f"samples in region: {len(df)}")

    rows = []
    for _, s in df.iterrows():
        pt = Point(s._lon, s._lat)
        base = None
        for j in tree.query(pt):
            if geoms[ids[j]].contains(pt):
                base = ids[j]; break
        if base is None:
            continue
        cat = upstream(base)
        tot = sum(info[h][1] for h in cat)
        ind = sum(info[h][1] for h in cat if h in in_dom)
        rec = {"lat": round(float(s._lat), 4), "lon": round(float(s._lon), 4),
               "catchment_km2": round(tot, 1), "pct_in_domain": round(100 * ind / tot, 1) if tot else 0.0}
        if args.idcol and args.idcol in df.columns:
            rec["id"] = s[args.idcol]
        rows.append(rec)
    out = pd.DataFrame(rows).sort_values("pct_in_domain", ascending=False)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False)
    for thr in (90, 75, 50, 25):
        print(f"  pct_in_domain >= {thr}: {int((out.pct_in_domain >= thr).sum())} samples")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
