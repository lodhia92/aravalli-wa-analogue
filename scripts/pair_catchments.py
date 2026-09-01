"""Build upstream catchment polygons for the 20 robust analogue-pair samples (one side per run).
Traces HydroBASINS L12 upstream (NEXT_DOWN) from each sample's basin and unions the basin polygons.

    python pair_catchments.py --side india
    python pair_catchments.py --side australia

Output: results/pair_catchments_india.geojson, results/pair_catchments_australia.geojson

Author: Bhavik Harish Lodhia, Curtin University
"""

import argparse
import json
import os

import pandas as pd
import shapefile
from shapely.geometry import Point, mapping
from shapely.geometry import shape as shp
from shapely.ops import unary_union
from shapely.strtree import STRtree

import paths

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(HERE)
RES = f"{HERE}/results"


def load_basins(path, region):
    r = shapefile.Reader(path)
    fl = [f[0] for f in r.fields[1:]]
    # This loader returns basin topology only. SUB_AREA is located but not read; the
    # loader in drainage_containment.py is the one that sums it to weight catchment area.
    iId, iND, iSub = fl.index("HYBAS_ID"), fl.index("NEXT_DOWN"), fl.index("SUB_AREA")  # noqa: F841
    info, geoms = {}, {}
    x0, y0, x1, y1 = region
    for sr in r.iterShapeRecords():
        b = sr.shape.bbox
        if b[2] < x0 or b[0] > x1 or b[3] < y0 or b[1] > y1:
            continue
        rec = sr.record
        hid = int(rec[iId])
        info[hid] = int(rec[iND])
        geoms[hid] = shp(sr.shape.__geo_interface__).buffer(0)
    return info, geoms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--side", required=True)
    a = ap.parse_args()
    rob = pd.read_csv(f"{RES}/analogue_pairs.csv")
    # pair numbering per domain by ascending distance
    rob["pair_no"] = 0
    for dom in rob.domain.unique():
        idx = rob[rob.domain == dom].sort_values("dist").index
        for n, i in enumerate(idx, 1):
            rob.loc[i, "pair_no"] = n
    if a.side == "india":
        basins = paths.HYBAS_AS
        region = (72, 23, 77, 28)
        pts = []
        for _, r in rob.iterrows():
            lat, lon = float(r.india_sid.split("_")[1]), float(r.india_sid.split("_")[2])
            pts.append((r.domain, r.pair_no, r.india_sid, lon, lat))
        out = f"{RES}/pair_catchments_india.geojson"
    else:
        basins = paths.HYBAS_AU
        region = (113, -34, 130, -14)
        ng = pd.read_csv(f"{RES}/ngsa_points.csv")
        ng["SITEID"] = pd.to_numeric(ng.SITEID, errors="coerce")
        ngi = ng.set_index("SITEID")
        pts = []
        for _, r in rob.iterrows():
            sid = float(r.aus_sid.split("_")[-1])
            if sid not in ngi.index:
                continue
            row = ngi.loc[sid]
            pts.append((r.domain, r.pair_no, r.aus_sid, float(row.LON), float(row.LAT)))
        out = f"{RES}/pair_catchments_australia.geojson"
    info, geoms = load_basins(basins, region)
    print("basins:", len(info))
    up = {}
    for hid, nd in info.items():
        up.setdefault(nd, []).append(hid)

    def upstream(start):
        seen, st = {start}, [start]
        while st:
            c = st.pop()
            for u in up.get(c, []):
                if u not in seen:
                    seen.add(u)
                    st.append(u)
        return seen

    ids = list(geoms)
    tree = STRtree([geoms[i] for i in ids])
    feats = []
    for dom, no, sid, lon, lat in pts:
        pt = Point(lon, lat)
        base = None
        for j in tree.query(pt):
            if geoms[ids[j]].contains(pt):
                base = ids[j]
                break
        if base is None:
            print("no basin for", sid)
            continue
        cat = unary_union([geoms[h] for h in upstream(base) if h in geoms]).buffer(0)
        feats.append(
            {
                "type": "Feature",
                "properties": {
                    "domain": dom,
                    "pair_no": int(no),
                    "sid": sid,
                    "lon": lon,
                    "lat": lat,
                    "area_km2": round(cat.area * 12321, 0),
                },
                "geometry": mapping(cat),
            }
        )
    json.dump({"type": "FeatureCollection", "features": feats}, open(out, "w"))
    print(f"wrote {out} ({len(feats)} catchments)")


if __name__ == "__main__":
    main()
