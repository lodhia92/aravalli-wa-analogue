"""
The Palaeoproterozoic Australian analogue pool combines the Halls Creek, Lamboo, Capricorn
and Gascoyne terranes. This script attributes every
drainage-selected Australian Palaeoproterozoic NGSA site, and every Australian member of
the twenty robust analogue pairs, to a named tectonic unit of the GSWA 1:500 000 tectonic
units layer, so the composition of the pool can be reported rather than asserted.

Method: point in polygon of the sample site against 500k_tectonicp.shp. Where a site falls
inside more than one polygon the smallest (most specific) polygon is used. Attributes kept
are the ones needed for the four-domain comparison table: unit name, orogen, province,
lithology, tectonic setting and age range.

Output: results/terrane_sites.csv
Run:    python scripts/terrane_attribution.py

Author: Bhavik Harish Lodhia, Curtin University
"""
import os
import numpy as np
import pandas as pd
import shapefile
import paths
from aravalli_wa.geometry import shape_to_path

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(os.path.dirname(HERE))
RES  = os.path.join(HERE, "results")
DR   = os.path.join(RES, "drainage")
SHP  = paths.TECTONIC
AUS_THR = 25.0   # same threshold as analogue_pairing.py

KEEP = ["TECTNAME", "OROGEN", "PROVINCE", "CRATON", "DOMAIN_", "LITHOLOGY",
        "TECTSETTIN", "TSETT_QUAL", "ERA_FROM", "ERA_TO", "MAX_AGE_MA", "MIN_AGE_MA"]


def main():
    wapp = pd.read_csv(os.path.join(DR, "wa_palaeoprot_contained.csv"))
    wapp = wapp[wapp.pct_in_domain >= AUS_THR].copy()
    wapp["pool"] = "WA_PP"
    yil = pd.read_csv(os.path.join(DR, "yilgarn_contained.csv"))
    yil = yil[yil.pct_in_domain >= AUS_THR].copy()
    yil["pool"] = "Y+N"
    sites = pd.concat([wapp, yil], ignore_index=True)
    sites["SITEID"] = pd.to_numeric(sites["id"], errors="coerce")

    # which sites are members of the twenty robust pairs
    pairs = pd.read_csv(os.path.join(RES, "analogue_pairs.csv"))
    pairs = pairs[pairs.mnn.astype(str).str.lower() == "true"].copy()
    pair_id = {}
    for n, row in enumerate(pairs.itertuples(), start=1):
        sid = float(str(row.aus_sid).split("_")[-1])
        pair_id[sid] = (n, row.domain, row.india_sid, row.dist)
    sites["pair_no"] = sites.SITEID.map(lambda s: pair_id.get(s, (None,))[0])
    sites["pair_domain"] = sites.SITEID.map(lambda s: pair_id.get(s, (None, None))[1] if s in pair_id else None)
    sites["india_sid"] = sites.SITEID.map(lambda s: pair_id[s][2] if s in pair_id else None)
    sites["match_dist"] = sites.SITEID.map(lambda s: pair_id[s][3] if s in pair_id else None)

    print("sites to attribute: %d (WA_PP %d, Y+N %d); of these %d are pair members"
          % (len(sites), (sites.pool == "WA_PP").sum(), (sites.pool == "Y+N").sum(),
             sites.pair_no.notna().sum()))
    r = shapefile.Reader(SHP)
    fields = [f[0] for f in r.fields[1:]]
    idx = {k: fields.index(k) for k in KEEP if k in fields}
    area_i = fields.index("SHAPE_STAr")
    n = len(r)
    bboxes = np.zeros((n, 4))
    for i, s in enumerate(r.iterShapes()):
        bboxes[i] = s.bbox
    print("tectonic polygons: %d" % n)

    lons = sites.lon.values
    lats = sites.lat.values
    best = [None] * len(sites)          # (area, record)
    for i in range(n):
        x0, y0, x1, y1 = bboxes[i]
        cand = np.where((lons >= x0) & (lons <= x1) & (lats >= y0) & (lats <= y1))[0]
        if cand.size == 0:
            continue
        p = shape_to_path(r.shape(i))
        if p is None:
            continue
        inside = p.contains_points(np.column_stack([lons[cand], lats[cand]]))
        if not inside.any():
            continue
        rec = r.record(i)
        a = float(rec[area_i])
        for j in cand[inside]:
            if best[j] is None or a < best[j][0]:
                best[j] = (a, rec)

    for k in KEEP:
        if k in idx:
            sites[k] = [b[1][idx[k]] if b else "" for b in best]
    sites["unit_area_deg2"] = [b[0] if b else np.nan for b in best]
    sites["attributed"] = [b is not None for b in best]

    out = os.path.join(RES, "terrane_sites.csv")
    cols = ["SITEID", "pool", "lat", "lon", "catchment_km2", "pct_in_domain",
            "pair_no", "pair_domain", "india_sid", "match_dist", "attributed",
            "unit_area_deg2"] + [k for k in KEEP if k in sites.columns]
    sites[cols].to_csv(out, index=False)
    print("wrote", out)
    print("unattributed (site falls outside every polygon): %d"
          % (~sites.attributed).sum())
    print()
    print("--- WA Palaeoproterozoic pool, by tectonic unit ---")
    print(sites[sites.pool == "WA_PP"].TECTNAME.value_counts().to_string())
    print()
    print("--- Australian members of the twenty pairs ---")
    pm = sites[sites.pair_no.notna()].sort_values("pair_no")
    print(pm[["pair_no", "pair_domain", "SITEID", "TECTNAME", "OROGEN", "PROVINCE"]].to_string(index=False))


if __name__ == "__main__":
    main()
