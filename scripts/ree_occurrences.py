"""
Which documented rare-earth occurrences lie INSIDE the twenty analogue catchments?

This supplies the rare-earth mineralisation row of the four-domain comparison table.

Motivation: naming a deposit that sits in the analogue REGION but outside the catchments the
analogue population actually drains is not evidence about the analogue. Mount Weld, for example,
lies at 122.5 E in the Eastern Goldfields Superterrane, which contributes 0 per cent of the
Archaean analogue catchment area. This script replaces recall with a containment test.

Source: Geological Survey of Western Australia MINEDEX database (CC-BY-4.0),
resolved by paths.MINEDEX. Sites are flagged as rare-earth by the
string "RARE EARTH" appearing in Commodities or TargetCommodityGroups.

Catchments: results/pair_catchments_australia.geojson, the same polygons used by
catchment_geology.py, tested with matplotlib.path containment.

Outputs
  results/ree_occurrences.csv          one row per (site, containing catchment)
  results/ree_occurrence_summary.csv   per-domain counts for the table cell
Run
  python scripts/ree_occurrences.py

Author: Bhavik Harish Lodhia, Curtin University
"""
import csv, json, os
import numpy as np
import paths
from aravalli_wa.geometry import geom_to_path

csv.field_size_limit(10 ** 7)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJ = os.path.dirname(os.path.dirname(HERE))
RES = os.path.join(HERE, "results")
MINEDEX = paths.MINEDEX

# columns carried through to the output
KEEP = ["SiteCode", "ShortTitle", "Title", "Type", "SubType", "Stage",
        "Commodities", "TargetCommodityGroups", "MineralizationStyle",
        "Latitude", "Longitude"]


def is_ree(rec):
    s = ((rec.get("Commodities") or "") + " " +
         (rec.get("TargetCommodityGroups") or "")).upper()
    return "RARE EARTH" in s


def main():
    feats = json.load(open(os.path.join(RES, "pair_catchments_australia.geojson")))["features"]
    print("catchments: %d" % len(feats))

    sites, no_position = [], 0
    with open(MINEDEX, encoding="utf-8-sig", errors="replace") as fh:
        for rec in csv.DictReader(fh):
            if not is_ree(rec):
                continue
            try:
                lat = float(rec["Latitude"]); lon = float(rec["Longitude"])
            except (TypeError, ValueError):
                # A MINEDEX record with no usable position. Counted, not dropped in silence.
                no_position += 1
                continue
            sites.append((lon, lat, rec))
    print("MINEDEX rare-earth sites with coordinates: %d" % len(sites))
    if no_position:
        print("  rare-earth records skipped for a missing or non-numeric position: %d" % no_position)

    # de-duplicate on (ShortTitle, rounded position): MINEDEX repeats group/infrastructure rows
    seen, uniq = set(), []
    for lon, lat, rec in sites:
        k = ((rec["ShortTitle"] or "").strip().upper(), round(lon, 3), round(lat, 3))
        if k in seen:
            continue
        seen.add(k)
        uniq.append((lon, lat, rec))
    print("after de-duplication: %d" % len(uniq))

    pts = np.array([[lon, lat] for lon, lat, _ in uniq])

    rows = []
    for f in feats:
        p = f["properties"]
        cpath, cbox = geom_to_path(f["geometry"])
        if cpath is None:
            continue
        x0, y0, x1, y1 = cbox
        near = np.where((pts[:, 0] >= x0) & (pts[:, 0] <= x1) &
                        (pts[:, 1] >= y0) & (pts[:, 1] <= y1))[0]
        if near.size == 0:
            continue
        ins = cpath.contains_points(pts[near])
        for j in near[ins]:
            rec = uniq[j][2]
            row = {k: rec.get(k, "") for k in KEEP}
            row.update(pair_no=p["pair_no"], domain=p["domain"], sid=p["sid"],
                       catchment_km2=p["area_km2"])
            rows.append(row)

    out = os.path.join(RES, "ree_occurrences.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["pair_no", "domain", "sid", "catchment_km2"] + KEEP)
        w.writeheader()
        w.writerows(rows)
    print("wrote %s (%d rows)" % (out, len(rows)))

    # per-domain summary
    summ = {}
    for f in feats:
        d = f["properties"]["domain"]
        summ.setdefault(d, {"catchments": 0, "catchments_with_ree": 0, "sites": 0})
        summ[d]["catchments"] += 1
    hit_cat = {}
    for r in rows:
        summ[r["domain"]]["sites"] += 1
        hit_cat.setdefault(r["domain"], set()).add(r["sid"])
    for d in summ:
        summ[d]["catchments_with_ree"] = len(hit_cat.get(d, ()))

    outs = os.path.join(RES, "ree_occurrence_summary.csv")
    with open(outs, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["domain", "n_catchments", "n_catchments_with_ree_site", "n_ree_sites_inside"])
        for d, v in sorted(summ.items()):
            w.writerow([d, v["catchments"], v["catchments_with_ree"], v["sites"]])
            print("  %-18s %d catchments, %d contain a rare-earth site, %d sites inside"
                  % (d, v["catchments"], v["catchments_with_ree"], v["sites"]))
    print("wrote %s" % outs)

    # where do the named deposits of interest actually sit?
    print("\nreference check, position of named sites:")
    for lon, lat, rec in uniq:
        t = (rec["ShortTitle"] or "").upper()
        if "WELD" in t or "YANGIBANA" in t or "GIFFORD" in t or "BALD HILL" in t:
            inside = [f["properties"]["pair_no"] for f in feats
                      if (lambda cp: cp is not None and cp.contains_point((lon, lat)))(
                          geom_to_path(f["geometry"])[0])]
            print("  %-34s %8.3f %8.3f  inside pairs: %s"
                  % (rec["ShortTitle"][:34], lon, lat, inside or "none"))


if __name__ == "__main__":
    main()
