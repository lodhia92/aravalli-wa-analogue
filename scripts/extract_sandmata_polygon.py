"""
Extract the Sandmata Complex polygon from Ghosh et al. (2026) Fig. 1a.

Why this exists: manual digitising of this figure kept distorting because the
figure is a generalised, Lambert-projected schematic and (in the cropped copy
first used) the E75 meridian at the very top was missing. This script instead
georeferences the FULL figure from its graticule labels (read by OCR) and
colour-segments the Sandmata unit, so the result is reproducible rather than
hand-traced.

Pipeline:
  1. Use the full-figure render digitising/ghosh2026_fig1a_FULL.png
     (pdftoppm -png -r 450 of Ghosh 2026 p.2, cropped to include the top edge).
  2. OCR the graticule labels (N24, N25, N26, E74, E75) to get their pixel
     positions; build an affine pixel->(lon,lat) transform. Latitude plane is
     fit from the three N labels; longitude from the two E labels plus an
     orthogonality constraint (parallels perpendicular to meridians).
  3. Colour-segment the warm tan Sandmata unit (R>G>B, moderate R-B, bright),
     restricted to the Fig.1a map area (inset and legend boxes excluded).
  4. Vectorise (union of grid cells), simplify, transform vertices to lon/lat.
  5. Write results/sandmata_complex.geojson (EPSG:4326, MultiPolygon).

Caveat: generalised source map -> ~10 km positional uncertainty. Adequate for
masking HydroBASINS level-12 catchments; not a survey-grade boundary. If a GSI
vector becomes available, replace the GeoJSON and rerun the drainage step.

Control points (pixel x, y in the FULL render; degree value), measured by OCR:
  N25 (190,1049)=25, N24 (190,1540)=24, N26 (1175,495)=26   [latitude]
  E74 (363,1761)=74, E75 (909,284)=75                         [longitude]
Re-measure these if the figure is re-rendered at a different size/crop.

Output: results/sandmata_complex.geojson

Author: Bhavik Harish Lodhia, Curtin University
"""
import os, json

import numpy as np

from PIL import Image, ImageFilter

from shapely.geometry import box as shbox, Polygon, MultiPolygon, Point

from shapely.ops import unary_union

def to_lonlat(x, y):
    return (d * x + e * y + f, a * x + b * y + c)

def main():
    global A, B, FULL, G, H, HERE, LAT, LON, OUT, R, W, _, a, area, b, boxes, c, coords, d, e, f, gj, i, j, ll, m, mask, mp, p, parts, poly, r, rgb, rings, s, tan, v, x, xs, y, ys
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    FULL = os.path.join(HERE, "digitising", "ghosh2026_fig1a_FULL.png")
    OUT = os.path.join(HERE, "results", "sandmata_complex.geojson")
    LAT = np.array([[190, 1049, 25.0], [190, 1540, 24.0], [1175, 495, 26.0]])
    LON = np.array([[363, 1761, 74.0], [909, 284, 75.0]])
    a, b, c = np.linalg.solve(np.c_[LAT[:, 0], LAT[:, 1], np.ones(3)], LAT[:, 2])
    A = np.array([[a, b, 0], [LON[0, 0], LON[0, 1], 1], [LON[1, 0], LON[1, 1], 1]])
    d, e, f = np.linalg.solve(A, np.array([0.0, 74.0, 75.0]))
    rgb = np.asarray(Image.open(FULL).convert("RGB")).astype(int)
    H, W, _ = rgb.shape
    R, G, B = rgb[:, :, 0], rgb[:, :, 1], rgb[:, :, 2]
    xs = np.arange(W)[None, :].repeat(H, 0)
    ys = np.arange(H)[:, None].repeat(W, 1)
    area = (xs > 125) & (xs < 1205) & (ys > 220) & (ys < 1825)
    area &= ~((xs < 560) & (ys < 660))      # India inset (top-left)
    area &= ~((xs > 800) & (ys > 1555))     # legend box (bottom)
    tan = (R > G) & (G >= B - 3) & ((R - B) >= 22) & ((R - B) <= 85) & (R >= 175) & (G >= 150)
    mask = tan & area
    m = Image.fromarray((mask * 255).astype("uint8"))
    m = m.filter(ImageFilter.MinFilter(3)).filter(ImageFilter.MaxFilter(7)).filter(ImageFilter.MinFilter(5))
    mask = np.asarray(m) > 127
    s = 6
    boxes = [shbox(i, j, i + s, j + s)
             for j in range(0, H - s, s) for i in range(0, W - s, s)
             if mask[j:j + s, i:i + s].mean() > 0.5]
    poly = unary_union(boxes).simplify(4)
    parts = list(poly.geoms) if poly.geom_type == "MultiPolygon" else [poly]
    parts = [p for p in parts if p.area > 1500]
    coords, rings = [], []
    for p in parts:
        ll = [list(map(float, to_lonlat(x, y))) for x, y in p.exterior.coords]
        coords.append([ll]); rings.append(ll)
    gj = {"type": "FeatureCollection",
          "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}},
          "features": [{"type": "Feature",
                        "properties": {"name": "Sandmata Complex",
                                       "source": "Auto colour-segmented + graticule-georeferenced from Ghosh et al. 2026 Fig 1a",
                                       "uncertainty_km": 10},
                        "geometry": {"type": "MultiPolygon", "coordinates": coords}}]}
    json.dump(gj, open(OUT, "w"))
    mp = MultiPolygon([Polygon(r) for r in rings]).buffer(0)
    print("wrote", OUT)
    print("parts:", len(parts), "| bbox lon/lat:", [round(v, 3) for v in mp.bounds])
    print("contains Bhilwara:", mp.contains(Point(74.636, 25.346)))

if __name__ == "__main__":
    main()
