"""
Converting shapefile and GeoJSON geometries to matplotlib paths for point-in-polygon tests.

matplotlib.path is used only as a geometry primitive here; nothing in this package plots.

Author: Bhavik Harish Lodhia, Curtin University
"""

import numpy as np
from matplotlib.path import Path


def shape_to_path(shp):
    """Path covering every ring of a pyshp shape; None when it has no usable ring."""
    pts = np.asarray(shp.points, dtype=float)
    parts = list(shp.parts) + [len(pts)]
    verts, codes = [], []
    for i in range(len(parts) - 1):
        seg = pts[parts[i] : parts[i + 1]]
        if len(seg) < 3:
            continue
        verts.append(seg)
        c = np.full(len(seg), Path.LINETO, dtype=np.uint8)
        c[0] = Path.MOVETO
        codes.append(c)
    if not verts:
        return None
    return Path(np.vstack(verts), np.concatenate(codes))


def geom_to_path(geom):
    """Path and bounding box for a GeoJSON Polygon or MultiPolygon; (None, None) if empty."""
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    verts, codes = [], []
    for poly in polys:
        for ring in poly:
            seg = np.asarray(ring, dtype=float)[:, :2]
            if len(seg) < 3:
                continue
            verts.append(seg)
            c = np.full(len(seg), Path.LINETO, dtype=np.uint8)
            c[0] = Path.MOVETO
            codes.append(c)
    if not verts:
        return None, None
    v = np.vstack(verts)
    return Path(v, np.concatenate(codes)), (
        v[:, 0].min(),
        v[:, 1].min(),
        v[:, 0].max(),
        v[:, 1].max(),
    )
