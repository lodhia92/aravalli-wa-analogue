"""
Reading the national survey tables and selecting the drainage-defined sample pools.

Author: Bhavik Harish Lodhia, Curtin University
"""

import pandas as pd

from aravalli_wa.constants import AUS_THR_PCT, IN_THR_PCT


def findcol(header, el):
    """Index of an element's column in the NGSA header, preferring ICP-MS over XRF.

    The NGSA header carries the element, the analytical method, the unit and the lower limit of
    detection in one string, and an element may be reported by more than one method. Returns
    None when the element is not present, which the caller must handle.
    """
    for meth in ("ICP-MS", "XRF"):
        for i, c in enumerate(header):
            if c.strip().startswith(f"{el} {meth}"):
                return i
    return None


def locality_key(frame):
    """The rounded lat/lon key that joins the Indian survey table to the drainage tables.

    The national geochemical mapping table carries no sample identifier that the drainage
    tables share, so a sample is matched on its coordinates rounded to four decimal places,
    about ten metres. Both sides must round identically or the join silently returns nothing.
    """
    return frame.lat.round(4).astype(str) + "_" + frame.lon.round(4).astype(str)


def indian_pool(ngcm, contained, name, threshold_pct=IN_THR_PCT):
    """Indian samples whose upstream catchment lies mostly inside one basement domain.

    ngcm is the national geochemical mapping table carrying a `k` locality key; contained is
    the drainage-containment table for that domain. A sample is kept when at least
    threshold_pct of its upstream catchment area falls inside the domain polygon.

    Returns the matching rows with a `sid` of the form "<name>_<locality key>".
    """
    kept = contained[contained.pct_in_domain >= threshold_pct]
    pool = ngcm[ngcm.k.isin(set(locality_key(kept)))].copy()
    pool["sid"] = name + "_" + pool.k
    return pool


def australian_pool(ngsa, contained, label, threshold_pct=AUS_THR_PCT):
    """Australian samples whose upstream catchment lies mostly inside one tectonic domain.

    ngsa must be INDEXED BY SITEID, which is how every caller builds it, with
    groupby("SITEID").median(). A frame that keeps SITEID as a column matches nothing here
    and yields an empty pool without raising, so check the index before calling.

    contained is the drainage-containment table for that domain, whose `id` column holds the
    site identifiers. Returns the matching rows, index reset, with a `sid` of the form
    "<label>_<SITEID>".
    """
    kept = contained[contained.pct_in_domain >= threshold_pct]
    ids = set(pd.to_numeric(kept["id"], errors="coerce").dropna())
    pool = ngsa[ngsa.index.isin(ids)].copy()
    pool["sid"] = [f"{label}_{i}" for i in pool.index]
    return pool.reset_index()
