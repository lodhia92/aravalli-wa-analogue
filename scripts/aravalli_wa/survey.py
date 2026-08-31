"""
Reading the national survey tables.

Author: Bhavik Harish Lodhia, Curtin University
"""


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
