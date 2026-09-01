"""
The weathering-robust provenance ratios used to match Indian and Australian catchments.

Author: Bhavik Harish Lodhia, Curtin University
"""

import numpy as np
import pandas as pd

from .constants import CI_CHONDRITE_PPM as CI


def ratios(df):
    """The six matching ratios, one column each, indexed as the input frame."""
    o = pd.DataFrame(index=df.index)
    o["Th/Sc"] = df.Th / df.Sc
    o["La/Sc"] = df.La / df.Sc
    o["Th/Co"] = df.Th / df.Co
    o["EuEu"] = (df.Eu / CI["Eu"]) / np.sqrt((df.Sm / CI["Sm"]) * (df.Gd / CI["Gd"]))
    o["La/Yb_n"] = (df.La / df.Yb) / (CI["La"] / CI["Yb"])
    o["Nb/Y"] = df.Nb / df.Y
    return o


def logr(df):
    """Log of the matching ratios; non-positive and non-finite ratios become NaN."""
    r = ratios(df).replace([np.inf, -np.inf], np.nan)
    return np.log(r.where(r > 0)).replace([np.inf, -np.inf], np.nan)
