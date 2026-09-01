"""
Resolve the location of every external dataset used by this workflow.

None of the source datasets are redistributed with this code. Download each one from the
provider listed in data/README.md and place it at the path given below, or set the
environment variable ARAVALLI_WA_DATA to a directory holding the same layout.

    data/
      ngsa/NGSA_data.csv
      hmma/HM_appendixA.XLSX
      ngcm/EarthBank_qcd_raw_data/<state>/<toposheet>/*.xlsx
      ngcm/ngcm_aravalli.csv                     (produced by prepare_ngcm.py)
      hydrosheds/hybas_as_lev12/hybas_as_india.shp
      hydrosheds/hybas_au_lev12/hybas_au_lev12_v1c.shp
      gswa/GEOLOGY_500k_Tectonics_GDA2020_SHP/ESRI/SHAPEFILES/500k_tectonicp.shp
      gswa/Minedex_GDA2020_CSV/Sites.csv

Author: Bhavik Harish Lodhia, Curtin University
"""

import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("ARAVALLI_WA_DATA", os.path.join(REPO, "data"))
RESULTS = os.path.join(REPO, "results")
DRAINAGE = os.path.join(RESULTS, "drainage")

NGSA = os.path.join(DATA, "ngsa", "NGSA_data.csv")
HMMA = os.path.join(DATA, "hmma", "HM_appendixA.XLSX")
NGCM_RAW = os.path.join(DATA, "ngcm", "EarthBank_qcd_raw_data")
NGCM_TABLE = os.path.join(DATA, "ngcm", "ngcm_aravalli.csv")
HYBAS_AS = os.path.join(DATA, "hydrosheds", "hybas_as_lev12", "hybas_as_india.shp")
HYBAS_AU = os.path.join(DATA, "hydrosheds", "hybas_au_lev12", "hybas_au_lev12_v1c.shp")
TECTONIC = os.path.join(
    DATA, "gswa", "GEOLOGY_500k_Tectonics_GDA2020_SHP", "ESRI", "SHAPEFILES", "500k_tectonicp"
)
MINEDEX = os.path.join(DATA, "gswa", "Minedex_GDA2020_CSV", "Sites.csv")

_REQUIRED = {
    "NGSA": NGSA,
    "HMMA": HMMA,
    "NGCM_RAW": NGCM_RAW,
    "NGCM_TABLE": NGCM_TABLE,
    "HYBAS_AS": HYBAS_AS,
    "HYBAS_AU": HYBAS_AU,
    "MINEDEX": MINEDEX,
}


def require(*names):
    """Fail early, with the provider named, rather than deep inside a read."""
    missing = [
        (n, _REQUIRED[n]) for n in names if n in _REQUIRED and not os.path.exists(_REQUIRED[n])
    ]
    if missing:
        lines = "\n".join(f"  {n}: expected at {p}" for n, p in missing)
        raise SystemExit("Missing input data. See data/README.md for sources.\n" + lines)
