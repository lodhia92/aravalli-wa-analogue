# Input data

None of the datasets below are redistributed with this code. Each is obtained from its provider
under that provider's licence. Download them and place them at the paths shown, or set the
environment variable `ARAVALLI_WA_DATA` to a directory holding the same layout.

```
data/
  ngsa/NGSA_data.csv
  hmma/HM_appendixA.XLSX
  ngcm/EarthBank_qcd_raw_data/<state>/<toposheet>/*.xlsx
  ngcm/ngcm_aravalli.csv
  hydrosheds/hybas_as_lev12/hybas_as_india.shp
  hydrosheds/hybas_au_lev12/hybas_au_lev12_v1c.shp
  gswa/GEOLOGY_500k_Tectonics_GDA2020_SHP/ESRI/SHAPEFILES/500k_tectonicp.shp
  gswa/Minedex_GDA2020_CSV/Sites.csv
```

`scripts/paths.py` resolves every one of these. Running a script without its inputs fails
immediately with the missing path named, rather than part way through a read.

## Sources

**National Geochemical Survey of Australia (NGSA)** — stream-sediment geochemistry, Australia.
Geoscience Australia. Released under Creative Commons Attribution 3.0 Australia.
de Caritat, P. and Cooper, M. (2016) and the Geochemical Atlas of Australia dataset,
https://www.ga.gov.au (GEOCAT 70478). File: `NGSA_data.csv`.
The column headers carry the element, the analytical method, the unit and the lower limit of
detection. Values below detection are recorded with a less-than operator.

**Heavy Mineral Map of Australia (HMMA)** — automated heavy-mineral grain counts made on the
NGSA samples themselves. Geoscience Australia. de Caritat et al. (2022); Walker et al. (2024).
File: `HM_appendixA.XLSX`, sheet `HMMA Dataset v.1.0`.

**National Geochemical Mapping programme (NGCM)** — stream-sediment geochemistry, India.
Geological Survey of India. India (2014); Bagchi et al. (2021). Supplied as survey packages by
analytical method: an X-ray fluorescence package, an inductively coupled plasma mass spectrometry
package, and a package covering the remaining elements, with a sample metadata workbook giving
coordinates. Obtain from the Geological Survey of India; the terms of use are the Survey's.

**HydroSHEDS HydroBASINS level 12** — drainage sub-basins and their downstream topology, used to
build true upstream catchments. Lehner, B. and Grill, G. (2013), Hydrological Processes 27,
2171-2186, doi:10.1002/hyp.9740. https://www.hydrosheds.org.

**GSWA 1:500 000 tectonic units** — Western Australian basement domains, used to attribute
catchments to mapped tectonic units. Geological Survey of Western Australia; Cutten et al. (2022).

**GSWA MINEDEX** — Western Australian mineral occurrence database, used to test which documented
rare-earth occurrences lie inside the analogue catchments. Geological Survey of Western Australia,
Creative Commons Attribution 4.0.

## The Indian geochemistry table

`ngcm/ngcm_aravalli.csv` is a single table of the north-west Indian NGCM samples, one row per
sample, with `LAT`, `LON` and one column per element or oxide. It is produced from the survey
packages by `scripts/prepare_ngcm.py`, which merges the three analytical packages per sample and
attaches the coordinates from the metadata workbook. Spatial selection is not performed here: the
drainage containment step selects the samples that belong to each basement domain.

## Derived files carried in this repository

Four files under `results/` are inputs to the workflow rather than outputs of it: the pipeline
reads them and no script here writes them. They are derived products of this study rather than
third-party data, they total about 64 kB, and they are carried because a clone cannot run
without them.

`sandmata_complex.geojson` is the Sandmata Complex boundary, used as the `--domain` polygon for
the Indian drainage step and to count samples per domain. `scripts/extract_sandmata_polygon.py`
produces it by georeferencing Fig. 1a of Ghosh et al. (2026) from its graticule labels and
colour-segmenting the Sandmata unit. That figure belongs to its publisher and is not
redistributed here, so the script cannot be rerun from a clone alone and the polygon it writes
is supplied instead. The source map is a generalised schematic, so the boundary carries a
positional uncertainty of order 10 km. That is adequate for masking level-12 catchments and is
not a survey-grade boundary; a Geological Survey of India vector, if one becomes available,
should replace it and the drainage step rerun.

`archaean_mangalwar_domain.geojson` is the Archaean Mangalwar basement, the second Indian domain
polygon and the Yilgarn analogue. It is Mangalwar Complex unit 22, georeferenced from the same
published map, with the Sandmata polygon subtracted so that the two domains do not overlap. Both
polygons record their own provenance in the `source` property of their single feature.

`ngsa_points.csv` carries the site identifier and coordinates of the 1315 NGSA sites. It is read
by `pair_catchments.py` and is the `--samples` input to the Australian drainage step. It is
derived from `NGSA_data.csv`.

`ngsa_draining_hallscreek.csv` gives the containment of each NGSA catchment within the Halls
Creek Orogen and is read by `provenance_ratios.py`. It is a `drainage_containment.py` output for
a domain that is not otherwise carried in `results/drainage/`.

The four `results/drainage/*_contained.csv` files are committed for the same reason. They are
`drainage_containment.py` outputs, but the `--domain`, `--region` and `--samples` arguments used
to produce each of them are not recorded here, so they are supplied rather than regenerated.
Everything downstream of the drainage step runs from them.
