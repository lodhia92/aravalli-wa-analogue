"""Does lithological mixing within the analogue catchments affect the transfer?

Stream sediments integrate the whole upstream drainage, so a catchment labelled by its dominant
basement domain may carry substantial contributions from other lithologies. The lithological
proportions within each selected catchment therefore matter more than the percentage of the
catchment lying inside a basement polygon.

The proportions themselves are computed for the Australian catchments by catchment_geology.py: the area fraction of every mapped GSWA tectonic unit within each of the
twenty pair catchments, point-counted on a cos-latitude-weighted grid. This script turns those
proportions into per-catchment measures of mixing and then asks whether the mixing matters.

Measures, per catchment:
    dominant_frac  area fraction of the single largest mapped unit; 1.0 means one unit only
    n_units        number of mapped units present
    shannon        Shannon diversity of the unit area fractions, in nats
    basin_frac     area fraction whose tectonic setting is basin rather than craton or orogen,
                   which is the specific dilution of concern, sedimentary cover
                   contributing to a catchment classified by its basement

Two tests, both following grain_size_morphometry.py and containment_sensitivity.py so that the treatment of a
candidate confound is the same throughout the study:
    --part scores    each measure against each fingerprint score over the selected pair members
    --part confound  each measure regressed out of both score vectors, transfer recomputed on the
                     residuals, permutation test on the residuals

A third check is reported alongside: each measure against catchment area. Australian catchment area
already correlates with the fingerprint scores at 0.73 (grain_size_morphometry.py), so a mixing measure that
merely tracks area is not new information and is said to be so.

Scope limit, stated rather than worked around: this covers the AUSTRALIAN catchments only. No
vector geological map of the Aravalli is held in this study. The Indian unit assignments come from
a georeferenced published map figure carrying about ten kilometres of positional uncertainty with a
partly decoded legend, which is not a basis for area fractions quoted to two decimals.

Inputs : results/catchment_geology.csv
         results/pair_catchments_australia.geojson   (areas, for the tracking check)
         the scores and pair selection from containment_sensitivity.build_D at the published thresholds
Outputs: results/lithological_mixing.csv
         results/lithological_mixing_scores.csv
         results/lithological_mixing_partial.csv

    python lithological_mixing.py --part scores
    python lithological_mixing.py --part confound

Author: Bhavik Harish Lodhia, Curtin University
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matching_sensitivity as core  # matching, rank statistics, BH correction
from aravalli_wa.stats import avg_rank, bh, corr_rows
import containment_sensitivity as cont   # loader with the published thresholds

RES = os.path.join(HERE, "results")
DOMS, LAB = core.DOMS, core.LAB
MEASURES = ["dominant_frac", "n_units", "shannon", "basin_frac"]
NICE = {"dominant_frac": "dominant unit fraction", "n_units": "number of mapped units",
        "shannon": "Shannon diversity of unit fractions", "basin_frac": "basin-hosted fraction"}


def mixing_table():
    """Per-catchment mixing measures from the A1 area fractions."""
    g = pd.read_csv(os.path.join(RES, "catchment_geology.csv"))
    tot = g.groupby("sid").area_frac.sum()
    assert np.allclose(tot.values, 1.0, atol=2e-3), "area fractions do not sum to one"
    rows = []
    for sid, sub in g.groupby("sid"):
        p = sub.area_frac.values
        p = p / p.sum()
        basin = sub.loc[sub.TECTSETTIN.astype(str).str.lower() == "basin", "area_frac"].sum()
        rows.append(dict(
            sid=sid, domain=sub.domain.iloc[0], pair_no=int(sub.pair_no.iloc[0]),
            n_units=len(p),
            dominant_frac=round(float(p.max()), 4),
            shannon=round(float(-(p * np.log(p)).sum()), 4),
            basin_frac=round(float(basin), 4),
            dominant_unit=sub.loc[sub.area_frac.idxmax(), "TECTNAME"]))
    M = pd.DataFrame(rows).sort_values(["domain", "pair_no"]).reset_index(drop=True)
    M.to_csv(os.path.join(RES, "lithological_mixing.csv"), index=False)
    return M


def selected(M):
    """Align the mixing measures to the twenty selected Australian pair members."""
    D = cont.build_D(cont.PUB_IN, cont.PUB_AU)
    sel = core.select(D, "published")
    asid, ipos_all, apos_all = [], [], []
    for dom in DOMS:
        ipos, apos = sel[dom]
        asid.extend(D[dom]["asid"][apos])
    lut = M.set_index("sid")
    missing = [s for s in asid if s not in lut.index]
    assert not missing, f"no geology for {missing}"
    X = lut.loc[asid, MEASURES].astype(float)
    return D, sel, np.array(asid), X


def areas(asid):
    """Australian catchment areas, for the check that a measure is not just area again."""
    geo = {}
    for ft in json.load(open(os.path.join(RES, "pair_catchments_australia.geojson")))["features"]:
        geo[ft["properties"]["sid"]] = float(ft["properties"]["area_km2"])
    return np.array([geo[s] for s in asid], float)


def scorevecs(D, sel, lab):
    """Indian and Australian score vectors over the selected pairs, in pair order."""
    a, b = [], []
    for dom in DOMS:
        ipos, apos = sel[dom]
        a.append(D[dom]["SI"][lab][ipos])
        b.append(D[dom]["SA"][lab][apos])
    return np.concatenate(a), np.concatenate(b)


def run_scores():
    M = mixing_table()
    D, sel, asid, X = selected(M)
    print("mixing measures over the twenty Australian pair catchments:")
    print(X.describe().loc[["min", "50%", "max"]].round(3).to_string())
    ar = areas(asid)

    rows = []
    for m in MEASURES:
        cov = X[m].values
        # the score test, one BH family over the five fingerprints
        ps, keep = [], {}
        for lab in LAB:
            _, b = scorevecs(D, sel, lab)
            ok = np.isfinite(b) & np.isfinite(cov)
            r, p = core.spearman_perm(b[ok], cov[ok])
            keep[lab] = (r, p, int(ok.sum()))
            ps.append(p)
        qs = bh(ps)
        # is this measure just catchment area again?
        ra, _ = core.spearman_perm(cov, ar)
        for j, lab in enumerate(LAB):
            r, p, n = keep[lab]
            rows.append(dict(measure=NICE[m], mineral=lab, n=n,
                             rho_australian_score_vs_measure=round(r, 3),
                             perm_p=round(p, 4), q=round(float(qs[j]), 4),
                             significant="Yes" if qs[j] < 0.05 else "No",
                             rho_measure_vs_catchment_area=round(ra, 3)))
        print("  %-34s monazite %+.3f (q %.4f)   this measure vs catchment area %+.3f"
              % (NICE[m], keep[LAB[0]][0], qs[0], ra), flush=True)
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(RES, "lithological_mixing_scores.csv"), index=False)
    print("\nwrote results/lithological_mixing.csv and results/lithological_mixing_scores.csv")


def run_confound():
    M = mixing_table()
    D, sel, asid, X = selected(M)
    ar = areas(asid)[:, None]
    controls = [(NICE[m], X[m].values[:, None]) for m in MEASURES]
    controls.append(("all four measures together", X[MEASURES].values))
    # Does mixing add anything beyond catchment size? The mixing measures correlate 0.70 to 0.92
    # with catchment area, and grain_size_morphometry.py already partialled area out alone, leaving 0.629.
    # Anything close to 0.629 below means the mixing measure carries no information of its own.
    controls.append(("Australian catchment area alone (morphometry anchor)", ar))
    for m in MEASURES:
        controls.append(("catchment area plus " + NICE[m],
                         np.column_stack([ar[:, 0], X[m].values])))
    controls.append(("catchment area plus all four measures",
                     np.column_stack([ar[:, 0], X[MEASURES].values])))
    rows = []
    for cname, cov in controls:
        C = np.column_stack([avg_rank(cov[:, j])[0] for j in range(cov.shape[1])])
        ps, keep = [], {}
        for lab in LAB:
            a, b = scorevecs(D, sel, lab)
            mask = np.isfinite(a) & np.isfinite(b)
            ra, rb = avg_rank(a[mask])[0], avg_rank(b[mask])[0]
            Xd = np.column_stack([np.ones(int(mask.sum())), C[mask]])
            resa = ra - Xd @ np.linalg.lstsq(Xd, ra, rcond=None)[0]
            resb = rb - Xd @ np.linalg.lstsq(Xd, rb, rcond=None)[0]
            r = float(corr_rows(resa[None, :], resb)[0])
            n = len(resb)
            null = np.abs(corr_rows(resb[core.bank(n)], resa))
            p = (np.sum(null >= abs(r)) + 1) / (core.NPERM + 1)
            keep[lab] = (r, p, n)
            ps.append(p)
        qs = bh(ps)
        for j, lab in enumerate(LAB):
            r, p, n = keep[lab]
            rows.append(dict(controlling_for=cname, mineral=lab, n=n, partial_rho=round(r, 3),
                             perm_p=round(p, 4), q=round(float(qs[j]), 4),
                             transfers="Yes" if qs[j] < 0.05 else "No"))
        print("  controlling for %-34s monazite %.3f (q %.4f)  xenotime %.3f"
              % (cname, keep[LAB[0]][0], qs[0], keep[LAB[1]][0]), flush=True)
    C2 = pd.DataFrame(rows)
    C2.to_csv(os.path.join(RES, "lithological_mixing_partial.csv"), index=False)
    pd.set_option("display.width", 260)
    print("\n" + C2.to_string(index=False))
    print("\nwrote results/lithological_mixing_partial.csv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True, choices=["scores", "confound"])
    A = ap.parse_args()
    if A.part == "scores":
        run_scores()
    else:
        run_confound()
