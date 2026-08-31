"""
Self-test: confirm that this repository reproduces the published headline results.

Checks, in order:
  1. the twenty analogue pairs are the published set;
  2. the light rare-earth transfer correlation is 0.696 and the heavy rare-earth 0.506.

Rank correlations are deterministic, so this runs in seconds and does not need the
hundred-thousand-permutation probabilities. Run it after placing the input data and before
trusting any other output.

    python scripts/verify_reproduction.py

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matching_sensitivity as core
import paths

EXPECT = {"monazite (Ce,Nd,Pr)": 0.696, "xenotime (Dy)": 0.506}
TOL = 0.001


def main():
    paths.require("NGSA", "NGCM_TABLE")
    D = core.prepare()
    sel = core.select(D, "published")
    npairs = len(core.pairkeys(D, sel))
    print(f"analogue pairs: {npairs}")
    ok = npairs == 20
    if not ok:
        print("  FAIL expected 20 pairs")
    for lab in core.LAB:
        a, b = core.vectors(D, sel, lab)[:2]
        ra, rb = core.avg_rank(a)[0], core.avg_rank(b)[0]
        rho = float(core.corr_rows(ra[None, :], rb)[0])
        if lab in EXPECT:
            hit = abs(rho - EXPECT[lab]) <= TOL
            ok &= hit
            print(f"  {lab:24s} rho={rho:.3f}  expected {EXPECT[lab]:.3f}  "
                  f"{'OK' if hit else 'FAIL'}")
        else:
            print(f"  {lab:24s} rho={rho:.3f}")
    print("\nREPRODUCTION", "CONFIRMED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
