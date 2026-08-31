"""Build the Table 9 cell content (mineralogical validation) from the mineralogical_validation.py and
alternative_fingerprints.py result files,
then verify every emitted cell back against them. Writes results/table09_cells.tsv.

Output: results/table09_cells.tsv

Author: Bhavik Harish Lodhia, Curtin University
"""
import os

import pandas as pd

def get(src, key, scale):
    d = A15 if src == "a15" else A8
    col = "fingerprint"
    r = d[(d[col] == key) & (d.scale == scale)]
    assert len(r) == 1, (src, key, scale, len(r))
    r = r.iloc[0]
    return int(r.n), float(r.rho), float(r.perm_p)

def fp(p):
    return "<0.001" if p < 0.001 else ("%.3f" % p)

def main():
    global A15, A8, HERE, RES, ROWS, SCALE, bad, f, hdr, j, key, label, mineral, n, n1, n2, out, p, p1, p2, pair_scale, pool_scale, r, r1, r2, row, scale, src, want
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    RES = os.path.join(HERE, "results")
    A8 = pd.read_csv(os.path.join(RES, "mineralogical_validation.csv"))
    A15 = pd.read_csv(os.path.join(RES, "alternative_fingerprints.csv"))
    ROWS = [
        ("Ce, Nd, Pr (as published)", "a15", "Monazite: Ce+Nd+Pr (as published)", "Monazite"),
        ("Ce", "a15", "Monazite: Ce", "Monazite"),
        ("La, Ce, Pr, Nd, Sm", "a15", "Monazite: La+Ce+Pr+Nd+Sm", "Monazite"),
        ("Ce, Nd, Pr, Th", "a15", "Monazite: Ce+Nd+Pr+Th", "Monazite"),
        ("Ce, Nd, Pr, Th, P", "a15", "Monazite: Ce+Nd+Pr+Th+P", "Monazite"),
        ("La, Ce, Ca", "a15", "Allanite: La+Ce+Ca", "Allanite"),
        ("Dy (as published)", "a15", "Xenotime: Dy (as published)", "Xenotime"),
        ("Y", "a15", "Xenotime: Y", "Xenotime"),
        ("Y, Dy", "a15", "Xenotime: Y+Dy", "Xenotime"),
        ("Dy, Ho, Er, Yb, Lu", "a15", "Xenotime: Dy+Ho+Er+Yb+Lu", "Xenotime"),
        ("Y, Dy, Ho, Er, Yb, Lu", "a15", "Xenotime: Y+Dy+Ho+Er+Yb+Lu", "Xenotime"),
        ("Zr, Hf", "a8", "Zircon", "Zircon"),
        ("Ti", "a8", "Ti oxides", "Ilmenite"),
        ("P", "a8", "Apatite", "Apatite"),
    ]
    SCALE = {"a15": ("pool (n=85)", "pairs (n=20)"),
             "a8": ("drainage-selected pool", "analogue pairs")}
    out = []
    for label, src, key, mineral in ROWS:
        pool_scale, pair_scale = SCALE[src]
        n1, r1, p1 = get(src, key, pool_scale)
        n2, r2, p2 = get(src, key, pair_scale)
        out.append([label, mineral, str(n1), "%.3f" % r1, fp(p1), str(n2), "%.3f" % r2, fp(p2)])
    hdr = ["Fingerprint elements", "Measured mineral", "Sites in pool", "rho (pool)", "p (pool)",
           "Analogue pairs", "rho (pairs)", "p (pairs)"]
    with open(os.path.join(RES, "table09_cells.tsv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\t".join(hdr) + "\n")
        for r in out:
            f.write("\t".join(r) + "\n")
    bad = 0
    for (label, src, key, mineral), row in zip(ROWS, out):
        pool_scale, pair_scale = SCALE[src]
        for j, scale in ((2, pool_scale), (5, pair_scale)):
            n, r, p = get(src, key, scale)
            want = [str(n), "%.3f" % r, fp(p)]
            if row[j:j + 3] != want:
                bad += 1
                print("MISMATCH", label, scale, row[j:j + 3], want)
        if row[1] != mineral:
            bad += 1
    print("rows: %d, cell mismatches against source: %d" % (len(out), bad))
    for r in [hdr] + out:
        print(" | ".join(r))

if __name__ == "__main__":
    main()
