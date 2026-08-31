"""Build the Table 8 cell content, the matching-sensitivity table, straight from the result files of
matching_sensitivity.py,
then verify every emitted cell back against those files. Writes results/table08_cells.tsv.

An asterisk marks a fingerprint whose Benjamini-Hochberg q-value is below 0.05 within that variant.
The two control rows carry the mean and the 95th percentile of 10 000 random-pairing draws.

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "results")
V = pd.read_csv(os.path.join(RES, "matching_variants.csv"))
C = pd.read_csv(os.path.join(RES, "matching_random_control.csv"))

COLS = ["monazite (Ce,Nd,Pr)", "xenotime (Dy)", "zircon (Zr,Hf)", "Ti-oxide (Ti)", "apatite (P)"]
ROWS = [
    ("Mutual nearest neighbour, ten pairs per domain (published)",
     "mutual nearest neighbour, top 10 per domain (published)"),
    ("One-directional nearest neighbour, mutuality not required",
     "one-directional nearest neighbour, mutuality not required"),
    ("Greedy unique assignment, each Indian catchment used once",
     "greedy unique, each Indian site used once"),
    ("Globally optimal one-to-one assignment",
     "globally optimal one-to-one assignment"),
    ("Mahalanobis distance in place of Euclidean",
     "Mahalanobis distance in place of Euclidean"),
    ("All mutual pairs, no distance ranking cap",
     "mutual pairs only, no distance ranking cap"),
    ("Five pairs per domain", "top 5 pairs per domain"),
    ("Eight pairs per domain", "top 8 pairs per domain"),
    ("Twelve pairs per domain", "top 12 pairs per domain"),
    ("Fifteen pairs per domain", "top 15 pairs per domain"),
    ("Twenty pairs per domain", "top 20 pairs per domain"),
    ("Thirty pairs per domain", "top 30 pairs per domain"),
    ("Mutual pairs within the 25th distance percentile",
     "mutual pairs within the 25th distance percentile"),
    ("Mutual pairs within the 50th distance percentile",
     "mutual pairs within the 50th distance percentile"),
    ("Mutual pairs within the 75th distance percentile",
     "mutual pairs within the 75th distance percentile"),
]


def fmt(rho, q):
    return ("%.2f" % rho) + ("*" if q < 0.05 else "")


out = []
for label, key in ROWS:
    sub = V[V.variant == key]
    assert len(sub) == 5, (key, len(sub))
    n = int(sub.n.iloc[0])
    shared = int(sub.shared_with_published.iloc[0])
    cells = []
    for m in COLS:
        r = sub[sub.mineral == m].iloc[0]
        cells.append(fmt(r.rho, r.q))
    out.append([label, str(n), str(shared)] + cells)

for stat, col in [("Random analogue control, mean of 10 000 draws", "rho_random_mean"),
                  ("Random analogue control, 95th percentile", "rho_random_p95")]:
    cells = []
    for m in COLS:
        r = C[C.mineral == m].iloc[0]
        cells.append("%.2f" % r[col])
    out.append([stat, "20", ""] + cells)

hdr = ["Matching variant", "Pairs", "Shared with published set",
       "Light rare earths (Ce, Nd, Pr)", "Heavy rare earths (Dy)", "Zircon (Zr, Hf)",
       "Titanium oxides (Ti)", "Apatite (P)"]
with open(os.path.join(RES, "table08_cells.tsv"), "w", encoding="utf-8", newline="\n") as f:
    f.write("\t".join(hdr) + "\n")
    for r in out:
        f.write("\t".join(r) + "\n")

# verification: every emitted cell re-derived from the source files independently of the loop above
bad = 0
for row in out:
    label = row[0]
    if label.startswith("Random"):
        col = "rho_random_mean" if "mean" in label else "rho_random_p95"
        for j, m in enumerate(COLS):
            want = "%.2f" % float(C.loc[C.mineral == m, col].iloc[0])
            if row[3 + j] != want:
                bad += 1
                print("MISMATCH", label, m, row[3 + j], want)
        continue
    key = dict(ROWS)[label]
    sub = V[V.variant == key]
    if row[1] != str(int(sub.n.iloc[0])) or row[2] != str(int(sub.shared_with_published.iloc[0])):
        bad += 1
        print("MISMATCH n/shared", label)
    for j, m in enumerate(COLS):
        r = sub[sub.mineral == m].iloc[0]
        want = ("%.2f" % r.rho) + ("*" if r.q < 0.05 else "")
        if row[3 + j] != want:
            bad += 1
            print("MISMATCH", label, m, row[3 + j], want)
print("rows: %d, cells verified against source: %d mismatches" % (len(out), bad))
for r in [hdr] + out:
    print(" | ".join(r))
