"""Does a broader heavy rare-earth fingerprint transfer between India and Australia?

The published heavy rare-earth fingerprint rests on Dy alone and gives a corrected q of 0.059, so
it is not counted as transferring. Whether a fingerprint built on additional heavy rare-earth and
yttrium elements transfers is a separate question from whether those elements predict measured
xenotime, which mineralogical_validation.py and alternative_fingerprints.py address.

Four variants are tested, using a fixed element list decided in advance, so that this is not a
search for a combination that passes:
    Y
    Y + Dy
    Dy + Ho + Er + Yb + Lu
    Y + Dy + Ho + Er + Yb + Lu
Tm is absent from both surveys, and Tb is excluded to keep the list fixed in advance.

Multiple comparisons. The four variants are corrected in ONE Benjamini-Hochberg family together
with the five published fingerprints, nine tests in total. This is the conservative treatment and
is deliberate: testing heavy rare-earth combinations until one passes is precisely the selection
concern addressed by matching_sensitivity.py. The q-values reported for the published five stand on
their own declared family of five and are printed here alongside for comparison, not replaced.

The analysis is otherwise identical to the published one: same drainage containment thresholds,
same matching, same twenty pairs, same within-survey standardisation, same hundred thousand
permutations and seed. Extending the element list cannot move the published scores, and the script
asserts that it does not: monazite must return 0.696 and the published Dy fingerprint 0.506.

Inputs : as containment_sensitivity (Indian and Australian geochemistry, drainage tables)
Outputs: results/hree_transfer.csv

Author: Bhavik Harish Lodhia, Curtin University, bhavik.lodhia@curtin.edu.au
Repository: aravalli-wa-analogue. Run order is given in README.md; data sources and
expected file locations are given in data/README.md.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import matching_sensitivity as core  # 
import containment_sensitivity as cont   # 

RES = os.path.join(HERE, "results")
DOMS = core.DOMS

NEW = {
    "Heavy rare earths, Y": ["Y"],
    "Heavy rare earths, Y+Dy": ["Y", "Dy"],
    "Heavy rare earths, Dy+Ho+Er+Yb+Lu": ["Dy", "Ho", "Er", "Yb", "Lu"],
    "Heavy rare earths, Y+Dy+Ho+Er+Yb+Lu": ["Y", "Dy", "Ho", "Er", "Yb", "Lu"],
}
EXTRA = ["Y", "Ho", "Er", "Yb", "Lu"]
PUB_RHO = {"monazite (Ce,Nd,Pr)": 0.696, "xenotime (Dy)": 0.506}


def main():
    # extend the element list and the fingerprint definitions, then rebuild from scratch
    cont.PATH = list(dict.fromkeys(list(core.PATH) + EXTRA))
    cont.MIN = dict(core.MIN)
    cont.MIN.update(NEW)
    cont._D.clear()
    D = cont.build_D(cont.PUB_IN, cont.PUB_AU)
    sel = core.select(D, "published")

    pub_labels = list(core.MIN)
    new_labels = list(NEW)
    all_labels = pub_labels + new_labels

    res = {}
    for lab in all_labels:
        a, b, _ = core.vectors(D, sel, lab)
        res[lab] = core.spearman_perm(a, b) + (len(a),)

    for lab, want in PUB_RHO.items():
        got = round(res[lab][0], 3)
        assert got == want, f"published anchor moved: {lab} {got} != {want}"
    print("published anchor holds after extending the element list: monazite %.3f, Dy %.3f"
          % (res[pub_labels[0]][0], res[pub_labels[1]][0]), flush=True)

    q_pub_family = dict(zip(pub_labels, core.bh([res[l][1] for l in pub_labels])))
    q_all_family = dict(zip(all_labels, core.bh([res[l][1] for l in all_labels])))
    q_new_only = dict(zip(new_labels, core.bh([res[l][1] for l in new_labels])))

    rows = []
    for lab in all_labels:
        r, p, n = res[lab]
        is_new = lab in NEW
        rows.append(dict(
            fingerprint=lab,
            elements="+".join(NEW[lab]) if is_new else "+".join(core.MIN[lab]),
            status="broadened heavy rare-earth variant" if is_new else "published fingerprint",
            n_pairs=n, rho=round(r, 3), perm_p=round(p, 4),
            q_published_family_of_five=("" if is_new else round(float(q_pub_family[lab]), 4)),
            q_new_variants_alone=(round(float(q_new_only[lab]), 4) if is_new else ""),
            q_combined_family_of_nine=round(float(q_all_family[lab]), 4),
            transfers_combined_family="Yes" if q_all_family[lab] < 0.05 else "No"))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "hree_transfer.csv"), index=False)
    pd.set_option("display.width", 300)
    print()
    print(out.to_string(index=False))

    print("\nsummary")
    best = max(new_labels, key=lambda l: res[l][0])
    print("  strongest new heavy rare-earth variant: %s, rho %.3f, q %.4f in the family of nine"
          % (best, res[best][0], q_all_family[best]))
    print("  any new variant clearing q<0.05 in the family of nine: %s"
          % (", ".join(l for l in new_labels if q_all_family[l] < 0.05) or "none"))
    print("  any new variant clearing q<0.05 even in a family of four: %s"
          % (", ".join(l for l in new_labels if q_new_only[l] < 0.05) or "none"))
    print("  monazite in the family of nine: rho %.3f, q %.4f"
          % (res[pub_labels[0]][0], q_all_family[pub_labels[0]]))
    print("\nwrote results/hree_transfer.csv")


if __name__ == "__main__":
    main()
