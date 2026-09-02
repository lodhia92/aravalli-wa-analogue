"""Are the single-element fingerprints too narrow to represent their host mineral?

Three of the five published fingerprints are single-element scores: the heavy rare earths on Dy,
the titanium oxides on Ti and apatite on P. Multivariate versions of the rare-earth fingerprints
already exist: alternative_fingerprints.py tests broader monazite and xenotime combinations
against measured Heavy Mineral Map of Australia grain counts, and hree_transfer.py tests
broader heavy rare-earth combinations for transfer between the matched catchments. What has not
been tested is apatite, and the remaining two fingerprints.

This script broadens each of them with the elements that actually substitute into or accompany the
host, and asks whether transfer improves:
    apatite            P, and P with Ca, with Sr, and with both; Ca and Sr are the major
                       substituting cations in the apatite structure
    titanium oxides    Ti, and Ti with Nb and Ta, which substitute into rutile and ilmenite
    zircon             Zr with Hf as published, and with U and Th added, which sit in zircon
All are present in both surveys.

Multiple comparisons are handled as in hree_transfer.py: the new variants are corrected in one
Benjamini-Hochberg family together with the five published fingerprints, which is the conservative
treatment. A per-fingerprint family is reported alongside so that the choice of family can be seen
not to drive the result.

The published scores cannot move when the element list is extended and the script asserts that
they do not, reproducing 0.696 for the light rare earths and 0.506 for the published heavy
rare-earth fingerprint before anything else is computed.

Outputs: results/broadened_fingerprints.csv

Author: Bhavik Harish Lodhia, Curtin University
"""

import os

import containment_sensitivity as cont
import matching_sensitivity as core
import pandas as pd

from aravalli_wa.stats import bh

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(HERE, "results")

NEW = {
    "Apatite, P+Ca": ["P", "Ca"],
    "Apatite, P+Sr": ["P", "Sr"],
    "Apatite, P+Ca+Sr": ["P", "Ca", "Sr"],
    "Titanium oxides, Ti+Ta": ["Ti", "Ta"],
    "Zircon, Zr+Hf+U": ["Zr", "Hf", "U"],
}

# Diagnostics, deliberately NOT corrected with the others and NOT offered as results. The six
# matching ratios are Th/Sc, La/Sc, Th/Co, Eu/Eu*, La/Yb and Nb/Y, so Th, La, Yb, Nb, Y, Sc, Co and
# Eu are matching variables. A fingerprint containing one of them is correlated across the pairs by
# construction, because the pairs were built on it. These two rows demonstrate that directly.
DEMO = {
    "Diagnostic: Th alone (a matching variable)": ["Th"],
    "Diagnostic: Zircon Zr+Hf+U+Th (contains Th)": ["Zr", "Hf", "U", "Th"],
}
GROUP = {
    "Apatite, P+Ca": "apatite",
    "Apatite, P+Sr": "apatite",
    "Apatite, P+Ca+Sr": "apatite",
    "Titanium oxides, Ti+Ta": "titanium oxides",
    "Zircon, Zr+Hf+U": "zircon",
}
BASE = {"apatite": "apatite (P)", "titanium oxides": "Ti-oxide (Ti)", "zircon": "zircon (Zr,Hf)"}
EXTRA = [
    "Ca",
    "Sr",
    "Nb",
    "Ta",
    "U",
    "Th",
]  # Th is used by the zircon variant and is not in the published PATH
PUB_RHO = {"monazite (Ce,Nd,Pr)": 0.696, "xenotime (Dy)": 0.506}


def main():
    # The Indian table reports calcium as the oxide and does not coerce the trace elements this
    # script needs, so both are prepared here exactly as the loader already does for P and Ti.
    # Units need not match between surveys: every score is a z-score of logs taken within one
    # survey, so a constant factor cancels.
    # PATH must be extended BEFORE the raw load, because the loader chooses which Australian
    # columns to read from it. Clear both caches so nothing built on the short list survives.
    cont.PATH = list(dict.fromkeys(list(core.PATH) + EXTRA))
    cont._RAW.clear()
    cont._D.clear()
    R = cont.load_raw()
    ar = R["ar"]
    if "Ca" not in ar.columns:
        ar["Ca"] = pd.to_numeric(ar["CaO"], errors="coerce") * 0.7147
    for c in EXTRA:
        if c in ar.columns:
            ar[c] = pd.to_numeric(ar[c], errors="coerce")
    missing = [c for c in EXTRA if c not in ar.columns]
    assert not missing, f"Indian table lacks {missing}"

    cont.MIN = dict(core.MIN)
    cont.MIN.update(NEW)
    cont.MIN.update(DEMO)
    D = cont.build_D(cont.IN_THR_PCT, cont.AUS_THR_PCT)
    sel = core.select(D, "published")

    pub, new, demo = list(core.MIN), list(NEW), list(DEMO)
    res = {}
    for lab in pub + new + demo:
        a, b, _ = core.vectors(D, sel, lab)
        res[lab] = core.spearman_perm(a, b) + (len(a),)
    for lab, want in PUB_RHO.items():
        got = round(res[lab][0], 3)
        assert got == want, f"published anchor moved: {lab} {got} != {want}"
    print("anchor holds: monazite %.3f, Dy %.3f" % (res[pub[0]][0], res[pub[1]][0]), flush=True)

    q_all = dict(zip(pub + new, bh([res[l][1] for l in pub + new])))
    q_pub = dict(zip(pub, bh([res[l][1] for l in pub])))
    q_grp = {}
    for g in set(GROUP.values()):
        members = [BASE[g]] + [l for l in new if GROUP[l] == g]
        q_grp.update(dict(zip(members, bh([res[l][1] for l in members]))))

    rows = []
    for lab in pub + new + demo:
        r, p, n = res[lab]
        is_new = lab in NEW
        is_demo = lab in DEMO
        g = GROUP.get(lab) or {v: k for k, v in BASE.items()}.get(lab, "")
        rows.append(
            dict(
                fingerprint=lab,
                elements="+".join((DEMO if is_demo else NEW if is_new else core.MIN)[lab]),
                group=g,
                status=(
                    "matching-element diagnostic, excluded from correction"
                    if is_demo
                    else "broadened variant"
                    if is_new
                    else "published fingerprint"
                ),
                n_pairs=n,
                rho=round(r, 3),
                perm_p=round(p, 4),
                q_published_family_of_five=(
                    "" if is_new or is_demo else round(float(q_pub[lab]), 4)
                ),
                q_within_its_own_mineral=(round(float(q_grp[lab]), 4) if lab in q_grp else ""),
                q_combined_family=("" if is_demo else round(float(q_all[lab]), 4)),
                transfers=(
                    "not corrected, see status" if is_demo else "Yes" if q_all[lab] < 0.05 else "No"
                ),
            )
        )
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "broadened_fingerprints.csv"), index=False)
    pd.set_option("display.width", 300)
    print()
    print(out.to_string(index=False))
    print("\nsummary")
    print("  matching-element diagnostics, NOT results:")
    for lab in demo:
        print("    %-46s rho %+.3f, uncorrected p %.4f" % (lab, res[lab][0], res[lab][1]))
    for g in ("apatite", "titanium oxides", "zircon"):
        base = BASE[g]
        members = [l for l in new if GROUP[l] == g]
        print(
            "  %-16s published %+.3f -> broadened %s"
            % (
                g,
                res[base][0],
                ", ".join(
                    "%s %+.3f (q %.3f)" % (l.split(", ")[1], res[l][0], q_all[l]) for l in members
                ),
            )
        )
    print(
        "  any broadened variant clearing q<0.05: %s"
        % (", ".join(l for l in new if q_all[l] < 0.05) or "none")
    )
    print("\nwrote results/broadened_fingerprints.csv")


if __name__ == "__main__":
    main()
