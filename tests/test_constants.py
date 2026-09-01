"""
Tests for the analytical constants.

These do not check arithmetic; they guard the values a published number depends on, so that a
change to any of them has to be deliberate rather than incidental.

Author: Bhavik Harish Lodhia, Curtin University
"""

import pytest

from aravalli_wa import constants as k


class TestPublishedConventions:
    def test_permutation_convention(self):
        assert k.NPERM == 100000
        assert k.SEED == 20260827

    def test_containment_cut_offs(self):
        assert k.IN_THR_PCT == 50.0
        assert k.AUS_THR_PCT == 25.0

    def test_the_five_fingerprints_and_their_elements(self):
        assert k.FINGERPRINTS == {
            "monazite (Ce,Nd,Pr)": ["Ce", "Nd", "Pr"],
            "xenotime (Dy)": ["Dy"],
            "zircon (Zr,Hf)": ["Zr", "Hf"],
            "Ti-oxide (Ti)": ["Ti"],
            "apatite (P)": ["P"],
        }

    def test_fingerprint_order_is_the_published_one(self):
        """The order is load-bearing: it decides which permutation block each fingerprint gets."""
        assert k.FINGERPRINT_LABELS == [
            "monazite (Ce,Nd,Pr)",
            "xenotime (Dy)",
            "zircon (Zr,Hf)",
            "Ti-oxide (Ti)",
            "apatite (P)",
        ]

    def test_every_fingerprint_element_is_a_pathfinder_element(self):
        for elements in k.FINGERPRINTS.values():
            for element in elements:
                assert element in k.PATHFINDER_ELEMENTS

    def test_the_eight_pathfinder_elements(self):
        assert k.PATHFINDER_ELEMENTS == ["Zr", "Hf", "Ti", "Ce", "Nd", "Pr", "Dy", "P"]

    def test_the_six_matching_ratios(self):
        assert k.MATCH_RATIOS == ["Th/Sc", "La/Sc", "Th/Co", "EuEu", "La/Yb_n", "Nb/Y"]

    def test_the_two_domains(self):
        assert k.DOMAINS == ["Palaeoproterozoic", "Archaean"]


class TestOxideConversions:
    """The factors are the element's mass fraction of the oxide, rounded to four decimals.

    The rounding is immaterial to every published number: the analysis takes logs, standardises,
    and then ranks, so a constant multiplier is an additive shift that the standardisation removes
    and the ranking ignores. test_a_constant_factor_cancels below holds that property in place.
    """

    def test_values_in_use(self):
        assert k.TI_MASS_FRACTION_OF_TIO2 == 0.5995
        assert k.P_MASS_FRACTION_OF_P2O5 == 0.4364
        assert k.AL_MASS_FRACTION_OF_AL2O3 == 0.5293

    def test_each_is_the_element_mass_fraction_of_its_oxide_to_four_decimals(self):
        assert k.TI_MASS_FRACTION_OF_TIO2 == pytest.approx(47.867 / 79.866, abs=2e-4)
        assert k.P_MASS_FRACTION_OF_P2O5 == pytest.approx(2 * 30.974 / 141.944, abs=2e-4)
        assert k.AL_MASS_FRACTION_OF_AL2O3 == pytest.approx(2 * 26.982 / 101.961, abs=2e-4)

    def test_one_weight_per_cent_is_ten_thousand_mg_per_kg(self):
        assert k.WT_PCT_TO_MG_KG == 1e4

    def test_a_constant_factor_cancels(self):
        """Why the four-decimal rounding cannot move a result."""
        import numpy as np

        from aravalli_wa.stats import avg_rank

        oxide = np.random.default_rng(0).lognormal(0.0, 1.0, size=40)

        def standardised(factor):
            x = np.log(oxide * k.WT_PCT_TO_MG_KG * factor)
            return (x - x.mean()) / x.std(ddof=0)

        rounded = standardised(k.TI_MASS_FRACTION_OF_TIO2)
        exact = standardised(47.867 / 79.866)
        assert np.allclose(rounded, exact, atol=1e-12)
        assert (avg_rank(rounded)[0] == avg_rank(exact)[0]).all()


class TestChondriteReference:
    def test_carries_every_element_the_ratios_normalise_by(self):
        assert set(k.CI_CHONDRITE_PPM) == {"La", "Yb", "Sm", "Eu", "Gd"}

    def test_values(self):
        assert k.CI_CHONDRITE_PPM == dict(La=0.237, Yb=0.170, Sm=0.148, Eu=0.0580, Gd=0.199)
