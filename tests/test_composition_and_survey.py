"""
Tests for the provenance ratios, the survey column lookup and the treatment of censored and
missing values, on a tiny synthetic frame.

The censored and missing handling is the behaviour most likely to break without anyone noticing,
because a substituted value produces a plausible number rather than an error.

Author: Bhavik Harish Lodhia, Curtin University
"""
import numpy as np
import pandas as pd
import pytest

from aravalli_wa.composition import logr, ratios
from aravalli_wa.constants import CI_CHONDRITE_PPM as CI
from aravalli_wa.constants import MATCH_RATIOS
from aravalli_wa.stats import zscore_elem
from aravalli_wa.survey import findcol


def frame(**overrides):
    """One sample carrying every element the ratios need, with optional overrides."""
    row = dict(Th=10.0, Sc=20.0, La=30.0, Co=5.0, Eu=1.0, Sm=4.0, Gd=5.0, Yb=2.0, Nb=8.0, Y=16.0)
    row.update(overrides)
    return pd.DataFrame([row])


class TestRatios:
    def test_produces_exactly_the_six_matching_ratios(self):
        assert list(ratios(frame()).columns) == MATCH_RATIOS

    def test_simple_ratios_are_the_plain_quotients(self):
        r = ratios(frame()).iloc[0]
        assert r["Th/Sc"] == pytest.approx(0.5)
        assert r["La/Sc"] == pytest.approx(1.5)
        assert r["Th/Co"] == pytest.approx(2.0)
        assert r["Nb/Y"] == pytest.approx(0.5)

    def test_europium_anomaly_uses_the_chondrite_normalisation(self):
        expected = (1.0 / CI["Eu"]) / np.sqrt((4.0 / CI["Sm"]) * (5.0 / CI["Gd"]))
        assert ratios(frame()).iloc[0]["EuEu"] == pytest.approx(expected)

    def test_normalised_lanthanum_ytterbium(self):
        expected = (30.0 / 2.0) / (CI["La"] / CI["Yb"])
        assert ratios(frame()).iloc[0]["La/Yb_n"] == pytest.approx(expected)

    def test_a_missing_element_gives_nan_not_a_guess(self):
        r = ratios(frame(Sc=np.nan)).iloc[0]
        assert np.isnan(r["Th/Sc"]) and np.isnan(r["La/Sc"])
        assert not np.isnan(r["Th/Co"])          # unaffected ratios still compute

    def test_a_zero_denominator_gives_infinity_which_logr_then_drops(self):
        assert np.isinf(ratios(frame(Sc=0.0)).iloc[0]["Th/Sc"])
        assert np.isnan(logr(frame(Sc=0.0)).iloc[0]["Th/Sc"])

    def test_the_index_is_preserved_so_rows_cannot_slip(self):
        df = frame()
        df.index = [17]
        assert list(ratios(df).index) == [17]


class TestLogRatios:
    def test_is_the_log_of_the_ratio(self):
        assert logr(frame()).iloc[0]["Th/Sc"] == pytest.approx(np.log(0.5))

    def test_a_negative_ratio_becomes_nan_rather_than_a_substitute(self):
        assert np.isnan(logr(frame(Th=-10.0)).iloc[0]["Th/Sc"])

    def test_a_zero_measurement_becomes_nan_rather_than_a_substitute(self):
        assert np.isnan(logr(frame(Th=0.0)).iloc[0]["Th/Sc"])

    def test_missing_stays_missing(self):
        assert np.isnan(logr(frame(Th=np.nan)).iloc[0]["Th/Sc"])


class TestCensoredAndMissingValues:
    """Values below detection and values not reported are dropped, never replaced."""

    def test_non_positive_values_do_not_enter_the_mean(self):
        pool = pd.DataFrame({"Ce": [10.0, 100.0, 0.0, -5.0]})
        mu, sd = zscore_elem(pool, ["Ce"])
        # only the two positive values contribute
        assert mu["Ce"] == pytest.approx(np.log([10.0, 100.0]).mean())

    def test_missing_values_do_not_enter_the_mean(self):
        pool = pd.DataFrame({"Ce": [10.0, 100.0, np.nan]})
        mu, _ = zscore_elem(pool, ["Ce"])
        assert mu["Ce"] == pytest.approx(np.log([10.0, 100.0]).mean())

    def test_a_substituted_floor_would_have_changed_the_answer(self):
        """Guards the choice: clipping at a small floor moves the mean a long way."""
        pool = pd.DataFrame({"Ce": [10.0, 100.0, 0.0]})
        dropped = zscore_elem(pool, ["Ce"])[0]["Ce"]
        clipped = np.log(pool["Ce"].clip(lower=1e-9)).mean()
        assert dropped == pytest.approx(np.log([10.0, 100.0]).mean())
        assert abs(clipped - dropped) > 5

    def test_standard_deviation_is_the_population_form(self):
        pool = pd.DataFrame({"Ce": [10.0, 100.0, 1000.0]})
        _, sd = zscore_elem(pool, ["Ce"])
        assert sd["Ce"] == pytest.approx(np.log([10.0, 100.0, 1000.0]).std(ddof=0))

    def test_a_column_with_no_usable_value_is_nan_not_an_exception(self):
        pool = pd.DataFrame({"Ce": [0.0, -1.0, np.nan]})
        mu, _ = zscore_elem(pool, ["Ce"])
        assert np.isnan(mu["Ce"])


class TestFindCol:
    HEADER = ["SITEID", "Th XRF ppm 1", "Sc ICP-MS ppm 0.1", "Th ICP-MS ppm 0.05",
              "La ICP-OES ppm 0.5", "Thorium total ppm"]

    def test_prefers_icp_ms_over_xrf_when_both_are_reported(self):
        assert findcol(self.HEADER, "Th") == 3

    def test_falls_back_to_xrf_when_that_is_all_there_is(self):
        assert findcol(["Nb XRF ppm 1"], "Nb") == 0

    def test_returns_none_when_the_element_is_absent(self):
        assert findcol(self.HEADER, "Dy") is None

    def test_returns_none_for_a_method_the_lookup_does_not_cover(self):
        """La is present but only by ICP-OES, which the lookup does not read."""
        assert findcol(self.HEADER, "La") is None

    def test_does_not_match_a_longer_element_name_by_prefix(self):
        """'Th' must not pick up a column whose name merely starts with those letters."""
        assert findcol(["Thorium ICP-MS ppm 1"], "Th") is None

    def test_tolerates_surrounding_whitespace_in_the_header(self):
        assert findcol(["  Zr ICP-MS ppm 2  "], "Zr") == 0

    def test_empty_header_returns_none(self):
        assert findcol([], "Zr") is None
