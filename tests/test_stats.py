"""
Tests for the rank statistics, the false-discovery-rate correction and the permutation
probabilities in aravalli_wa.stats.

These run without any of the licensed input datasets.

Author: Bhavik Harish Lodhia, Curtin University
"""

import numpy as np
import pytest
from scipy.stats import spearmanr

from aravalli_wa.stats import avg_rank, bh, corr_rows, perm_p, spearman_perm


def spearman(a, b):
    """Spearman correlation via the package's own rank and correlation primitives."""
    ra, rb = avg_rank(a)[0], avg_rank(b)[0]
    return float(corr_rows(ra[None, :], rb)[0])


class TestAvgRank:
    def test_hand_worked_example_without_ties(self):
        assert avg_rank([10.0, 30.0, 20.0, 40.0])[0].tolist() == [1.0, 3.0, 2.0, 4.0]

    def test_ties_take_the_average_of_the_ranks_they_span(self):
        # values 5,5 occupy ranks 2 and 3, so both take 2.5; 9,9,9 occupy 4,5,6 -> 5.0
        assert avg_rank([1.0, 5.0, 5.0, 9.0, 9.0, 9.0])[0].tolist() == [
            1.0,
            2.5,
            2.5,
            5.0,
            5.0,
            5.0,
        ]

    def test_all_values_equal_gives_every_element_the_mid_rank(self):
        assert avg_rank([7.0, 7.0, 7.0, 7.0])[0].tolist() == [2.5, 2.5, 2.5, 2.5]

    def test_single_element(self):
        assert avg_rank([42.0])[0].tolist() == [1.0]

    def test_matches_scipy_on_random_input_including_ties(self):
        rng = np.random.default_rng(12345)
        for _ in range(50):
            n = int(rng.integers(4, 40))
            # integers so ties occur often, which is where a rank implementation goes wrong
            x = rng.integers(0, 8, size=n).astype(float)
            y = rng.integers(0, 8, size=n).astype(float)
            if np.ptp(x) == 0 or np.ptp(y) == 0:
                continue
            assert spearman(x, y) == pytest.approx(spearmanr(x, y).statistic, abs=1e-12)


class TestCorrRows:
    def test_perfect_positive_and_negative(self):
        z = np.array([1.0, 2.0, 3.0, 4.0])
        assert corr_rows(z[None, :], z)[0] == pytest.approx(1.0)
        assert corr_rows((-z)[None, :], z)[0] == pytest.approx(-1.0)

    def test_a_row_with_no_spread_is_nan_not_a_division_error(self):
        flat = np.ones((1, 5))
        assert np.isnan(corr_rows(flat, np.arange(5.0))[0])


class TestBenjaminiHochberg:
    # The fifteen probabilities of Benjamini & Hochberg (1995), Table 1, from the multiple
    # endpoints of Needleman et al. At alpha = 0.05 the procedure rejects the first four.
    BH1995 = [
        0.0001,
        0.0004,
        0.0019,
        0.0095,
        0.0201,
        0.0278,
        0.0298,
        0.0344,
        0.0459,
        0.3240,
        0.4262,
        0.5719,
        0.6528,
        0.7590,
        1.0000,
    ]

    def test_rejects_exactly_four_at_five_per_cent(self):
        assert int((bh(self.BH1995) <= 0.05).sum()) == 4

    def test_worked_q_values(self):
        q = bh(self.BH1995)
        assert q[0] == pytest.approx(0.0015, abs=1e-9)
        assert q[1] == pytest.approx(0.0030, abs=1e-9)
        assert q[2] == pytest.approx(0.0095, abs=1e-9)
        assert q[3] == pytest.approx(0.035625, abs=1e-9)
        assert q[4] == pytest.approx(0.0603, abs=1e-9)
        assert q[14] == pytest.approx(1.0, abs=1e-9)

    def test_returned_in_input_order_not_sorted_order(self):
        shuffled = [0.0095, 0.0001, 1.0000, 0.0004]
        q = bh(shuffled)
        assert np.argmin(q) == 1  # the smallest p is still at index 1
        assert q[2] == pytest.approx(max(q))

    def test_monotone_non_decreasing_in_the_sorted_probabilities(self):
        q = bh(self.BH1995)
        assert np.all(np.diff(q) >= -1e-12)

    def test_never_exceeds_one(self):
        assert np.all(bh([0.9, 0.95, 0.99]) <= 1.0)

    def test_single_hypothesis_is_its_own_probability(self):
        assert bh([0.031])[0] == pytest.approx(0.031)


class TestPermutationProbabilities:
    def test_perm_p_is_deterministic_for_a_fixed_seed(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=25)
        y = x * 0.7 + rng.normal(size=25) * 0.7
        rho = spearman(x, y)
        first = perm_p(x, y, rho, n=2000, seed=20260827)
        second = perm_p(x, y, rho, n=2000, seed=20260827)
        assert first == second

    def test_perm_p_does_not_depend_on_earlier_draws(self):
        """perm_p seeds per call, so an unrelated draw beforehand must not move it."""
        rng = np.random.default_rng(0)
        x = rng.normal(size=20)
        y = x + rng.normal(size=20)
        rho = spearman(x, y)
        before = perm_p(x, y, rho, n=2000, seed=20260827)
        np.random.default_rng(1).random(1000)
        after = perm_p(x, y, rho, n=2000, seed=20260827)
        assert before == after

    def test_spearman_perm_is_reproducible_from_an_equally_seeded_generator(self):
        x = np.arange(20.0)
        y = np.array([3.0, 1, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18, 20])
        a = spearman_perm(x, y, np.random.default_rng(20260827), nperm=2000)
        b = spearman_perm(x, y, np.random.default_rng(20260827), nperm=2000)
        assert a == b

    def test_spearman_perm_rho_agrees_with_scipy(self):
        rng = np.random.default_rng(7)
        x = rng.normal(size=30)
        y = x * 0.5 + rng.normal(size=30)
        rho, _ = spearman_perm(x, y, np.random.default_rng(1), nperm=200)
        assert rho == pytest.approx(spearmanr(x, y).statistic, abs=1e-12)

    def test_probability_can_never_be_zero(self):
        """The (hits + 1) / (nperm + 1) form keeps a probability strictly positive."""
        x = np.arange(30.0)
        p = perm_p(x, x, 1.0, n=1000, seed=1)
        assert p > 0
        assert p == pytest.approx(1.0 / 1001, abs=1e-12)
