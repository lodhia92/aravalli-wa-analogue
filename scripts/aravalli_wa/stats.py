"""
Rank statistics, permutation probabilities and the false-discovery-rate correction.

Every reported correlation, probability and corrected probability in the study comes from the
functions below, so that no two analyses can use different implementations of the same statistic.

Author: Bhavik Harish Lodhia, Curtin University
"""

import numpy as np
import pandas as pd

from .constants import NPERM, SEED


def avg_rank(A):
    """Row-wise ranks with ties given their average rank, as Spearman requires."""
    A = np.atleast_2d(A)
    m, n = A.shape
    order = np.argsort(A, axis=1, kind="mergesort")
    ranks = np.empty((m, n), float)
    ranks[np.arange(m)[:, None], order] = np.arange(1, n + 1)[None, :]
    S = np.take_along_axis(A, order, axis=1)
    for i in range(m):
        j = 0
        while j < n:
            k = j
            while k + 1 < n and S[i, k + 1] == S[i, j]:
                k += 1
            if k > j:
                ranks[i, order[i, j : k + 1]] = (j + k + 2) / 2.0
            j = k + 1
    return ranks


def corr_rows(Z, z):
    """Pearson correlation of every row of Z against z; NaN where a row has no spread."""
    Zc = Z - Z.mean(1, keepdims=True)
    zc = z - z.mean()
    den = np.sqrt((Zc**2).sum(1) * (zc**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (Zc @ zc) / den, np.nan)


def bh(ps):
    """Benjamini-Hochberg corrected probabilities, in the order the inputs were given.

    The family is whatever the caller passes. The paper corrects across the five fingerprints;
    passing a different set corrects across a different family and is a different claim.
    """
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    m = len(ps)
    q = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(o[::-1]):
        prev = min(prev, ps[idx] * m / (m - rank))
        q[idx] = prev
    return q


def spearman_perm(a, b, rng, nperm=NPERM):
    """Spearman correlation and its two-sided permutation probability.

    The generator is passed in rather than created here. Callers that draw repeatedly from one
    module-level generator depend on their own draw order for the exact probability, so moving
    this function must not change when each script draws. Use perm_p instead where an
    order-independent probability is wanted.
    """
    ra = pd.Series(a).rank().values
    rb = pd.Series(b).rank().values
    rho = float(np.corrcoef(ra, rb)[0, 1])
    B = rb[rng.random((nperm, len(rb))).argsort(axis=1)]
    Bc = B - B.mean(1, keepdims=True)
    ac = ra - ra.mean()
    den = np.sqrt((Bc**2).sum(1) * (ac**2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        null = np.abs(np.where(den > 0, (Bc @ ac) / den, np.nan))
    return rho, (np.sum(null >= abs(rho)) + 1) / (nperm + 1)


def perm_p(x, y, rho, n=NPERM, seed=SEED):
    """Two-sided permutation probability from a generator seeded per call, in blocks.

    Seeding per call makes the result independent of how many draws the calling script made
    earlier, unlike spearman_perm.
    """
    rng = np.random.default_rng(seed)
    xr = pd.Series(x).rank().values
    yr = pd.Series(y).rank().values
    xr = xr - xr.mean()
    sx = np.sqrt((xr**2).sum())
    hits, done = 0, 0
    while done < n:
        m = min(20000, n - done)
        Y = yr[rng.random((m, len(yr))).argsort(axis=1)]
        Y = Y - Y.mean(1, keepdims=True)
        den = sx * np.sqrt((Y**2).sum(1))
        with np.errstate(invalid="ignore", divide="ignore"):
            null = np.abs(np.where(den > 0, (Y @ xr) / den, np.nan))
        hits += int(np.sum(null >= abs(rho)))
        done += m
    return (hits + 1) / (n + 1)


def zscore_elem(pool, elems):
    """Mean and population standard deviation of log concentrations, ignoring non-positives.

    Values at or below zero become NaN rather than being replaced by a substitute, which is the
    study's treatment of values below detection and values not reported.
    """
    x = np.log(pool[elems].where(pool[elems] > 0))
    return x.mean(), x.std(ddof=0)
