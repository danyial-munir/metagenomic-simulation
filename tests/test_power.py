"""Sanity checks for the power-study code. Run: pytest -q"""
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova
import pandas as pd

from src.power import balanced_months, cosinor_p, driver, permanova_p, pseudo_F, gower_centre, simulate_seasonal


def test_pseudo_F_matches_skbio():
    rng = np.random.default_rng(0)
    Y = rng.normal(size=(40, 5)); Y[:10] += 0.8
    D = squareform(pdist(Y))
    lab = np.repeat(np.arange(4), 10)
    F = pseudo_F(gower_centre(D), lab[None, :], 4)[0]
    ids = [str(i) for i in range(40)]
    ref = permanova(DistanceMatrix(D, ids), pd.DataFrame({"g": lab}, index=ids), column="g", permutations=0)
    assert np.isclose(F, float(ref["test statistic"]))


def test_permanova_type_I_error():
    """Under H0 the permutation test rejects ~5 % (MC tolerance: 3 x sqrt(.05 .95 / 400) = 0.033)."""
    rng = np.random.default_rng(1)
    lab = np.repeat(np.arange(4), 12)
    p = [permanova_p(squareform(pdist(rng.lognormal(size=(48, 6)), "braycurtis")), lab, 99, rng)[0] for _ in range(400)]
    assert abs(np.mean(np.array(p) <= 0.05) - 0.05) < 0.033


def test_cosinor_null_and_alternative():
    rng = np.random.default_rng(2)
    m = balanced_months(96, rng)
    p0 = cosinor_p(rng.normal(size=(96, 5000)), m)             # 5000 null series at once
    assert abs((p0 <= 0.05).mean() - 0.05) < 0.01
    p1 = cosinor_p(rng.normal(size=(96, 2000)) + 1.0 * driver(m)[:, None], m)
    assert (p1 <= 0.05).mean() > 0.95


def test_seasonal_simulator_imposes_amplitude():
    """With almost no biological noise and huge depth, the log-ratio of an affected taxon to Other follows mu + A*s(m)."""
    rng = np.random.default_rng(3)
    mu, Sigma = np.full(21, -1.0), np.eye(21) * 1e-6
    X, months, taxa = simulate_seasonal(mu, Sigma, 1200, A=0.7, K=2, rng=rng, depth=10 ** 7)
    for j, sign in zip(taxa, [1, -1]):
        slope = np.polyfit(driver(months), np.log(X[:, j] / X[:, 21]), 1)[0]
        assert abs(slope - sign * 0.7) < 0.02
    assert np.bincount(months).tolist()[1:] == [100] * 12
