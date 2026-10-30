"""Sanity checks for the generative models. Run: pytest -q"""
import itertools

import numpy as np
from scipy.optimize import check_grad

from src.models import (calibrate_zilnm, dm_moments, fit_dm, sample_zilnm, fit_dm_reg, fit_lnm, fit_multinomial, loglik_dm, loglik_dm_reg,
                        loglik_lnm, loglik_lnm_is, loglik_multinomial, sample_dm, sample_lnm)


def _all_compositions(N, D):
    """Every count vector with D parts summing to N."""
    return np.array([c for c in itertools.product(range(N + 1), repeat=D) if sum(c) == N])


def test_likelihoods_are_probabilities():
    """Summed over all possible count vectors, each model's probabilities add to 1."""
    X = _all_compositions(5, 3)
    p, alpha = np.array([.2, .3, .5]), np.array([.4, 1.5, .7])
    assert np.isclose(np.exp(loglik_multinomial(X, p)).sum(), 1)
    assert np.isclose(np.exp(loglik_dm(X, alpha)).sum(), 1)


def test_dm_tends_to_multinomial():
    """As precision A -> infinity the DM loses its overdispersion and becomes the multinomial."""
    X = _all_compositions(6, 3)
    p = np.array([.2, .3, .5])
    assert np.allclose(loglik_dm(X, p * 1e8), loglik_multinomial(X, p), atol=1e-5)


def test_dm_simulator_matches_theory():
    rng = np.random.default_rng(0)
    alpha, N, n = np.array([0.5, 2.0, 1.0, 6.5]), 200, 40000
    X = sample_dm(alpha, N, n, rng)
    m, v = dm_moments(alpha, N)
    assert np.allclose(X.mean(0), m, rtol=0.02)              # Monte Carlo tolerance (n = 40000)
    assert np.allclose(X.var(0), v, rtol=0.05)


def test_fits_recover_parameters():
    rng = np.random.default_rng(1)
    alpha = np.array([0.3, 1.2, 0.5, 3.0])
    X = sample_dm(alpha, 500, 3000, rng)
    assert np.allclose(fit_multinomial(X), alpha / alpha.sum(), atol=0.02)
    assert np.allclose(fit_dm(X), alpha, rtol=0.1)


def test_dm_reg_gradient_and_fit():
    rng = np.random.default_rng(2)
    n, D, top = 1500, 5, np.array([0, 1])
    Z = np.c_[rng.standard_normal(n), rng.integers(0, 2, n)]
    b, B = np.log([1.0, 2.0, 0.5, 0.8, 1.5]), np.array([[0.5, -0.3], [0.0, 0.6]])
    eta = np.tile(b, (n, 1)); eta[:, top] += Z @ B
    X = sample_dm(np.exp(eta), 300, n, rng)
    theta = fit_dm_reg(X, Z, top)
    assert np.allclose(theta[D:].reshape(2, 2), B, atol=0.15)
    # analytic gradient used by L-BFGS agrees with finite differences
    f = lambda t: -loglik_dm_reg(X[:50], Z[:50], t, top).sum()
    from src.models import _dm_grad_alpha, _m2_alpha
    def g(t):
        a = _m2_alpha(t, Z[:50], D, top); G = _dm_grad_alpha(X[:50].astype(float), a) * a
        return -np.r_[G.sum(0), (Z[:50].T @ G[:, top]).ravel()]
    assert check_grad(f, g, theta) < 1e-3 * np.abs(g(theta)).max() + 1e-4


def test_lnm_laplace_close_to_importance_sampling_and_fit():
    rng = np.random.default_rng(3)
    mu, Sigma = np.array([0.5, -0.5, 1.0]), np.array([[1.0, .3, 0], [.3, .8, .2], [0, .2, .5]])
    X = sample_lnm(mu, Sigma, 400, 2000, rng)
    for x in X[:5]:
        la = loglik_lnm(x[None, :], mu, Sigma)[0]
        is_, ess = loglik_lnm_is(x, mu, Sigma, 20000, rng)
        assert abs(la - is_) < 0.05 and ess > 1000
    mu_h, S_h, trace = fit_lnm(X, shrink=0.0)
    assert np.allclose(mu_h, mu, atol=0.1) and np.allclose(S_h, Sigma, atol=0.15)


def test_zilnm_calibration_recovers_zero_inflation():
    """Data from a zero-inflated LNM with known pi: calibrating from the true Sigma recovers pi and mu."""
    rng = np.random.default_rng(4)
    mu, Sigma, pi = np.array([0.0, -1.0, -2.0]), np.diag([1.0, 1.5, 1.0]), np.array([0.3, 0.0, 0.5])
    X = sample_zilnm(mu, Sigma, pi, 500, 40000, rng)
    m4, pi4, e_prev, e_mean = calibrate_zilnm(X, mu, Sigma, rng, n_sim=40000)
    assert e_prev < 1e-3 and e_mean < 1e-3
    assert np.allclose(pi4, pi, atol=0.03) and np.allclose(m4, mu, atol=0.1)
