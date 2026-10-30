"""Generative models for compositional count data (samples x categories, every row sums to N).

M0  multinomial, one common probability vector p
M1  Dirichlet-multinomial (DM): p_i ~ Dirichlet(alpha), x_i ~ Multinomial(N, p_i)
M2  DM regression: log alpha_ij = b_j + z_i . beta_j for the top-K taxa, beta_j = 0 for the rest
M3  logistic-normal multinomial (LNM): z_i ~ Normal(mu, Sigma) in ALR coordinates, p_i = softmax([z_i, 0])

Every model has fit_*, loglik_* (one value per sample, the full log-probability including the multinomial
coefficient, so values are comparable across models) and sample_*.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import digamma, gammaln, logsumexp


def log_mult_coef(X):
    """log N! - sum_j log x_j!  (identical for every model; kept so log-likelihoods are true probabilities)."""
    X = np.asarray(X, float)
    return gammaln(X.sum(axis=1) + 1) - gammaln(X + 1).sum(axis=1)


# ---------------------------------------------------------------- M0: multinomial
def fit_multinomial(X):
    """MLE: p_j = total reads of j / total reads."""
    X = np.asarray(X, float)
    return X.sum(axis=0) / X.sum()


def loglik_multinomial(X, p):
    X = np.asarray(X, float)
    assert np.all(p[X.sum(axis=0) > 0] > 0), "a category unseen in training appears in test data"
    with np.errstate(divide="ignore"):
        return log_mult_coef(X) + np.where(X > 0, X * np.log(p), 0.0).sum(axis=1)


def sample_multinomial(p, N, n, rng):
    return rng.multinomial(N, p, size=n)


# ---------------------------------------------------------------- M1: Dirichlet-multinomial
def loglik_dm(X, alpha):
    """alpha: (D,) shared or (n, D) per sample.
    log P = coef + lgamma(A) - lgamma(N + A) + sum_j [lgamma(x_j + a_j) - lgamma(a_j)],  A = sum_j a_j."""
    X = np.asarray(X, float)
    alpha = np.broadcast_to(alpha, X.shape)
    A, N = alpha.sum(axis=1), X.sum(axis=1)
    return log_mult_coef(X) + gammaln(A) - gammaln(N + A) + (gammaln(X + alpha) - gammaln(alpha)).sum(axis=1)


def _dm_grad_alpha(X, alpha):
    """d loglik_i / d alpha_ij = psi(A_i) - psi(N_i + A_i) + psi(x_ij + a_ij) - psi(a_ij)."""
    A, N = alpha.sum(axis=1, keepdims=True), X.sum(axis=1, keepdims=True)
    return digamma(A) - digamma(N + A) + digamma(X + alpha) - digamma(alpha)


def fit_dm(X):
    """MLE of alpha by L-BFGS on eta = log(alpha) (keeps alpha > 0). Returns alpha (D,)."""
    X = np.asarray(X, float)
    p0 = (X.sum(axis=0) + 0.5) / (X.sum() + 0.5 * X.shape[1])
    eta0 = np.log(p0 * 1.0)                                  # start: mean = observed proportions, precision A = 1

    def nll(eta):
        a = np.broadcast_to(np.exp(eta), X.shape)
        g = (_dm_grad_alpha(X, a) * a).sum(axis=0)           # chain rule: d/d eta = alpha * d/d alpha
        return -loglik_dm(X, a[0]).sum(), -g

    res = minimize(nll, eta0, jac=True, method="L-BFGS-B", options={"maxiter": 5000})
    assert res.success, res.message
    return np.exp(res.x)


def dm_moments(alpha, N):
    """Theory: E x_j = N p_j ; Var x_j = N p_j (1 - p_j) (N + A) / (1 + A), with p = alpha / A."""
    A = alpha.sum()
    p = alpha / A
    return N * p, N * p * (1 - p) * (N + A) / (1 + A)


def sample_dm(alpha, N, n, rng):
    """alpha (D,) or (n, D). Draw p_i ~ Dirichlet(alpha_i) via normalised Gammas, then x_i ~ Multinomial(N, p_i)."""
    alpha = np.broadcast_to(alpha, (n, len(alpha) if np.ndim(alpha) == 1 else alpha.shape[1]))
    G = rng.gamma(alpha)
    P = G / G.sum(axis=1, keepdims=True)
    return np.array([rng.multinomial(N, p) for p in P])


# ---------------------------------------------------------------- M2: DM regression
def _m2_alpha(theta, Z, D, top):
    k, q = len(top), Z.shape[1]
    b, B = theta[:D], theta[D:].reshape(q, k)
    eta = np.tile(b, (Z.shape[0], 1))
    eta[:, top] += Z @ B
    return np.exp(eta)


def fit_dm_reg(X, Z, top, b0=None):
    """MLE of DM regression. Z: (n, q) covariates (standardised / 0-1), top: indices of the K taxa that get slopes.

    simplification: only the top-K taxa get covariate effects; the remaining ones keep a constant alpha, so the number of
    parameters is D + K*q instead of D*(q+1). Returns theta = [b (D,), beta (q*K,)].
    """
    X, Z = np.asarray(X, float), np.asarray(Z, float)
    D, q, k = X.shape[1], Z.shape[1], len(top)
    theta0 = np.r_[np.log(fit_dm(X)) if b0 is None else b0, np.zeros(q * k)]

    def nll(theta):
        a = _m2_alpha(theta, Z, D, top)
        G = _dm_grad_alpha(X, a) * a                         # d loglik_i / d eta_ij
        return -loglik_dm(X, a).sum(), -np.r_[G.sum(axis=0), (Z.T @ G[:, top]).ravel()]

    res = minimize(nll, theta0, jac=True, method="L-BFGS-B", options={"maxiter": 10000})
    assert res.success, res.message
    return res.x


def loglik_dm_reg(X, Z, theta, top):
    return loglik_dm(X, _m2_alpha(theta, np.asarray(Z, float), np.asarray(X).shape[1], top))


def sample_dm_reg(theta, Z, top, N, rng):
    a = _m2_alpha(theta, np.asarray(Z, float), len(theta) - Z.shape[1] * len(top), top)
    return sample_dm(a, N, len(a), rng)


# ---------------------------------------------------------------- M3: logistic-normal multinomial
# Coordinates: z = ALR = log(p_j / p_ref) for the D-1 non-reference categories; the reference column is the LAST column
# of X (the caller reorders). p = softmax([z, 0]).

def _laplace_modes(X, mu, P, Z0, iters=50, tol=1e-8):
    """Per-sample mode of  f(z) = x.z - N log(1 + sum exp z) - 0.5 (z-mu)' P (z-mu)  by damped Newton (f is concave).
    Returns modes Z (n, d), negative Hessians H (n, d, d), and f at the modes."""
    Xn, N = X[:, :-1], X.sum(axis=1)
    Z = Z0.copy()

    def f(Z):
        lse = logsumexp(np.c_[Z, np.zeros(len(Z))], axis=1)
        r = Z - mu
        return (Xn * Z).sum(1) - N * lse - 0.5 * np.einsum("ij,jk,ik->i", r, P, r)

    fz = f(Z)
    for _ in range(iters):
        p = np.exp(Z - logsumexp(np.c_[Z, np.zeros(len(Z))], axis=1, keepdims=True))
        g = Xn - N[:, None] * p - (Z - mu) @ P
        H = N[:, None, None] * (np.einsum("ij,jk->ijk", p, np.eye(p.shape[1])) - p[:, :, None] * p[:, None, :]) + P
        step = np.linalg.solve(H, g[..., None])[..., 0]
        t = np.ones(len(Z))
        for _ in range(30):                                  # backtracking: accept only steps that increase f
            Znew = Z + t[:, None] * step
            fnew = f(Znew)
            bad = fnew < fz - 1e-12
            if not bad.any():
                break
            t[bad] *= 0.5
        done = np.abs(fnew - fz).max() < tol
        Z, fz = Znew, fnew
        if done:
            break
    p = np.exp(Z - logsumexp(np.c_[Z, np.zeros(len(Z))], axis=1, keepdims=True))
    H = N[:, None, None] * (np.einsum("ij,jk->ijk", p, np.eye(p.shape[1])) - p[:, :, None] * p[:, None, :]) + P
    return Z, H, fz


def loglik_lnm(X, mu, Sigma, Z0=None, return_modes=False):
    """Laplace approximation of log p(x) = log integral Mult(x | softmax z) N(z | mu, Sigma) dz:
       log p(x) ~ coef + f(z_hat) - 0.5 log det Sigma - 0.5 log det H(z_hat)   (the (2 pi)^(d/2) factors cancel).
    """
    X = np.asarray(X, float)
    P = np.linalg.inv(Sigma)
    if Z0 is None:
        Z0 = np.log((X[:, :-1] + 0.5) / (X[:, -1:] + 0.5))
    Z, H, fz = _laplace_modes(X, mu, P, Z0)
    ll = log_mult_coef(X) + fz - 0.5 * np.linalg.slogdet(Sigma)[1] - 0.5 * np.linalg.slogdet(H)[1]
    return (ll, Z, H) if return_modes else ll


def fit_lnm(X, shrink=0.5, max_iter=300, tol=1e-4):
    """Approximate EM (Laplace E-step) for mu, Sigma.

    E-step: per-sample Gaussian approximation N(z_hat_i, H_i^-1) of the posterior of z_i.
    M-step: mu = mean z_hat ; S = mean[(z_hat - mu)(z_hat - mu)' + H_i^-1] ; Sigma = (1-shrink) S + shrink diag(S).
    simplification: Laplace instead of exact (Monte-Carlo) EM; shrinkage towards the diagonal instead of a low-rank
    factor model. Shrinkage makes Sigma well conditioned with 113 dimensions and 273 samples; chosen by cross-validation.
    Returns mu, Sigma, trace of mean approx log-lik per iteration.
    """
    X = np.asarray(X, float)
    Zh = np.log((X[:, :-1] + 0.5) / (X[:, -1:] + 0.5))       # start from the pseudocount ALR
    mu = Zh.mean(0)
    S = np.cov(Zh.T)
    Sigma = (1 - shrink) * S + shrink * np.diag(np.diag(S))
    trace = []
    for _ in range(max_iter):
        ll, Zh, H = loglik_lnm(X, mu, Sigma, Z0=Zh, return_modes=True)
        trace.append(ll.mean())
        V = np.linalg.inv(H)
        mu = Zh.mean(0)
        R = Zh - mu
        S = (R.T @ R) / len(X) + V.mean(0)
        Sigma = (1 - shrink) * S + shrink * np.diag(np.diag(S))
        if len(trace) > 1 and abs(trace[-1] - trace[-2]) < tol:
            break
    return mu, Sigma, np.array(trace)


def loglik_lnm_is(x, mu, Sigma, n_draws, rng, scale=1.2):
    """Importance-sampling check of the Laplace value for ONE sample. Proposal = Gaussian at the Laplace mode with
    covariance scale^2 * H^-1. Returns (log p(x) estimate, effective sample size)."""
    x = np.asarray(x, float)[None, :]
    _, Zh, H = loglik_lnm(x, mu, Sigma, return_modes=True)
    C = np.linalg.inv(H[0]) * scale ** 2
    L = np.linalg.cholesky(C)
    d = len(mu)
    E = rng.standard_normal((n_draws, d))
    Zs = Zh[0] + E @ L.T
    log_q = -0.5 * (E ** 2).sum(1) - np.log(np.diag(L)).sum() - 0.5 * d * np.log(2 * np.pi)
    lse = logsumexp(np.c_[Zs, np.zeros(n_draws)], axis=1)
    loglik = log_mult_coef(x)[0] + (x[0, :-1] * Zs).sum(1) - x[0].sum() * lse
    Pinv = np.linalg.inv(Sigma)
    R = Zs - mu
    log_prior = -0.5 * np.einsum("ij,jk,ik->i", R, Pinv, R) - 0.5 * np.linalg.slogdet(Sigma)[1] - 0.5 * d * np.log(2 * np.pi)
    w = loglik + log_prior - log_q
    ess = np.exp(2 * logsumexp(w) - logsumexp(2 * w))
    return logsumexp(w) - np.log(n_draws), ess


def sample_lnm(mu, Sigma, N, n, rng):
    Z = rng.multivariate_normal(mu, Sigma, size=n)
    P = np.exp(np.c_[Z, np.zeros(n)] - logsumexp(np.c_[Z, np.zeros(n)], axis=1, keepdims=True))
    return np.array([rng.multinomial(N, p) for p in P])
