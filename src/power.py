"""Phase 4: seasonal scenario simulator and detection tools (frozen plan: report/analysis_plan_phase4.md).

Column convention of the D = 22 matrix (counts_top20.csv): columns 0-19 = the 20 named taxa, 20 = Unassigned,
21 = Other (ALR reference). Distances and CLR use the 21 parts without Unassigned (as in 02_eda / 04_validation).
"""
import numpy as np
from scipy import stats
from scipy.spatial.distance import pdist, squareform

from .eda import clr
from .models import _zilnm_probs

PHI = 7                         # peak month of the driver (July)
NAMED = np.arange(20)           # the 20 named taxa (candidates for a seasonal effect)
UNASSIGNED = 20
PARTS = np.r_[NAMED, 21]        # parts used for distances / CLR: named taxa + Other
SEASON = {12: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2, 9: 3, 10: 3, 11: 3}   # DJF, MAM, JJA, SON


def driver(months, phi=PHI):
    """s(t) = cos(2 pi (month - phi) / 12): +1 at the peak month, -1 six months later."""
    return np.cos(2 * np.pi * (np.asarray(months) - phi) / 12)


def balanced_months(n, rng):
    """Every month equally often (n must be a multiple of 12), in random order."""
    assert n % 12 == 0, "n must be a multiple of 12"
    m = np.tile(np.arange(1, 13), n // 12)
    rng.shuffle(m)
    return m


def simulate_seasonal(mu, Sigma, n, A, K, rng, depth=1000, od=1.0, pi=None):
    """One simulated study.

    z_i ~ N(mu, od * Sigma) (ALR coordinates, reference Other);  for K randomly chosen named taxa
    z_ij += A * s(month_i) * sign_j  with signs +1, -1, +1, ... (some taxa peak in summer, some in winter);
    pi given -> zero-inflated M4 (presence is NOT seasonal: simplification); reads ~ Multinomial(depth, p_i).
    Returns counts (n, 22), months (n,), indices of the affected taxa (K,).
    """
    months = balanced_months(n, rng)
    taxa = rng.choice(NAMED, K, replace=False)
    signs = np.resize([1.0, -1.0], K)
    Z = rng.multivariate_normal(mu, od * Sigma, size=n)
    Z[:, taxa] += A * driver(months)[:, None] * signs
    present = np.ones(Z.shape, bool) if pi is None else rng.random(Z.shape) >= pi
    return rng.multinomial(depth, _zilnm_probs(Z, present)), months, taxa


# ---------------------------------------------------------------- PERMANOVA (fast, many permutations at once)
def gower_centre(D):
    """G = -1/2 J D^2 J (J = centring matrix). tr(G) = total sum of squares of PERMANOVA."""
    A = -0.5 * D ** 2
    return A - A.mean(0) - A.mean(1)[:, None] + A.mean()


def pseudo_F(G, labels, g):
    """labels: (P, n) integer group labels (one row per permutation). Returns (P,) pseudo-F.
    SS_between = tr(H G) with H the projection on group means = sum_g (1/n_g) 1_g' G 1_g."""
    M = np.eye(g)[labels]                                      # (P, n, g) one-hot
    GM = np.matmul(G, M)                                       # (P, n, g)
    ss_b = ((M * GM).sum(1) / M.sum(1)).sum(1)
    ss_t = np.trace(G)
    n = G.shape[0]
    return (ss_b / (g - 1)) / ((ss_t - ss_b) / (n - g))


def permanova_p(D, labels, n_perm, rng):
    """Permutation p-value (1 + #{F_perm >= F_obs}) / (1 + n_perm). Exact level for alpha when alpha*(n_perm+1) is an integer."""
    labels = np.asarray(labels)
    g = labels.max() + 1
    L = np.vstack([labels, np.array([rng.permutation(labels) for _ in range(n_perm)])])
    F = pseudo_F(gower_centre(D), L, g)
    return (1 + (F[1:] >= F[0]).sum()) / (1 + n_perm), F[0]


# ---------------------------------------------------------------- cosinor regression on CLR
def cosinor_p(Y, months):
    """Per column of Y (n, q): OLS  y = b0 + b1 cos(2 pi m/12) + b2 sin(2 pi m/12);  F-test of b1 = b2 = 0 (2 df).
    The phase is not assumed known: cos + sin together fit any phase."""
    t = 2 * np.pi * np.asarray(months) / 12
    Xd = np.c_[np.ones_like(t), np.cos(t), np.sin(t)]
    B = np.linalg.lstsq(Xd, Y, rcond=None)[0]
    rss1 = ((Y - Xd @ B) ** 2).sum(0)
    rss0 = ((Y - Y.mean(0)) ** 2).sum(0)
    n = len(Y)
    F = ((rss0 - rss1) / 2) / (rss1 / (n - 3))
    return stats.f.sf(F, 2, n - 3)


def run_point(spec):
    """All R replicates of one grid point (top-level function so it can run in a worker process).
    spec: dict with point_id, seed, R, mu, Sigma, pi (None for M3), n, A, K, depth, od, shuffle."""
    out = []
    for r in range(spec["R"]):
        rng = np.random.default_rng([spec["seed"], 5, spec["point_id"], r])
        X, m, taxa = simulate_seasonal(spec["mu"], spec["Sigma"], spec["n"], spec["A"], spec["K"], rng,
                                       depth=spec["depth"], od=spec["od"], pi=spec["pi"])
        out.append({"point_id": spec["point_id"], "rep": r, **run_tests(X, m, taxa, rng, shuffle=spec["shuffle"])})
    return out


def run_tests(X, months, taxa, rng, n_perm=199, q=0.05, pseudo=0.5, shuffle=False):
    """All pre-specified tests on one simulated study. shuffle=True -> negative control (labels permuted first)."""
    months = rng.permutation(months) if shuffle else np.asarray(months)
    seasons = np.array([SEASON[m] for m in months])
    T = X[:, PARTS].astype(float)
    C = clr(T, pseudo)
    p_bc, _ = permanova_p(squareform(pdist(T, "braycurtis")), seasons, n_perm, rng)
    p_ait, _ = permanova_p(squareform(pdist(C, "euclidean")), seasons, n_perm, rng)
    p_tax = cosinor_p(C[:, :len(NAMED)], months)              # CLR of the 20 named taxa
    rej = stats.false_discovery_control(p_tax) <= q            # Benjamini-Hochberg at q
    truth = np.isin(NAMED, taxa)
    return {"p_permanova_bc": p_bc, "p_permanova_ait": p_ait, "cosinor_any": bool(rej.any()),
            "tp": int((rej & truth).sum()), "fp": int((rej & ~truth).sum()), "n_rej": int(rej.sum())}
