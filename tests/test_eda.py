import numpy as np
from skbio.diversity import alpha_diversity

from src.eda import aitchison, clr, richness, shannon, simpson


def _X():
    rng = np.random.default_rng(1)
    return rng.multinomial(1000, rng.dirichlet(np.ones(20) * 0.3), size=15)


def test_alpha_matches_skbio():
    X = _X()
    assert np.allclose(shannon(X), alpha_diversity("shannon", X).values)
    assert np.allclose(simpson(X), alpha_diversity("simpson", X).values)
    assert (richness(X) == alpha_diversity("sobs", X).values).all()


def test_clr_rows_sum_to_zero_and_scale_invariance():
    X = _X()
    assert np.allclose(clr(X).sum(axis=1), 0)
    # Aitchison geometry is scale-invariant: multiplying a composition by a constant changes nothing (no zeros here)
    Y = X + 1
    ids = range(len(Y))
    assert np.allclose(aitchison(Y, ids, pseudo=0).data, aitchison(3 * Y, ids, pseudo=0).data)
