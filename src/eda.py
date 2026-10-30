"""Small, explainable helpers for diversity, compositional transforms and distances.

Counts matrices are samples x taxa (rows = samples) everywhere in src/.
"""
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix


def richness(X):
    """Observed taxa = number of non-zero entries per sample."""
    return (np.asarray(X) > 0).sum(axis=1)


def shannon(X):
    """H = -sum p ln p over non-zero p (natural log)."""
    X = np.asarray(X, float)
    p = X / X.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        return -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)


def simpson(X):
    """Gini-Simpson 1 - sum p^2: probability that two random reads belong to different taxa."""
    X = np.asarray(X, float)
    p = X / X.sum(axis=1, keepdims=True)
    return 1 - (p ** 2).sum(axis=1)


def clr(X, pseudo=0.5):
    """Centred log-ratio: log(x+pseudo) minus the row mean of the logs. Rows then sum to 0.

    simplification: one additive pseudocount (default 0.5 reads) replaces zeros; sensitivity is checked in 02_eda.
    """
    L = np.log(np.asarray(X, float) + pseudo)
    return L - L.mean(axis=1, keepdims=True)


def bray_curtis(X, ids):
    """Bray-Curtis on raw counts (depth is a constant 1000, so no renormalisation is needed)."""
    return DistanceMatrix(squareform(pdist(np.asarray(X, float), "braycurtis")), ids=list(ids))


def aitchison(X, ids, pseudo=0.5):
    """Aitchison distance = Euclidean distance between CLR vectors."""
    return DistanceMatrix(squareform(pdist(clr(X, pseudo), "euclidean")), ids=list(ids))
