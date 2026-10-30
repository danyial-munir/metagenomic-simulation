"""Global settings: one seed, fixed paths. Every script/notebook imports from here."""
from pathlib import Path

SEED = 20261030  # single global seed; derive child RNGs with np.random.default_rng([SEED, k])
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
READS = 1000  # rarefaction depth (environment.csv: Reads == 1000 for every subject)
