"""Parsing and tidying of taxonomy.csv (taxa x samples) and environment.csv.

Facts established by the audit (see notebooks/01_audit.ipynb):
  * real table = rows up to the first row with an empty Kingdom cell; everything after is junk
    (richness row, repeated header, a side table, a stray note)
  * values are decimal-comma fractions of READS=1000 reads -> exact integer counts
  * 'Other' and an empty 'x__' both mean "not resolved at this rank"
"""
import csv

import numpy as np
import pandas as pd

from .config import READS

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
N_TAXCOLS = 7  # 6 ranks + Family_Genus


def parse_decimal(s: str) -> float:
    """'0,003' -> 0.003 ; '' -> 0.0 (empty cell = taxon absent)."""
    s = s.strip()
    return float(s.replace(",", ".")) if s else 0.0


def read_taxonomy(path):
    """Return (frac, tax, junk) from taxonomy.csv.

    frac : DataFrame taxa x samples, fractions (decimal comma parsed)
    tax  : DataFrame, one row per taxon, rank columns (index aligned with frac)
    junk : dict with the trailing rows, kept only as audit evidence
    """
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    header, body = rows[0], rows[1:]
    # real table ends at the first row with an empty Kingdom cell (the richness row follows directly, no blank line)
    end = next(i for i, r in enumerate(body) if not r[0].strip())
    table, rest = body[:end], body[end:]
    ids = header[N_TAXCOLS:]
    assert len(set(ids)) == len(ids), "duplicate sample IDs"
    assert all(len(r) == len(header) for r in table), "ragged rows in real table"
    assert all(r[0].startswith(("k__", "Unclassified")) for r in table), "non-taxon row inside table"

    frac = pd.DataFrame([[parse_decimal(c) for c in r[N_TAXCOLS:]] for r in table], columns=ids)
    tax = pd.DataFrame([r[:N_TAXCOLS] for r in table], columns=RANKS + ["Family_Genus"])
    nonblank = [r for r in rest if any(c.strip() for c in r)]
    junk = {"n_rows_after_table": len(rest), "nonblank_rows": nonblank, "table_end_line": end + 1}
    return frac, tax, junk


def to_counts(frac: pd.DataFrame, reads: int = READS) -> pd.DataFrame:
    """Fractions -> integer counts; asserts that the data really are k/reads."""
    c = frac * reads
    assert np.allclose(c, c.round(), atol=1e-6), "values are not multiples of 1/reads"
    return c.round().astype(int)


def resolved(v: str) -> bool:
    """A rank value carries information unless it is 'Other', empty or a bare 'x__'."""
    v = v.strip()
    return v not in ("", "Other", "Unclassified") and not (len(v) == 3 and v.endswith("__"))


def taxon_label(r: pd.Series) -> str:
    """Deepest resolved rank -> label. Genus resolved -> 'Family;Genus' (genus names repeat across families)."""
    if resolved(r["Genus"]):
        g = r["Genus"][3:]
        return f"{r['Family'][3:]};{g}" if resolved(r["Family"]) else g
    for rank in ["Family", "Order", "Class", "Phylum", "Kingdom"]:
        if resolved(r[rank]):
            return f"{r[rank][3:]} (unresolved below {rank.lower()})"
    return "Unclassified"


def aggregate(counts: pd.DataFrame, tax: pd.DataFrame):
    """Sum rows sharing the same label. Returns (counts_agg, label_table)."""
    labels = tax.apply(taxon_label, axis=1)
    agg = counts.groupby(labels.values).sum()
    # deepest resolved rank per label, for the audit
    lab_tab = pd.DataFrame({"label": labels.values, "genus_resolved": tax["Genus"].map(resolved).values})
    return agg, lab_tab
