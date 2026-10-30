"""Sanity checks for the parser. Run: pytest -q"""
import numpy as np
import pandas as pd

from src.config import RAW, READS
from src.parse import aggregate, parse_decimal, read_taxonomy, resolved, taxon_label, to_counts


def test_parse_decimal():
    assert parse_decimal("0,003") == 0.003
    assert parse_decimal("") == 0.0
    assert parse_decimal(" 1 ") == 1.0


def test_resolved_and_label():
    assert not resolved("Other") and not resolved("g__") and not resolved("")
    assert resolved("g__Staphylococcus")
    r = pd.Series(dict(Kingdom="k__Bacteria", Phylum="p__Firmicutes", Class="c__Bacilli", Order="o__Bacillales",
                       Family="f__Staphylococcaceae", Genus="g__Staphylococcus"))
    assert taxon_label(r) == "Staphylococcaceae;Staphylococcus"
    r["Genus"] = "g__"
    assert taxon_label(r) == "Staphylococcaceae (unresolved below family)"
    assert taxon_label(pd.Series(dict(Kingdom="Unclassified", Phylum="Other", Class="Other",
                                      Order="Other", Family="Other", Genus="Other"))) == "Unclassified"


def test_real_file():
    frac, tax, junk = read_taxonomy(RAW / "taxonomy.csv")
    assert frac.shape == (649, 273) and len(tax) == 649           # table ends after 649 taxa
    assert junk["n_rows_after_table"] > 0                          # junk was found and cut off
    counts = to_counts(frac)                                       # also asserts values are k/1000
    assert (counts.values >= 0).all() and counts.sum().max() <= READS
    # richness row in the junk block equals observed taxa per sample (explains the 'median 31' column sum)
    rich_row = next(r for r in junk["nonblank_rows"] if r[7:10] != [""] * 3 and r[0] == "")
    assert np.array_equal(np.array([int(x) for x in rich_row[7:]]), (counts > 0).sum().values)
    # aggregation conserves total counts and removes duplicate labels
    agg, _ = aggregate(counts, tax)
    assert agg.values.sum() == counts.values.sum() and agg.index.is_unique and len(agg) < 649
