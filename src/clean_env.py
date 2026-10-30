"""Clean environment.csv. Every change is logged as (decision, reason, rows_affected)."""
import numpy as np
import pandas as pd

UNUSABLE_FRAC = 0.90  # drop a column if >= 90% of its values are Unknown/empty


class Log(list):
    def add(self, decision, reason, n):
        self.append({"decision": decision, "reason": reason, "rows_affected": int(n)})

    def table(self):
        return pd.DataFrame(self)


def _missing(s: pd.Series) -> pd.Series:
    return s.isna() | s.astype(str).str.strip().str.lower().isin(["", "unknown"])


def clean_env(path):
    """Return (clean DataFrame indexed by sample ID as str, decision Log, missingness Series of RAW file)."""
    raw = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    log = Log()
    raw = raw.rename(columns={"Ethnicity_Catagory": "Ethnicity_Category",
                              "Wash_Freq of BB (# times per week)": "Wash_Freq"})
    df = raw.copy()
    df["ID"] = df["ID"].str.strip()
    assert df["ID"].is_unique

    # raw missingness per column (before any cleaning), reported as-is
    miss_raw = pd.Series({c: _missing(raw[c]).mean() for c in raw.columns if c != "ID"})

    # --- trim whitespace everywhere ---
    n_ws = sum((df[c] != df[c].str.strip()).sum() for c in df.columns)
    df = df.apply(lambda s: s.str.strip())
    log.add("strip leading/trailing whitespace in all columns", "e.g. 'Innie ' vs 'Innie', 'NC ' vs 'NC'", n_ws)

    # --- Unknown -> NaN (kept as NaN, never imputed) ---
    for c in df.columns:
        m = _missing(df[c])
        df.loc[m, c] = np.nan
    log.add("'Unknown'/empty -> NaN, no imputation", "missing stays missing; models use complete cases and report n",
            int(sum(_missing(raw[c]).sum() for c in raw.columns if c != "ID")))

    # --- navel type ---
    before = df["Innie_or_outie"].copy()
    df["Navel"] = df["Innie_or_outie"].str.capitalize()
    log.add("Innie_or_outie: case-fold to {Innie, Outie}", "'innie'/'Innie '/'Innie' are one label",
            (before.dropna() != df["Navel"].dropna()).sum())

    # --- numeric columns ---
    df["Age"] = pd.to_numeric(df["Age"], errors="raise")
    df["Wash_Freq"] = pd.to_numeric(df["Wash_Freq"].str.replace(",", ".", regex=False), errors="raise")
    log.add("Wash_Freq: decimal comma -> point, to float", "'3,5' -> 3.5", raw["Wash_Freq"].str.contains(",").sum())

    # Wash_Range disagrees with Wash_Freq in a few rows: Wash_Freq (numeric) is the reference.
    rng = df["Wash_Range"].str.lower().replace({"none": "none"})
    expected = pd.cut(df["Wash_Freq"], [-0.1, 0.0, 6.99, 7.0, np.inf],
                      labels=["none", "less than 7", "7 days a week", "more than 7"]).astype(object)
    bad = (rng.notna() & expected.notna() & (rng != expected))
    log.add("Wash_Range inconsistent with Wash_Freq: keep Wash_Freq, drop Wash_Range from models",
            "e.g. range 'none' with freq 0.25; 'None'/'none' also differ only by case (0 = never washes, a real value, not missing)",
            bad.sum())
    df["Wash_Group"] = expected  # derived from the numeric value only

    # --- ethnicity: merge tiny groups (n<10) so a model can use it ---
    eth = df["Ethnicity_Category"].replace({"Latino": "Hispanic"})
    eth = eth.replace({"Black": "Other", "non-Hispanic Black": "Other", "Asian": "Other", "Indian": "Other"})
    df["Ethnicity"] = eth
    log.add("Ethnicity: Latino->Hispanic; Black/Asian/Indian -> 'Other'",
            "groups with n<10 cannot be estimated; Ethnicity_Quest (free text, 22 labels) not used",
            (df["Ethnicity_Category"] != eth).sum())

    # --- geography ---
    state = df["Current State"]
    df["In_NC"] = np.where(state.isna(), np.nan, (state == "NC").astype(float))
    log.add("Geography: keep only In_NC (is current state NC) from Current State; drop ZIP/city/country-by-age",
            "ZIPs are >=70% missing and partly corrupt ('1130 ??'); 'Austria'/'France' sit in the State column",
            state.notna().sum())

    drop = [c for c in raw.columns if c != "ID" and miss_raw[c] >= UNUSABLE_FRAC]
    log.add(f"Drop columns with >= {UNUSABLE_FRAC:.0%} missing: {drop}", "no information to model", len(df))
    log.add("Drop Reads", "constant 1000 for all subjects (rarefied)", len(df))

    keep = ["Gender", "Age", "Wash_Freq", "Wash_Group", "Navel", "Ethnicity", "In_NC"]
    clean = df.set_index("ID")[keep]
    return clean, log, miss_raw
