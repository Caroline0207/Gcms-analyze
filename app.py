import streamlit as st
import pandas as pd
import numpy as np
import io
import re

st.set_page_config(layout="wide", page_title="GC-MS Multi-Trial Analysis")
st.title("GC-MS Multi-Trial Analysis 🧪")
st.caption(
    "Paste cleaned trial tables (TSV copied from Excel). The app matches compounds across trials "
    "using alias-aware Name/Species logic, keeps only compounds present in ALL trials, and outputs "
    "AVG AREA and CUMULATIVE % with validation."
)

# ------------------ helpers ------------------
def normalize_text(s: str) -> str:
    """Aggressive-but-safe normalization for matching."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ""
    s = str(s).lower().strip()

    # unify unicode dashes
    s = s.replace("–", "-").replace("—", "-")

    # remove common punctuation that causes mismatches
    s = re.sub(r"[(),;]", " ", s)

    # keep hyphens as separators (turn into spaces)
    s = s.replace("-", " ")

    # collapse whitespace
    s = re.sub(r"\s+", " ", s).strip()
    return s


def extract_aliases(name, species):
    """
    Build alias candidates from Name and Species.
    Handles NaN safely and splits "or" phrases.
    """
    aliases = []
    for x in [name, species]:
        if x is None or pd.isna(x):
            continue
        x = str(x).strip()
        if not x or x.lower() == "nan":
            continue

        x = x.lower()

        # Split only as a word-ish separator.
        # Example: "or o-Cymene" / "X or Y" / "X, or Y"
        x = re.sub(r"\s+or\s+", "|", x)
        x = re.sub(r"^\s*or\s+", "", x)

        parts = [p.strip() for p in x.split("|") if p.strip()]
        aliases.extend(parts)

    return aliases


def is_trivial_name(s: str) -> bool:
    """
    Heuristic: short/common names (e.g., o-cymene) preferred over long structural descriptors.
    """
    if not s:
        return False
    # After normalization, commas/parens removed; still guard against long IUPAC-ish descriptors.
    bad_tokens = {"benzene", "methyl", "ethyl", "propyl", "butyl", "hexenyl", "dimethyl", "trimethyl"}
    toks = set(s.split())
    return (len(s) <= 20) and (len(toks & bad_tokens) == 0)


def canonical_key(name, species) -> str:
    aliases = extract_aliases(name, species)
    aliases = [normalize_text(a) for a in aliases if a]
    aliases = [a for a in aliases if a]

    if not aliases:
        return ""  # will be excluded naturally

    trivial = [a for a in aliases if is_trivial_name(a)]
    if trivial:
        return trivial[0]

    # fallback: choose the shortest normalized alias
    return min(aliases, key=len)


def parse_trial(text: str) -> pd.DataFrame:
    """
    Parse TSV copied from Excel. Requires columns: RT, Area, Area %, Score, Name
    Species optional.
    """
    df = pd.read_csv(io.StringIO(text.strip()), sep="\t", dtype=str, engine="python")
    df.columns = [c.strip() for c in df.columns]

    required = ["RT", "Area", "Area %", "Score", "Name"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in a trial: {missing}. Found: {list(df.columns)}")

    # Ensure optional cols exist
    if "Species" not in df.columns:
        df["Species"] = ""

    # Clean strings (avoid NaN -> float issues)
    df["Name"] = df["Name"].fillna("").astype(str).str.strip()
    df["Species"] = df["Species"].fillna("").astype(str).str.strip()
    if "Formula" in df.columns:
        df["Formula"] = df["Formula"].fillna("").astype(str).str.strip()
    else:
        df["Formula"] = ""

    # Numeric conversions
    def to_num(s):
        return pd.to_numeric(s.astype(str).str.replace(",", "", regex=False).str.strip(), errors="coerce")

    df["RT"] = to_num(df["RT"])
    df["Area"] = to_num(df["Area"])
    df["Area %"] = to_num(df["Area %"])
    df["Score"] = to_num(df["Score"])

    # Key
    df["key"] = df.apply(lambda r: canonical_key(r["Name"], r["Species"]), axis=1)

    # Drop rows with no key or no Area %
    df = df[df["key"].astype(str).str.strip() != ""].copy()

    return df


# ------------------ UI ------------------
n_trials = st.selectbox("Number of trials", [2, 3, 4, 5, 6], index=2)
st.write("Paste each trial table below (TSV copied from Excel).")

trial_dfs = {}
for i in range(n_trials):
    st.subheader(f"T{i+1}")
    txt = st.text_area(f"Paste T{i+1} table", height=180, key=f"trial_{i+1}")
    if txt.strip():
        try:
            trial_dfs[f"T{i+1}"] = parse_trial(txt)
            st.success(f"T{i+1} parsed ✅")
        except Exception as e:
            st.error(f"T{i+1} error: {e}")

if len(trial_dfs) != n_trials:
    st.info("Please paste valid data for ALL trials.")
    st.stop()

# ------------------ analysis ------------------
# Intersection of keys across all trials
common_keys = set.intersection(*[set(df["key"]) for df in trial_dfs.values()])

if not common_keys:
    st.warning("No common compounds found across all trials (after alias matching).")
    st.stop()

rows = []
for key in sorted(common_keys):
    # Use T1 as the display anchor (Name/Formula/Species); you can switch to "best name" later if needed
    base = trial_dfs["T1"].loc[trial_dfs["T1"]["key"] == key].iloc[0]

    row = {
        "RT": base["RT"],
        "Name": base["Name"],
        "Formula": base.get("Formula", ""),
        "Species": base.get("Species", ""),
    }

    avg = 0.0
    for t, df in trial_dfs.items():
        r = df[df["key"] == key].iloc[0]
        row[f"{t} RT"] = r["RT"]
        row[f"{t} Area"] = r["Area"]
        row[f"{t} Area %"] = r["Area %"]
        row[f"{t} Score"] = r["Score"]
        avg += (r["Area %"] if pd.notna(r["Area %"]) else 0.0)

    row["AVG AREA"] = avg / n_trials
    rows.append(row)

out = pd.DataFrame(rows)

# Sort like your example: increasing AVG AREA, then cumulative
out = out.sort_values("AVG AREA", ascending=True, na_position="last").reset_index(drop=True)
out["CUMULATIVE %"] = out["AVG AREA"].cumsum()

# ------------------ validation ------------------
st.subheader("Validation (Area % sums)")
val_cols = st.columns(n_trials + 1)
for idx, (t, df) in enumerate(trial_dfs.items()):
    val_cols[idx].metric(f"{t} Area % sum", f"{df['Area %'].sum(skipna=True):.2f}%")

val_cols[-1].metric("AVG AREA sum", f"{out['AVG AREA'].sum(skipna=True):.2f}%")

# ------------------ output ------------------
st.subheader("Final Output Table")
st.dataframe(out, use_container_width=True)

csv = out.to_csv(index=False).encode("utf-8")
st.download_button("Download CSV", csv, "gcms_multi_trial_output.csv", mime="text/csv")
