import streamlit as st
import pandas as pd
import numpy as np
import io
import re

st.set_page_config(layout="wide", page_title="GC-MS Multi-Trial Analysis")
st.title("GC-MS Multi-Trial Analysis 🧪")
st.caption(
    "Paste trial tables (TSV copied from Excel). The app matches compounds across trials using "
    "alias-aware Name/Species logic, keeps only compounds present in ALL trials, then recalculates "
    "Area % within the common-compound set so each trial sums to 100."
)

# ------------------ helpers ------------------
def normalize_text(s: str) -> str:
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ""
    s = str(s).lower().strip()
    s = s.replace("–", "-").replace("—", "-")
    s = re.sub(r"[(),;]", " ", s)
    s = s.replace("-", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s

def extract_aliases(name, species):
    aliases = []
    for x in [name, species]:
        if x is None or pd.isna(x):
            continue
        x = str(x).strip()
        if not x or x.lower() == "nan":
            continue
        x = x.lower()
        # split "or" as separator
        x = re.sub(r"\s+or\s+", "|", x)
        x = re.sub(r"^\s*or\s+", "", x)
        parts = [p.strip() for p in x.split("|") if p.strip()]
        aliases.extend(parts)
    return aliases

def is_trivial_name(s: str) -> bool:
    if not s:
        return False
    bad_tokens = {"benzene", "methyl", "ethyl", "propyl", "butyl", "hexenyl", "dimethyl", "trimethyl"}
    toks = set(s.split())
    return (len(s) <= 20) and (len(toks & bad_tokens) == 0)

def canonical_key(name, species) -> str:
    aliases = extract_aliases(name, species)
    aliases = [normalize_text(a) for a in aliases if a]
    aliases = [a for a in aliases if a]
    if not aliases:
        return ""
    trivial = [a for a in aliases if is_trivial_name(a)]
    return trivial[0] if trivial else min(aliases, key=len)

def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series.astype(str).str.replace(",", "", regex=False).str.strip(), errors="coerce")

def parse_trial(text: str) -> pd.DataFrame:
    df = pd.read_csv(io.StringIO(text.strip()), sep="\t", dtype=str, engine="python")
    df.columns = [c.strip() for c in df.columns]

    # 최소 요구 컬럼: RT, Area, Name (Species optional)
    required = ["RT", "Area", "Name"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Found: {list(df.columns)}")

    if "Species" not in df.columns:
        df["Species"] = ""
    if "Formula" not in df.columns:
        df["Formula"] = ""
    if "Score" not in df.columns:
        df["Score"] = np.nan

    df["Name"] = df["Name"].fillna("").astype(str).str.strip()
    df["Species"] = df["Species"].fillna("").astype(str).str.strip()
    df["Formula"] = df["Formula"].fillna("").astype(str).str.strip()

    df["RT"] = to_num(df["RT"])
    df["Area"] = to_num(df["Area"])
    df["Score"] = to_num(df["Score"])

    # canonical key
    df["key"] = df.apply(lambda r: canonical_key(r["Name"], r["Species"]), axis=1)

    # drop unusable rows
    df = df[df["key"].astype(str).str.strip() != ""].copy()
    df = df[pd.notna(df["Area"]) & (df["Area"] > 0)].copy()

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
common_keys = set.intersection(*[set(df["key"]) for df in trial_dfs.values()])

if not common_keys:
    st.warning("No common compounds found across all trials (after alias matching).")
    st.stop()

# Build a per-trial lookup (key -> row)
trial_maps = {t: df.set_index("key") for t, df in trial_dfs.items()}

# For display fields (Name/Formula/Species/RT): use T1 as anchor
anchor = trial_maps["T1"].loc[list(common_keys)].copy()

rows = []
for key in sorted(common_keys):
    base = trial_maps["T1"].loc[key]

    row = {
        "RT": base["RT"],
        "Name": base["Name"],
        "Formula": base.get("Formula", ""),
        "Species": base.get("Species", ""),
    }

    # compute trial-wise Area% AFTER restricting to common_keys
    # first collect areas for this key across trials and also totals per trial
    for t in trial_maps.keys():
        row[f"{t} RT"] = trial_maps[t].loc[key]["RT"]
        row[f"{t} Area"] = trial_maps[t].loc[key]["Area"]
        row[f"{t} Score"] = trial_maps[t].loc[key].get("Score", np.nan)

    rows.append(row)

out = pd.DataFrame(rows)

# Now compute trial totals on common set and Area%
for t in trial_maps.keys():
    # sum of Area for common compounds in this trial
    areas_common = out[f"{t} Area"].astype(float)
    total_area = areas_common.sum(skipna=True)

    if total_area <= 0 or np.isnan(total_area):
        out[f"{t} Area %"] = np.nan
    else:
        out[f"{t} Area %"] = (areas_common / total_area) * 100

# AVG AREA = mean of trial Area% columns
area_pct_cols = [f"{t} Area %" for t in trial_maps.keys()]
out["AVG AREA"] = out[area_pct_cols].mean(axis=1, skipna=True)

# Sort + cumulative like your example
out = out.sort_values("AVG AREA", ascending=True, na_position="last").reset_index(drop=True)
out["CUMULATIVE %"] = out["AVG AREA"].cumsum()

# ------------------ validation ------------------
st.subheader("Validation (must be ~100%)")
cols = st.columns(n_trials + 1)
for idx, t in enumerate(trial_maps.keys()):
    s = out[f"{t} Area %"].sum(skipna=True)
    cols[idx].metric(f"{t} Area % sum", f"{s:.2f}%")

avg_sum = out["AVG AREA"].sum(skipna=True)
cols[-1].metric("AVG AREA sum", f"{avg_sum:.2f}%")

# ------------------ input overview ------------------
st.subheader("Input Overview (RT-sorted by trial)")

show_overview = st.toggle(
    "Show RT-sorted input tables for each trial",
    value=False,
    help="Displays each pasted trial table sorted by RT for quick sanity checks."
)

if show_overview:
    for t, df in trial_dfs.items():
        st.markdown(f"### {t}")

        # Columns to show (존재하는 것만)
        cols = ["RT", "Name", "Formula", "Species", "Area", "Score"]
        cols = [c for c in cols if c in df.columns]

        df_view = (
            df[cols]
            .sort_values("RT", ascending=True, na_position="last")
            .reset_index(drop=True)
        )

        st.dataframe(
            df_view,
            use_container_width=True,
            hide_index=True
        )
# ------------------ output ------------------
st.subheader("Final Output Table")
st.dataframe(out, use_container_width=True)

csv = out.to_csv(index=False).encode("utf-8")
st.download_button("Download CSV", csv, "gcms_multi_trial_output.csv", mime="text/csv")
