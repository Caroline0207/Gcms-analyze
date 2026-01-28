import streamlit as st
import pandas as pd
import numpy as np
import io
import re

st.set_page_config(layout="wide", page_title="GC-MS Multi-Trial Analysis")
st.title("GC-MS Multi-Trial Analysis 🧪")
st.caption(
    "Paste trial tables (TSV copied from Excel). The app matches compounds across trials using "
    "alias-aware Name/Species logic, keeps only compounds present in ALL trials, recalculates Area % "
    "within the common-compound set (so each trial sums to 100), and outputs AVG AREA and CUMULATIVE %."
)

# =========================
# Helpers: matching / parsing
# =========================
def normalize_text(s: str) -> str:
    """Normalization for matching keys (not for display)."""
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return ""
    s = str(s).lower().strip()
    s = s.replace("–", "-").replace("—", "-")
    s = re.sub(r"[(),;]", " ", s)
    s = s.replace("-", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s

def extract_aliases(name, species):
    """Extract aliases from Name + Species, split 'or' phrases, handle NaN safely."""
    aliases = []
    for x in [name, species]:
        if x is None or pd.isna(x):
            continue
        x = str(x).strip()
        if not x or x.lower() == "nan":
            continue
        x = x.lower()

        # split "or" only as a separator between names
        x = re.sub(r"\s+or\s+", "|", x)
        x = re.sub(r"^\s*or\s+", "", x)

        parts = [p.strip() for p in x.split("|") if p.strip()]
        aliases.extend(parts)
    return aliases

def is_trivial_name(s: str) -> bool:
    """Heuristic for choosing short/common names (e.g., o-cymene) over long descriptors."""
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
    return pd.to_numeric(
        series.astype(str).str.replace(",", "", regex=False).str.strip(),
        errors="coerce",
    )

def parse_trial(text: str) -> pd.DataFrame:
    """Parse TSV (Excel paste). Requires RT, Area, Name. Species/Formula/Score optional."""
    df = pd.read_csv(io.StringIO(text.strip()), sep="\t", dtype=str, engine="python")
    df.columns = [c.strip() for c in df.columns]

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

    df["key"] = df.apply(lambda r: canonical_key(r["Name"], r["Species"]), axis=1)

    # usable rows only
    df = df[df["key"].astype(str).str.strip() != ""].copy()
    df = df[pd.notna(df["Area"]) & (df["Area"] > 0)].copy()

    return df

def first_row_for_key(df: pd.DataFrame, key: str) -> pd.Series:
    """Return the first matching row for a key (no collapsing)."""
    return df.loc[df["key"] == key].iloc[0]

# =========================
# UI: trial inputs
# =========================
n_trials = st.selectbox("Number of trials", [2, 3, 4, 5, 6], index=2)

trial_dfs = {}
dup_warnings = {}

for i in range(n_trials):
    t = f"T{i+1}"
    st.subheader(t)
    txt = st.text_area(f"Paste {t} table (TSV)", height=180, key=f"trial_{t}")
    if txt.strip():
        try:
            df = parse_trial(txt)
            trial_dfs[t] = df

            # detect duplicates (same key appears multiple times)
            dup_counts = df["key"].value_counts()
            dups = dup_counts[dup_counts > 1]
            if len(dups) > 0:
                dup_warnings[t] = dups

            st.success(f"{t} parsed ✅  ({len(df)} rows)")
        except Exception as e:
            st.error(f"{t} error: {e}")

if len(trial_dfs) != n_trials:
    st.info("Please paste valid data for ALL trials.")
    st.stop()

# duplicates warning (no collapse; just warn)
if dup_warnings:
    with st.expander("⚠️ Duplicate key warnings (no collapsing applied)"):
        st.write(
            "Some trials contain multiple rows mapped to the same compound key. "
            "Per your request, the app does NOT collapse them; it uses the FIRST occurrence only. "
            "If this is not desired, fix upstream cleaning or adjust matching rules."
        )
        for t, dups in dup_warnings.items():
            st.write(f"**{t}** duplicates (key → count):")
            st.dataframe(dups.rename("count").reset_index().rename(columns={"index": "key"}), hide_index=True)

# =========================
# Compute common compounds (intersection)
# =========================
common_keys = set.intersection(*[set(df["key"]) for df in trial_dfs.values()])
if not common_keys:
    st.warning("No common compounds found across all trials (after alias matching).")
    st.stop()

# =========================
# RT Overview (All trials combined, RT ascending)
# =========================
st.subheader("RT Overview (All trials combined, RT ascending)")
show_rt_overview = st.toggle(
    "Show combined RT-sorted table (all trials)",
    value=False,
    help="All pasted rows from all trials combined and sorted by RT."
)

if show_rt_overview:
    dfs = []
    for t, df in trial_dfs.items():
        tmp = df.copy()
        tmp["Trial"] = t
        dfs.append(tmp)

    rt_all = pd.concat(dfs, ignore_index=True)
    cols = [c for c in ["Trial", "RT", "Name", "Formula", "Species", "Area", "Score", "key"] if c in rt_all.columns]
    rt_all = rt_all[cols].sort_values("RT", ascending=True, na_position="last").reset_index(drop=True)
    st.dataframe(rt_all, use_container_width=True, hide_index=True)

st.divider()

# =========================
# Build Final Output Table (common keys only)
# =========================
# For display name/formula/species, anchor to T1 first occurrence
anchor_t = "T1"

rows = []
for key in sorted(common_keys):
    # per trial: grab FIRST row for this key
    per_trial = {t: first_row_for_key(df, key) for t, df in trial_dfs.items()}

    # Display RT (leftmost): mean RT across trials (robust + good for sorting)
    rts = [per_trial[t]["RT"] for t in trial_dfs.keys()]
    rt_display = pd.to_numeric(pd.Series(rts), errors="coerce").mean(skipna=True)

    base = per_trial[anchor_t]
    row = {
        "RT": rt_display,
        "Name": base.get("Name", ""),
        "Formula": base.get("Formula", ""),
        "Species": base.get("Species", ""),
    }

    # fill per-trial RT/Area/Score
    for t in trial_dfs.keys():
        row[f"{t} RT"] = per_trial[t]["RT"]
        row[f"{t} Area"] = per_trial[t]["Area"]
        row[f"{t} Score"] = per_trial[t].get("Score", np.nan)

    rows.append(row)

out = pd.DataFrame(rows)

# =========================
# Recompute Area % within common set so each trial sums 100
# =========================
for t in trial_dfs.keys():
    areas = pd.to_numeric(out[f"{t} Area"], errors="coerce")
    total_area = areas.sum(skipna=True)
    out[f"{t} Area %"] = (areas / total_area) * 100 if pd.notna(total_area) and total_area > 0 else np.nan

# AVG AREA = mean of trial Area% columns
area_pct_cols = [f"{t} Area %" for t in trial_dfs.keys()]
out["AVG AREA"] = out[area_pct_cols].mean(axis=1, skipna=True)

# Sort Final Output by RT ascending (your request)
out = out.sort_values("RT", ascending=True, na_position="last").reset_index(drop=True)

# CUMULATIVE % in RT order (based on AVG AREA)
out["CUMULATIVE %"] = out["AVG AREA"].cumsum()

# =========================
# Format columns exactly like your desired table layout
# =========================
final_cols = ["RT", "Name", "Formula", "Species", "AVG AREA", "CUMULATIVE %"]
for t in trial_dfs.keys():
    final_cols += [f"{t} RT", f"{t} Area", f"{t} Area %", f"{t} Score"]

out_final = out[final_cols].copy()

# =========================
# Validation
# =========================
st.subheader("Validation (each trial Area % sum must be 100; AVG AREA sum must be 100)")
cols = st.columns(n_trials + 1)

for idx, t in enumerate(trial_dfs.keys()):
    s = pd.to_numeric(out_final[f"{t} Area %"], errors="coerce").sum(skipna=True)
    cols[idx].metric(f"{t} Area % sum", f"{s:.2f}%")

avg_sum = pd.to_numeric(out_final["AVG AREA"], errors="coerce").sum(skipna=True)
cols[-1].metric("AVG AREA sum", f"{avg_sum:.2f}%")

# =========================
# Final Output Table
# =========================
st.subheader("Final Output Table (common compounds only, RT ascending)")
st.dataframe(out_final, use_container_width=True)

csv = out_final.to_csv(index=False).encode("utf-8")
st.download_button("Download CSV", csv, "gcms_multi_trial_final_output.csv", mime="text/csv")
