import streamlit as st
import pandas as pd
import numpy as np
import io
import re

# ------------------ Page ------------------
st.set_page_config(layout="wide", page_title="GC-MS Multi-Trial Analysis")
st.title("GC-MS Multi-Trial Analysis 🧪")
st.caption(
    "Paste trial tables (TSV copied from Excel). The app matches compounds across trials using "
    "alias-aware Name/Species logic, provides RT overview, overlap summaries (in 2/3/... trials), "
    "and a final common-compounds table with recomputed Area % and AVG AREA."
)

# ------------------ Matching helpers ------------------
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
        x = re.sub(r"\s+or\s+", "|", x)
        x = re.sub(r"^\s*or\s+", "", x)
        aliases.extend([p.strip() for p in x.split("|") if p.strip()])
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

# ------------------ Parsing helpers ------------------
def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(
        series.astype(str).str.replace(",", "", regex=False).str.strip(),
        errors="coerce",
    )

def parse_trial(text: str) -> pd.DataFrame:
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

    df = df[df["key"].astype(str).str.strip() != ""].copy()
    df = df[pd.notna(df["Area"]) & (df["Area"] > 0)].copy()
    return df

def first_row_for_key(df: pd.DataFrame, key: str) -> pd.Series:
    return df.loc[df["key"] == key].iloc[0]

# ------------------ UI: Trials ------------------
n_trials = st.selectbox("Number of trials", [2, 3, 4, 5, 6], index=2)

trial_dfs = {}
dup_warnings = {}

for i in range(n_trials):
    t = f"T{i+1}"
    st.subheader(t)
    txt = st.text_area(f"Paste {t} table (TSV)", height=170, key=f"trial_{t}")
    if txt.strip():
        try:
            df = parse_trial(txt)
            trial_dfs[t] = df

            dups = df["key"].value_counts()
            dups = dups[dups > 1]
            if len(dups) > 0:
                dup_warnings[t] = dups

            st.success(f"{t} parsed ✅  ({len(df)} rows)")
        except Exception as e:
            st.error(f"{t} error: {e}")

if len(trial_dfs) != n_trials:
    st.info("Please paste valid data for ALL trials.")
    st.stop()

if dup_warnings:
    with st.expander("⚠️ Duplicate key warnings (no collapsing; first match is used)"):
        for t, s in dup_warnings.items():
            st.write(f"**{t}** duplicates (key → count)")
            st.dataframe(s.rename("count").reset_index().rename(columns={"index": "key"}), hide_index=True)

st.divider()

# ------------------ RT Overview (All trials combined) ------------------
st.subheader("RT Overview (All trials combined, RT ascending)")

all_rows = []
for t, df in trial_dfs.items():
    tmp = df.copy()
    tmp["Trial"] = t
    all_rows.append(tmp)

rt_all = pd.concat(all_rows, ignore_index=True)
rt_all = rt_all.sort_values("RT", ascending=True, na_position="last").reset_index(drop=True)

overview_cols = ["Trial", "RT", "Name", "Formula", "Species", "Area", "Score", "key"]
overview_cols = [c for c in overview_cols if c in rt_all.columns]

st.dataframe(rt_all[overview_cols], use_container_width=True, hide_index=True)

# ------------------ Repeated Compounds ------------------
st.subheader("Repeated Compounds")
st.caption("Compounds that appear in two or more trials.")

# Build key -> set(trials) locally (so no NameError)
key_trials = {}
for t, df in trial_dfs.items():
    for k in set(df["key"].tolist()):
        key_trials.setdefault(k, set()).add(t)

rows = []
for key, ts in key_trials.items():
    if len(ts) < 2:
        continue

    # Representative info: use first trial where it appears (alphabetical T1, T2, ...)
    first_trial = sorted(ts)[0]
    rep = trial_dfs[first_trial].loc[trial_dfs[first_trial]["key"] == key].iloc[0]

    rows.append({
        "Name": rep.get("Name", ""),
        "Formula": rep.get("Formula", ""),
        "Trials": ", ".join(sorted(ts)),
    })

if not rows:
    st.info("No repeated compounds found.")
else:
    repeated_df = pd.DataFrame(rows)

    # Optional: sort by how many trials (desc), then name (asc)
    repeated_df["#Trials"] = repeated_df["Trials"].apply(lambda x: len([p for p in x.split(",") if p.strip()]))
    repeated_df = repeated_df.sort_values(["#Trials", "Name"], ascending=[False, True]).drop(columns=["#Trials"])

    st.dataframe(repeated_df, use_container_width=True, hide_index=True)

# ------------------ Final Output: common compounds only (intersection of all trials) ------------------
common_keys = set.intersection(*[set(df["key"]) for df in trial_dfs.values()])
if not common_keys:
    st.warning("No common compounds found across all trials.")
    st.stop()

rows = []
for key in sorted(common_keys):
    per_trial = {t: first_row_for_key(df, key) for t, df in trial_dfs.items()}

    # Use mean RT across trials for stable RT sorting
    rt_vals = pd.to_numeric(pd.Series([per_trial[t]["RT"] for t in trial_dfs.keys()]), errors="coerce")
    rt_display = rt_vals.mean(skipna=True)

    base = per_trial["T1"]
    row = {
        "RT": rt_display,
        "Name": base.get("Name", ""),
        "Formula": base.get("Formula", ""),
        "Species": base.get("Species", ""),
    }

    for t in trial_dfs.keys():
        row[f"{t} RT"] = per_trial[t]["RT"]
        row[f"{t} Area"] = per_trial[t]["Area"]
        row[f"{t} Score"] = per_trial[t].get("Score", np.nan)

    rows.append(row)

out = pd.DataFrame(rows)

# recompute Area % within common set so each trial sums to 100
for t in trial_dfs.keys():
    areas = pd.to_numeric(out[f"{t} Area"], errors="coerce")
    total = areas.sum(skipna=True)
    out[f"{t} Area %"] = (areas / total) * 100 if pd.notna(total) and total > 0 else np.nan

area_pct_cols = [f"{t} Area %" for t in trial_dfs.keys()]
out["AVG AREA"] = out[area_pct_cols].mean(axis=1, skipna=True)

# sort RT ascending + cumulative
out = out.sort_values("RT", ascending=True, na_position="last").reset_index(drop=True)
out["CUMULATIVE %"] = out["AVG AREA"].cumsum()

# ------------------ Top 10 by AVG AREA ------------------
st.subheader("Top 10 Compounds (by AVG AREA)")
top10 = out.sort_values("AVG AREA", ascending=False, na_position="last").head(10).copy()

top_cols = ["RT", "Name", "Formula", "Species", "AVG AREA"]
for t in trial_dfs.keys():
    top_cols += [f"{t} Area %", f"{t} Score"]
st.dataframe(top10[top_cols], use_container_width=True, hide_index=True)

# ------------------ Validation ------------------
st.subheader("Validation (each trial Area % sum must be 100; AVG AREA sum must be 100)")
cols = st.columns(n_trials + 1)
for idx, t in enumerate(trial_dfs.keys()):
    s = pd.to_numeric(out[f"{t} Area %"], errors="coerce").sum(skipna=True)
    cols[idx].metric(f"{t} Area % sum", f"{s:.2f}%")
cols[-1].metric("AVG AREA sum", f"{pd.to_numeric(out['AVG AREA'], errors='coerce').sum(skipna=True):.2f}%")

# ------------------ Final Output Table (exact layout) ------------------
st.subheader("Final Output Table (common compounds only, RT ascending)")

final_cols = ["RT", "Name", "Formula", "Species", "AVG AREA", "CUMULATIVE %"]
for t in trial_dfs.keys():
    final_cols += [f"{t} RT", f"{t} Area", f"{t} Area %", f"{t} Score"]

st.dataframe(out[final_cols], use_container_width=True, hide_index=True)

csv = out[final_cols].to_csv(index=False).encode("utf-8")
st.download_button("Download CSV", csv, "gcms_multi_trial_final_output.csv", mime="text/csv")
