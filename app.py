import streamlit as st
import pandas as pd
import numpy as np
import io
import re

st.set_page_config(layout="wide")
st.title("GC-MS Multi-Trial Analysis 🧪")

# ------------------ helpers ------------------
def normalize_text(s):
    if not s or pd.isna(s):
        return ""
    s = s.lower()
    s = re.sub(r"[\(\),\-]", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s

def extract_aliases(name, species):
    aliases = []
    for x in [name, species]:
        if not x:
            continue
        x = x.lower().replace("or", "|")
        aliases += [i.strip() for i in x.split("|")]
    return aliases

def canonical_key(name, species):
    aliases = extract_aliases(name, species)
    aliases = [normalize_text(a) for a in aliases if a]

    trivial = [
        a for a in aliases
        if len(a) <= 20 and "benzene" not in a and "methyl" not in a
    ]
    return trivial[0] if trivial else min(aliases, key=len)

def parse_trial(text):
    df = pd.read_csv(io.StringIO(text), sep="\t")
    df["Area %"] = pd.to_numeric(df["Area %"], errors="coerce")
    df["Area"] = pd.to_numeric(df["Area"], errors="coerce")
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    df["RT"] = pd.to_numeric(df["RT"], errors="coerce")

    df["key"] = df.apply(
        lambda r: canonical_key(r["Name"], r.get("Species", "")),
        axis=1,
    )
    return df

# ------------------ UI ------------------
n_trials = st.selectbox("Number of trials", [2, 3, 4, 5, 6], index=2)

trial_dfs = {}
for i in range(n_trials):
    st.subheader(f"T{i+1}")
    txt = st.text_area(f"Paste T{i+1} table", height=180)
    if txt.strip():
        trial_dfs[f"T{i+1}"] = parse_trial(txt)

if len(trial_dfs) != n_trials:
    st.info("Please paste data for all trials.")
    st.stop()

# ------------------ analysis ------------------
common_keys = set.intersection(*[set(df["key"]) for df in trial_dfs.values()])

rows = []
for key in sorted(common_keys):
    base = trial_dfs["T1"].loc[trial_dfs["T1"]["key"] == key].iloc[0]

    row = {
        "RT": base["RT"],
        "Name": base["Name"],
        "Formula": base["Formula"],
        "Species": base.get("Species", ""),
    }

    avg = 0
    for t, df in trial_dfs.items():
        r = df[df["key"] == key].iloc[0]
        row[f"{t} RT"] = r["RT"]
        row[f"{t} Area"] = r["Area"]
        row[f"{t} Area %"] = r["Area %"]
        row[f"{t} Score"] = r["Score"]
        avg += r["Area %"]

    row["AVG AREA"] = avg / n_trials
    rows.append(row)

out = pd.DataFrame(rows)
out = out.sort_values("AVG AREA", ascending=True)
out["CUMULATIVE %"] = out["AVG AREA"].cumsum()

# ------------------ validation ------------------
st.subheader("Validation")
for t, df in trial_dfs.items():
    st.write(f"{t} Area % sum:", round(df["Area %"].sum(), 2))

st.write("AVG AREA sum:", round(out["AVG AREA"].sum(), 2))

# ------------------ output ------------------
st.subheader("Final Output Table")
st.dataframe(out, use_container_width=True)

csv = out.to_csv(index=False).encode()
st.download_button("Download CSV", csv, "gcms_multi_trial_output.csv")
