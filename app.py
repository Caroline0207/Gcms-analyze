import streamlit as st
import pandas as pd
import numpy as np
import io
import re
st.markdown("""
<style>
/* ---------- REMOVE STREAMLIT DIVIDER / EMPTY BARS ---------- */
hr { display: none !important; }

div[data-testid="stVerticalBlock"] > div:has(> hr),
div[data-testid="stVerticalBlock"] > div:empty {
  display: none !important;
}

div[data-testid="stSpacer"] {
  display: none !important;
}

/* ---------- Layout ---------- */
.block-container {
  max-width: 1180px;
  padding-top: 1.2rem;
  padding-bottom: 2.2rem;
  background: #fafbfc;
}

/* ---------- Typography ---------- */
h1, h2, h3 { letter-spacing: -0.3px; }
h1 { font-weight: 750; margin-bottom: 0.2rem; }
h2 { margin-top: 1.2rem; }

/* ---------- Simple section divider ---------- */
.section-divider {
  border-top: 1px solid rgba(0,0,0,0.08);
  margin: 1.4rem 0 1.2rem 0;
}

/* ---------- Soft badge ---------- */
.badge {
  : inline-block;
  padding: 0.22rem 0.6rem;
  border-radius: 999px;
  font-size: 0.78rem;
  font-weight: 600;
  background: linear-gradient(135deg, #e6f3ff, #f2ecff);
  color: #4a4a68;
  border: 1px solid rgba(120,120,180,0.25);
  margin-left: 0.45rem;
}

/* ---------- Section emoji header ---------- */
.section-title {
  : flex;
  align-items: center;
  gap: 0.45rem;
}

/* ---------- Subtle caption ---------- */
.smallcap {
  font-size: 0.9rem;
  opacity: 0.7;
  margin-top: -0.25rem;
}

/* ---------- Dataframe polish ---------- */
div[data-testid="stDataFrame"] {
  border-radius: 12px;
  overflow: hidden;
}

/* ---------- Expander ---------- */
details > summary {
  font-weight: 600;
}
</style>
""", unsafe_allow_html=True)

st.markdown(
    """
    <div class="section-title">
      <h1>GC-MS Multi-Trial Analysis</h1>
      <span class="badge">🧪 lab tool</span>
      <span class="badge">🫧 clean data</span>
    </div>
    <div class="smallcap">
      Paste cleaned GC-MS tables, inspect overlaps, and export a publication-ready summary.
    </div>
    """,
    unsafe_allow_html=True
)


# =====================================================
# Helper functions
# =====================================================
def normalize_text(s):
    if s is None or pd.isna(s):
        return ""
    s = str(s).lower()
    s = re.sub(r"[(),;]", " ", s)
    s = s.replace("-", " ")
    s = re.sub(r"\s+", " ", s).strip()
    return s

def canonical_key(name, species):
    parts = []
    for x in [name, species]:
        if x and str(x).strip():
            x = str(x).lower()
            x = re.sub(r"\s+or\s+", "|", x)
            parts += [p.strip() for p in x.split("|") if p.strip()]
    parts = [normalize_text(p) for p in parts if p]
    return min(parts, key=len) if parts else ""

def to_num(s):
    return pd.to_numeric(
        s.astype(str).str.replace(",", "", regex=False),
        errors="coerce"
    )

def parse_trial(text):
    df = pd.read_csv(io.StringIO(text.strip()), sep="\t", dtype=str, engine="python")
    df.columns = [c.strip() for c in df.columns]

    for c in ["RT", "Area", "Score"]:
        if c not in df.columns:
            df[c] = np.nan
    for c in ["Name", "Formula", "Species"]:
        if c not in df.columns:
            df[c] = ""

    df["RT"] = to_num(df["RT"])
    df["Area"] = to_num(df["Area"])
    df["Score"] = to_num(df["Score"])

    df["Name"] = df["Name"].fillna("").str.strip()
    df["Species"] = df["Species"].fillna("").str.strip()
    df["Formula"] = df["Formula"].fillna("").str.strip()

    df["key"] = df.apply(lambda r: canonical_key(r["Name"], r["Species"]), axis=1)

    df = df[df["key"] != ""].copy()
    df = df[df["Area"] > 0].copy()
    return df
ATOMIC_WEIGHTS = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.06,
    "P": 30.974,
    "F": 18.998,
    "Cl": 35.45,
    "Br": 79.904,
    "I": 126.904,
    "Si": 28.085
}

def parse_formula_counts(formula):
    """
    Parse simple molecular formulas like C10H18O, C12H22O11, C7H8O2.
    Returns dict of element counts.
    """
    if formula is None or pd.isna(formula):
        return {}

    formula = str(formula).strip()
    if formula == "":
        return {}

    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    if not tokens:
        return {}

    counts = {}
    for elem, num in tokens:
        count = int(num) if num else 1
        counts[elem] = counts.get(elem, 0) + count
    return counts

def molecular_weight(formula):
    counts = parse_formula_counts(formula)
    if not counts:
        return np.nan

    mw = 0.0
    for elem, count in counts.items():
        if elem not in ATOMIC_WEIGHTS:
            return np.nan
        mw += ATOMIC_WEIGHTS[elem] * count
    return mw

def carbon_fraction(formula):
    """
    Returns carbon mass fraction in the compound:
    (mass of carbon in molecule) / (molecular weight)
    """
    counts = parse_formula_counts(formula)
    if not counts:
        return np.nan

    if "C" not in counts or counts["C"] == 0:
        return 0.0

    mw = molecular_weight(formula)
    if pd.isna(mw) or mw == 0:
        return np.nan

    carbon_mass = counts["C"] * ATOMIC_WEIGHTS["C"]
    return carbon_mass / mw
# =====================================================
# Trial input
# =====================================================
n_trials = st.selectbox("Number of trials", [2, 3, 4, 5], index=2)

trial_dfs = {}
st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)

for i in range(n_trials):
    t = f"T{i+1}"
    txt = st.text_area(f"{t} (paste TSV from Cleaning web)", height=160, key=t)
    if txt.strip():
        try:
            trial_dfs[t] = parse_trial(txt)
            st.caption(f"{t}: {len(trial_dfs[t])} rows parsed")
        except Exception as e:
            st.error(f"{t}: {e}")

st.markdown("</div>", unsafe_allow_html=True)

if len(trial_dfs) != n_trials:
    st.stop()

# =====================================================
# RT Overview (all trials combined)
# =====================================================
st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("🧭 RT Overview")
st.caption("All trials combined, sorted by retention time.")

rt_all = []
for t, df in trial_dfs.items():
    tmp = df.copy()
    tmp["Trial"] = t
    rt_all.append(tmp)

rt_all = pd.concat(rt_all).sort_values("RT").reset_index(drop=True)

st.dataframe(
    rt_all[["Trial", "RT", "Name", "Formula", "Species", "Area", "Score"]],
    use_container_width=True,
    height=320
)
st.markdown("</div>", unsafe_allow_html=True)

# =====================================================
# Repeated Compounds (≥2 trials)
# =====================================================
st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("🔁 Repeated Compounds")
st.caption("Compounds appearing in two or more trials.")

key_trials = {}
for t, df in trial_dfs.items():
    for k in df["key"].unique():
        key_trials.setdefault(k, set()).add(t)

rows = []
for k, ts in key_trials.items():
    if len(ts) >= 2:
        t0 = sorted(ts)[0]
        rep = trial_dfs[t0].loc[trial_dfs[t0]["key"] == k].iloc[0]
        rows.append({
            "Name": rep["Name"],
            "Formula": rep["Formula"],
            "Trials": ", ".join(sorted(ts))
        })

if rows:
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
else:
    st.info("No repeated compounds.")

st.markdown("</div>", unsafe_allow_html=True)

# =====================================================
# Final Output Table (common to ALL trials)
# =====================================================
common_keys = set.intersection(*[set(df["key"]) for df in trial_dfs.values()])


# ---------- Non-common compound summary ----------
summary_rows = []
for t, df in trial_dfs.items():
    total_area = df["Area"].sum()
    common_area = df.loc[df["key"].isin(common_keys), "Area"].sum()
    non_common_area = total_area - common_area

    common_pct_total = (common_area / total_area * 100) if total_area > 0 else np.nan
    non_common_pct_total = (non_common_area / total_area * 100) if total_area > 0 else np.nan

    summary_rows.append({
        "Trial": t,
        "Total Area": total_area,
        "Common Area": common_area,
        "Common % of Total": common_pct_total,
        "Non-common Area": non_common_area,
        "Non-common % of Total": non_common_pct_total
    })

non_common_summary = pd.DataFrame(summary_rows)

# ---------- Build final common-compound table ----------
rows = []
for k in sorted(common_keys):
    per = {}
    for t in trial_dfs:
        match = trial_dfs[t].loc[trial_dfs[t]["key"] == k]
        if not match.empty:
            per[t] = match.iloc[0]

    # safety check: must exist in all trials
    if len(per) != len(trial_dfs):
        continue

    rt_mean = np.mean([per[t]["RT"] for t in per])

    base = per[list(per.keys())[0]]
    row = {
        "RT": rt_mean,
        "Name": base["Name"],
        "Formula": base["Formula"],
        "Species": base["Species"],
    }

    for t in per:
        row[f"{t} RT"] = per[t]["RT"]
        row[f"{t} Area"] = per[t]["Area"]
        row[f"{t} Score"] = per[t]["Score"]

    rows.append(row)

out = pd.DataFrame(rows)
final_cols = []

if out.empty:
    st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
    st.subheader("📊 Final Output Table")
    st.warning("No compounds were common to all trials.")
else:
    out = out.sort_values("RT").reset_index(drop=True)

    for t in trial_dfs:
        trial_total_area = out[f"{t} Area"].sum()
        out[f"{t} Area %"] = (
            out[f"{t} Area"] / trial_total_area * 100
            if trial_total_area > 0 else np.nan
        )

    area_pct_cols = [f"{t} Area %" for t in trial_dfs]
    out["AVG AREA"] = out[area_pct_cols].mean(axis=1)
    out["CUMULATIVE %"] = out["AVG AREA"].cumsum()
    out["C Fraction"] = out["Formula"].apply(carbon_fraction)
    out["Organic Carbon %"] = out["C Fraction"] * 100

    final_cols = ["RT", "Name", "Formula", "Species", "Organic Carbon %", "AVG AREA", "CUMULATIVE %"]
    for t in trial_dfs:
        final_cols += [f"{t} RT", f"{t} Area", f"{t} Area %", f"{t} Score"]

    st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
    st.subheader("📊 Final Output Table")
    st.caption("Compounds common to all trials.")

    st.dataframe(
        out[final_cols],
        use_container_width=True,
    )

st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("🧩 Non-common Compound Summary")
st.caption("Compounds not shared by all trials, shown as percentage of each trial's total area.")

st.dataframe(
    non_common_summary,
    use_container_width=True,
    hide_index=True
)
# ---------- Copy Final Output Table (toggle) ----------
if "show_copy" not in st.session_state:
    st.session_state.show_copy = False

col1, col2 = st.columns([1, 6])

with col1:
    if st.button("📋 Copy"):
        st.session_state.show_copy = not st.session_state.show_copy

if st.session_state.show_copy:
    if out.empty or not final_cols:
        st.info("No final output table available to copy.")
    else:
        tsv_text = out[final_cols].to_csv(sep="\t", index=False)
        st.text_area(
            "Copy & paste (TSV, header included):",
            tsv_text,
            height=180
        )
# =====================================================
# TOC-based Compound Concentration Estimation
# =====================================================
st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("🧪 TOC-based Compound Concentration")
st.caption("Estimate each compound concentration from TOC (ppm = mg C/L), AVG AREA, and formula-based carbon fraction.")

if out.empty or "AVG AREA" not in out.columns:
    st.info("TOC calculation is unavailable because there are no compounds common to all trials.")
else:
    toc_value = st.number_input(
        "Enter TOC value (ppm = mg C/L)",
        min_value=0.0,
        value=12.4,
        step=0.1
    )

    toc_df = out.copy()
    toc_df["C Fraction"] = toc_df["Formula"].apply(carbon_fraction)
    toc_df["Organic Carbon %"] = toc_df["C Fraction"] * 100
    toc_df["Molecular Weight"] = toc_df["Formula"].apply(molecular_weight)
    toc_df["Carbon Count"] = toc_df["Formula"].apply(
        lambda f: parse_formula_counts(f).get("C", np.nan) if pd.notna(f) and str(f).strip() else np.nan
    )
    toc_df["Area Fraction"] = toc_df["AVG AREA"] / 100

    # remove rows with invalid formula interpretation
    toc_df = toc_df.dropna(subset=["C Fraction"]).copy()

    if toc_df.empty:
        st.warning("Could not calculate concentrations because no valid molecular formulas were available.")
    else:
        toc_df["Weighted C Term"] = toc_df["Area Fraction"] * toc_df["C Fraction"]
        denom = toc_df["Weighted C Term"].sum()

        if denom > 0:
            toc_df["Estimated Compound Conc. (mg/L)"] = (
                toc_value * toc_df["Area Fraction"] / denom
            )

            toc_df["Estimated Carbon Conc. (mg C/L)"] = (
                toc_df["Estimated Compound Conc. (mg/L)"] * toc_df["C Fraction"]
            )

            toc_df["% of Total TOC"] = (
                toc_df["Estimated Carbon Conc. (mg C/L)"] / toc_value * 100
            )

            toc_df["Rank"] = (
                toc_df["Estimated Compound Conc. (mg/L)"]
                .rank(method="dense", ascending=False)
                .astype(int)
            )

            toc_df = toc_df.sort_values(
                "Estimated Compound Conc. (mg/L)", ascending=False
            ).reset_index(drop=True)

            _cols = [
                "Rank",
                "RT",
                "Name",
                "Species",
                "Formula",
                "Carbon Count",
                "Molecular Weight",
                "AVG AREA",
                "Area Fraction",
                "C Fraction",
                "Organic Carbon %",
                "Weighted C Term",
                "Estimated Compound Conc. (mg/L)",
                "Estimated Carbon Conc. (mg C/L)",
                "% of Total TOC"
            ]
          toc_df["Molecular Weight"] = toc_df["Molecular Weight"].round(2)
      toc_df["AVG AREA"] = toc_df["AVG AREA"].round(2)
toc_df["Area Fraction"] = toc_df["Area Fraction"].round(4)
toc_df["C Fraction"] = toc_df["C Fraction"].round(4)
toc_df["Organic Carbon %"] = toc_df["Organic Carbon %"].round(2)
toc_df["Weighted C Term"] = toc_df["Weighted C Term"].round(4)
toc_df["Estimated Compound Conc. (mg/L)"] = toc_df["Estimated Compound Conc. (mg/L)"].round(3)
toc_df["Estimated Carbon Conc. (mg C/L)"] = toc_df["Estimated Carbon Conc. (mg C/L)"].round(3)
toc_df["% of Total TOC"] = toc_df["% of Total TOC"].round(2)
toc_df["RT"] = toc_df["RT"].round(3)

            st.dataframe(
                toc_df[display_cols],
                use_container_width=True,
                hide_index=True,
                height=420
            )

            st.caption(
                f"Check: Estimated carbon sum = "
                f"{toc_df['Estimated Carbon Conc. (mg C/L)'].sum():.3f} mg C/L "
                f"(target TOC = {toc_value:.3f} mg C/L)"
            )

            if "show_toc_copy" not in st.session_state:
                st.session_state.show_toc_copy = False

            c1, c2 = st.columns([1, 6])
            with c1:
                if st.button("📋 Copy TOC Table"):
                    st.session_state.show_toc_copy = not st.session_state.show_toc_copy

            if st.session_state.show_toc_copy:
                toc_tsv = toc_df[display_cols].to_csv(sep="\t", index=False)
                st.text_area(
                    "Copy & paste TOC result table (TSV, header included):",
                    toc_tsv,
                    height=180
                )
        else:
            st.warning("Could not calculate concentrations. Check AVG AREA values and molecular formulas.")
# =====================================================
# Top 10 (after Final Output)
# =====================================================
st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("✨ Top 10 Compounds")
st.caption("Top contributors based on average area percentage.")

if out.empty or "AVG AREA" not in out.columns:
    st.info("No Top 10 table available.")
else:
    top10 = out.sort_values("AVG AREA", ascending=False).head(10)[["Name", "Formula", "AVG AREA"]]
    st.dataframe(top10, use_container_width=True, hide_index=True, height=260)

st.markdown("</div>", unsafe_allow_html=True)


# =====================================================
# Validation (quiet)
# =====================================================
with st.expander("Validation (Area % sums)", expanded=False):
    if out.empty or "AVG AREA" not in out.columns:
        st.info("Validation unavailable because no common compounds were found.")
    else:
        cols = st.columns(n_trials + 1)
        for i, t in enumerate(trial_dfs):
            cols[i].caption(f"{t}: {out[f'{t} Area %'].sum():.2f}%")
        cols[-1].caption(f"AVG AREA: {out['AVG AREA'].sum():.2f}%")
