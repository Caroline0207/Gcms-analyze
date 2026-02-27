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
  display: inline-block;
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
  display: flex;
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

# common_keys 계산 (NaN 방어)
key_sets = [set(df["key"].dropna()) for df in trial_dfs.values()]
common_keys = set.intersection(*key_sets) if key_sets else set()

rows = []
for k in sorted(common_keys):
    per = {}
    ok = True
    for t, df in trial_dfs.items():
        m = df.loc[df["key"] == k]
        if m.empty:
            ok = False
            break
        per[t] = m.iloc[0]
    if not ok:
        continue

    rt_vals = [per[t]["RT"] for t in per if pd.notna(per[t]["RT"])]
    rt_mean = float(np.mean(rt_vals)) if rt_vals else np.nan

    base = per[next(iter(per))]
    row = {
        "RT": rt_mean,
        "Name": base.get("Name", ""),
        "Formula": base.get("Formula", ""),
        "Species": base.get("Species", ""),
    }

    for t in per:
        row[f"{t} RT"] = per[t].get("RT", np.nan)
        row[f"{t} Area"] = per[t].get("Area", np.nan)
        row[f"{t} Score"] = per[t].get("Score", np.nan)

    rows.append(row)

# out 만들기 (빈 경우에도 스키마 유지)
base_cols = ["RT", "Name", "Formula", "Species"]
trial_cols = []
for t in trial_dfs:
    trial_cols += [f"{t} RT", f"{t} Area", f"{t} Score"]
out = pd.DataFrame(rows)
if out.empty:
    out = pd.DataFrame(columns=base_cols + trial_cols)
else:
    out = out.sort_values("RT", kind="mergesort").reset_index(drop=True)

# Area% + AVG + CUMULATIVE 계산 (sum=0 / 빈 out 방어)
for t in trial_dfs:
    area_col = f"{t} Area"
    pct_col = f"{t} Area %"

    if area_col in out.columns and not out.empty:
        denom = pd.to_numeric(out[area_col], errors="coerce").sum()
        if denom and denom != 0:
            out[pct_col] = pd.to_numeric(out[area_col], errors="coerce") / denom * 100
        else:
            out[pct_col] = np.nan
    else:
        out[pct_col] = np.nan

area_pct_cols = [f"{t} Area %" for t in trial_dfs]
if not out.empty:
    out["AVG AREA"] = pd.to_numeric(out[area_pct_cols], errors="coerce").mean(axis=1)
    out["CUMULATIVE %"] = out["AVG AREA"].cumsum()
else:
    out["AVG AREA"] = pd.Series(dtype=float)
    out["CUMULATIVE %"] = pd.Series(dtype=float)

final_cols = ["RT", "Name", "Formula", "Species", "AVG AREA", "CUMULATIVE %"]
for t in trial_dfs:
    final_cols += [f"{t} RT", f"{t} Area", f"{t} Area %", f"{t} Score"]

st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("📊 Final Output Table")
st.caption("Compounds common to all trials.")

if out.empty:
    st.info("No common compounds across all trials.")
else:
    st.dataframe(out[final_cols], use_container_width=True)

# ---------- Copy Final Output Table (toggle) ----------
if "show_copy" not in st.session_state:
    st.session_state.show_copy = False

col1, col2 = st.columns([1, 6])
with col1:
    if st.button("📋 Copy"):
        st.session_state.show_copy = not st.session_state.show_copy

if st.session_state.show_copy:
    if out.empty:
        st.info("Nothing to copy (final table is empty).")
    else:
        tsv_text = out[final_cols].to_csv(sep="\t", index=False)
        st.text_area("Copy & paste (TSV, header included):", tsv_text, height=180)

# =====================================================
# Top 10 (after Final Output)
# =====================================================
st.markdown("<div class='section-divider'></div>", unsafe_allow_html=True)
st.subheader("✨ Top 10 Compounds")
st.caption("Top contributors based on average area percentage.")

if out.empty or "AVG AREA" not in out.columns:
    st.info("Top 10 not available (final table is empty).")
else:
    top10 = out.sort_values("AVG AREA", ascending=False).head(10)[["Name", "Formula", "AVG AREA"]]
    st.dataframe(top10, use_container_width=True, hide_index=True, height=260)

st.markdown("</div>", unsafe_allow_html=True)

# =====================================================
# Validation (quiet)
# =====================================================
with st.expander("Validation (Area % sums)", expanded=False):
    if out.empty:
        st.info("No data to validate.")
    else:
        cols = st.columns(n_trials + 1)
        for i, t in enumerate(trial_dfs):
            cols[i].caption(f"{t}: {out[f'{t} Area %'].sum(skipna=True):.2f}%")
        cols[-1].caption(f"AVG AREA: {out['AVG AREA'].sum(skipna=True):.2f}%")
