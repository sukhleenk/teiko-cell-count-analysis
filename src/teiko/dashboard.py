"""Interactive Streamlit dashboard over the cell-count database.

Run with `make dashboard` (or `streamlit run src/teiko/dashboard.py`).
"""
from __future__ import annotations

import sys
from pathlib import Path

if __package__ in (None, ""):  # `streamlit run src/teiko/dashboard.py`
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from teiko import analysis, config, db

# Validated categorical slots 1 & 2 (see plots.py).
BLUE, ORANGE = "#2a78d6", "#eb6834"
RESPONSE_COLORS = {"yes": BLUE, "no": ORANGE}
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e5e4e0"

st.set_page_config(page_title="Loblaw Bio — Cell Population Dashboard", layout="wide")


# --------------------------------------------------------------------------
# Data access (cached so filter changes don't re-query the whole database)
# --------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_frequencies() -> pd.DataFrame:
    return analysis.frequency_with_metadata()


@st.cache_data(show_spinner=False)
def load_cohort(condition: str, treatment: str, sample_type: str) -> pd.DataFrame:
    return analysis.responder_cohort(
        condition=condition, treatment=treatment, sample_type=sample_type
    )


def _base_layout(fig: go.Figure, height: int = 460) -> go.Figure:
    fig.update_layout(
        height=height,
        paper_bgcolor=SURFACE,
        plot_bgcolor=SURFACE,
        font=dict(color=INK, size=13),
        margin=dict(l=10, r=10, t=50, b=10),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0, title=None),
        hovermode="closest",
    )
    fig.update_xaxes(showgrid=False, linecolor=GRID, tickcolor=GRID, color=MUTED)
    fig.update_yaxes(gridcolor=GRID, zeroline=False, linecolor=GRID, color=MUTED)
    return fig


def _label(population: str) -> str:
    return config.POPULATION_LABELS.get(population, population)


# --------------------------------------------------------------------------
if not config.DB_PATH.exists():
    st.error(
        f"`{config.DB_PATH.name}` not found. Run `python load_data.py` "
        "(or `make pipeline`) first."
    )
    st.stop()

freq = load_frequencies()

st.title("Immune cell populations — miraclib trial")
st.caption(
    "Relative frequency of five immune cell populations across "
    f"{freq['sample'].nunique():,} samples from {freq['subject_id'].nunique():,} subjects "
    f"in {freq['project_id'].nunique()} projects."
)

tab_overview, tab_stats, tab_baseline, tab_schema = st.tabs(
    ["Overview (Part 2)", "Responders (Part 3)", "Baseline subset (Part 4)", "Database"]
)

# ------------------------------------------------------- Part 2: overview --
with tab_overview:
    st.subheader("Relative frequency of each population in each sample")

    c1, c2, c3, c4 = st.columns(4)
    projects = c1.multiselect("Project", sorted(freq["project_id"].unique()))
    conditions = c2.multiselect("Condition", sorted(freq["condition"].unique()))
    sample_types = c3.multiselect("Sample type", sorted(freq["sample_type"].unique()))
    search = c4.text_input("Sample or subject id contains", "")

    view = freq
    if projects:
        view = view[view["project_id"].isin(projects)]
    if conditions:
        view = view[view["condition"].isin(conditions)]
    if sample_types:
        view = view[view["sample_type"].isin(sample_types)]
    if search:
        needle = search.strip().lower()
        view = view[
            view["sample"].str.lower().str.contains(needle)
            | view["subject_id"].str.lower().str.contains(needle)
        ]

    m1, m2, m3 = st.columns(3)
    m1.metric("Samples", f"{view['sample'].nunique():,}")
    m2.metric("Rows in summary", f"{len(view):,}")
    m3.metric(
        "Median total count",
        f"{view.drop_duplicates('sample')['total_count'].median():,.0f}" if len(view) else "—",
    )

    summary = view[["sample", "total_count", "population", "count", "percentage"]]
    st.dataframe(summary, use_container_width=True, hide_index=True, height=380)
    st.download_button(
        "Download summary table (CSV)",
        summary.to_csv(index=False).encode(),
        file_name="cell_frequency_summary.csv",
        mime="text/csv",
    )

    if len(view):
        means = (
            view.groupby("population")["percentage"].mean().sort_values(ascending=True).reset_index()
        )
        fig = go.Figure(
            go.Bar(
                x=means["percentage"],
                y=[_label(p) for p in means["population"]],
                orientation="h",
                marker=dict(color=BLUE, line=dict(color=SURFACE, width=2)),
                text=[f"{v:.1f}%" for v in means["percentage"]],
                textposition="outside",
                hovertemplate="%{y}: %{x:.2f}%<extra></extra>",
            )
        )
        fig.update_layout(title="Mean relative frequency across the filtered samples")
        fig.update_xaxes(title="% of sample")
        st.plotly_chart(_base_layout(fig, 360), use_container_width=True)

# -------------------------------------------------------- Part 3: stats ----
with tab_stats:
    st.subheader("Responders vs non-responders")

    c1, c2, c3, c4 = st.columns(4)
    condition = c1.selectbox(
        "Condition", sorted(freq["condition"].unique()), index=sorted(freq["condition"].unique()).index(config.MELANOMA)
    )
    treatments = sorted(t for t in freq["treatment"].dropna().unique())
    treatment = c2.selectbox("Treatment", treatments, index=treatments.index(config.MIRACLIB))
    stypes = sorted(freq["sample_type"].unique())
    sample_type = c3.selectbox("Sample type", stypes, index=stypes.index(config.PBMC))
    unit = c4.selectbox(
        "Observation unit",
        ["sample", "subject"],
        help="'subject' averages a patient's repeat visits first, removing pseudo-replication.",
    )

    cohort = load_cohort(condition, treatment, sample_type)
    if cohort.empty:
        st.warning("No samples match that combination.")
    else:
        results = analysis.compare_responders(cohort, unit=unit)

        m1, m2, m3 = st.columns(3)
        m1.metric("Samples", f"{cohort['sample'].nunique():,}")
        m2.metric("Subjects", f"{cohort['subject_id'].nunique():,}")
        m3.metric("Significant populations (FDR 5%)", int(results["significant"].sum()))

        # ---- boxplot -----------------------------------------------------
        populations = sorted(cohort["population"].unique())
        plot_df = cohort
        if unit == "subject":
            plot_df = cohort.groupby(
                ["subject_id", "population", "response"], as_index=False
            )["percentage"].mean()

        fig = go.Figure()
        for response in ("yes", "no"):
            sub = plot_df[plot_df["response"] == response]
            fig.add_trace(
                go.Box(
                    x=[_label(p) for p in sub["population"]],
                    y=sub["percentage"],
                    name=RESPONSE_LABELS[response],
                    marker_color=RESPONSE_COLORS[response],
                    line=dict(width=1.5),
                    fillcolor=RESPONSE_COLORS[response],
                    boxpoints=False,
                    hovertemplate="%{x}<br>%{y:.2f}%<extra>"
                    + RESPONSE_LABELS[response]
                    + "</extra>",
                )
            )
        fig.update_layout(
            boxmode="group",
            boxgap=0.35,
            boxgroupgap=0.12,
            title=f"Relative frequency by population ({condition}, {treatment}, {sample_type})",
        )
        fig.update_yaxes(title="% of sample")
        st.plotly_chart(_base_layout(fig, 500), use_container_width=True)

        # ---- statistics --------------------------------------------------
        st.markdown("**Statistical test**")
        st.write(analysis.significance_narrative(results))

        display = results[
            [
                "population",
                "n_responder",
                "n_non_responder",
                "mean_responder",
                "mean_non_responder",
                "mean_difference",
                "rank_biserial",
                "p_value",
                "p_adjusted",
                "significant",
            ]
        ].copy()
        display["population"] = display["population"].map(_label)
        st.dataframe(
            display.style.format(
                {
                    "mean_responder": "{:.2f}",
                    "mean_non_responder": "{:.2f}",
                    "mean_difference": "{:+.2f}",
                    "rank_biserial": "{:+.3f}",
                    "p_value": "{:.2e}",
                    "p_adjusted": "{:.2e}",
                }
            ),
            use_container_width=True,
            hide_index=True,
        )
        st.caption(
            "Two-sided Mann-Whitney U per population (rank-based, no normality "
            "assumption), Benjamini-Hochberg FDR correction across the five "
            "populations. `rank_biserial` is the effect size: +1 means responders "
            "are always higher, -1 always lower. `mean_difference` is in "
            "percentage points."
        )
        st.download_button(
            "Download statistics (CSV)",
            results.to_csv(index=False).encode(),
            file_name=f"responder_statistics_{unit}.csv",
            mime="text/csv",
        )

# ------------------------------------------------------ Part 4: baseline ---
with tab_baseline:
    st.subheader("Baseline melanoma PBMC samples from miraclib-treated patients")

    c1, c2, c3 = st.columns(3)
    b_condition = c1.selectbox(
        "Condition ", sorted(freq["condition"].unique()),
        index=sorted(freq["condition"].unique()).index(config.MELANOMA), key="b_cond",
    )
    b_treatment = c2.selectbox(
        "Treatment ", treatments, index=treatments.index(config.MIRACLIB), key="b_trt"
    )
    times = sorted(int(t) for t in freq["time_from_treatment_start"].dropna().unique())
    b_time = c3.selectbox(
        "Time from treatment start (days)", times, index=times.index(config.BASELINE_DAY)
    )

    subset = analysis.baseline_subset(
        condition=b_condition, treatment=b_treatment, time_point=b_time
    )
    if subset.empty:
        st.warning("No samples match that combination.")
    else:
        breakdowns = analysis.baseline_breakdowns(subset)

        m1, m2 = st.columns(2)
        m1.metric("Samples", f"{len(subset):,}")
        m2.metric("Subjects", f"{subset['subject_id'].nunique():,}")

        cols = st.columns(3)
        panels = [
            ("samples_per_project", "project", "samples", "Samples per project"),
            ("subjects_by_response", "response", "subjects", "Subjects by response"),
            ("subjects_by_sex", "sex", "subjects", "Subjects by sex"),
        ]
        for col, (key, label_col, value_col, title) in zip(cols, panels):
            table = breakdowns[key]
            colors = (
                [RESPONSE_COLORS.get(v, BLUE) for v in table[label_col]]
                if key == "subjects_by_response"
                else [BLUE] * len(table)
            )
            fig = go.Figure(
                go.Bar(
                    x=table[label_col].astype(str),
                    y=table[value_col],
                    marker=dict(color=colors, line=dict(color=SURFACE, width=2)),
                    text=[f"{v:,}" for v in table[value_col]],
                    textposition="outside",
                    hovertemplate="%{x}: %{y:,}<extra></extra>",
                )
            )
            fig.update_layout(title=title, showlegend=False)
            col.plotly_chart(_base_layout(fig, 320), use_container_width=True)
            col.dataframe(table, use_container_width=True, hide_index=True)

        st.markdown("**Matching samples**")
        st.dataframe(subset, use_container_width=True, hide_index=True, height=320)
        st.download_button(
            "Download subset (CSV)",
            subset.to_csv(index=False).encode(),
            file_name="baseline_subset.csv",
            mime="text/csv",
        )

# ----------------------------------------------------------- schema tab ----
with tab_schema:
    st.subheader("Database")
    tables = db.query(
        "SELECT type, name FROM sqlite_master WHERE name NOT LIKE 'sqlite_%' "
        "AND type IN ('table','view') ORDER BY type, name"
    )
    rows = []
    for r in tables.itertuples():
        n = db.query(f"SELECT COUNT(*) AS n FROM {r.name}").iloc[0]["n"]
        rows.append({"object": r.name, "type": r.type, "rows": int(n)})
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    st.code(config.SCHEMA_PATH.read_text(), language="sql")
