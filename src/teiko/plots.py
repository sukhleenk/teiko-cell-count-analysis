"""Static figures for the pipeline (matplotlib). The dashboard renders its own
interactive versions of the same charts with Plotly.

Colors come from the validated two-slot categorical palette: slot 1 blue for
responders, slot 2 orange for non-responders. Identity is never color-alone --
every figure carries a legend and axis labels.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from . import config  # noqa: E402

RESPONSE_COLORS = {"yes": "#2a78d6", "no": "#eb6834"}
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
GRID = "#e5e4e0"


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, length=0)


def responder_boxplot(
    cohort: pd.DataFrame,
    results: pd.DataFrame | None = None,
    path: Path | str | None = None,
) -> Path:
    """Grouped boxplot of relative frequency by population and response."""
    populations = sorted(cohort["population"].unique())
    fig, ax = plt.subplots(figsize=(11, 6), facecolor=SURFACE)

    width, offset = 0.34, 0.19
    for i, response in enumerate(("yes", "no")):
        positions = [x - offset + i * 2 * offset for x in range(len(populations))]
        data = [
            cohort.loc[
                (cohort["population"] == p) & (cohort["response"] == response), "percentage"
            ].to_numpy()
            for p in populations
        ]
        color = RESPONSE_COLORS[response]
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=width,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color=SURFACE, linewidth=2),
            whiskerprops=dict(color=color, linewidth=1.5),
            capprops=dict(color=color, linewidth=1.5),
            boxprops=dict(facecolor=color, edgecolor=SURFACE, linewidth=2),
        )
        bp["boxes"][0].set_label(RESPONSE_LABELS[response])

    ax.set_xticks(range(len(populations)))
    ax.set_xticklabels([config.POPULATION_LABELS.get(p, p) for p in populations], color=INK)
    ax.set_ylabel("Relative frequency (% of sample)", color=INK_MUTED)
    ax.set_title(
        "Immune cell population frequencies, responders vs non-responders\n"
        "Melanoma patients on miraclib, PBMC samples",
        color=INK,
        loc="left",
        fontsize=13,
        pad=14,
    )
    _style(ax)

    # Significance markers: identity is carried by the annotation, not color.
    if results is not None and not results.empty:
        top = cohort["percentage"].max()
        marks = results.set_index("population")
        for i, population in enumerate(populations):
            if population in marks.index and bool(marks.loc[population, "significant"]):
                q = marks.loc[population, "p_adjusted"]
                ax.text(
                    i,
                    top * 1.02,
                    f"q={q:.1e} *",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color=INK_MUTED,
                )
        ax.set_ylim(top=top * 1.12)

    ax.legend(frameon=False, loc="upper right", labelcolor=INK)
    fig.tight_layout()

    path = Path(path or config.FIGURE_DIR / "responder_vs_nonresponder_boxplot.png")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return path


def baseline_breakdown_figure(
    breakdowns: dict[str, pd.DataFrame], path: Path | str | None = None
) -> Path:
    """Three small bar charts summarising the Part 4 baseline subset."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), facecolor=SURFACE)

    panels = [
        ("samples_per_project", "project", "samples", "Samples per project"),
        ("subjects_by_response", "response", "subjects", "Subjects by response"),
        ("subjects_by_sex", "sex", "subjects", "Subjects by sex"),
    ]
    for ax, (key, label_col, value_col, title) in zip(axes, panels):
        df = breakdowns[key]
        colors = (
            [RESPONSE_COLORS.get(v, "#2a78d6") for v in df[label_col]]
            if key == "subjects_by_response"
            else ["#2a78d6"] * len(df)
        )
        bars = ax.bar(df[label_col].astype(str), df[value_col], color=colors, width=0.55)
        for rect, value in zip(bars, df[value_col]):
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                rect.get_height(),
                f"{value:,}",
                ha="center",
                va="bottom",
                fontsize=10,
                color=INK,
            )
        ax.set_title(title, color=INK, loc="left", fontsize=11)
        ax.set_ylabel(value_col, color=INK_MUTED)
        ax.margins(y=0.18)
        _style(ax)

    fig.suptitle(
        "Baseline melanoma PBMC samples from miraclib-treated patients",
        color=INK,
        x=0.01,
        ha="left",
        fontsize=13,
    )
    fig.tight_layout()

    path = Path(path or config.FIGURE_DIR / "baseline_subset_breakdown.png")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, facecolor=SURFACE)
    plt.close(fig)
    return path
