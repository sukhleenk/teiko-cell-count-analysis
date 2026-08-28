"""Analytical queries for Parts 2-4.

Every function reads from the SQLite database so that the database -- not a
pandas script -- remains the single source of truth for the data.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from . import config
from .db import query


# --------------------------------------------------------------------------
# Part 2: data overview
# --------------------------------------------------------------------------
def frequency_summary(db_path: Path | str = config.DB_PATH) -> pd.DataFrame:
    """Relative frequency of each population in each sample.

    Columns: sample, total_count, population, count, percentage.
    """
    return query(
        """
        SELECT sample, total_count, population, count, percentage
        FROM sample_population_frequency
        ORDER BY sample, population
        """,
        db_path=db_path,
    )


def frequency_with_metadata(db_path: Path | str = config.DB_PATH) -> pd.DataFrame:
    """Part 2 summary joined to sample/subject metadata (used by Parts 3-4)."""
    return query(
        """
        SELECT f.sample, f.total_count, f.population, f.count, f.percentage,
               d.subject_id, d.project_id, d.condition, d.age, d.sex,
               d.treatment, d.response, d.sample_type, d.time_from_treatment_start
        FROM sample_population_frequency f
        JOIN sample_detail d ON d.sample_id = f.sample
        ORDER BY f.sample, f.population
        """,
        db_path=db_path,
    )


# --------------------------------------------------------------------------
# Part 3: responders vs non-responders
# --------------------------------------------------------------------------
def responder_cohort(
    db_path: Path | str = config.DB_PATH,
    condition: str = config.MELANOMA,
    treatment: str = config.MIRACLIB,
    sample_type: str = config.PBMC,
) -> pd.DataFrame:
    """Per-sample population frequencies for the responder-comparison cohort.

    Scoped to `condition` patients on `treatment` with a known response,
    `sample_type` samples only.
    """
    return query(
        """
        SELECT f.sample, f.population, f.percentage, f.count, f.total_count,
               d.subject_id, d.project_id, d.response, d.sex, d.age,
               d.time_from_treatment_start
        FROM sample_population_frequency f
        JOIN sample_detail d ON d.sample_id = f.sample
        WHERE d.condition   = ?
          AND d.treatment   = ?
          AND d.sample_type = ?
          AND d.response IN ('yes', 'no')
        ORDER BY f.sample, f.population
        """,
        params=(condition, treatment, sample_type),
        db_path=db_path,
    )


def _describe(values: np.ndarray) -> dict[str, float]:
    return {
        "n": len(values),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
    }


def compare_responders(
    cohort: pd.DataFrame, alpha: float = config.ALPHA, unit: str = "sample"
) -> pd.DataFrame:
    """Test responders vs non-responders for each cell population.

    Two-sided Mann-Whitney U (rank based, so it makes no normality assumption
    about the percentages) per population, with Benjamini-Hochberg control of
    the false discovery rate across the five populations. Effect size is
    reported as the rank-biserial correlation plus the raw difference in means.

    `unit`:
        "sample"  -- every sample is an observation (as specified in Part 3).
        "subject" -- percentages are averaged per subject first, which removes
                     the pseudo-replication from repeat visits of one patient.
    """
    df = cohort.copy()
    if unit == "subject":
        df = (
            df.groupby(["subject_id", "population", "response"], as_index=False)["percentage"]
            .mean()
        )

    rows = []
    for population in sorted(df["population"].unique()):
        sub = df[df["population"] == population]
        yes = sub.loc[sub["response"] == "yes", "percentage"].to_numpy()
        no = sub.loc[sub["response"] == "no", "percentage"].to_numpy()
        if len(yes) < 2 or len(no) < 2:
            continue

        u_stat, p_value = stats.mannwhitneyu(yes, no, alternative="two-sided")
        # Rank-biserial r: +1 => responders always higher, -1 => always lower.
        rank_biserial = 2.0 * u_stat / (len(yes) * len(no)) - 1.0

        y, n = _describe(yes), _describe(no)
        rows.append(
            {
                "population": population,
                "n_responder": y["n"],
                "n_non_responder": n["n"],
                "mean_responder": y["mean"],
                "mean_non_responder": n["mean"],
                "median_responder": y["median"],
                "median_non_responder": n["median"],
                "std_responder": y["std"],
                "std_non_responder": n["std"],
                "mean_difference": y["mean"] - n["mean"],
                "rank_biserial": rank_biserial,
                "u_statistic": float(u_stat),
                "p_value": float(p_value),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    out["p_adjusted"] = multipletests(out["p_value"], alpha=alpha, method="fdr_bh")[1]
    out["significant"] = out["p_adjusted"] < alpha
    out["direction"] = np.where(
        out["mean_difference"] > 0, "higher in responders", "higher in non-responders"
    )
    out["unit"] = unit
    return out.sort_values("p_adjusted").reset_index(drop=True)


def significance_narrative(results: pd.DataFrame, alpha: float = config.ALPHA) -> str:
    """One-paragraph plain-language readout of `compare_responders`."""
    if results.empty:
        return "No populations had enough observations to test."

    hits = results[results["significant"]]
    lines = [
        f"Mann-Whitney U with Benjamini-Hochberg FDR control at alpha={alpha:g} "
        f"across {len(results)} cell populations "
        f"({int(results['n_responder'].iloc[0])} responder vs "
        f"{int(results['n_non_responder'].iloc[0])} non-responder observations)."
    ]
    if hits.empty:
        lines.append("No population showed a significant difference after correction.")
    else:
        names = ", ".join(
            f"{config.POPULATION_LABELS.get(r.population, r.population)} "
            f"({r.direction}, {r.mean_difference:+.2f} pp, q={r.p_adjusted:.2e})"
            for r in hits.itertuples()
        )
        lines.append(f"Significant: {names}.")
    misses = results[~results["significant"]]
    if not misses.empty:
        lines.append(
            "Not significant: "
            + ", ".join(
                f"{config.POPULATION_LABELS.get(r.population, r.population)} (q={r.p_adjusted:.2f})"
                for r in misses.itertuples()
            )
            + "."
        )
    return " ".join(lines)


# --------------------------------------------------------------------------
# Part 4: baseline subset
# --------------------------------------------------------------------------
def baseline_subset(
    db_path: Path | str = config.DB_PATH,
    condition: str = config.MELANOMA,
    treatment: str = config.MIRACLIB,
    sample_type: str = config.PBMC,
    time_point: int = config.BASELINE_DAY,
) -> pd.DataFrame:
    """Melanoma PBMC samples at baseline from miraclib-treated patients."""
    return query(
        """
        SELECT sample_id AS sample, subject_id, project_id, condition, age, sex,
               treatment, response, sample_type, time_from_treatment_start
        FROM sample_detail
        WHERE condition   = ?
          AND treatment   = ?
          AND sample_type = ?
          AND time_from_treatment_start = ?
        ORDER BY sample_id
        """,
        params=(condition, treatment, sample_type, time_point),
        db_path=db_path,
    )


def _subject_breakdown(subset: pd.DataFrame, column: str, label: str) -> pd.DataFrame:
    """Count distinct subjects per value of `column` (not samples)."""
    subjects = subset.drop_duplicates("subject_id")
    out = (
        subjects.groupby(column, dropna=False)["subject_id"]
        .nunique()
        .reset_index(name="subjects")
        .rename(columns={column: label})
    )
    out[label] = out[label].fillna("unknown")
    return out.sort_values(label).reset_index(drop=True)


def baseline_breakdowns(subset: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """The three Part 4 roll-ups: samples per project, and subjects by
    response and by sex."""
    per_project = (
        subset.groupby("project_id")
        .agg(samples=("sample", "nunique"), subjects=("subject_id", "nunique"))
        .reset_index()
        .rename(columns={"project_id": "project"})
        .sort_values("project")
        .reset_index(drop=True)
    )
    return {
        "samples_per_project": per_project,
        "subjects_by_response": _subject_breakdown(subset, "response", "response"),
        "subjects_by_sex": _subject_breakdown(subset, "sex", "sex"),
    }
