"""Correctness checks for the loader and the three analyses."""
from pathlib import Path

import pandas as pd
import pytest

from teiko import analysis, config, db

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def database(tmp_path_factory):
    path = tmp_path_factory.mktemp("db") / "test.db"
    db.initialize(path, config.CSV_PATH)
    return path


def test_row_counts_match_csv(database):
    csv = pd.read_csv(config.CSV_PATH)
    assert db.query("SELECT COUNT(*) n FROM sample", db_path=database).iloc[0]["n"] == len(csv)
    assert db.query("SELECT COUNT(*) n FROM subject", db_path=database).iloc[0]["n"] == csv[
        "subject"
    ].nunique()
    assert db.query("SELECT COUNT(*) n FROM cell_count", db_path=database).iloc[0]["n"] == len(
        csv
    ) * len(config.POPULATIONS)


def test_percentages_sum_to_100(database):
    summary = analysis.frequency_summary(database)
    totals = summary.groupby("sample")["percentage"].sum()
    assert totals.between(99.99, 100.01).all()


def test_total_count_matches_raw_sum(database):
    csv = pd.read_csv(config.CSV_PATH).set_index("sample")
    summary = analysis.frequency_summary(database).drop_duplicates("sample").set_index("sample")
    expected = csv[config.POPULATIONS].sum(axis=1)
    assert (summary["total_count"] == expected.loc[summary.index]).all()


def test_responder_cohort_is_scoped(database):
    cohort = analysis.responder_cohort(database)
    assert set(cohort["response"]) == {"yes", "no"}
    detail = db.query(
        "SELECT DISTINCT condition, treatment, sample_type FROM sample_detail d "
        "JOIN sample s ON s.sample_id = d.sample_id",
        db_path=database,
    )
    assert not cohort.empty and len(detail) > 1  # cohort is a strict subset


def test_statistics_have_expected_shape(database):
    results = analysis.compare_responders(analysis.responder_cohort(database))
    assert len(results) == len(config.POPULATIONS)
    assert (results["p_adjusted"] >= results["p_value"] - 1e-12).all()
    assert results["rank_biserial"].between(-1, 1).all()


def test_baseline_subset_filters(database):
    subset = analysis.baseline_subset(database)
    assert (subset["condition"] == config.MELANOMA).all()
    assert (subset["treatment"] == config.MIRACLIB).all()
    assert (subset["sample_type"] == config.PBMC).all()
    assert (subset["time_from_treatment_start"] == 0).all()

    breakdowns = analysis.baseline_breakdowns(subset)
    assert breakdowns["samples_per_project"]["samples"].sum() == len(subset)
    assert (
        breakdowns["subjects_by_sex"]["subjects"].sum()
        == breakdowns["subjects_by_response"]["subjects"].sum()
        == subset["subject_id"].nunique()
    )
