"""End-to-end pipeline: load the database, then write every Part 2-4 output.

Run with `python -m teiko.pipeline` or, from the repo root, `make pipeline`.
"""
from __future__ import annotations

import json
import sys

from . import analysis, config, db, plots


def _write(df, name: str) -> None:
    path = config.TABLE_DIR / f"{name}.csv"
    df.to_csv(path, index=False)
    print(f"  wrote {path.relative_to(config.ROOT)}  ({len(df):,} rows)")


def run() -> dict:
    config.TABLE_DIR.mkdir(parents=True, exist_ok=True)
    config.FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    if not config.DB_PATH.exists():
        raise SystemExit(f"{config.DB_PATH} not found -- run `python load_data.py` first.")

    # ---- Part 2 ---------------------------------------------------------
    print("Part 2: relative frequency summary")
    summary = analysis.frequency_summary()
    _write(summary, "part2_cell_frequency_summary")

    # ---- Part 3 ---------------------------------------------------------
    print("Part 3: responders vs non-responders (melanoma, miraclib, PBMC)")
    cohort = analysis.responder_cohort()
    _write(cohort, "part3_responder_cohort_frequencies")

    by_sample = analysis.compare_responders(cohort, unit="sample")
    by_subject = analysis.compare_responders(cohort, unit="subject")
    _write(by_sample, "part3_statistics_by_sample")
    _write(by_subject, "part3_statistics_by_subject")

    narrative = analysis.significance_narrative(by_sample)
    print(f"  {narrative}")

    fig1 = plots.responder_boxplot(cohort, by_sample)
    print(f"  wrote {fig1.relative_to(config.ROOT)}")

    # ---- Part 4 ---------------------------------------------------------
    print("Part 4: baseline subset")
    subset = analysis.baseline_subset()
    _write(subset, "part4_baseline_samples")

    breakdowns = analysis.baseline_breakdowns(subset)
    for name, table in breakdowns.items():
        _write(table, f"part4_{name}")

    fig2 = plots.baseline_breakdown_figure(breakdowns)
    print(f"  wrote {fig2.relative_to(config.ROOT)}")

    # ---- machine-readable summary of the run ----------------------------
    report = {
        "samples": int(summary["sample"].nunique()),
        "populations": int(summary["population"].nunique()),
        "responder_cohort_samples": int(cohort["sample"].nunique()),
        "responder_cohort_subjects": int(cohort["subject_id"].nunique()),
        "significant_populations_by_sample": by_sample.loc[
            by_sample["significant"], "population"
        ].tolist(),
        "significant_populations_by_subject": by_subject.loc[
            by_subject["significant"], "population"
        ].tolist(),
        "narrative": narrative,
        "baseline_samples": int(len(subset)),
        "baseline_subjects": int(subset["subject_id"].nunique()),
        "baseline_breakdowns": {k: v.to_dict("records") for k, v in breakdowns.items()},
    }
    report_path = config.OUTPUT_DIR / "run_report.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"  wrote {report_path.relative_to(config.ROOT)}")
    return report


def main() -> int:
    if not config.DB_PATH.exists():
        print("Database missing; initializing it first.")
        db.initialize()
    run()
    print("Pipeline complete.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
