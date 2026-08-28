# Loblaw Bio — Immune Cell Population Analysis

Analysis of `cell-count.csv`: a SQLite-backed pipeline that loads the trial data,
computes per-sample relative frequencies of five immune cell populations, tests
whether those frequencies differ between miraclib responders and non-responders,
and profiles the baseline melanoma cohort. Results are browsable in an
interactive dashboard.

**Dashboard (hosted, read-only):**
<https://claude.ai/code/artifact/529456ad-eacd-4644-b541-1464a17e5bf0>

**Dashboard (interactive, local):** `make dashboard` &rarr; <http://localhost:8501>

---

## Quick start (GitHub Codespaces or any machine with Python 3.10+)

```bash
make setup      # install dependencies from requirements.txt
make pipeline   # build the database, then write every table and figure
make dashboard  # start the interactive dashboard on port 8501
```

`make pipeline` is fully non-interactive and runs:

1. `python load_data.py` — creates `cell_counts.db` in the repo root and loads all
   10,500 rows of `cell-count.csv`.
2. `python -m teiko.pipeline` — writes every Part 2–4 table to `outputs/tables/`,
   every figure to `outputs/figures/`, and a machine-readable `outputs/run_report.json`.

Optional: `make test` runs the test suite; `make clean` deletes generated artifacts.

In Codespaces, `make dashboard` binds to `0.0.0.0:8501`; Codespaces will offer to
forward the port — open the forwarded URL. Override the port with
`make dashboard PORT=8080`.

---

## Part 1 — Database schema

`load_data.py` builds this schema (full DDL in `src/teiko/schema.sql`):

```
project ──┐
          │
     subject (subject_id, project_id, condition_id, age, sex, treatment_id, response)
          │
       sample (sample_id, subject_id, sample_type_id, time_from_treatment_start)
          │
    cell_count (sample_id, population_id, count)          <- fact table, one row
          │                                                  per (sample, population)
   cell_population (population_id, name)

lookup tables: condition, treatment, sample_type, cell_population
views:         sample_detail, sample_population_frequency
```

### Rationale

**The counts are stored long, not wide.** The CSV has one column per cell
population. The database instead stores one row per `(sample, population)` pair.
This is the single most important choice in the schema:

- Adding a sixth or twentieth population is an `INSERT`, never an `ALTER TABLE`,
  and no downstream query or dashboard code changes.
- "Frequency of every population in every sample" becomes a `GROUP BY`, which is
  what the view `sample_population_frequency` does — rather than five hand-written
  column expressions that must be edited whenever the panel changes.
- Sparse panels (an assay that only measured three populations) cost nothing;
  wide tables would store NULLs.

**Subject attributes live on `subject`, not on every sample row.** In the source
CSV, condition/age/sex/treatment/response are repeated on all three visits of a
patient — 10,500 copies of 3,500 facts, and three chances to disagree. Normalising
them onto `subject` makes a contradiction unrepresentable, and the loader asserts
that every subject's attributes are consistent before inserting (a hard error if
the CSV ever violates it).

**Categoricals are lookup tables.** `condition`, `treatment`, `sample_type` and
`cell_population` keep the vocabularies closed, keep the fact rows narrow
(integers, not repeated strings), and give a natural place to hang future
attributes (a drug's mechanism, an assay's panel version) without touching the
fact table.

**Views encode the definitions.** `sample_population_frequency` is *the*
definition of "relative frequency" — the pipeline, the tests and the dashboard all
read it, so there is no way for two consumers to disagree about the denominator.
`sample_detail` flattens the joins for ad-hoc querying.

**Treatment and response are modelled on the subject** because in this dataset a
subject has exactly one arm and one outcome (verified during loading). For a
crossover or multi-line trial, this is the one place the schema would move: a
`subject_treatment(subject_id, treatment_id, line, start_date, response)` table,
with `sample` gaining a foreign key to it. Nothing else in the design changes.

### Scaling to hundreds of projects and thousands of samples

- **Row counts stay modest.** 5 populations × 10,500 samples = 52,500 fact rows
  today. A hundred projects at the same density is a few million rows — still
  comfortable for SQLite, and the exact same DDL ports to Postgres.
- **Indexes match the access patterns.** The queries this analysis actually runs
  are "filter a cohort, then aggregate": `idx_subject_cohort`
  `(condition_id, treatment_id, response)` serves the cohort filter,
  `idx_sample_type_time (sample_type_id, time_from_treatment_start)` serves the
  baseline/timepoint filter, and `idx_cell_count_pop` serves per-population
  aggregation.
- **The natural growth axis is rows, not columns.** New projects, new subjects,
  new timepoints, new populations, and new assay types are all inserts. That is
  what keeps the schema stable as the data grows.
- **Beyond SQLite**, the migration path is mechanical: move to Postgres,
  range-partition `cell_count` by project or load date, and — for analytics that
  scan the whole fact table — materialise `sample_population_frequency` and
  refresh it on load, since it is a pure function of `cell_count`.
- **Other analytics fit without redesign.** Longitudinal trajectories, per-project
  QC, batch effects, and cross-treatment comparisons are all filters and group-bys
  over the same three tables. Additional assay families (e.g. cytokines) attach as
  sibling fact tables keyed on `sample_id`, reusing the whole subject/sample spine.

---

## Part 2 — Data overview

`outputs/tables/part2_cell_frequency_summary.csv` — 52,500 rows, one per
(sample, population), with columns `sample`, `total_count`, `population`,
`count`, `percentage`. The total count is the sum across all five populations for
that sample; the percentage is `100 × count / total_count`. A test asserts the
five percentages sum to 100 for every sample.

The dashboard's **Overview** tab shows the same table with project / condition /
sample-type / id filters and a CSV download.

---

## Part 3 — Responders vs non-responders

Cohort: **melanoma** patients on **miraclib**, **PBMC** samples only, with a
recorded response — 1,968 samples from 656 subjects (993 responder / 975
non-responder samples).

**Method.** For each of the five populations, a two-sided **Mann-Whitney U test**
compares the relative frequencies of responders against non-responders. It is
rank-based, so it assumes nothing about the shape of the percentage distributions
(which are bounded at 0 and 100 and mildly skewed). Because five populations are
tested, p-values are corrected with **Benjamini-Hochberg FDR control** at
α = 0.05. Effect size is reported two ways: the difference in means in percentage
points, and the **rank-biserial correlation** (+1 = responders always higher,
−1 = always lower).

**Result (`outputs/tables/part3_statistics_by_sample.csv`):**

| population | mean responder | mean non-responder | difference (pp) | p | q (BH) | significant |
|---|---|---|---|---|---|---|
| CD4+ T cell | 30.54 | 29.90 | **+0.64** | 0.013 | 0.067 | no |
| B cell      |  9.80 |  9.996 | −0.20 | 0.056 | 0.139 | no |
| NK cell     | 14.84 | 15.07 | −0.23 | 0.121 | 0.202 | no |
| Monocyte    | 19.94 | 20.08 | −0.14 | 0.163 | 0.204 | no |
| CD8+ T cell | 24.88 | 24.94 | −0.06 | 0.639 | 0.639 | no |

**What to tell Yah:** *no* cell population shows a statistically significant
difference in relative frequency between responders and non-responders after
correcting for testing five populations. CD4+ T cells are the only suggestive
signal — responders average 0.64 percentage points higher, nominally p = 0.013,
but q = 0.067, above the 5% FDR threshold. The effect sizes are small in absolute
terms (rank-biserial +0.064 for CD4+, i.e. barely above chance separation), so
even if a larger trial confirmed it, CD4+ frequency alone would be a weak
predictor of miraclib response. This dataset does not support a claim that
baseline population frequencies predict response.

**A caveat worth stating up front, before Yah raises it.** Each patient
contributes up to three samples (days 0, 7, 14), so treating samples as
independent observations overstates the effective sample size. The pipeline
therefore also reports a subject-level sensitivity analysis
(`part3_statistics_by_subject.csv`), averaging each patient's visits first: the
conclusion is unchanged — CD4+ T cell p = 0.012, q = 0.062, nothing significant.
The dashboard exposes both via the *Observation unit* selector.

**Figure:** `outputs/figures/responder_vs_nonresponder_boxplot.png` — grouped
boxplot of relative frequency per population, responders vs non-responders. The
dashboard renders an interactive version.

---

## Part 4 — Baseline subset

Melanoma **PBMC** samples at **`time_from_treatment_start = 0`** from
**miraclib**-treated patients: **656 samples from 656 subjects** (one baseline
sample per subject). Full listing in `outputs/tables/part4_baseline_samples.csv`.

| Samples per project | | Subjects by response | | Subjects by sex | |
|---|---|---|---|---|---|
| prj1 | 384 | responders (yes) | 331 | male | 344 |
| prj3 | 272 | non-responders (no) | 325 | female | 312 |

(prj2 contributes no samples to this subset.) Response and sex are counted per
**subject**, as asked; the project roll-up reports both samples and subjects.
The dashboard's **Baseline subset** tab lets Bob change condition, treatment and
timepoint and recomputes all three breakdowns live.

---

## Dashboard

Hosted read-only copy of the results:
<https://claude.ai/code/artifact/529456ad-eacd-4644-b541-1464a17e5bf0>
(the same page is checked in at `docs/hosted_dashboard.html`).

The full interactive dashboard, which queries the live database:

```bash
make dashboard      # http://localhost:8501  (PORT=… to change)
```

Four tabs:

- **Overview (Part 2)** — the full frequency summary, filterable by project,
  condition, sample type and id, with a mean-frequency bar chart and CSV download.
- **Responders (Part 3)** — cohort selectors (condition / treatment / sample type /
  observation unit), the grouped boxplot, the statistics table, and a
  plain-language readout of what is and isn't significant.
- **Baseline subset (Part 4)** — the filtered sample list plus the three
  breakdowns, all recomputed from the database as the filters change.
- **Database** — object/row inventory and the live schema DDL.

Charts use a two-colour categorical palette (blue = responder, orange =
non-responder) validated for colour-vision deficiency separation; every chart also
carries a legend and direct labels, so identity is never conveyed by colour alone.

---

## Code structure

```
.
├── load_data.py                  # Part 1 entry point (required name/location)
├── Makefile                      # setup / pipeline / dashboard / test / clean
├── requirements.txt
├── cell-count.csv                # input
├── cell_counts.db                # generated by load_data.py
├── outputs/
│   ├── tables/                   # every Part 2-4 table as CSV
│   ├── figures/                  # static PNGs
│   └── run_report.json           # machine-readable summary of the run
├── docs/hosted_dashboard.html    # self-contained hosted copy of the results
├── src/teiko/
│   ├── config.py                 # paths and constants -- no magic strings elsewhere
│   ├── schema.sql                # the schema, as SQL
│   ├── db.py                     # connect / create schema / load CSV / query
│   ├── analysis.py               # Parts 2-4 as pure, parameterised functions
│   ├── plots.py                  # static matplotlib figures
│   ├── pipeline.py               # orchestrates: analysis -> tables + figures
│   └── dashboard.py              # Streamlit UI
└── tests/test_pipeline.py        # loader and analysis correctness checks
```

**Why this shape:**

- **`load_data.py` is a thin entry point.** It satisfies the required
  name-and-location contract while the real logic lives in `src/teiko/db.py`,
  importable by the pipeline, the tests and the dashboard.
- **`analysis.py` returns DataFrames and takes no side effects.** Every function
  is parameterised on cohort (`condition`, `treatment`, `sample_type`, timepoint)
  rather than hard-coding "melanoma / miraclib / PBMC". That is exactly why the
  dashboard can offer those as dropdowns without duplicating a line of analysis
  code — the CLI pipeline and the UI call the same functions with different
  arguments.
- **The database is the source of truth, not a pandas script.** Filtering and the
  frequency calculation happen in SQL (`sample_population_frequency`,
  `sample_detail`), so the pipeline, the tests and the dashboard cannot drift
  apart on what "relative frequency" or "the baseline cohort" means. It also means
  the analysis scales past what fits in memory.
- **Presentation is separated from computation.** `plots.py` (static, for the
  pipeline) and `dashboard.py` (interactive) both consume `analysis.py` output;
  neither contains statistics.
- **Tests guard the invariants that matter** — row counts match the CSV,
  percentages sum to 100, `total_count` equals the raw row sum, the cohort filters
  are actually applied, adjusted p-values dominate raw ones, and the Part 4
  breakdowns are internally consistent.

## Requirements

Python 3.10+ and the packages in `requirements.txt` (pandas, numpy, scipy,
statsmodels, matplotlib, seaborn, streamlit, plotly, pytest). SQLite ships with
Python — no database server needed.
