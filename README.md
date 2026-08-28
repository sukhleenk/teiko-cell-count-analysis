# Loblaw Bio: Immune Cell Population Analysis

A SQLite-backed pipeline over `cell-count.csv`. It loads the trial data, computes
per-sample relative frequencies of five immune cell populations, tests whether
those frequencies differ between miraclib responders and non-responders, and
profiles the baseline melanoma cohort. Results are browsable in an interactive
dashboard.

**Dashboard link:** <http://localhost:8501>, served by `make dashboard`.
In Codespaces, open the forwarded URL that appears for port 8501 instead.

## Quick start

Requires Python 3.10+. Everything else is installed by `make setup`.

```bash
make setup      # install dependencies from requirements.txt
make pipeline   # build the database, then write every table and figure
make dashboard  # start the dashboard on port 8501
```

`make pipeline` is fully non-interactive and runs two steps:

1. `python load_data.py` creates `cell_counts.db` in the repo root and loads all
   10,500 rows of `cell-count.csv`.
2. `python -m teiko.pipeline` writes every Part 2-4 table to `outputs/tables/`
   and both figures to `outputs/figures/`.

`make test` runs the test suite. `make clean` deletes the generated database,
tables and figures.

In Codespaces, `make dashboard` binds to `0.0.0.0:8501` and Codespaces offers to
forward the port; open the forwarded URL. Use `make dashboard PORT=8080` to
change it.

## Part 1: database schema

`load_data.py` builds this schema. The full DDL is in `src/teiko/schema.sql`.

```
project
   |
subject ......... subject_id, project_id, condition_id, age, sex,
   |              treatment_id, response
sample .......... sample_id, subject_id, sample_type_id,
   |              time_from_treatment_start
cell_count ...... sample_id, population_id, count
                  (fact table: one row per sample and population)

lookups: condition, treatment, sample_type, cell_population
views:   sample_detail, sample_population_frequency
```

### Why this design

*Counts are stored long, not wide.* The CSV has one column per cell population;
the database stores one row per `(sample, population)` pair instead. This is the
most consequential choice in the schema. Adding a sixth or twentieth population
becomes an `INSERT` rather than an `ALTER TABLE`, and no downstream query or
dashboard code changes. "Frequency of every population in every sample" becomes a
`GROUP BY`, which is what the `sample_population_frequency` view does, instead of
five hand-written column expressions that must be edited whenever the assay panel
changes. Sparse panels, where a run only measured three populations, cost nothing;
a wide table would store NULLs.

*Subject attributes live on `subject`, not on every sample row.* In the source
CSV, condition, age, sex, treatment and response are repeated across all three
visits of a patient: 10,500 copies of 3,500 facts, and three chances to disagree.
Normalising them onto `subject` makes a contradiction unrepresentable, and the
loader raises a hard error if the CSV ever violates it.

*Categoricals are lookup tables.* `condition`, `treatment`, `sample_type` and
`cell_population` keep the vocabularies closed, keep fact rows narrow (integers
rather than repeated strings), and give a place to hang future attributes such as
a drug's mechanism or an assay's panel version without touching the fact table.

*Views encode the definitions.* `sample_population_frequency` is the definition of
"relative frequency". The pipeline, the tests and the dashboard all read it, so
two consumers cannot drift apart on what the denominator is. `sample_detail`
flattens the joins for ad-hoc querying.

*Treatment and response are modelled on the subject* because in this dataset a
subject has exactly one arm and one outcome, which the loader verifies. For a
crossover or multi-line trial this is the one place the schema would move: a
`subject_treatment(subject_id, treatment_id, line, start_date, response)` table,
with `sample` gaining a foreign key to it. Nothing else in the design changes.

### Scaling

Five populations across 10,500 samples is 52,500 fact rows today. A hundred
projects at the same density is a few million rows, still comfortable for SQLite,
and the same DDL ports to Postgres unchanged.

The indexes match the access pattern this analysis actually uses, which is "filter
a cohort, then aggregate". `idx_subject_cohort (condition_id, treatment_id,
response)` serves the cohort filter, `idx_sample_type_time (sample_type_id,
time_from_treatment_start)` serves the timepoint filter, and `idx_cell_count_pop`
serves per-population aggregation.

Growth happens in rows, not columns. New projects, subjects, timepoints,
populations and assay types are all inserts, which is what keeps the schema stable
as the data grows. Beyond SQLite the migration is mechanical: move to Postgres,
range-partition `cell_count` by project or load date, and materialise
`sample_population_frequency` for analytics that scan the whole fact table, since
it is a pure function of `cell_count`. Longitudinal trajectories, per-project QC,
batch effects and cross-treatment comparisons are all filters and group-bys over
the same three tables. Additional assay families such as cytokines attach as
sibling fact tables keyed on `sample_id`, reusing the whole subject/sample spine.

## Part 2: data overview

`outputs/tables/part2_cell_frequency_summary.csv` has 52,500 rows, one per
(sample, population), with columns `sample`, `total_count`, `population`, `count`
and `percentage`. The total count is the sum across all five populations for that
sample; the percentage is `100 * count / total_count`. A test asserts that the
five percentages sum to 100 for every sample.

The dashboard's Overview tab shows the same table with project, condition, sample
type and id filters, plus a CSV download.

## Part 3: responders vs non-responders

Cohort: melanoma patients on miraclib, PBMC samples only, with a recorded
response. That is 1,968 samples from 656 subjects (993 responder and 975
non-responder samples).

### Method

For each of the five populations, a two-sided Mann-Whitney U test compares the
relative frequencies of responders against non-responders. It is rank-based, so
it assumes nothing about the shape of the percentage distributions, which are
bounded at 0 and 100 and mildly skewed. Because five populations are tested,
p-values are corrected with Benjamini-Hochberg FDR control at alpha = 0.05. Effect
size is reported two ways: the difference in means in percentage points, and the
rank-biserial correlation, where +1 means responders are always higher and -1
always lower.

### Result

From `outputs/tables/part3_statistics_by_sample.csv`:

| population | mean responder | mean non-responder | difference (pp) | p | q (BH) | significant |
|---|---|---|---|---|---|---|
| CD4+ T cell | 30.54 | 29.90 | +0.64 | 0.013 | 0.067 | no |
| B cell | 9.80 | 10.00 | -0.20 | 0.056 | 0.139 | no |
| NK cell | 14.84 | 15.07 | -0.23 | 0.121 | 0.202 | no |
| Monocyte | 19.94 | 20.08 | -0.14 | 0.163 | 0.204 | no |
| CD8+ T cell | 24.88 | 24.94 | -0.06 | 0.639 | 0.639 | no |

### Interpretation

No cell population shows a statistically significant difference in relative
frequency between responders and non-responders after correcting for testing five
populations. CD4+ T cells are the only suggestive signal: responders average 0.64
percentage points higher, nominally p = 0.013, but q = 0.067, above the 5% FDR
threshold. The effect sizes are small in absolute terms; a rank-biserial of +0.064
for CD4+ is barely above chance separation. Even if a larger trial confirmed it,
CD4+ frequency alone would be a weak predictor of miraclib response. This dataset
does not support a claim that population frequencies predict response.

### Caveat: repeated measures

Each patient contributes up to three samples (days 0, 7 and 14), so treating
samples as independent observations overstates the effective sample size. The
pipeline therefore also reports a subject-level sensitivity analysis in
`part3_statistics_by_subject.csv`, averaging each patient's visits first. The
conclusion is unchanged: CD4+ T cell p = 0.012, q = 0.062, nothing significant.
The dashboard exposes both through the "Observation unit" selector.

The figure is `outputs/figures/responder_vs_nonresponder_boxplot.png`, a grouped
boxplot of relative frequency per population. The dashboard renders an interactive
version.

## Part 4: baseline subset

Melanoma PBMC samples at `time_from_treatment_start = 0` from miraclib-treated
patients: 656 samples from 656 subjects, one baseline sample per subject. Full
listing in `outputs/tables/part4_baseline_samples.csv`.

| samples per project | | subjects by response | | subjects by sex | |
|---|---|---|---|---|---|
| prj1 | 384 | yes | 331 | M | 344 |
| prj3 | 272 | no | 325 | F | 312 |

prj2 contributes nothing to this subset because it is whole-blood only; it has
melanoma patients on miraclib, but no PBMC samples. Response and sex are counted
per subject, as asked; the project roll-up reports both samples and subjects.

The dashboard's Baseline subset tab changes condition, treatment and timepoint and
recomputes all three breakdowns live.

## Dashboard

```bash
make dashboard
```

Serves the dashboard at <http://localhost:8501>. In Codespaces, port 8501 is
forwarded automatically; open the forwarded URL from the Ports panel. Use
`make dashboard PORT=8080` to serve it elsewhere. Run `make pipeline` first, or
at least `python load_data.py`, so the database exists.

Three tabs, one per analysis:

- **Overview (Part 2)** is the full frequency summary, filterable by project,
  condition, sample type and id, with a mean-frequency bar chart and CSV download.
- **Responders (Part 3)** has cohort selectors (condition, treatment, sample type,
  observation unit), the grouped boxplot, the statistics table, and a plain-language
  readout of what is and is not significant.
- **Baseline subset (Part 4)** has the filtered sample list and the three
  breakdowns, all recomputed from the database as the filters change.

Charts use two colours, blue for responders and orange for non-responders,
checked for colour-vision-deficiency separation. Every chart also carries a legend
and labels, so identity is never conveyed by colour alone.

## Code structure

```
.
├── load_data.py                  # Part 1 entry point (required name and location)
├── Makefile                      # setup / pipeline / dashboard / test / clean
├── requirements.txt
├── cell-count.csv                # input
├── cell_counts.db                # generated by load_data.py
├── outputs/
│   ├── tables/                   # every Part 2-4 table as CSV
│   └── figures/                  # static PNGs
├── src/teiko/
│   ├── config.py                 # paths and constants
│   ├── schema.sql                # the schema, as SQL
│   ├── db.py                     # connect, create schema, load CSV, query
│   ├── analysis.py               # Parts 2-4 as pure functions
│   ├── plots.py                  # static matplotlib figures
│   ├── pipeline.py               # analysis to tables and figures
│   └── dashboard.py              # Streamlit UI
└── tests/test_pipeline.py        # loader and analysis correctness checks
```

`load_data.py` is a thin entry point. It satisfies the required name-and-location
contract while the real logic lives in `src/teiko/db.py`, importable by the
pipeline, the tests and the dashboard.

`analysis.py` returns DataFrames and has no side effects. Every function is
parameterised on cohort (condition, treatment, sample type, timepoint) rather than
hard-coding melanoma, miraclib and PBMC. That is why the dashboard can offer those
as dropdowns without duplicating a line of analysis code: the CLI pipeline and the
UI call the same functions with different arguments.

The database is the source of truth, not a pandas script. Filtering and the
frequency calculation happen in SQL, via `sample_population_frequency` and
`sample_detail`, so the pipeline, the tests and the dashboard cannot disagree
about what "relative frequency" or "the baseline cohort" means. It also means the
analysis is not bounded by what fits in memory.

Presentation is separated from computation. `plots.py` (static, for the pipeline)
and `dashboard.py` (interactive) both consume `analysis.py` output, and neither
contains statistics.

The tests guard the invariants that matter: row counts match the CSV, percentages
sum to 100, `total_count` equals the raw row sum, the cohort filters are actually
applied, adjusted p-values dominate raw ones, and the Part 4 breakdowns are
internally consistent.
