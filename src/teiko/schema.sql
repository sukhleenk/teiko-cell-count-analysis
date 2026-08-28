-- Relational schema for Loblaw Bio clinical trial cell-count data.
--
-- Design notes:
--   * One row per (sample, population) in cell_count -- a long/tidy fact table.
--     Adding a sixth immune population is a data insert, never a schema change.
--   * Subject-level attributes (condition, demographics, trial arm, outcome) are
--     stored once on `subject` instead of being repeated on every sample row.
--   * Lookup tables keep the categorical vocabularies (condition, treatment,
--     sample type, population) closed and cheap to join/index.

PRAGMA foreign_keys = ON;

DROP VIEW  IF EXISTS sample_population_frequency;
DROP VIEW  IF EXISTS sample_detail;
DROP TABLE IF EXISTS cell_count;
DROP TABLE IF EXISTS sample;
DROP TABLE IF EXISTS subject;
DROP TABLE IF EXISTS project;
DROP TABLE IF EXISTS cell_population;
DROP TABLE IF EXISTS condition;
DROP TABLE IF EXISTS treatment;
DROP TABLE IF EXISTS sample_type;

-- ---------------------------------------------------------------- dimensions
CREATE TABLE project (
    project_id   TEXT PRIMARY KEY               -- e.g. 'prj1'
);

CREATE TABLE condition (
    condition_id INTEGER PRIMARY KEY,
    name         TEXT NOT NULL UNIQUE           -- melanoma, carcinoma, healthy
);

CREATE TABLE treatment (
    treatment_id INTEGER PRIMARY KEY,
    name         TEXT NOT NULL UNIQUE           -- miraclib, phauximab, none
);

CREATE TABLE sample_type (
    sample_type_id INTEGER PRIMARY KEY,
    name           TEXT NOT NULL UNIQUE         -- PBMC, WB
);

CREATE TABLE cell_population (
    population_id INTEGER PRIMARY KEY,
    name          TEXT NOT NULL UNIQUE          -- b_cell, cd8_t_cell, ...
);

-- ------------------------------------------------------------------ entities
CREATE TABLE subject (
    subject_id   TEXT PRIMARY KEY,              -- e.g. 'sbj000'
    project_id   TEXT    NOT NULL REFERENCES project(project_id),
    condition_id INTEGER NOT NULL REFERENCES condition(condition_id),
    age          INTEGER CHECK (age IS NULL OR age >= 0),
    sex          TEXT    CHECK (sex IN ('M', 'F') OR sex IS NULL),
    treatment_id INTEGER REFERENCES treatment(treatment_id),
    response     TEXT    CHECK (response IN ('yes', 'no') OR response IS NULL)
);

CREATE TABLE sample (
    sample_id                 TEXT PRIMARY KEY, -- e.g. 'sample00000'
    subject_id                TEXT    NOT NULL REFERENCES subject(subject_id),
    sample_type_id            INTEGER NOT NULL REFERENCES sample_type(sample_type_id),
    time_from_treatment_start INTEGER            -- days; 0 == baseline
);

-- --------------------------------------------------------------------- facts
CREATE TABLE cell_count (
    sample_id     TEXT    NOT NULL REFERENCES sample(sample_id) ON DELETE CASCADE,
    population_id INTEGER NOT NULL REFERENCES cell_population(population_id),
    count         INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

-- ------------------------------------------------------------------- indexes
CREATE INDEX idx_subject_project   ON subject (project_id);
CREATE INDEX idx_subject_cohort    ON subject (condition_id, treatment_id, response);
CREATE INDEX idx_sample_subject    ON sample  (subject_id);
CREATE INDEX idx_sample_type_time  ON sample  (sample_type_id, time_from_treatment_start);
CREATE INDEX idx_cell_count_pop    ON cell_count (population_id);

-- --------------------------------------------------------------------- views
-- Flat, analysis-ready row per sample with all descriptive attributes resolved.
CREATE VIEW sample_detail AS
SELECT s.sample_id,
       sub.subject_id,
       sub.project_id,
       c.name  AS condition,
       sub.age,
       sub.sex,
       t.name  AS treatment,
       sub.response,
       st.name AS sample_type,
       s.time_from_treatment_start
FROM sample s
JOIN subject     sub ON sub.subject_id     = s.subject_id
JOIN condition   c   ON c.condition_id     = sub.condition_id
JOIN sample_type st  ON st.sample_type_id  = s.sample_type_id
LEFT JOIN treatment t ON t.treatment_id    = sub.treatment_id;

-- Part 2 summary table, computed in SQL so every consumer sees one definition
-- of "relative frequency".
CREATE VIEW sample_population_frequency AS
SELECT cc.sample_id                                        AS sample,
       tot.total_count                                     AS total_count,
       p.name                                              AS population,
       cc.count                                            AS count,
       ROUND(100.0 * cc.count / tot.total_count, 6)        AS percentage
FROM cell_count cc
JOIN cell_population p ON p.population_id = cc.population_id
JOIN (SELECT sample_id, SUM(count) AS total_count
      FROM cell_count
      GROUP BY sample_id) tot ON tot.sample_id = cc.sample_id
WHERE tot.total_count > 0;
