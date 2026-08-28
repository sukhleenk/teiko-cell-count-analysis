"""SQLite connection helpers and schema/loading logic."""
from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

from . import config


def connect(db_path: Path | str = config.DB_PATH) -> sqlite3.Connection:
    """Open a connection with foreign keys enforced and rows as tuples."""
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def create_schema(conn: sqlite3.Connection) -> None:
    """(Re)create all tables, indexes and views from schema.sql."""
    conn.executescript(config.SCHEMA_PATH.read_text())
    conn.commit()


def read_csv(csv_path: Path | str = config.CSV_PATH) -> pd.DataFrame:
    """Read the wide source CSV, normalising empty strings to NULL."""
    df = pd.read_csv(csv_path)
    missing = {"project", "subject", "sample"} - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path} is missing required columns: {sorted(missing)}")
    for col in ("response", "treatment", "sex", "condition"):
        if col in df.columns:
            df[col] = df[col].astype("object").where(df[col].notna(), None)
    return df


def _lookup(conn: sqlite3.Connection, table: str, values) -> dict[str, int]:
    """Insert distinct `values` into a lookup table and return name -> id."""
    names = sorted({v for v in values if v is not None and pd.notna(v)})
    conn.executemany(f"INSERT INTO {table} (name) VALUES (?)", [(n,) for n in names])
    return dict(conn.execute(f"SELECT name, {table}_id FROM {table}"))


def load_dataframe(conn: sqlite3.Connection, df: pd.DataFrame) -> dict[str, int]:
    """Populate every table from the wide dataframe. Returns row counts."""
    projects = sorted(df["project"].dropna().unique())
    conn.executemany("INSERT INTO project (project_id) VALUES (?)", [(p,) for p in projects])

    conditions = _lookup(conn, "condition", df["condition"])
    treatments = _lookup(conn, "treatment", df["treatment"])
    sample_types = _lookup(conn, "sample_type", df["sample_type"])

    populations = [p for p in config.POPULATIONS if p in df.columns]
    conn.executemany(
        "INSERT INTO cell_population (population_id, name) VALUES (?, ?)",
        list(enumerate(populations, start=1)),
    )
    pop_ids = {name: i for i, name in enumerate(populations, start=1)}

    subj_cols = ["subject", "project", "condition", "age", "sex", "treatment", "response"]
    subjects = df[subj_cols].drop_duplicates(subset="subject")
    conflicting = df[subj_cols].drop_duplicates()
    if len(conflicting) != len(subjects):
        dupes = conflicting[conflicting.duplicated("subject", keep=False)]["subject"].unique()
        raise ValueError(
            "Subjects carry conflicting attributes across rows: " f"{sorted(dupes)[:5]}"
        )
    conn.executemany(
        "INSERT INTO subject (subject_id, project_id, condition_id, age, sex, "
        "treatment_id, response) VALUES (?, ?, ?, ?, ?, ?, ?)",
        [
            (
                r.subject,
                r.project,
                conditions[r.condition],
                None if pd.isna(r.age) else int(r.age),
                r.sex,
                treatments.get(r.treatment),
                r.response,
            )
            for r in subjects.itertuples()
        ],
    )

    conn.executemany(
        "INSERT INTO sample (sample_id, subject_id, sample_type_id, "
        "time_from_treatment_start) VALUES (?, ?, ?, ?)",
        [
            (
                r.sample,
                r.subject,
                sample_types[r.sample_type],
                None if pd.isna(r.time_from_treatment_start) else int(r.time_from_treatment_start),
            )
            for r in df.itertuples()
        ],
    )

    long = df.melt(
        id_vars="sample", value_vars=populations, var_name="population", value_name="count"
    ).dropna(subset=["count"])
    conn.executemany(
        "INSERT INTO cell_count (sample_id, population_id, count) VALUES (?, ?, ?)",
        [(r.sample, pop_ids[r.population], int(r.count)) for r in long.itertuples()],
    )

    conn.commit()
    return {
        "project": len(projects),
        "subject": len(subjects),
        "sample": len(df),
        "cell_population": len(populations),
        "cell_count": len(long),
    }


def initialize(
    db_path: Path | str = config.DB_PATH, csv_path: Path | str = config.CSV_PATH
) -> dict[str, int]:
    """Create the schema and load the CSV. Idempotent: safe to re-run."""
    with connect(db_path) as conn:
        create_schema(conn)
        return load_dataframe(conn, read_csv(csv_path))


def query(sql: str, params: tuple | dict = (), db_path: Path | str = config.DB_PATH):
    """Run a read-only query and return a DataFrame."""
    with connect(db_path) as conn:
        return pd.read_sql_query(sql, conn, params=params)
