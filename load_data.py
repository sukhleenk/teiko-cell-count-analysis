#!/usr/bin/env python3
"""Initialize the SQLite database and load cell-count.csv into it.

Usage:  python load_data.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from teiko import config, db  # noqa: E402


def main() -> int:
    if not config.CSV_PATH.exists():
        print(f"ERROR: {config.CSV_PATH} not found.", file=sys.stderr)
        return 1

    print(f"Loading {config.CSV_PATH.name} -> {config.DB_PATH.name}")
    counts = db.initialize(config.DB_PATH, config.CSV_PATH)
    for table, n in counts.items():
        print(f"  {table:<16} {n:>8,} rows")

    total = db.query("SELECT COUNT(*) AS n FROM sample_population_frequency").iloc[0]["n"]
    print(f"  {'(freq view)':<16} {total:>8,} rows")
    print(f"Done. Database written to {config.DB_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
