"""Project-wide paths and constants."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
SCHEMA_PATH = Path(__file__).resolve().parent / "schema.sql"

OUTPUT_DIR = ROOT / "outputs"
TABLE_DIR = OUTPUT_DIR / "tables"
FIGURE_DIR = OUTPUT_DIR / "figures"

# The five immune cell populations, in the order they appear in the CSV.
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

# Pretty labels for plots and dashboard tables.
POPULATION_LABELS = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8+ T cell",
    "cd4_t_cell": "CD4+ T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}

# Cohort that Part 3 and Part 4 are scoped to.
MELANOMA = "melanoma"
MIRACLIB = "miraclib"
PBMC = "PBMC"
BASELINE_DAY = 0

ALPHA = 0.05
