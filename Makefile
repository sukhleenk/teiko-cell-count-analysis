# Loblaw Bio cell-count analysis
PYTHON ?= python3
PORT   ?= 8501
export PYTHONPATH := src

.PHONY: setup pipeline dashboard test clean

## Install all dependencies.
setup:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

## Initialize + load the database (Part 1), then generate every table and figure (Parts 2-4).
pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) -m teiko.pipeline

## Start the interactive dashboard.
dashboard:
	$(PYTHON) -m streamlit run src/teiko/dashboard.py \
		--server.port=$(PORT) --server.address=0.0.0.0 --server.headless=true

## Run the test suite.
test:
	$(PYTHON) -m pytest -q

## Remove generated artifacts.
clean:
	rm -f cell_counts.db
	rm -rf outputs/tables outputs/figures
