# OmicsIntegration-Immune: Multi-Omics Integration & Immune Discovery APIs

Consolidates 4 immune-profiling modalities - **scRNA-seq**, **spatial
(CODEX/Xenium-style)**, **flow cytometry**, and **TCR/BCR repertoire** - from
**55 simulated vaccine studies** into one unified SQLite store, runs real
statistical biomarker discovery against a known ground truth, and exposes
everything through a **researcher-facing Flask REST API** - a reusable
access layer, not just a one-off analysis notebook.

> **Disclaimer:** This is a portfolio/demo project. All study, sample,
> omics, and response data is **synthetically simulated** by
> `simulate_multiomics_data.py` with a fixed random seed. No real patient,
> subject, or sponsor data is used. Feature names (genes, cell-surface
> markers, spatial niches) are drawn from general immunology/vaccinology
> knowledge for realism, not copied from any real dataset or publication.

## Differentiation from sibling portfolio projects

This repository sits alongside other multi-omics / immune-profiling
portfolio projects (e.g. `VaccineMultiOmics`, `ImmuneProfiler-CellTherapy`,
`MultiOmicsIntegration`, `FlowCytometry-ImmuneCell-Phenotyping`). Where those
focus on the *integration pipeline itself* (harmonizing raw modality data
into one schema) or on *QA/validation frameworks* for a single modality,
**this project's distinguishing angle is the access layer**: once data is
consolidated, what do researchers actually query day-to-day? This repo
answers that with a real, running Flask API (`api/app.py`) - study/sample
lookup, per-modality retrieval, a biomarker ranking/lookup service, and a
cross-modality join endpoint - the kind of internal "data API" a
bioinformatics platform team would build on top of a data lake so other
scientists and tools don't need direct SQL/file access. The omics simulation
and biomarker discovery exist to give that API something real to serve, not
as the end deliverable themselves.

Conventions (project layout, simulate -> analyze -> `run_all.py` ->
`pytest`, honest-numbers README) follow this portfolio's
`CharacterizationPipeline` project.

## What it does

1. **`simulate_multiomics_data.py`** simulates **55 vaccine studies**
   (mRNA, protein-subunit, viral-vector, inactivated, live-attenuated; 8-16
   subjects each) across 4 modalities per subject:
   - **scRNA-seq**: pseudobulk expression of 40 immune genes
   - **Spatial (CODEX/Xenium-style)**: 15 lymph-node niche /
     cell-colocalization scores (germinal center density, Tfh-GC
     colocalization, etc.)
   - **Flow cytometry**: 20 immune cell subset frequencies (%)
   - **TCR/BCR repertoire**: 15 clonality/diversity/clonotype-frequency
     features

   A ground-truth set of **exactly 18 features, spread across all 4
   modalities**, is wired with known signed weights into a simulated
   continuous vaccine-response score (a neutralizing-titer-percentile
   proxy), with realistic noise on top. Everything else is pure noise. The
   ground truth is saved separately and used only to *validate* discovery,
   never to perform it. Everything is consolidated into one unified SQLite
   store, `data/omics.db` (`studies`, `samples`, `modality_scrna`,
   `modality_spatial`, `modality_flow`, `modality_tcr_bcr` tables).

2. **`discover_biomarkers.py`** runs real statistics: Spearman correlation
   of every one of the **90 candidate features** against the response
   score, Benjamini-Hochberg FDR correction across all of them, and ranks
   the significant hits. Validates the discovered set against the injected
   18-biomarker ground truth (precision/recall/F1), and writes a
   `biomarkers` table back into the SQLite store for the API to serve.

3. **`api/app.py`** is the researcher-facing Flask REST API:
   - `GET /api/studies`, `GET /api/studies/<study_id>` - study/sample query
   - `GET /api/samples`, `GET /api/samples/<sample_id>` - sample query,
     filterable by `study_id` / `responder`
   - `GET /api/modality/<modality>/<sample_id>` - per-modality raw data
     retrieval (`scrna`, `spatial`, `flow`, `tcr_bcr`)
   - `GET /api/biomarkers`, `GET /api/biomarkers/<feature>` - ranked
     biomarker lookup, filterable by modality
   - `GET /api/query/cross_modality/<sample_id>` - all 4 modalities joined
     for one sample in a single call
   - `GET /api/health`

4. **`run_all.py`** runs the whole pipeline end-to-end, including a live
   smoke test of the API via Flask's test client.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
pytest -q
python api/app.py   # serves the API at http://localhost:5050
```

## Results

_All numbers below are copied directly from an actual run of `python
run_all.py` followed by `pytest -q` - nothing is hardcoded. Fixed random
seed makes them reproducible._

- **55 vaccine studies, 667 samples** consolidated across 4 modalities into
  one unified SQLite store (`data/omics.db`).
- **90 candidate features** tested (40 scRNA genes + 15 spatial + 20 flow +
  15 TCR/BCR).
- **21 statistically significant biomarkers discovered** (Spearman
  correlation vs. vaccine-response score, Benjamini-Hochberg FDR q < 0.05).
- **Validated against 18 injected ground-truth biomarkers:
  precision = 0.86, recall = 1.00, F1 = 0.92.** All 18 true biomarkers were
  recovered (0 false negatives); 3 of 21 discovered hits (`GENE_CD38`,
  `GENE_STAT1`, `BCR_Vgene_skew`) were not part of the injected signal -
  expected false-discovery behavior at FDR 5% when testing 90 features, and
  reported honestly rather than rounded down to exactly 18.
- Discovered biomarkers span **all 4 modalities**: scRNA (e.g. `GENE_IFNG`,
  r = 0.332), spatial (`GC_density_LN`, r = 0.293), flow
  (`pct_plasmablast`, r = 0.280), and TCR/BCR
  (`vaccine_specific_clonotype_freq`, r = 0.272) - consistent with the
  known multi-compartment biology of a germinal-center vaccine response
  (T follicular helper activity, germinal center formation, plasmablast
  differentiation, and antigen-specific clonal expansion).
- **17/17 pytest tests pass**, including 10 that exercise the Flask API via
  a real test client (study/sample CRUD-style lookups, per-modality
  retrieval for all 4 modalities, biomarker ranking/lookup, 404 handling,
  and the cross-modality join endpoint).

Full ranked list and stats: `output/biomarker_discovery_results.csv` (all
90 features) and `output/discovered_biomarkers.json` (significant subset +
validation summary).

## Project structure

```
simulate_multiomics_data.py   # 55 studies x 4 modalities x ground-truth-driven response -> data/omics.db
discover_biomarkers.py        # Spearman + BH-FDR discovery, validated vs. ground truth -> output/
api/app.py                    # Flask REST API over the consolidated SQLite store
run_all.py                    # orchestrates simulate -> discover -> API smoke test
tests/test_pipeline.py        # pytest: data consolidation, biomarker validation, full API test-client coverage
data/ground_truth_biomarkers.json   # injected ground truth (used only for validation)
output/                       # discovery results + validation summary
```

## Tech stack

- Python 3
- pandas / numpy / scipy.stats for simulation and Spearman correlation /
  Benjamini-Hochberg FDR
- Flask for the researcher-facing REST API
- SQLite as the unified consolidated data store
- pytest + Flask test client for end-to-end testing (no live server
  required for CI)
