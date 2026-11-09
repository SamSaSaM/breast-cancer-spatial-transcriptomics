# Wu et al. (2021) — Breast Cancer Spatial Transcriptomics Reanalysis

An exploratory, reproducible analysis of publicly released 10x Visium data from [Wu et al. (2021)](https://www.nature.com/articles/s41588-021-00911-1). **CID44971** is the completed pilot; five other samples are available for additional analysis. This is **not** a validated reproduction of the original paper’s figures.

## Contents

- [`SPATIAL_ANALYSIS_REPORT.md`](SPATIAL_ANALYSIS_REPORT.md) — scientific questions, results, interpretations and caveats.
- [`figures/`](figures/) — 20 actual pilot PNG figures (S1, S2, S3, S4 and S6).
- [`results/CID44971/leiden_top20_markers.csv`](results/CID44971/leiden_top20_markers.csv) — 120 top marker rows across six clusters.
- [`scripts/pilot/`](scripts/pilot/) — original S1–S6 exploratory plotting workflow, with gzip magic-byte handling.
- [`scripts/advanced/run_analysis.py`](scripts/advanced/run_analysis.py) — advanced workflow for cluster markers, gene signatures, exploratory neighborhoods, metadata checks, and six-sample descriptive comparison. Includes fixes for single-column features and pandas index names.

## Data (not included)

Download `filtered_count_matrices.tar.gz`, `spatial.tar.gz` and `metadata.tar.gz` from [Zenodo](https://zenodo.org/records/4739739), then extract to:

```text
~/Downloads/wu2021_spatial/data/
  filtered_count_matrices/CID44971_filtered_count_matrix/
  spatial/CID44971_spatial/
  metadata/CID44971_metadata.csv
```

## Run on macOS

```bash
cd ~/Downloads/wu2021_spatial_github_complete
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/advanced/run_analysis.py --sample CID44971
# Once the pilot succeeds:
python scripts/advanced/run_analysis.py --sample all
```

By default, advanced results are written to `~/Downloads/wu2021_spatial/advanced_results/`. To run the pilot S1–S6 pipeline instead, use `cd scripts/pilot && python run_all.py` after installing its requirements. **Both pipelines are included as distinct workflows; they may produce different Leiden labels and are not interchangeable without checking preprocessing settings.**

## Interpretation and reproducibility

- 10x Visium spots contain multiple cells. Gene signatures are **not** cell-type proportions, and epithelial expression does not establish malignancy.
- Leiden groups are expression-defined tissue regions; marker annotations are provisional.
- Spot-wise differential-expression p-values are exploratory because spots are spatially correlated.
- The advanced pipeline does **not** perform formal scRNA-seq deconvolution, validate macrophage–tumor contact, or establish patient-level reproducibility.
- Outputs in this repository are from the supplied CID44971 pilot; other five samples have **not** been reported as successfully analyzed here.
- Spatial input archives, raw matrices and intermediate h5ad files are excluded due to size; obtain the original data from Zenodo.
