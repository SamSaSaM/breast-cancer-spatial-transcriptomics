# CID44971 Spatial Transcriptomics — Results Supplement (Part 4)

Add these files to your **existing spatial-transcriptomics repository**; do not create a new repository.

## Contents

- `CID44971_qc_summary.csv`: summary of UMI counts, detected genes and mitochondrial transcript percentage across 1,160 tissue spots.
- `CID44971_cluster_counts.csv`: Leiden cluster spot counts (0–5).
- `CID44971_signature_scores.csv`: per-spot expression signature scores; scores are not cell counts.
- `CID44971_signature_gene_coverage.csv`: genes included/missing in each signature.
- `CID44971_metadata_inspection.txt`: inspection of sample metadata and barcode matching.

## Upload

1. Unzip this archive locally.
2. In your existing GitHub repository, open `results/CID44971/` (create it if needed).
3. Upload the five data files from this folder and commit.
4. Optionally place this README in `results/CID44971/` as `RESULTS_README.md`.

## Notes

This archive deliberately excludes `.h5ad` files (large, and two exceed GitHub's browser upload limit), original Visium matrices, the virtual environment, and previously uploaded figures/scripts. No existing results have been altered.

These results represent exploratory analysis of **CID44971** only, not a six-sample validation. QC metrics do not by themselves establish that mitochondrial filtering was applied.
