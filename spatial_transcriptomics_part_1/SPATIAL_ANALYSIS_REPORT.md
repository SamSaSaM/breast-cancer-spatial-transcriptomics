# Spatial Transcriptomics Reanalysis of Human Breast Cancer

**Dataset:** Wu et al. (2021), *A single-cell and spatially resolved atlas of human breast cancers*  
**Pilot sample:** CID44971  
**Analysis type:** Exploratory reanalysis of publicly available Visium spatial transcriptomics data  
**Status:** Single-sample exploratory analysis completed; multi-sample validation and cell-type deconvolution pending

## 1. Research questions

1. Does a breast cancer tissue section exhibit spatial transcriptional heterogeneity?
2. Which tissue regions show epithelial-, stromal-, and immune-associated expression programs?
3. What transcriptional programs characterize the six Leiden clusters?
4. Where are macrophage-associated signals enriched, and how do they relate to epithelial-associated regions?
5. Are these patterns reproducible across the other five spatial samples, and do they agree with pathology annotations?

The current pilot addresses questions 1–3 descriptively and provides exploratory evidence relevant to question 4. Question 5, and quantitative macrophage–epithelial spatial relationships, remain open.

## 2. Data and workflow

The original spatial dataset contains six tissue sections: CID44971, CID4535, 1160920F, CID4290, 1142243F, and CID4465. The pilot focuses on **CID44971**. Data were obtained from the [public Zenodo record](https://zenodo.org/records/4739739), including `filtered_count_matrices`, `spatial`, and `metadata` archives.

The analysis workflow includes:

- Loading sparse 10x Matrix Market counts and Visium spot coordinates.
- Overlaying tissue spots on the H&E image.
- Calculating and visualizing per-spot QC metrics: total counts, detected genes, and mitochondrial fraction.
- Plotting selected epithelial, stromal, and immune-associated genes.
- Performing expression-based dimensionality reduction and Leiden clustering.
- Visualizing exploratory myeloid/macrophage-related gene-set scores.
- Ranking marker genes for each Leiden cluster.

**Observed pilot output:** **1,160 in-tissue spots**, **six Leiden clusters**, and **120 reported marker rows** (top 20 genes per cluster) in `leiden_top20_markers.csv`.

The exact normalization, feature-selection, clustering, differential-expression, and gene-set-scoring parameters should be documented from the executed scripts before claiming methodological equivalence to the original paper.

## 3. Figure inventory and interpretation

| Series | Output files | Purpose / interpretation |
|---|---|---|
| S1 | `S1_CID44971_tissue_spots.png` | H&E image with Visium spot overlay; checks tissue–coordinate alignment. |
| S2 | `S2_CID44971_total_counts_spatial.png`, `S2_CID44971_n_genes_by_counts_spatial.png`, `S2_CID44971_pct_counts_mt_spatial.png`, `S2_CID44971_qc_distributions.png` | Describes heterogeneity in sequencing depth, gene detection, and mitochondrial fraction; does **not** establish that stringent QC filtering was performed. |
| S3 | `S3_CID44971_EPCAM.png`, `KRT19`, `COL1A1`, `DCN`, `PTPRC`, `CD3D`, `LYZ`, `C1QA` (same filename prefix) | Maps selected epithelial-, matrix-, leukocyte-, T-cell-, and myeloid-associated genes. |
| S4 | `S4_CID44971_spatial_leiden.png`, `S4_CID44971_umap_leiden.png` | Displays six expression-derived Leiden clusters in tissue space and UMAP space. |
| S6 | `S6_CID44971_macrophage_core.png`, `antigen_presentation`, `lipid_associated_myeloid`, `inflammatory_myeloid`, `immunoregulatory_myeloid` (same filename prefix) | Exploratory spatial gene-set scores; not macrophage counts or definitive cell-state assignments. |

*The 20 figures were produced in the local analysis workflow; this report does not embed them. To display them on GitHub, upload the PNG files and update their relative paths.*

### 3.1 Spatial heterogeneity

The Leiden analysis identifies six transcriptionally distinguishable spot clusters, with spatial structure in the tissue section. This supports **within-section spatial transcriptional heterogeneity**. It does not by itself establish which clusters represent malignant epithelium, stroma, or immune compartments, and some separation may reflect distinct tissue fragments, local cell composition, or technical differences.

### 3.2 Marker-expression patterns

EPCAM and KRT19 are epithelial-associated markers; COL1A1 and DCN are extracellular-matrix/stromal-associated; PTPRC and CD3D are immune/T-cell-associated; LYZ and C1QA provide myeloid/macrophage-related clues. Spatial expression of these genes suggests heterogeneous transcriptional programs. **A Visium spot generally contains multiple cells**, so a single marker cannot identify all constituent cells or establish malignancy.

### 3.3 Exploratory macrophage-associated patterns

The previously inspected spatial maps suggested relatively elevated C1QA and selected macrophage-related signatures in part of the lower tissue region. This is **a qualitative observation requiring quantitative confirmation**. Signature scores are not cell abundances; antigen-presentation genes are also expressed by multiple immune and nonimmune cell types. Cross-signature color scales should not be compared as absolute expression measurements.

## 4. Leiden cluster annotation from marker-gene ranking

The annotations below are **provisional transcriptional-region labels**, not verified cell types or pathology categories. Marker lists were checked against the supplied `leiden_top20_markers.csv`.

| Cluster | Selected ranked markers | Provisional interpretation | Key caveat |
|---|---|---|---|
| **0** | MTRNR2L12, MALAT1, KRT15, MT-ND3, RPS6 | KRT15-associated epithelial signal with strong mitochondrial/ribosomal contributions | Potential technical or cellular-stress influence; insufficient for a definitive epithelial subtype |
| **1** | IGFBP7, IGKC, IGHG4, IGLC2, CCL19, VWF, PLVAP | Mixed vascular, immunoglobulin, and immune-associated region | Multiple transcriptional sources likely |
| **2** | XBP1, TFF3, AZGP1, TFF1, AGR2 | Secretory/luminal epithelial-associated program | Cannot infer malignant status from markers alone |
| **3** | EEF1A1, TMSB10, STMN1, CENPW, FRZB, CRABP1 | Possible proliferation/cell-cycle-related program with mixed features | Proliferation is tentative; requires cell-cycle and broader DEG validation |
| **4** | PHLDA2, KRT19, CCND1, KRT5, KRT16, CLDN3 | Epithelial-enriched program with mixed keratin expression | Not proof of a basal-like tumor subtype |
| **5** | CD74, HLA-DRA, FN1, C1QA, COL1A1, HLA-DPA1 | Antigen-presentation/myeloid- and stromal-associated mixed region | Not a pure macrophage cluster; may include several cell types |

### 4.1 Cluster 5: mixed immune/stromal program

Selected differential-expression statistics from the uploaded marker table:

| Gene | Reported log fold change | Adjusted p-value | Interpretation |
|---|---:|---:|---|
| C1QA | 2.469 | 2.59 × 10^-63 | Complement/macrophage-associated expression |
| HLA-DRA | 2.107 | 5.44 × 10^-74 | MHC class II / antigen-presentation-associated expression |
| FN1 | 3.042 | 2.09 × 10^-73 | Extracellular-matrix-associated expression |
| COL1A1 | 2.385 | 1.37 × 10^-56 | Collagen/stromal-associated expression |

**Interpretation:** Cluster 5 is compatible with a **mixed myeloid/antigen-presentation and stromal transcriptional region**. It is not valid to label it a macrophage-only region or to infer immunosuppressive function without additional evidence.

### 4.2 Cluster 2: secretory/luminal epithelial-associated program

| Gene | Reported log fold change | Adjusted p-value |
|---|---:|---:|
| TFF3 | 4.513 | 2.04 × 10^-76 |
| TFF1 | 4.535 | 5.47 × 10^-72 |
| AGR2 | 3.953 | 7.03 × 10^-67 |
| AZGP1 | 3.673 | 1.29 × 10^-75 |

These markers support a **secretory/luminal epithelial-associated expression program**, distinct from the keratin-enriched program of Cluster 4. They do not independently establish malignancy or receptor subtype.

**Statistical caution:** These adjusted p-values describe the current spot-level cluster comparison, not independent validation across patients. Spatial autocorrelation and data-driven clustering can make naive differential-expression significance overly optimistic. Reported fold changes retain the pipeline's own definition and should not be interpreted as absolute fold differences without checking its implementation.

## 5. Answers and unresolved questions

| Question | Current evidence-based answer | Needed for stronger inference |
|---|---|---|
| Is the tissue transcriptionally heterogeneous? | **Yes, descriptively:** six expression-derived clusters with spatial structure. | Robustness to QC and clustering parameters; pathology context. |
| Which regions have epithelial/stromal/immune programs? | Marker maps and cluster marker rankings support heterogeneous epithelial, stromal, and immune-associated programs. | Pathology labels and multi-marker / reference-based annotation. |
| What are the six clusters? | Provisional labels are given in Section 4. | Broader DEG analysis, marker co-expression, H&E review, and reference integration. |
| Where are macrophages? | C1QA and macrophage-associated scores identify **candidate enriched spots**, not confirmed macrophage counts. | scRNA-seq-based deconvolution and validation with multiple markers. |
| Are macrophage-associated spots near tumor epithelial regions? | **Not established** from the supplied marker CSV and images alone. | Quantitative neighborhood/distance analysis using validated tumor-region definitions, with appropriate null models. |
| Are patterns reproducible across six samples? | **Not yet established:** only CID44971 results are verified here. | Run and evaluate all six sections separately; account for patient/sample effects. |
| Do clusters match pathology? | **Not yet established.** | Verify barcode alignment and actual pathology metadata fields, then compare annotations. |
| Does this reproduce a specific original-paper panel? | **Not established.** This is a reanalysis of public data. | Match exact samples, processing parameters, and original figure panels. |

## 6. Suggested next analyses

1. **Validate cluster annotations:** inspect additional marker genes and the full DEG table; review image morphology and any available pathology labels.
2. **Reference-based deconvolution:** integrate the Wu et al. scRNA-seq atlas to estimate cell-type mixtures per Visium spot, with appropriate reference preprocessing and uncertainty checks.
3. **Spatial proximity:** define epithelial/tumor-associated regions using validated annotations, then test distances or neighborhood enrichment relative to macrophage-associated spots; use within-section spatial permutations rather than treating spots as independent patients.
4. **Cross-sample replication:** repeat analyses on the five remaining sections and report both shared and sample-specific findings.
5. **Paper reproduction:** compare a selected original-paper spatial panel against a matched sample and author-provided analysis settings before claiming a figure has been reproduced.

## 7. Reproducibility and limitations

- **Platform:** 10x Genomics Visium; spots capture mixtures of cells.
- **Pilot sample size:** one tissue section (CID44971), 1,160 in-tissue spots; no patient-level replication in this report.
- **Clustering:** six Leiden clusters; expression-derived clustering does not necessarily incorporate spatial adjacency.
- **QC:** metrics and visualizations available; QC exclusion thresholds and filtering effects require explicit reporting.
- **Annotation:** based on selected marker genes and top-20 DEG rankings; not confirmed cell identities.
- **Macrophage scores:** exploratory gene-set measures, not deconvolution or functional assays.
- **Statistical inference:** spot-level p-values should not be generalized to patients or breast cancer subtypes.
- **Original study:** do not describe these visualizations as exact reproductions of published figure panels without direct figure/method matching.
- **Data availability:** source data are public; do not commit large raw count matrices, images with restrictive reuse terms, or private files without checking licenses and repository size constraints.

## 8. Files to include in the GitHub repository

Recommended layout (upload only files you have actually generated):

```text
spatial-analysis/
├── SPATIAL_ANALYSIS_REPORT.md
├── leiden_top20_markers.csv
├── scripts/
│   └── ... (executed analysis scripts)
└── figures/
    └── ... (S1–S6 PNG figures)
```

To include the images inline in this Markdown report, add relative links after uploading them, e.g.:

```markdown
![CID44971 spatial Leiden clusters](figures/S4_CID44971_spatial_leiden.png)
```

## 9. References

- Wu, S. Z. et al. (2021). *A single-cell and spatially resolved atlas of human breast cancers*. **Nature Genetics**. https://doi.org/10.1038/s41588-021-00911-1
- Public spatial dataset: https://zenodo.org/records/4739739
- Original analysis code: https://github.com/Swarbricklab-code/BrCa_cell_atlas

---

**Summary:** The CID44971 pilot demonstrates spatially structured transcriptional heterogeneity and identifies provisional epithelial, vascular/immune, and myeloid/stromal-associated programs. Cluster 5 is particularly notable for co-occurring C1QA, HLA-DRA, FN1, and COL1A1 markers. The data support **hypothesis generation**, not definitive macrophage localization, tumor–immune interaction, or cross-patient generalization.
