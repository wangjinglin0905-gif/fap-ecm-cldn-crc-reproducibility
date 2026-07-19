# FAP–ECM–CLDN colorectal cancer reproducibility package

Version: v1.1.0 candidate, 2026-07-19, R 4.6.1 audit release

## v1.1.0 revision scope

This candidate release supersedes v1.0.0 for the audited manuscript. It adds:

- member-specific TCGA correlations for CLDN1, CLDN2, CLDN4 and the prespecified CLDN core;
- a 380-primary-tumour sensitivity analysis;
- a 32-pair matched TCGA primary-tumour versus adjacent-normal sensitivity analysis;
- corrected sparse-matrix dimension-name handling and barcode identity checks in both conventional CellChat scripts;
- regenerated Figures 1 and 6 and updated figure source data;
- a consolidated validation manifest that distinguishes the historical full run from the 2026-07-19 reruns.

This package supports the manuscript:

> Multiscale analyses distinguish a FAP-associated extracellular matrix programme from epithelial claudin expression in colorectal cancer

## Scientific scope

The reproduced analyses support a FAP-associated stromal and extracellular-matrix programme in colorectal cancer. Single-cell and spatial CellChat analyses identify expression-compatible collagen/FN1–SDC4/CD44 interactions. The analysed public data do not establish direct FAP–CLDN regulation, receptor activation, downstream claudin remodelling or a causal therapeutic mechanism.

No simulated dataset is used as evidence. UALCAN and the direct CPTAC analysis use the same CPTAC source and are not counted as independent proteomic cohorts.

The internal result-column name `FAP_CAF_de_ligand` is retained for traceability. In the manuscript and figures it is called the **reduced FAP–CAF score**. The score removes only `TGFB1`, `INHBA`, `WNT2` and `WNT5A` from the original FAP–CAF gene set.

## Package contents

- `scripts/`: R-first analysis scripts and the Python utilities used for public-data download or memory-efficient extraction.
- `results/`: machine-readable analysis outputs used in the manuscript.
- `cellchat_results/`: conventional CellChat and SMC20-excluded sensitivity outputs.
- `figures/`: six main and five supplementary figures in PNG and TIFF, plus source data and the figure manifest.
- `manifests/`: public-data source, checksum, download and run manifests.
- `logs/`: the completed R 4.6.1 stage logs and session information.
- `audits/`: evidence ledger, result-comparison, reference-verification and document-structure audit outputs.
- `FILE_MANIFEST.csv`: SHA-256, size and relative path for every packaged file.

Large public raw files are retained locally but are not redistributed in this ZIP. Their filenames, sizes and checksums are listed in `manifests/raw_input_manifest_r461.csv` (408 files; 5.678 GiB in the retained audit set).

## Required R environment

All statistical inference, omics analyses, spatial analyses, CellChat workflows and figure generation in the final audit were run with:

```text
R version 4.6.1 (2026-06-24 ucrt)
```

The default Windows executable expected by the portable runner is:

```text
C:/Program Files/R/R-4.6.1/bin/Rscript.exe
```

Set `FAP_R461` if R 4.6.1 is installed elsewhere. Set `FAP_R_LIBRARY` or `R_LIBS_USER` to a compatible package library when packages are not installed in R's default library.

Key package versions in the completed audit were:

- Seurat 5.5.1
- Harmony 2.0.5
- CellChat 2.2.0.9001
- NicheNet 2.2.1.1
- edgeR 4.10.1
- PROGENy 1.34.0
- DoRothEA 1.23.0
- MCPcounter 1.2.0
- EPIC 1.1.7

Exact package and platform details are recorded in `logs/r461_final_sessionInfo.txt` and the analysis-specific session files.

Python did not perform inferential statistics or draw manuscript figures. It was used only for segmented downloads, extraction of selected genes from dense matrices, and document or archive processing. The Python scripts are retained.

## External inputs

The source repositories and accessions are listed in `manifests/public_data_sources_manifest.csv`. Place downloaded files under `work/reproducibility/inputs/`, or set the relevant environment variables:

| Variable | Input |
|---|---|
| `FAP_TCGA_EXPRESSION` | TCGA-COAD/READ expression matrix |
| `FAP_GSE132465_MATRIX` | GSE132465 raw UMI matrix |
| `FAP_GSE132465_METADATA` | GSE132465 cell annotation file |
| `FAP_QI_MATRIX` | Qi et al. single-cell count object |
| `FAP_QI_METADATA` | Qi et al. cell metadata workbook |
| `FAP_NICHENET_RESOURCES` | NicheNet ligand-target and ligand-receptor resource directory |
| `FAP_LEGACY_CMS` | Retained CMS labels used for sensitivity analysis |
| `FAP_LEGACY_MSI` | Retained MSI labels used for cross-checking |
| `FAP_LEGACY_HABITAT` | Legacy Qi habitat assignments used only for agreement checks |
| `FAP_LEGACY_RPPA` | Optional legacy RPPA table used only for comparison |
| `FAP_R_LIBRARY` | R package library containing required packages |

The Valdeolivas spatial download manifest records the source MD5 checksums. The audit-level raw-input manifest records SHA-256 checksums for retained inputs.

## Preparation steps that use Python

Run these only when the corresponding prepared inputs are absent:

1. `scripts/02_singlecell_l1.py`
2. `scripts/07b_extract_gse_harmony_matrix.py`
3. `scripts/11_download_valdeolivas_spatial.py`
4. `scripts/12_gse166555_validation.py`
5. `scripts/download_segmented.py` as called by the download scripts

These utilities require public raw files or network access. They do not replace the R analysis stages.

## R 4.6.1 run order

From the project root, run:

```powershell
$env:FAP_R_LIBRARY='D:/path/to/R/library'
& 'work/reproducibility/run_r461_pipeline.ps1'
```

The runner verifies that the configured executable is R 4.6.1 and executes 19 R stages in this order:

1. TCGA bulk analysis
2. Qi spatial analysis
3. Conventional single-cell CellChat
4. SMC20-excluded CellChat sensitivity analysis
5. NicheNet
6. Immune deconvolution
7. Qi variable-gene selection
8. Harmony integration
9. RPPA validation
10. PROGENy and DoRothEA
11. CPTAC protein validation
12. GSE166555 patient-level statistics
13. Valdeolivas spatial validation
14. Spatial CellChat
15. Spatial CellChat post-processing
16. UALCAN CPTAC parsing
17. Primary-tumour-only TCGA sensitivity analysis
18. Individual-claudin and matched primary-tumour/normal sensitivity analysis
19. Final figure regeneration

The historical full-run manifest is retained in `manifests/run_manifest.csv`. The consolidated `manifests/current_validation_manifest.csv` records the historical run, the 2026-07-19 reruns, the added stage 18 sensitivity analysis and the corrected CellChat reruns. Reference verification is kept separate because it requires live Crossref or NCBI access.

`scripts/06_oncopredict_l4.R` is retained for audit completeness. Its output is excluded from the revised manuscript because it did not provide externally validated clinical prediction.

## CellChat result correction

The previous draft's value of 85 interactions was incorrect. The R 4.6.1 rerun produced:

- Full pooled analysis: 80 significant ligand–receptor pairs, including 29 COLLAGEN-pathway and 3 FN1-pathway pairs.
- SMC20-excluded sensitivity analysis: 93 significant pairs, including 35 COLLAGEN-pathway and 4 FN1-pathway pairs.

The conventional analysis uses `type = "triMean"`. The spatial CellChat analysis uses a 10% truncated mean. These outputs represent expression-compatible communication probabilities and do not prove receptor activation or causality.

Both conventional CellChat scripts now preserve matrix dimension names after sparse-matrix normalisation and assert exact barcode identity between the expression matrix and metadata before creating the CellChat object. The corrected reruns retained the same reported pair counts: 80 in the pooled analysis and 93 after excluding SMC20.

## Member-specific claudin and matched-pair sensitivity

`scripts/18_tcga_individual_claudin_and_matched_sensitivity.R` prevents opposing claudin directions from being hidden by the composite score.

- Across all 434 profiles, reduced FAP-CAF score correlations were positive for CLDN1 (rho=0.155, FDR=0.00255), inverse for CLDN4 (rho=-0.154, FDR=0.00255), and not significant for CLDN2 or the CLDN core.
- In 380 primary tumours, correlations were inverse for CLDN2, CLDN4 and the CLDN core; the CLDN1 estimate attenuated to rho=0.100 (FDR=0.0510).
- In 32 matched primary-tumour/adjacent-normal pairs, FAP, CLDN1 and CLDN2 were higher in tumours, whereas CLDN4 was not significantly different after correction.

These analyses support weak, member-specific associations rather than a consistent family-wide positive FAP-CLDN relationship.

## Analysis units and multiplicity

- TCGA analyses use tissue profiles; the sensitivity analysis restricts tumours to 380 primary profiles.
- GSE166555 inference uses patient pseudobulk rather than individual cells.
- Valdeolivas spatial metrics are calculated by section and averaged within patient before patient-level tests.
- Spatial CellChat results are corrected across six prespecified pairs within each section and summarised by patient-level replicate reproducibility.
- Benjamini–Hochberg correction is applied within analysis families defined by one biological question in one dataset.

## Figure regeneration

`scripts/20_regenerate_all_figures.R` writes six main and five supplementary figures as 300 dpi PNG and 600 dpi LZW TIFF files. It also copies the source tables used by the plots and writes a figure manifest and R session information.

All final drawing, export and visual correction were performed in R 4.6.1.

## Licence

Code and scripts are licensed under the **MIT License** (see `LICENSE`). Generated result tables and figures are licensed under **Creative Commons Attribution 4.0 International (CC-BY 4.0)**. Public datasets remain subject to their source repositories' terms; this package does not grant new rights over third-party data.

## Public release and unresolved fields

Before GitHub and Zenodo deposition:

1. ~~Select and add an explicit software licence.~~ **Done** — MIT (code) + CC-BY 4.0 (data and figures).
2. Create a versioned GitHub release and record the commit and release tag.
3. Archive the same release in Zenodo and obtain the version DOI.
4. Replace the repository and DOI placeholder in the manuscript.
5. Test the public ZIP in a clean directory and confirm the checksums.

The public datasets remain subject to their source repositories' terms. This package does not grant new rights over third-party data.
