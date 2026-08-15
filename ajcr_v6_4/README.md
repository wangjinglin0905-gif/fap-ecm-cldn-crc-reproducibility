# AJCR v6.4 frozen analysis bundle

Frozen: 2026-08-15  
Repository release: v1.4.0

This directory contains the scripts, derived result tables, quality-control records and independent audits supporting the AJCR manuscript:

> Fibroblast-enriched senescence-associated transcription and a distinct FAP-linked matrix state in colorectal cancer

The release supersedes the senescence interpretation in v1.3.0 without deleting the earlier provenance. Raw public datasets are intentionally not redistributed.

## Scientific claim boundary

- The primary GSE132465 analysis is restricted to 1,501 source-annotated tumour fibroblast-lineage cells and 17,469 tumour epithelial cells from 23 patients.
- SenMayo119 and SASP25 scores are higher in the tumour fibroblast compartment than in the paired epithelial compartment. These are lineage-enriched senescence-associated transcriptional features; they do not diagnose senescent cells.
- Within tumour fibroblasts, the FAP-detection coefficients for SenMayo119 and SASP25 become null after log-UMI depth, fibroblast-subtype, MKI67 and patient adjustment. The release therefore does not support a FAP-specific senescence claim.
- The fibroblast-versus-epithelial SenMayo direction recurs in the independent GSE166555 cohort, but effect sizes are not treated as directly interchangeable because the annotation and pseudobulk procedures differ.
- FAP-matrix coupling is reproducible. FAP-SenMayo and SenMayo-matrix coupling are not reproducible at patient level.
- No multiplicity-supported cognate SASP-immune receptor association was found. Continuous matrix and senescence scores were not independently prognostic in the two bulk cohorts.

## Frozen numerical anchors

- GSE132465 SenMayo119: median paired standardized difference 0.466, 95% bootstrap CI 0.410-0.553, exact P = 2.38e-7; 23/23 patients in the same direction.
- GSE132465 SASP25: mean paired difference 0.236, 95% bootstrap CI 0.192-0.277, exact P = 4.77e-7.
- GSE132465 MKI67: mean paired difference -0.104, 95% bootstrap CI -0.142 to -0.063, exact P = 1.81e-4.
- Adjusted tumour-fibroblast FAP-detection coefficient: SenMayo119 beta = 0.0238, P = 0.106; SASP25 beta = 0.0323, P = 0.400.
- GSE166555: SenMayo119 median paired difference 0.149, exact P = 0.0195; 9/10 eligible patients in the same direction.
- CPTAC: FAP protein versus the FAP-excluded COL1A1/COL1A2/FN1 score, rho = 0.812, 95% bootstrap CI 0.706-0.879, BH FDR = 1.33e-23, n = 97.
- Cox boundary: every continuous-score Cox model in `results/frozen_boundaries/continuous_score_cox_prognosis.csv` has P > 0.10.

## Gene-set and circularity audit

- The source SenMayo list contains 125 genes.
- Four genes overlap FAP13 and are removed before SenMayo scoring: `CXCL12`, `MMP2`, `MMP9` and `WNT2`.
- Two non-overlapping candidates are absent from GSE132465 (`BEX3`, `CCL3L1`), leaving the frozen 119-gene score.
- `SASP25 intersect FAP13 = MMP9`. MMP9 remains in SASP25 because SASP25 is a separate, prespecified descriptive/sensitivity signature. It is not used as the overlap-removed SenMayo specificity test, and no SASP-ECM causal claim is made.

## P-value and interval conventions

- Patient-paired GSE132465 and GSE166555 contrasts use two-sided exact Wilcoxon signed-rank P values.
- Reported paired-effect intervals use 10,000 patient-level percentile bootstrap resamples. SASP25 and MKI67 seeds are 20260818 and 20260819, respectively.
- Expression-matched null P values use `(1 + exceedances) / (5,000 + 1)`; the minimum attainable value is 2.00e-4.
- Mixed-model P values use `lmerTest` Satterthwaite degrees of freedom. Marker-detection models are logistic mixed models.
- Spearman correlation P values use the asymptotic `exact = FALSE` calculation. Correlation intervals use patient/sample bootstrap resampling where reported.
- CPTAC FDR is Benjamini-Hochberg correction across the two prespecified protein-score comparisons in the frozen summary.
- Cox P values are model-based Wald P values from the frozen continuous-score analyses.

## Directory map

```text
config/                         frozen SenMayo source list
figures/                        frozen PNG, TIFF and PDF manuscript figures
ledgers/                        input, claim and reference metadata
qa/                             manuscript, reproduction and 600-dpi figure QA records
reports/                        independent editor/domain/bioinformatics audits
results/bulk/                   composition-sensitivity and mixed-model anchors
results/GSE132465/              tumour-only source tables and model outputs
results/GSE166555/              external pseudobulk source tables and ledger
results/frozen_boundaries/      CPTAC and Cox boundary tables
scripts/                        analysis, figure and verification scripts
```

## Main execution order

Run R explicitly with R 4.6.1 and the R 4.6 user library. The commands below assume execution from the repository root.
Ensure the packages recorded in the session-information files are installed on the active `.libPaths()` before using `--vanilla`.

```text
"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/recompute_single_cell_senescence.R \
  <GSE132465_seurat.rds> ajcr_v6_4/config/senmayo_genes.txt <gse132465_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/domain_tumor_only_check.R <gse132465_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/domain_marker_check.R \
  <GSE132465_seurat.rds> <gse132465_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/test_signature_specificity.R \
  <GSE132465_seurat.rds> ajcr_v6_4/config/senmayo_genes.txt <gse132465_output>

python ajcr_v6_4/scripts/gse166555_senescence_validation.py \
  --input-dir <GSE166555_directory> \
  --gene-availability <gse132465_output>/gene_availability.csv \
  --output-dir <gse166555_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/gse166555_external_stats_and_plot.R \
  <gse166555_output>/patient_compartment_scores.csv \
  <gse166555_output>/analysis_ledger.json \
  <gse166555_stats_output> <figure_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/compute_minor_revision_ci.R \
  <gse132465_output> <ci_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/make_revised_figures.R \
  <gse132465_output> <figure_output>

"C:/Program Files/R/R-4.6.1/bin/Rscript.exe" --vanilla \
  ajcr_v6_4/scripts/recompute_cptac_frozen_summary.R \
  <cptac_protein_sdc4_cd44.csv> <boundary_output>
```

Every multi-panel figure produced by this release uses lowercase panel labels (`a`, `b`, `c`, ...). PNG and LZW-compressed TIFF exports are 600 dpi RGB.

## Input integrity and redistribution

See `ledgers/input_manifest.json` for accession-level source records and hashes. The 559 MB GSE132465 Seurat object, GSE166555 raw count archive and metadata, and the CPTAC source table are not included. Small, non-identifying, derived tables required to audit reported values are included.

## Integrity verification

`CHECKSUMS_SHA256.csv` contains the SHA-256 digest and byte size of every file in this directory except the checksum file itself. Recompute hashes after checkout and compare them before rerunning analyses.
The directory-level `.gitattributes` disables text conversion so the frozen bytes and hashes remain stable across Git clients and operating systems.

```text
python ajcr_v6_4/scripts/validate_release.py ajcr_v6_4
```

The frozen plotting code intentionally retains the ggplot2 calls used for the submission rasters. Under ggplot2 4.0 these calls emit deprecation warnings, but the regenerated PNG and TIFF files were byte-identical to the frozen 600-dpi files; see `qa/reproduction_validation_v1.4.0.json`.
