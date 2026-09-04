# PLOS ONE v7.1 reproducibility package

This directory contains the frozen analysis outputs, analysis and figure scripts, figure source data, supplementary tables, environment records and integrity checks supporting the manuscript:

> **Fibroblast-enriched senescence-associated transcription and measurement-context-dependent FAP–matrix covariation in colorectal cancer**

Release: **v2.0.0 (2026-09-04)**.

## Scientific scope

The package supports four bounded conclusions:

1. Senescence-associated transcription is higher in fibroblast-lineage than epithelial tumour compartments at the patient level in GSE132465 and GSE166555.
2. Patient-cluster analyses do not support additional SenMayo or SASP25 enrichment among FAP-detected GSE132465 tumour fibroblasts; this is not an equivalence result.
3. FAP and extracellular-matrix transcription covary within fibroblasts and across bulk and protein measurements, while the bulk effect size is sensitive to the composition-proxy definition.
4. Spatial generalizability is heterogeneous: the 111-gene common-core direction is positive in GSE280315 and Valdeolivas but mixed or negative in GSE334323.

Transcript scores do not establish durable cellular senescence, causal stromal signalling, trogocytosis, mitoxyperilysis or a distinct FAP-defined senescent cell state.

## Contents

```text
results/                 frozen single-cell, bulk, common-core and spatial-null outputs
derived_inputs/          small public-data-derived inputs needed to rebuild figures/tables
figures/main/            Fig 1–Fig 4 in PNG and 600-dpi LZW TIFF plus source data
figures/supporting/      S1–S6 in PNG and 600-dpi LZW TIFF plus source data
tables/                  S1–S7 tables and four companion source-data files
ledgers/                 numeric, claim, endpoint and reference-verification ledgers
scripts/                 analysis, table, figure, checksum and validation scripts
qa/                      frozen session information and release validation
CHECKSUMS_SHA256.csv     file-level release integrity manifest
```

All multi-panel figure labels are lowercase.

## Data sources and redistribution boundary

Raw public datasets are not redistributed. Repository identifiers, eligibility rules and hashes of the consumed artifacts are recorded in `tables/S1_Table_dataset_accessions_eligibility_and_checksums.csv`.

- TCGA COAD/READ: Genomic Data Commons and UCSC Xena
- GSE39582, GSE132465, GSE166555, GSE280315 and GSE334323: NCBI Gene Expression Omnibus
- Valdeolivas spatial cohort: Zenodo version DOI `10.5281/zenodo.7760264`
- CPTAC colon proteomics: cBioPortal

`scripts/download_gse132465_public_inputs.py` and `scripts/merge_verified_range_parts.py` document the resumable GSE132465 retrieval used for public-data reconstruction. Source-dataset licences and terms remain controlling.

## Runtime

The frozen run used R 4.6.1 on Windows 11. Principal packages were Matrix 1.7-5, lme4 2.0-6, lmerTest 3.2-1, ggplot2 4.0.3, patchwork 1.3.2, arrow 25.0.0, MCPcounter and EPIC. Exact session records are stored with their result modules and in `qa/R_sessionInfo_figures.txt`.

Python 3.11 or later is sufficient for the checksum and release validators. Rebuilding supplementary tables additionally requires pandas. The reference ledger was verified live against Crossref before this release.

## Verification

From this directory:

```text
python scripts/09_build_sha256_manifest.py
python scripts/validate_release.py
```

The expected release-validation result is zero failures. The checksum builder intentionally excludes `CHECKSUMS_SHA256.csv` and the generated `qa/release_validation.json` to avoid self-referential hashes.

## Rebuilding derived tables and figures

The frozen result files are already supplied. To recreate the claim/numeric locks and manuscript-associated tables from the included derived inputs:

```text
python scripts/08_build_rerun_lock_artifacts.py
python scripts/15_build_supplementary_tables.py . derived_inputs/v6_4_qa derived_inputs/route_b_results
```

To recreate Fig 1–Fig 4 and S1–S6:

```text
Rscript scripts/21_build_four_main_and_supplement_previews_v7_1.R results derived_inputs/v6_4_qa derived_inputs/route_b_results regenerated/main regenerated/supporting
```

The earlier-stage scripts `01`–`06` and `13a`–`13b` accept explicit input and output paths. They are retained so that the frozen analyses can be rerun after the corresponding public raw data are obtained. `06_rescore_spatial_common_core.R` requires the full Route-B working dataset because it reconstructs spatial scores from the downloaded matrices; the compact derived inputs in this release are sufficient for table and figure regeneration but do not redistribute the full public matrices.

## Provenance and interpretation controls

- `ledgers/numeric_lock_v01_2026-09-04.csv` links manuscript-level numerical claims to source outputs.
- `ledgers/claim_lock_v01_2026-09-04.csv` records the permitted interpretation ceiling and prohibited extensions.
- `ledgers/scrna_endpoint_estimand_lock_v01_2026-09-04.csv` fixes patient-level single-cell estimands.
- `ledgers/reference_verification_v7.1_2026-09-04.csv` records Crossref verification for all 31 DOI-bearing references.

See `RELEASE_NOTES.md` for the change from v1.4.0.
