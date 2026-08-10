# FAP-associated stromal matrix state in colorectal cancer

Repository-ready reproducibility archive for the Scientific Reports R1 manuscript:

> A reproducible FAP-associated stromal matrix state with senescence-related transcriptional covariation in colorectal cancer

Release date: 2026-08-10.

## What is in this release

This package preserves the full code lineage assembled for the manuscript and adds the audited R1 analyses. It contains 58 R analysis/provenance scripts, 14 Python acquisition/figure/manuscript utilities, one candidate R template and one candidate Python script retained only as audit evidence, derived result tables, figure source data, sample metadata, environment records and checksums.

The inferential R1 additions are:

- `scripts/locked_release/R1_audited/12_fibroblast_adjustment_audited.R`
  - selects MCP-counter `Fibroblasts` and EPIC `CAFs` by exact row name;
  - calculates standard partial Spearman coefficients as Pearson correlations of rank residuals;
  - uses 5,000 row bootstraps with seed 20260810;
  - creates the 380-sample primary manifest and the values in Supplementary Table S5.
- `scripts/locked_release/R1_audited/13_satterthwaite_mixed_models_audited.R`
  - normalises raw UMI counts with each cell's full-library total;
  - fits patient-random-intercept models with `lmerTest`;
  - reports Satterthwaite degrees of freedom;
  - retains the historical seven-target-gene-denominator calculation only as an audit comparison.
- `scripts/locked_release/R1_audited/14_verify_references_crossref_pubmed.R` and `15_finalize_reference_verification.R`
  - verify the 44 references against Crossref and PubMed;
  - record DOI, PMID, title/year concordance and retraction flags.

Machine-readable outputs are in `results/R1_audited/`.

## Scientific status of the R1 manuscript

- Primary TCGA score analysis: 380 unique primary tumours (COAD 286; READ 94).
- Stage-stratified TCGA subset: 89 tumours.
- GSE39582 replication: 566 tumours after excluding 19 non-tumour mucosa arrays.
- FAP13–matrix4: rho = 0.930 in TCGA and 0.913 in GSE39582.
- Fibroblast-abundance adjustment attenuates the TCGA association to rho = 0.113 with MCP-counter and rho = 0.288 with EPIC; the residual estimate is proxy-dependent.
- Full-library-normalised mixed models support matrix4–FAP coupling (1,501 fibroblast-lineage cells; Satterthwaite P = 9.17e-33) but not receptor2–FAP coupling (P = 0.302).
- The overlap-removed SenMayo–FAP13 relationship persists across fibroblast-abundance specifications, whereas the SenMayo–matrix4 relationship is largely attenuated.
- The results do not establish senescent CAFs, SASP causality, prognosis or a clinical decision tool.

## Repository layout

```text
scripts/
  locked_release/          frozen core scripts and R1_audited additions
  review_corrected/        corrected tumour-only senescence analysis
  legacy_exploratory/      historical provenance; not the inferential pipeline
  python_utils/            acquisition and metadata utilities
  manuscript_assembly/     manuscript and figure build utilities
config/
  gene_sets/               locked analysis gene lists
data/                       intentionally empty; public raw data are downloaded locally
derived_inputs/             permitted derived deconvolution inputs
derived_results/            corrected historical output tables
results/R1_audited/         audited R1 outputs, manifests and provenance
figure_source_data/         machine-readable figure source values
audit/                      candidate scripts retained only for comparison
docs/                       data-retrieval and lineage notes
checksums/                  release checksum records
```

## R1 execution order

1. Install R 4.6.1 and the packages in `R_PACKAGES.txt`.
2. Follow `docs/DATA_DOWNLOAD.md`; keep public raw data outside Git unless redistribution is explicitly permitted.
3. Run the established analysis scripts in `scripts/locked_release/` as required for the original figures and tables.
4. Run the R1 fibroblast sensitivity analysis:

   ```text
   Rscript scripts/locked_release/R1_audited/12_fibroblast_adjustment_audited.R INPUT_DIR OUTPUT_DIR
   ```

   `INPUT_DIR` must contain `TCGA.COAD.HiSeqV2.gz`, `TCGA.READ.HiSeqV2.gz`, `mcpcounter_tumor_scores.csv`, `epic_tumor_cell_fractions.csv` and `senmayo_genes.txt`.

5. Run the R1 mixed models:

   ```text
   Rscript scripts/locked_release/R1_audited/13_satterthwaite_mixed_models_audited.R R_LIBRARY CDS_RDS OUTPUT_DIR
   ```

6. Compare rerun outputs with `results/R1_audited/` and verify the source hashes recorded in `analysis_provenance.txt` and `satterthwaite_provenance.txt`.

## Data redistribution policy

Raw TCGA, GEO, cBioPortal, CPTAC and spatial-transcriptomic datasets are not bundled. Public availability does not automatically grant unrestricted redistribution, and GitHub is not an authoritative data repository. The package instead supplies accessions, URLs, access dates, retrieval logic, inclusion maps and source hashes. Small derived, non-restricted tables needed to audit the reported values are included.

See `DATA_MANIFEST.csv`, `DATA_USE_NOTICE.md` and `docs/DATA_DOWNLOAD.md`.

## Lineage warnings

- The cBioPortal TCGA PanCancer Atlas entry (n = 592) is an overlapping processing-route sensitivity analysis, not an independent cohort.
- Files under `scripts/legacy_exploratory/` are preserved for provenance and are not automatically valid for inference.
- The scripts in `audit/candidate_scripts_not_used/` did not generate the R1 manuscript results.
- The historical target-gene-denominator mixed model is retained only to explain earlier values; the manuscript uses full-library normalisation.
- MCP-counter and EPIC outputs must be selected by name, never by row position.

## Citation and licence

**GitHub release (v1.3.0):** https://github.com/wangjinglin0905-gif/fap-ecm-cldn-crc-reproducibility/releases/tag/v1.3.0
**Zenodo DOI:** https://doi.org/10.5281/zenodo.21873369

Please cite the associated manuscript, this repository release, and the source datasets used. Code and derived-data terms are described in `LICENSE`, `LICENSE-DATA.md` and `DATA_USE_NOTICE.md`; source-dataset terms continue to apply.

