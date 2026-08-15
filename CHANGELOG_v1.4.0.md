# v1.4.0 — AJCR v6.4 frozen analysis release

Release date: 2026-08-15

## Added

- `ajcr_v6_4/`, a self-contained release directory for the current AJCR analysis boundary.
- Tumour-only GSE132465 compartment analyses using 1,501 fibroblast-lineage and 17,469 epithelial cells from 23 patients.
- Independent GSE166555 tumour-compartment direction check.
- Expression-matched signature null analyses and FAP-overlap removal records.
- Patient-level bootstrap intervals added during the v6.4 minor revision.
- Frozen CPTAC protein and continuous-score Cox boundary tables.
- Publication figures in PNG, TIFF and PDF, with lowercase labels in every multi-panel figure.
- Independent bioinformatics, domain, editorial and reference audits.
- Input-source hashes, release checksums and environment records.

## Corrected interpretation

- The 3,462-fibroblast/18,539-epithelial whole-object comparison is not used as the primary tumour analysis; the current primary analysis is tumour-only.
- The overlap-removed GSE132465 SenMayo score contains 119 represented genes, not 120.
- Within tumour fibroblasts, FAP-detection associations with SenMayo119, SASP25 and individual marker detection are null after sequencing-depth and subtype adjustment. The release therefore does not claim FAP-specific senescent CAFs.
- The fibroblast-versus-epithelial result is described as senescence-associated transcription, not as a diagnostic demonstration of cellular senescence.
- The FAP-linked matrix state is treated as distinct from the compartment-level SenMayo pattern.
- CPTAC supports FAP–ECM protein covariation, not SASP–ECM coupling.
- The cognate SASP–immune screen and continuous-score Cox analyses are retained as negative boundaries.

## Reproducibility notes

- Primary R scripts were parsed and rerun with R 4.6.1 and the R 4.6 user library.
- Exact, bootstrap, empirical-null, mixed-model, correlation and Cox P-value conventions are stated separately in `ajcr_v6_4/README.md`.
- Public raw datasets are not redistributed. Accessions, expected hashes and retrieval information are provided instead.
