# Release v2.0.2 reporting correction

Date: 2026-09-07

This snapshot supports PLOS ONE manuscript v7.2 while retaining the stable package directory name. Five S4 source-note cells now describe the marker-detection model without a historical manuscript-version label. Sixteen S6 P-value cells have been restored from the unchanged frozen result CSVs and are reported in scientific notation rather than zero. The portable table builder preserves small P values through decimal-string parsing. The cell-level patch ledger provides repository-relative sources, and the current DOI ledger covers 34 references including three linked correction notices.

No model was rerun. All frozen scientific results, effect estimates, confidence intervals, sample counts and figure images are unchanged. Both manifests have been refreshed for the new reporting files. Prior release notes below describe historical snapshots, not the current manuscript.

## Prior release v2.0.1 integrity correction

Date: 2026-09-07

This patch repairs the integrity metadata for the PLOS ONE v7.1 reproducibility package.

- Freeze LF line endings for package text files through `.gitattributes`; calculate SHA-256 against the actual Git archive bytes.
- Refresh the supplementary-table manifest after final formatting and exclude its own entry. The manifest now covers exactly seven tables and four companion files.
- Add a table-manifest validator and a portable manifest-refresh script.
- Correct the README verification order: validate the supplied manifest before rebuilding any checksum baseline.

In v2.0.0, the package-level checksum manifest was calculated from LF-normalized working files while some committed files retained different line endings. The supplementary-table manifest also included its own earlier version. These issues affected integrity checks, not scientific values. Scientific inputs, results, evidence ledgers and PNG/TIFF figure files are unchanged. Prior releases and their DOIs remain available for provenance.

The manuscript remains v7.1. A local supporting-document media synchronization aligns the embedded S6 preview with the already published deterministic S6 PNG; that DOCX is not part of this public code package.

## Prior release history

# Release notes — v2.0.0

Date: 2026-09-04

Version 2.0.0 is a major evidence-architecture update relative to the AJCR v6.4 archive in v1.4.0. Earlier releases remain immutable and available through GitHub and Zenodo version history.

## Added

- A patient-level single-cell localization hierarchy using GSE132465 and the independent GSE166555 cohort.
- Patient-cluster bootstrap inference for FAP-detected SenMayo/SASP25 models and the continuous FAP–matrix4 model.
- Direct bootstrap contrasts between dependent patient-level correlations.
- Expression-, variability- and detection-matched null analyses.
- Target-purged and globally disjoint MCP-counter/EPIC composition-proxy analyses.
- A fixed 111-gene SenMayo common-core sensitivity across two single-cell and three spatial cohorts, without cross-platform pooling.
- A three-cohort spatial stress test that retains the mixed/negative GSE334323 result and its host-transcript denominator correction provenance.
- Four main figures and six supporting figures in PNG and 600-dpi LZW TIFF, all with lowercase panel labels.
- Seven supporting tables, four companion source-data files and machine-readable claim, numeric, reference and checksum ledgers.

## Interpretation changes

- The current manuscript no longer describes a distinct FAP-linked matrix state as a demonstrated cell state.
- FAP transcript detection is not used as evidence of a uniquely senescent fibroblast subgroup.
- The bulk FAP13–matrix4 result is described as composition-proxy sensitive rather than composition independent or composition driven.
- Spatial data are presented as heterogeneous generalizability evidence rather than universal validation.
- Organoid, trogocytosis and mitoxyperilysis analyses are excluded from the Results because the available data do not directly measure those mechanisms.

## Integrity

The release contains no raw TCGA, GEO, CPTAC/cBioPortal or spatial matrices. Public identifiers and consumed-artifact hashes are retained instead. Run `python scripts/validate_release.py` after downloading the release to verify its structure and SHA-256 manifest.
