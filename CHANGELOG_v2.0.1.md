# Release v2.0.1 integrity correction

Date: 2026-09-07

This patch repairs the integrity metadata for the PLOS ONE v7.1 reproducibility package.

- Freeze LF line endings for package text files through `.gitattributes`; calculate SHA-256 against the actual Git archive bytes.
- Refresh the supplementary-table manifest after final formatting and exclude its own entry. The manifest now covers exactly seven tables and four companion files.
- Add a table-manifest validator and a portable manifest-refresh script.
- Correct the README verification order: validate the supplied manifest before rebuilding any checksum baseline.

In v2.0.0, the package-level checksum manifest was calculated from LF-normalized working files while some committed files retained different line endings. The supplementary-table manifest also included its own earlier version. These issues affected integrity checks, not scientific values. Scientific inputs, results, evidence ledgers and PNG/TIFF figure files are unchanged. Prior releases and their DOIs remain available for provenance.

The manuscript remains v7.1. A local supporting-document media synchronization aligns the embedded S6 preview with the already published deterministic S6 PNG; that DOCX is not part of this public code package.
