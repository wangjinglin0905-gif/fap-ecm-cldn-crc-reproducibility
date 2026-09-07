# Release v2.0.2 — PLOS ONE v7.2 reporting corrections

Date: 2026-09-07

This release synchronizes the revised supplementary data with manuscript v7.2. It preserves all previous GitHub releases and Zenodo records.

- S4: replace five historical manuscript-version labels with an explicit description of the adjusted marker-detection model.
- S6: restore 16 small nonzero P values from the frozen bulk-result CSVs. The previous supplementary-table serialization showed zero; this was a CSV parsing/reporting error, not a statistical recomputation. Scientific notation now preserves the source values to 15 significant digits.
- Update the portable table builder to read P-value source strings before decimal formatting and to emit the current S4 descriptions.
- Add a row/cell-level correction ledger with repository-relative provenance and a 34-reference DOI ledger including three linked correction notices.
- Refresh both integrity manifests and release validation.

No statistical analyses were rerun. All source results, effect estimates, confidence intervals, sample counts, scientific conclusions, and 20 PNG/TIFF figure files remain unchanged from v2.0.1. The stable `plos_one_v7_1/` directory name is retained for script and data-path continuity; its v2.0.2 snapshot supports manuscript v7.2.

Raw public matrices, local manuscript Word files, private metadata-cleaning reports, machine paths and credentials are not part of this patch. The parent/concept DOI remains `10.5281/zenodo.21441731`; the release-specific DOI is assigned by Zenodo after GitHub publication.
