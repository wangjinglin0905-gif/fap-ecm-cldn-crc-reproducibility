# Changelog — v2.0.0

Release date: 2026-09-04

This major release replaces the AJCR v6.4 evidence architecture with the PLOS ONE v7.1 analysis package. Previous releases remain immutable and available through GitHub and Zenodo version history.

## Added

- Patient-paired localization of senescence-associated transcription in GSE132465 and independent GSE166555.
- Patient-cluster bootstrap inference for FAP-detected fibroblast analyses and direct contrasts between dependent correlations.
- Expression-, variability- and detection-matched null signatures.
- Target-purged and globally disjoint MCP-counter/EPIC composition-proxy sensitivity analyses.
- A fixed 111-gene SenMayo common-core sensitivity across two single-cell and three spatial cohorts without cross-platform pooling.
- A three-cohort spatial stress test retaining the mixed or negative GSE334323 result and its corrected provenance.
- Four main figures and six supporting figures as PNG and 600-dpi LZW TIFF with lowercase panel labels.
- Seven supporting tables, four companion source-data files, reference verification and machine-readable claim, numeric and checksum ledgers.

## Interpretation changes

- Senescence-associated transcription is localized preferentially to fibroblast compartments but is not presented as proof of durable cellular senescence.
- FAP transcript detection does not identify a demonstrated uniquely senescent fibroblast subgroup.
- FAP–matrix covariation is reported as measurement-context and composition-proxy sensitive, not as a distinct cell state.
- Spatial effects are treated as heterogeneous generalizability evidence rather than universal validation.
- Organoid, trogocytosis and mitoxyperilysis analyses remain outside the Results because the available datasets do not directly measure those mechanisms.

## Integrity

Raw public matrices are not redistributed. The package records accession identifiers, eligibility rules, consumed-artifact hashes, portable reconstruction scripts and a file-level SHA-256 manifest.
