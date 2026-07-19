# Data availability

All source datasets analysed in this study are publicly available from TCGA/UCSC Xena, GEO (GSE132465 and GSE166555), cBioPortal/CPTAC, the supplementary files accompanying Qi et al., and the Valdeolivas et al. Zenodo record (concept DOI 10.5281/zenodo.7551712).

This v1.1.0 candidate package contains the analysis scripts, processed result tables, corrected CellChat outputs, run logs, session information, figure source data and checksum manifests used for the final audited manuscript. Large third-party raw files are not redistributed in the release ZIP. The retained local audit set contains 408 files totalling 5.678 GiB; filenames, sizes and SHA-256 values are recorded in `manifests/raw_input_manifest_r461.csv`, and repository URLs/accessions are recorded in `manifests/public_data_sources_manifest.csv`.

The permanent v1.1.0 GitHub release URL and Zenodo version DOI must be inserted in the manuscript after deposition and before submission.

## Repository actions still required

- Publish the GitHub v1.1.0 release from the final ZIP.
- Archive the same release in Zenodo and record the version DOI.
- Replace `[FINAL REVISION DOI AND RELEASE URL]` in the manuscript.
- Test the public download in a clean directory and verify `FILE_MANIFEST.csv`.

## FAIR audit status

- Findable: public accessions and URLs are recorded; the final v1.1.0 DOI remains pending.
- Accessible: source repository links and accessions are documented; third-party terms still apply.
- Interoperable: processed tables use CSV/TSV and figures use PNG/TIFF.
- Reusable: code, package versions, logs, checksums and figure source data are present. Code is MIT licensed; generated tables and figures are CC-BY 4.0.
