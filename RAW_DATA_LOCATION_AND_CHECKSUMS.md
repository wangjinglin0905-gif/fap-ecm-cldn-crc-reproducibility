# Raw-data retention and checksums

The complete local audit set contains 408 public-data files totalling 6,096,798,339 bytes (5.678 GiB).

Raw files are intentionally excluded from the GitHub-ready code ZIP and the processed-data ZIP because several source repositories retain their own redistribution terms and the files are large.

Use:

- `manifests/public_data_sources_manifest.csv` for source repositories, accessions and URLs;
- `manifests/raw_input_manifest_r461.csv` for relative filenames, exact byte sizes and SHA-256 checksums;
- `manifests/valdeolivas_download_manifest.csv` for the source-provided Valdeolivas spatial checksums;
- `scripts/11_download_valdeolivas_spatial.py`, `scripts/12_gse166555_validation.py` and `scripts/download_segmented.py` for reproducible retrieval steps.

Before rerunning the pipeline, place the raw files below `work/reproducibility/inputs/` or set the environment variables documented in `README.md`.
