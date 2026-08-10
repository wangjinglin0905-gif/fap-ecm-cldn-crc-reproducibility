# Public-data download and placement

Raw public datasets are intentionally excluded. Download them from their authoritative repositories and keep their original filenames or document any change.

Expected layout:

```text
data/
  TCGA/
    TCGA_COADREAD_expression.txt.gz
    COADREAD_clinicalMatrix.txt
  GSE39582/
    GSE39582_series_matrix.txt.gz
    gene2probe.csv
  GSE132465/
    [source matrices or documented processed object]
  GSE166555/
    [source matrices or documented processed object]
  TCGA592/
    TCGA592_target_genes_expression.csv
```

## GSE39582

Download the series matrix from GEO accession GSE39582. Use `scripts/python_utils/inspect_gse39582_metadata.py` to rebuild the sample-to-dataset map and verify that the analysis set contains 566 tumours after exclusion of 19 samples labelled `Non Tumoral`.

## TCGA primary set

Download TCGA-COAD/READ HiSeqV2 expression from UCSC Xena and clinical data from the GDC/Xena clinical matrix. The primary set is defined by TCGA tumour barcodes followed by patient-level deduplication of the first 12 barcode characters. Preserve the downloaded filename and record its SHA-256 hash.

## Alternate TCGA route

Run `scripts/python_utils/download_tcga592_cbioportal.py` from the repository root after installing `requests`. This targeted file is used as a processing-route sensitivity analysis and must not be labelled an independent cohort.

## Single-cell and spatial datasets

Download GSE132465 and GSE166555 from GEO. Obtain the Qi and Valdeolivas spatial data from their cited article repositories or supplementary files. Preserve source annotations and record every transformation used to construct processed objects.

## CPTAC

Retrieve protein abundance through cBioPortal study `coad_cptac_2019`. The primary protein score is matrix3 (COL1A1, COL1A2 and FN1) and excludes FAP.

## Hashes and access log

For each downloaded file, add filename, source URL, access date and SHA-256 to a local `data_download_log.csv`. Do not commit raw files unless the source licence explicitly permits redistribution and the journal/repository size limits are satisfied.
