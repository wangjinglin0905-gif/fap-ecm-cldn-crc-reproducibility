# Exploratory CAF marker programme extension

This module supports manuscript v7.5, supplementary Fig S7, Table S8 and S1 Data. Core analyses and Fig 1-4/S1-S6 are unchanged. No programme met the positive-association criterion in both single-cell cohorts. Bulk associations do not establish a CAF subtype, durable senescence or organoid mechanism.

## Reproduction

Work on a **copy** of this module: scripts write into results/adapted_v02 and overwrite the frozen results in that copy. Python 3.12 with the numerical versions in requirements.txt was used. Do not replace system R libraries. Optional plotting uses R 4.6.1 and ggplot2/patchwork/ragg/svglite; FAP_R_LIBRARY may specify a compatible library.

1. The compact derived counts retain all represented genes in 1,435 GSE132465 and 1,315 GSE166555 eligible fibroblasts. Metadata, full-library checks, patient eligibility and both marker selections are included. Run `python scripts/03_single_cell_associations.py` to reconstruct patient correlations, sensitivity analyses, 5,000-bootstrap intervals, t-based tests, q values and 1,000 matched random programmes. This is compute intensive; preserve the frozen snapshot.
2. Run `python scripts/04b_bulk_from_frozen_scores.py` to reproduce every bulk estimate, interval and q value from the archived 380-patient scores. To regenerate scores from public UCSC Xena TCGA.COAD/READ HiSeqV2 matrices, put the two original .gz matrices in inputs/ or set FAP_TCGA_INPUT_DIR, then run scripts/04_bulk_associations.py. The consumed primary-patient manifest and locked score controls are in the parent core package.
3. Run `Rscript scripts/06_plot_CAF_extension.R .` for Fig S7; scripts/05_environment_and_statistical_QA.R independently checks patient correlations and t tests. The executed independent-check output is preserved in qa/.
4. Optional marker re-extraction: download the publisher Supplementary Tables workbook using inputs/source_download_manifest.json, verify SHA256 and place it in inputs/. Run scripts/02c_adapt_source_detectability.py. The full workbook and article are not redistributed. Frozen marker extracts allow reproduction without downloading the workbook.

## Provenance and limitations

Marker source: Buissant des Amorie et al., Nature 2026, doi:10.1038/s41586-026-10344-7, Supplementary Table 5. The unadapted top-100 sets failed the target detectability gate before testing target associations. Source-only detectability adaptation used pct.1 >=0.25 and pct.1-pct.2 >=0.10, adjusted P <0.05, positive log2FC; top 100 followed by overlap purging without backfilling. Both selections and full coverage are retained. The CAF-specific shared SenMayo set contains 114 genes; it is distinct from the 111-gene cross-platform core analysis.

Scripts 03 and all numerical outputs are copied unchanged from the executed analysis. Scripts 02c/04/05/06 have path portability edits only; 04b is the unchanged statistical block with a frozen-score reader. Full raw GEO/Xena matrices are obtainable from their original public repositories, not duplicated here. The compressed derived arrays are public-data subsets, not identifiable clinical records. Source data retain original source attribution/terms; MIT applies to the authors' original code. No full copyrighted article, private report, unpublished manuscript, authentication material or local home-directory path is included.
