options(stringsAsFactors = FALSE, warn = 1)
root <- normalizePath(getwd(), winslash = "/")
analysis_dir <- file.path(root, "results", "analysis")
publication_dir <- file.path(root, "results", "publication")
read_required <- function(name) {
  path <- file.path(analysis_dir, name); stopifnot(file.exists(path));
  read.csv(path, check.names = FALSE)
}
bulk <- read_required("Figure1_bulk_source_data.csv")
cptac <- read_required("Figure1_CPTAC_source_data.csv")
tcga <- bulk[bulk$cohort == "TCGA-COADREAD", ]
gse <- bulk[bulk$cohort == "GSE39582", ]
stopifnot(nrow(tcga) == 380L, nrow(gse) == 566L, nrow(cptac) == 97L)
stopifnot(abs(cor(tcga$FAP13, tcga$matrix4, method = "spearman") - 0.9296) < 1e-3)
stopifnot(abs(cor(gse$FAP13, gse$matrix4, method = "spearman") - 0.9134) < 1e-3)
stopifnot(abs(cor(cptac$FAP, cptac$matrix3, method = "spearman") - 0.8116) < 1e-3)
integrity <- read_required("extended_analysis_integrity_checks.csv")
stopifnot(nrow(integrity) == 10L, all(integrity$pass))
tables <- file.path(publication_dir, "tables", c("Table1_dataset_overview.csv",
  "Table2_primary_correlations.csv", "Table3_robustness_analyses.csv",
  "Table4_exploratory_pathways.csv"))
stopifnot(all(file.exists(tables)))
cat("Frozen release validation passed.
")
