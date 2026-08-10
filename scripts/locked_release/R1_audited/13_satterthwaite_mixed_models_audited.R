#!/usr/bin/env Rscript

# Rebuild the GSE132465 fibroblast-lineage mixed models from the archived
# Monocle3 cell_data_set and obtain Satterthwaite degrees of freedom.
#
# Two normalisation variants are emitted:
#   1. corrected_full_library: log1p(raw UMI / total cellular UMI * 10,000),
#      which matches the manuscript's stated method;
#   2. legacy_target_gene_denominator: the historical seven-gene denominator,
#      retained only to audit/reproduce previously reported slopes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript satterthwaite_mixed_models_audited.R R_LIB CDS_RDS OUTPUT_DIR")
}
.libPaths(c(normalizePath(args[[1]], mustWork = TRUE), .libPaths()))
cds_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(monocle3)
  library(SummarizedExperiment)
  library(lmerTest)
  library(Matrix)
})

FIB_LINEAGE <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")
MATRIX4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
RECEPTOR2 <- c("SDC4", "CD44")
TARGET <- c("FAP", MATRIX4, RECEPTOR2)

zscore <- function(v) {
  as.numeric((v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
}

build_cell_data <- function(raw_counts, metadata, denominator, label) {
  if (any(denominator <= 0 | is.na(denominator))) {
    keep <- denominator > 0 & !is.na(denominator)
    raw_counts <- raw_counts[, keep, drop = FALSE]
    metadata <- metadata[keep, , drop = FALSE]
    denominator <- denominator[keep]
  }
  norm <- log1p(t(t(as.matrix(raw_counts)) / denominator) * 10000)
  gene_z <- t(apply(norm, 1, zscore))
  dat <- data.frame(
    Patient = factor(metadata$Patient),
    FAP = as.numeric(gene_z["FAP", ]),
    matrix4 = colMeans(gene_z[MATRIX4, , drop = FALSE]),
    receptor2 = colMeans(gene_z[RECEPTOR2, , drop = FALSE]),
    normalisation = label,
    stringsAsFactors = FALSE
  )
  if (!all(complete.cases(dat))) stop("Non-complete model data after ", label)
  dat
}

fit_one <- function(dat, outcome) {
  formula <- as.formula(paste0(outcome, " ~ FAP + (1 | Patient)"))
  model <- lmerTest::lmer(formula, data = dat, REML = TRUE)
  coef_table <- summary(model)$coefficients
  row <- coef_table["FAP", ]
  vc <- as.data.frame(VarCorr(model))
  data.frame(
    normalisation = unique(dat$normalisation),
    model = paste0(outcome, " ~ FAP + (1|Patient)"),
    n_cells = nrow(dat),
    n_patients = nlevels(dat$Patient),
    estimate = unname(row[["Estimate"]]),
    std_error = unname(row[["Std. Error"]]),
    df_satterthwaite = unname(row[["df"]]),
    t_value = unname(row[["t value"]]),
    p_satterthwaite = unname(row[["Pr(>|t|)"]]),
    patient_random_intercept_sd = vc$sdcor[vc$grp == "Patient"][[1]],
    residual_sd = sigma(model),
    singular = lme4::isSingular(model, tol = 1e-4),
    stringsAsFactors = FALSE
  )
}

cds <- readRDS(cds_path)
metadata <- as.data.frame(colData(cds))
counts <- assay(cds, "counts")
if (!all(TARGET %in% rownames(counts))) {
  stop("Missing required genes: ", paste(setdiff(TARGET, rownames(counts)), collapse = ", "))
}

keep <- metadata$Class == "Tumor" & metadata$Cell_subtype %in% FIB_LINEAGE
metadata_keep <- metadata[keep, , drop = FALSE]
target_counts <- counts[TARGET, keep, drop = FALSE]
full_library_denominator <- Matrix::colSums(counts[, keep, drop = FALSE])
legacy_denominator <- Matrix::colSums(target_counts)

dat_corrected <- build_cell_data(target_counts, metadata_keep, full_library_denominator,
                                 "corrected_full_library")
dat_legacy <- build_cell_data(target_counts, metadata_keep, legacy_denominator,
                              "legacy_target_gene_denominator")

results <- do.call(rbind, list(
  fit_one(dat_corrected, "matrix4"),
  fit_one(dat_corrected, "receptor2"),
  fit_one(dat_legacy, "matrix4"),
  fit_one(dat_legacy, "receptor2")
))

write.csv(results, file.path(output_dir, "satterthwaite_mixed_model_results.csv"), row.names = FALSE)
write.csv(dat_corrected, file.path(output_dir, "mixed_model_cell_data_corrected.csv"), row.names = FALSE)
write.csv(dat_legacy, file.path(output_dir, "mixed_model_cell_data_legacy_audit.csv"), row.names = FALSE)

patient_summary <- aggregate(cbind(FAP, matrix4, receptor2) ~ Patient,
                             data = dat_corrected, FUN = mean)
patient_n <- as.data.frame(table(dat_corrected$Patient), stringsAsFactors = FALSE)
colnames(patient_n) <- c("Patient", "n_cells")
patient_summary <- merge(patient_summary, patient_n, by = "Patient", all.x = TRUE)
write.csv(patient_summary, file.path(output_dir, "mixed_model_patient_summary_corrected.csv"), row.names = FALSE)

# Rebuild the four prespecified patient-level comparisons with the same full-
# library normalisation. Epithelial receptor2 is z-scored within the tumour
# epithelial compartment before patient aggregation.
epi_keep <- metadata$Class == "Tumor" & grepl("Epithelial", metadata$Cell_subtype)
epi_dat <- data.frame()
patient_valid <- data.frame()
patient_results <- data.frame()
if (sum(epi_keep) > 0L) {
  epi_metadata <- metadata[epi_keep, , drop = FALSE]
  epi_counts <- counts[RECEPTOR2, epi_keep, drop = FALSE]
  epi_denominator <- Matrix::colSums(counts[, epi_keep, drop = FALSE])
  epi_norm <- log1p(t(t(as.matrix(epi_counts)) / epi_denominator) * 10000)
  epi_gene_z <- t(scale(t(epi_norm), center = TRUE, scale = TRUE))
  epi_dat <- data.frame(
    Patient = factor(epi_metadata$Patient),
    epithelial_receptor2 = colMeans(epi_gene_z[RECEPTOR2, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  epi_summary <- aggregate(epithelial_receptor2 ~ Patient, data = epi_dat, FUN = mean)
  epi_n <- as.data.frame(table(epi_dat$Patient), stringsAsFactors = FALSE)
  colnames(epi_n) <- c("Patient", "n_epi")
  epi_summary <- merge(epi_summary, epi_n, by = "Patient", all.x = TRUE)

  names(patient_summary)[names(patient_summary) == "n_cells"] <- "n_fib"
  patient_complete <- merge(patient_summary, epi_summary, by = "Patient", all = TRUE)
  patient_valid <- patient_complete[
    patient_complete$n_fib >= 20 & patient_complete$n_epi >= 20 &
      complete.cases(patient_complete[, c("FAP", "matrix4", "receptor2", "epithelial_receptor2")]),
    , drop = FALSE
  ]
  patient_pairs <- list(
    list("Stromal FAP vs stromal matrix4", patient_valid$FAP, patient_valid$matrix4),
    list("Stromal FAP vs epithelial receptor2", patient_valid$FAP, patient_valid$epithelial_receptor2),
    list("Stromal matrix4 vs epithelial receptor2", patient_valid$matrix4, patient_valid$epithelial_receptor2),
    list("Stromal FAP vs stromal receptor2", patient_valid$FAP, patient_valid$receptor2)
  )
  patient_results <- do.call(rbind, lapply(patient_pairs, function(item) {
    test <- suppressWarnings(cor.test(item[[2]], item[[3]], method = "spearman", exact = FALSE))
    data.frame(comparison = item[[1]], n_patients = length(item[[2]]),
               rho = unname(test$estimate), p_value = test$p.value,
               stringsAsFactors = FALSE)
  }))
  patient_results$fdr_bh <- p.adjust(patient_results$p_value, method = "BH")
  write.csv(patient_complete, file.path(output_dir, "patient_compartment_scores_corrected.csv"), row.names = FALSE)
  write.csv(patient_results, file.path(output_dir, "patient_level_correlations_corrected.csv"), row.names = FALSE)
} else {
  writeLines("The archived Monocle3 object contains stromal/endothelial cells but no epithelial cells; cross-compartment patient-level results require the public full UMI matrix.",
             file.path(output_dir, "patient_level_cross_compartment_note.txt"))
}

provenance <- c(
  paste0("R_version=", R.version.string),
  paste0("monocle3=", as.character(packageVersion("monocle3"))),
  paste0("lme4=", as.character(packageVersion("lme4"))),
  paste0("lmerTest=", as.character(packageVersion("lmerTest"))),
  paste0("SummarizedExperiment=", as.character(packageVersion("SummarizedExperiment"))),
  paste0("source_rds=", cds_path),
  paste0("source_rds_md5=", unname(tools::md5sum(cds_path))),
  paste0("cds_dimensions=", paste(dim(cds), collapse = "x")),
  paste0("tumour_fibroblast_lineage_cells=", nrow(dat_corrected)),
  paste0("tumour_epithelial_cells=", nrow(epi_dat)),
  paste0("patients=", nlevels(dat_corrected$Patient)),
  paste0("patients_with_at_least_20_cells_per_compartment=", nrow(patient_valid)),
  paste0("fibroblast_lineage_definition=", paste(FIB_LINEAGE, collapse = ";")),
  "primary_normalisation=log1p(raw UMI / full-cell UMI total * 10000)",
  "historical_normalisation_audit=log1p(raw target-gene UMI / seven-target-gene UMI total * 10000)"
)
writeLines(provenance, file.path(output_dir, "satterthwaite_provenance.txt"), useBytes = TRUE)

print(provenance)
print(results, digits = 8, row.names = FALSE)
if (nrow(patient_results)) print(patient_results, digits = 8, row.names = FALSE)
