#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript 06_rescore_spatial_common_core.R",
    "<Route-B root> <common-core genes.txt> <SenMayo source.txt> <output directory>"
  ))
}

suppressPackageStartupMessages({
  library(arrow)
  library(Matrix)
})

route_root <- normalizePath(args[[1]], mustWork = TRUE)
core_path <- normalizePath(args[[2]], mustWork = TRUE)
sen_path <- normalizePath(args[[3]], mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

core <- unique(trimws(readLines(core_path, warn = FALSE)))
core <- core[nzchar(core)]
sen_source <- unique(trimws(readLines(sen_path, warn = FALSE)))
sen_source <- sen_source[nzchar(sen_source)]
fap13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
           "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
sen_nonoverlap <- setdiff(sen_source, fap13)

one_sample_test <- function(x) {
  suppressWarnings(wilcox.test(x, mu = 0, alternative = "two.sided", exact = TRUE,
                               correct = FALSE)$p.value)
}

score_from_gene_effects <- function(stats, genes) {
  hit <- stats[stats$gene %in% genes & is.finite(stats$stromal_minus_tumour_standardized_mean), , drop = FALSE]
  list(effect = mean(hit$stromal_minus_tumour_standardized_mean), retained = nrow(hit))
}

section_rows <- list()
calibration_rows <- list()
counter <- 1L
cal_counter <- 1L

# GSE280315: the stored per-gene z contrasts allow exact reconstruction of
# group-mean score contrasts without reopening the 500k-bin matrices.
g280_patients <- c("P1", "P2", "P5")
g280_h1 <- read.csv(file.path(route_root, "results_full", "GSE280315", "patient_H1_compartment.csv"),
                    check.names = FALSE)
for (patient in g280_patients) {
  stats_path <- file.path(route_root, "results_full", "GSE280315", paste0(patient, "_gene_stats.parquet"))
  stats <- as.data.frame(read_parquet(stats_path))
  old <- score_from_gene_effects(stats, sen_nonoverlap)
  common <- score_from_gene_effects(stats, core)
  archived <- g280_h1$mean_difference[g280_h1$patient == patient]
  if (length(archived) != 1L) stop("Missing archived GSE280315 H1 row for ", patient)
  calibration_rows[[cal_counter]] <- data.frame(
    cohort = "GSE280315", patient = patient, archived_mean_difference = archived,
    reconstructed_mean_difference = old$effect,
    absolute_difference = abs(archived - old$effect), represented_nonoverlap_n = old$retained
  )
  cal_counter <- cal_counter + 1L
  section_rows[[counter]] <- data.frame(
    cohort = "GSE280315", patient = patient, sample = patient,
    definition = "author-deconvolution singlet; exact mean-score reconstruction from stored gene z contrasts",
    common_core_n = length(core), retained_nonzero_variance_n = common$retained,
    stromal_minus_tumour_mean_score = common$effect
  )
  counter <- counter + 1L
}

# GSE334323: rescore the compact Visium matrices directly and retain the frozen
# marker-margin compartment labels from the Route-B spot ledger.
g334_samples <- data.frame(
  patient = c("CRC03", "CRC07", "CRC08"),
  prefix = c("GSM9785429_CRC_03_Tumor", "GSM9785431_CRC_07_Tumor", "GSM9785433_CRC_08_Tumor")
)
g334_h1 <- read.csv(file.path(route_root, "results_full", "GSE334323", "patient_H1_compartment.csv"),
                    check.names = FALSE)

score_matrix <- function(counts, totals, genes) {
  present <- intersect(genes, rownames(counts))
  x <- counts[present, , drop = FALSE]
  lognorm <- t(t(x) / totals)
  lognorm@x <- log1p(lognorm@x * 10000)
  means <- as.numeric(rowMeans(lognorm))
  second <- as.numeric(rowMeans(lognorm ^ 2))
  n <- ncol(lognorm)
  sds <- sqrt(pmax(0, (second - means ^ 2) * n / pmax(1, n - 1)))
  keep <- is.finite(sds) & sds > 0
  dense <- as.matrix(lognorm[keep, , drop = FALSE])
  z <- sweep(dense, 1, means[keep], "-")
  z <- sweep(z, 1, sds[keep], "/")
  list(score = colMeans(z), represented = length(present), retained = sum(keep))
}

for (i in seq_len(nrow(g334_samples))) {
  patient <- g334_samples$patient[[i]]
  prefix <- g334_samples$prefix[[i]]
  raw_dir <- file.path(route_root, "data_raw", "GSE334323")
  matrix_path <- file.path(raw_dir, paste0(prefix, "_matrix.mtx.gz"))
  feature_path <- file.path(raw_dir, paste0(prefix, "_features.tsv.gz"))
  barcode_path <- file.path(raw_dir, paste0(prefix, "_barcodes.tsv.gz"))
  features <- read.delim(gzfile(feature_path), header = FALSE, sep = "\t", quote = "",
                         comment.char = "", stringsAsFactors = FALSE)
  barcodes <- readLines(gzfile(barcode_path), warn = FALSE)
  counts_all <- as(readMM(gzfile(matrix_path)), "dgCMatrix")
  symbols <- trimws(features[[2]])
  target_rows <- symbols %in% sen_nonoverlap
  target_symbols <- symbols[target_rows]
  target_counts <- counts_all[target_rows, , drop = FALSE]
  factor_symbols <- factor(target_symbols, levels = unique(target_symbols))
  aggregator <- sparseMatrix(
    i = as.integer(factor_symbols), j = seq_along(factor_symbols), x = 1,
    dims = c(nlevels(factor_symbols), length(factor_symbols))
  )
  counts <- aggregator %*% target_counts
  rownames(counts) <- levels(factor_symbols)
  colnames(counts) <- barcodes
  spot_path <- file.path(route_root, "data_processed", "GSE334323", paste0(patient, "_spot_scores.parquet"))
  spots <- as.data.frame(read_parquet(spot_path))
  matched <- match(spots$barcode, colnames(counts))
  if (anyNA(matched)) stop("Unmatched GSE334323 processed barcodes for ", patient)
  counts <- counts[, matched, drop = FALSE]
  old <- score_matrix(counts, spots$total_count, sen_nonoverlap)
  common <- score_matrix(counts, spots$total_count, core)
  fib <- spots$compartment == "Stromal"
  tumour <- spots$compartment == "Epithelial"
  old_effect <- mean(old$score[fib]) - mean(old$score[tumour])
  common_effect <- mean(common$score[fib]) - mean(common$score[tumour])
  archived <- g334_h1$mean_difference[g334_h1$patient == patient]
  if (length(archived) != 1L) stop("Missing archived GSE334323 H1 row for ", patient)
  calibration_rows[[cal_counter]] <- data.frame(
    cohort = "GSE334323", patient = patient, archived_mean_difference = archived,
    reconstructed_mean_difference = old_effect,
    absolute_difference = abs(archived - old_effect), represented_nonoverlap_n = old$represented
  )
  cal_counter <- cal_counter + 1L
  section_rows[[counter]] <- data.frame(
    cohort = "GSE334323", patient = patient, sample = patient,
    definition = "frozen stromal-vs-epithelial marker margin; direct common-core rescore",
    common_core_n = length(core), retained_nonzero_variance_n = common$retained,
    stromal_minus_tumour_mean_score = common_effect
  )
  counter <- counter + 1L
  rm(counts_all, target_counts, counts, old, common)
  gc(verbose = FALSE)
}

# Valdeolivas Visium: exact section-level mean-score reconstruction from stored
# per-gene z contrasts, followed by the frozen replicate-to-patient aggregation.
val_dir <- file.path(route_root, "results_full", "Valdeolivas_Visium")
val_stats <- read.csv(gzfile(file.path(val_dir, "section_gene_stats.csv.gz")), check.names = FALSE)
val_h1 <- read.csv(file.path(val_dir, "section_H1_compartment.csv"), check.names = FALSE)
eligible <- val_h1[val_h1$eligible_primary %in% c(TRUE, "TRUE", 1, "1"), c("patient", "sample")]
for (i in seq_len(nrow(eligible))) {
  patient <- eligible$patient[[i]]
  sample_id <- eligible$sample[[i]]
  stats <- val_stats[val_stats$patient == patient & val_stats$sample == sample_id, , drop = FALSE]
  old <- score_from_gene_effects(stats, sen_nonoverlap)
  common <- score_from_gene_effects(stats, core)
  calibration_rows[[cal_counter]] <- data.frame(
    cohort = "Valdeolivas_Visium", patient = patient, archived_mean_difference = NA_real_,
    reconstructed_mean_difference = old$effect, absolute_difference = NA_real_,
    represented_nonoverlap_n = old$retained
  )
  cal_counter <- cal_counter + 1L
  section_rows[[counter]] <- data.frame(
    cohort = "Valdeolivas_Visium", patient = patient, sample = sample_id,
    definition = "pathologist fibroblastic-stroma vs tumour; exact mean-score reconstruction from stored gene z contrasts",
    common_core_n = length(core), retained_nonzero_variance_n = common$retained,
    stromal_minus_tumour_mean_score = common$effect
  )
  counter <- counter + 1L
}

section_effects <- do.call(rbind, section_rows)
calibration <- do.call(rbind, calibration_rows)
if (any(calibration$absolute_difference[is.finite(calibration$absolute_difference)] > 1e-10)) {
  stop("Spatial reconstruction failed archived-score calibration")
}

patient_effects <- aggregate(
  stromal_minus_tumour_mean_score ~ cohort + patient,
  section_effects, mean
)
patient_counts <- aggregate(sample ~ cohort + patient, section_effects, length)
names(patient_counts)[names(patient_counts) == "sample"] <- "n_sections"
patient_effects <- merge(patient_effects, patient_counts, by = c("cohort", "patient"), sort = FALSE)

set.seed(2026090313L)
summaries <- do.call(rbind, lapply(unique(patient_effects$cohort), function(cohort) {
  x <- patient_effects$stromal_minus_tumour_mean_score[patient_effects$cohort == cohort]
  boot <- replicate(10000L, median(sample(x, length(x), replace = TRUE)))
  data.frame(
    cohort = cohort, common_core_n = length(core), patients = length(x),
    mean_patient_effect = mean(x), median_patient_effect = median(x),
    median_ci_low = unname(quantile(boot, 0.025, type = 6)),
    median_ci_high = unname(quantile(boot, 0.975, type = 6)),
    positive_patients = sum(x > 0), negative_patients = sum(x < 0),
    wilcoxon_signed_rank_p = one_sample_test(x),
    stringsAsFactors = FALSE
  )
}))

write.csv(section_effects, file.path(out_dir, "spatial_common_core_section_effects.csv"), row.names = FALSE)
write.csv(patient_effects, file.path(out_dir, "spatial_common_core_patient_effects.csv"), row.names = FALSE)
write.csv(summaries, file.path(out_dir, "spatial_common_core_summary.csv"), row.names = FALSE)
write.csv(calibration, file.path(out_dir, "spatial_common_core_reconstruction_calibration.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_spatial_common_core.txt"))

print(calibration, digits = 7, row.names = FALSE)
print(summaries, digits = 6, row.names = FALSE)
