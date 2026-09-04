#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "Usage: Rscript 05b_rescore_gse166_common_core.R",
    "<target_gene_pseudobulk.csv> <common-core genes.txt> <output directory>"
  ))
}

pseudo_path <- normalizePath(args[[1]], mustWork = TRUE)
core_path <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

core <- unique(trimws(readLines(core_path, warn = FALSE)))
core <- core[nzchar(core)]
pseudo <- read.csv(pseudo_path, check.names = FALSE)
cell_ledger <- unique(pseudo[, c("patient", "compartment", "cell_count")])
cell_wide <- reshape(cell_ledger, idvar = "patient", timevar = "compartment", direction = "wide")
eligible_patients <- cell_wide$patient[
  is.finite(cell_wide$cell_count.Fibroblast) & is.finite(cell_wide$cell_count.Epithelial) &
    cell_wide$cell_count.Fibroblast >= 20 & cell_wide$cell_count.Epithelial >= 100
]
pseudo <- pseudo[pseudo$patient %in% eligible_patients, , drop = FALSE]
pseudo <- pseudo[pseudo$gene %in% core, , drop = FALSE]
pseudo$key <- paste(pseudo$patient, pseudo$compartment, sep = "||")
keys <- unique(pseudo$key)
log_cpm <- matrix(NA_real_, nrow = length(keys), ncol = length(core), dimnames = list(keys, core))
for (i in seq_len(nrow(pseudo))) {
  log_cpm[pseudo$key[[i]], pseudo$gene[[i]]] <- log1p(
    as.numeric(pseudo$count_sum[[i]]) / as.numeric(pseudo$total_umi[[i]]) * 1e6
  )
}
if (anyNA(log_cpm)) stop("Incomplete GSE166555 common-core matrix")

means <- colMeans(log_cpm)
sds <- apply(log_cpm, 2, sd)
keep <- is.finite(sds) & sds > 0
z <- sweep(log_cpm[, keep, drop = FALSE], 2, means[keep], "-")
z <- sweep(z, 2, sds[keep], "/")
parts <- strsplit(rownames(z), "||", fixed = TRUE)
patient_scores <- data.frame(
  cohort = "GSE166555",
  patient = vapply(parts, `[[`, character(1), 1L),
  compartment = vapply(parts, `[[`, character(1), 2L),
  score = rowMeans(z),
  cell_count = vapply(keys, function(key) unique(pseudo$cell_count[pseudo$key == key])[[1]], numeric(1)),
  common_core_n = length(core),
  retained_nonzero_variance_n = sum(keep),
  stringsAsFactors = FALSE
)
wide <- reshape(patient_scores[, c("patient", "compartment", "score")],
                idvar = "patient", timevar = "compartment", direction = "wide")
wide <- wide[complete.cases(wide[, c("score.Fibroblast", "score.Epithelial")]), , drop = FALSE]
wide$difference <- wide$score.Fibroblast - wide$score.Epithelial
set.seed(2026090312L)
boot <- replicate(10000L, mean(sample(wide$difference, nrow(wide), replace = TRUE)))
test <- suppressWarnings(wilcox.test(wide$difference, mu = 0, alternative = "two.sided",
                                     exact = TRUE, correct = FALSE))
summary <- data.frame(
  cohort = "GSE166555",
  normalization = "patient-compartment pseudobulk log1p(CPM), gene-wise z across eligible compartments",
  common_core_n = length(core),
  retained_nonzero_variance_n = sum(keep),
  patients = nrow(wide),
  paired_mean_difference = mean(wide$difference),
  paired_median_difference = median(wide$difference),
  mean_difference_ci_low = unname(quantile(boot, 0.025, type = 6)),
  mean_difference_ci_high = unname(quantile(boot, 0.975, type = 6)),
  positive_differences = sum(wide$difference > 0),
  negative_differences = sum(wide$difference < 0),
  zero_differences = sum(wide$difference == 0),
  wilcoxon_signed_rank_statistic = unname(test$statistic),
  wilcoxon_signed_rank_p = test$p.value,
  stringsAsFactors = FALSE
)

paired <- data.frame(
  cohort = "GSE166555",
  patient = wide$patient,
  fibroblast_score = wide$score.Fibroblast,
  epithelial_score = wide$score.Epithelial,
  difference = wide$difference,
  stringsAsFactors = FALSE
)
write.csv(patient_scores, file.path(out_dir, "gse166_common_core_patient_compartment.csv"), row.names = FALSE)
write.csv(paired, file.path(out_dir, "gse166_common_core_paired_values.csv"), row.names = FALSE)
write.csv(summary, file.path(out_dir, "gse166_common_core_summary.csv"), row.names = FALSE)
write.csv(
  data.frame(
    cohort = "GSE166555", requested_core_n = length(core), represented_core_n = ncol(log_cpm),
    retained_nonzero_variance_n = sum(keep), zero_variance_genes = paste(colnames(log_cpm)[!keep], collapse = ";")
  ),
  file.path(out_dir, "gse166_common_core_scoring_coverage.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_gse166_common_core.txt"))

print(summary, digits = 6, row.names = FALSE)
