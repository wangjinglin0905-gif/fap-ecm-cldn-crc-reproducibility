#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript 05_rescore_scrna_common_core.R",
    "<GSE132 Seurat RDS> <common-core genes.txt>",
    "<GSE166 target_gene_pseudobulk.csv> <output directory>"
  ))
}

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

gse132_rds <- normalizePath(args[[1]], mustWork = TRUE)
core_path <- normalizePath(args[[2]], mustWork = TRUE)
gse166_path <- normalizePath(args[[3]], mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

core <- unique(trimws(readLines(core_path, warn = FALSE)))
core <- core[nzchar(core)]
if (length(core) < 100L) stop("Unexpectedly small common core")

score_gene_z <- function(mat, genes, cells) {
  present <- intersect(genes, rownames(mat))
  x <- as.matrix(mat[present, cells, drop = FALSE])
  mu <- rowMeans(x)
  sig <- apply(x, 1, sd)
  keep <- is.finite(sig) & sig > 0
  if (!any(keep)) stop("All common-core genes have zero variance")
  z <- sweep(x[keep, , drop = FALSE], 1, mu[keep], "-")
  z <- sweep(z, 1, sig[keep], "/")
  list(score = colMeans(z), represented = present, retained = present[keep], zero_variance = present[!keep])
}

paired_summary <- function(df, cohort, local_seed) {
  wide <- reshape(df[, c("patient", "compartment", "score")],
                  idvar = "patient", timevar = "compartment", direction = "wide")
  required <- c("score.Fibroblast", "score.Epithelial")
  if (!all(required %in% names(wide))) stop("Missing paired compartments for ", cohort)
  wide <- wide[complete.cases(wide[, required]), , drop = FALSE]
  wide$difference <- wide$score.Fibroblast - wide$score.Epithelial
  set.seed(local_seed)
  boot_mean <- replicate(10000L, mean(sample(wide$difference, nrow(wide), replace = TRUE)))
  wt <- suppressWarnings(wilcox.test(wide$difference, mu = 0, alternative = "two.sided",
                                     exact = TRUE, correct = FALSE))
  summary <- data.frame(
    cohort = cohort,
    score = "SenMayo common core",
    common_core_n = length(core),
    patients = nrow(wide),
    fibroblast_mean = mean(wide$score.Fibroblast),
    epithelial_mean = mean(wide$score.Epithelial),
    paired_mean_difference = mean(wide$difference),
    paired_median_difference = median(wide$difference),
    mean_difference_ci_low = unname(quantile(boot_mean, 0.025, type = 6)),
    mean_difference_ci_high = unname(quantile(boot_mean, 0.975, type = 6)),
    positive_differences = sum(wide$difference > 0),
    negative_differences = sum(wide$difference < 0),
    zero_differences = sum(wide$difference == 0),
    wilcoxon_signed_rank_p = wt$p.value,
    stringsAsFactors = FALSE
  )
  list(summary = summary, paired = wide)
}

message("Loading GSE132465 Seurat object")
obj <- readRDS(gse132_rds)
meta <- obj@meta.data
rna_data <- LayerData(obj[["RNA"]], layer = "data")
fib_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")
is_fib <- meta$Cell_subtype %in% fib_subtypes
is_epi <- meta$Cell_type == "Epithelial cells"
is_tumour <- meta$Class == "Tumor"
cells_132 <- rownames(meta)[is_tumour & (is_fib | is_epi)]
comp_132 <- ifelse(meta[cells_132, "Cell_subtype"] %in% fib_subtypes, "Fibroblast", "Epithelial")
scored_132 <- score_gene_z(rna_data, core, cells_132)
cell_132 <- data.frame(
  patient = as.character(meta[cells_132, "Patient"]),
  compartment = comp_132,
  score = scored_132$score,
  stringsAsFactors = FALSE
)
patient_132 <- aggregate(score ~ patient + compartment, cell_132, mean)
sum_132 <- paired_summary(patient_132, "GSE132465", 2026090311L)
rm(obj, rna_data, cell_132)
gc()

message("Reconstructing GSE166555 patient-compartment pseudobulk")
pseudo <- read.csv(gse166_path, check.names = FALSE)
pseudo <- pseudo[pseudo$gene %in% core, , drop = FALSE]
pseudo$key <- paste(pseudo$patient, pseudo$compartment, sep = "||")
keys <- unique(pseudo$key)
log_cpm <- matrix(NA_real_, nrow = length(keys), ncol = length(core),
                  dimnames = list(keys, core))
for (i in seq_len(nrow(pseudo))) {
  log_cpm[pseudo$key[[i]], pseudo$gene[[i]]] <- log1p(
    as.numeric(pseudo$count_sum[[i]]) / as.numeric(pseudo$total_umi[[i]]) * 1e6
  )
}
if (anyNA(log_cpm)) stop("Incomplete GSE166555 common-core pseudobulk matrix")
mu_166 <- colMeans(log_cpm)
sd_166 <- apply(log_cpm, 2, sd)
keep_166 <- is.finite(sd_166) & sd_166 > 0
z_166 <- sweep(log_cpm[, keep_166, drop = FALSE], 2, mu_166[keep_166], "-")
z_166 <- sweep(z_166, 2, sd_166[keep_166], "/")
parts <- strsplit(rownames(z_166), "||", fixed = TRUE)
patient_166 <- data.frame(
  patient = vapply(parts, `[[`, character(1), 1L),
  compartment = vapply(parts, `[[`, character(1), 2L),
  score = rowMeans(z_166),
  stringsAsFactors = FALSE
)
sum_166 <- paired_summary(patient_166, "GSE166555", 2026090312L)

patient_scores <- rbind(
  transform(patient_132, cohort = "GSE132465", common_core_n = length(core)),
  transform(patient_166, cohort = "GSE166555", common_core_n = length(core))
)
patient_scores <- patient_scores[, c("cohort", "patient", "compartment", "score", "common_core_n")]
paired_values <- rbind(
  transform(sum_132$paired, cohort = "GSE132465"),
  transform(sum_166$paired, cohort = "GSE166555")
)
paired_values <- paired_values[, c("cohort", "patient", "score.Fibroblast", "score.Epithelial", "difference")]
summaries <- rbind(sum_132$summary, sum_166$summary)
coverage <- data.frame(
  cohort = c("GSE132465", "GSE166555"),
  requested_core_n = length(core),
  represented_core_n = c(length(scored_132$represented), ncol(log_cpm)),
  retained_nonzero_variance_n = c(length(scored_132$retained), sum(keep_166)),
  zero_variance_genes = c(paste(scored_132$zero_variance, collapse = ";"),
                          paste(colnames(log_cpm)[!keep_166], collapse = ";")),
  stringsAsFactors = FALSE
)

write.csv(patient_scores, file.path(out_dir, "scrna_common_core_patient_compartment.csv"), row.names = FALSE)
write.csv(paired_values, file.path(out_dir, "scrna_common_core_paired_values.csv"), row.names = FALSE)
write.csv(summaries, file.path(out_dir, "scrna_common_core_summary.csv"), row.names = FALSE)
write.csv(coverage, file.path(out_dir, "scrna_common_core_scoring_coverage.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_scrna_common_core.txt"))

print(coverage, row.names = FALSE)
print(summaries, digits = 6, row.names = FALSE)
