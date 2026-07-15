#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 300)
set.seed(20260714)
user_library <- Sys.getenv("FAP_R_LIBRARY", unset = "")
if (nzchar(user_library)) .libPaths(c(user_library, .libPaths()))

suppressPackageStartupMessages(library(data.table))

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = TRUE)
output_dir <- file.path(task_root, "results", "L0_TCGA_primary_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- Sys.getenv("FAP_TCGA_EXPRESSION", file.path(task_root, "inputs", "TCGA_COADREAD_expression.txt.gz"))
cms_file <- file.path(task_root, "results", "L0_TCGA", "cms_classifier_results.csv")

stopifnot(file.exists(expression_file), file.exists(cms_file))

raw_table <- fread(expression_file)
gene_symbols <- raw_table[[1]]
expression_matrix <- as.matrix(raw_table[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_symbols

sample_ids <- colnames(expression_matrix)
sample_type <- substr(sample_ids, 14, 15)
keep <- sample_type %in% c("01", "11")
expression_matrix <- expression_matrix[, keep, drop = FALSE]
sample_ids <- sample_ids[keep]
sample_type <- sample_type[keep]
sample_scope <- ifelse(sample_type == "01", "Primary_tumour", "Adjacent_normal")

score_zmean <- function(matrix_input, genes) {
  gene_z <- t(scale(t(matrix_input[genes, , drop = FALSE]), center = TRUE, scale = TRUE))
  colMeans(gene_z, na.rm = TRUE)
}

fap_genes <- c("FAP", "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "THY1", "PDPN",
               "TAGLN", "ACTA2", "MMP2", "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
deligand_genes <- setdiff(fap_genes, c("TGFB1", "INHBA", "WNT2", "WNT5A"))
cldn_genes <- c("CLDN1", "CLDN2", "CLDN4")

scores <- data.table(
  sample_id = sample_ids,
  sample_scope = sample_scope,
  FAP_CAF = score_zmean(expression_matrix, fap_genes),
  FAP_CAF_de_ligand = score_zmean(expression_matrix, deligand_genes),
  CLDN_core = score_zmean(expression_matrix, cldn_genes)
)

correlation_result <- function(label, x, y) {
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.table(comparison = label, n = sum(complete.cases(x, y)), rho = unname(test$estimate), p_value = test$p.value)
}

score_correlations <- rbindlist(list(
  correlation_result("FAP_CAF_vs_de_ligand_primary_and_normal", scores$FAP_CAF, scores$FAP_CAF_de_ligand),
  correlation_result("de_ligand_vs_CLDN_core_primary_and_normal", scores$FAP_CAF_de_ligand, scores$CLDN_core),
  correlation_result("de_ligand_vs_CLDN_core_primary_only", scores[sample_scope == "Primary_tumour", FAP_CAF_de_ligand], scores[sample_scope == "Primary_tumour", CLDN_core])
))
score_correlations[, fdr_bh_three_tests := p.adjust(p_value, method = "BH")]

tumour_normal <- rbindlist(lapply(c("FAP", "CLDN1", "CLDN2", "CLDN4"), function(gene) {
  tumour <- expression_matrix[gene, sample_scope == "Primary_tumour"]
  normal <- expression_matrix[gene, sample_scope == "Adjacent_normal"]
  test <- wilcox.test(tumour, normal, exact = FALSE)
  data.table(
    gene = gene,
    n_primary_tumour = length(tumour),
    n_adjacent_normal = length(normal),
    median_primary_tumour = median(tumour),
    median_adjacent_normal = median(normal),
    median_difference = median(tumour) - median(normal),
    p_value = test$p.value
  )
}))
tumour_normal[, fdr_bh := p.adjust(p_value, method = "BH")]

cms <- fread(cms_file)
cms <- cms[substr(sample_id, 14, 15) == "01" & !is.na(SSP_predicted) & SSP_predicted != ""]
cms_test <- kruskal.test(FAP_CAF ~ SSP_predicted, data = cms)
cms_summary <- cms[, .(n = .N, median = median(FAP_CAF), q1 = quantile(FAP_CAF, 0.25), q3 = quantile(FAP_CAF, 0.75)), by = SSP_predicted]
cms_overall <- data.table(n = nrow(cms), statistic = unname(cms_test$statistic), df = unname(cms_test$parameter), p_value = cms_test$p.value)

fwrite(scores, file.path(output_dir, "tcga_primary_scores.csv"))
fwrite(score_correlations, file.path(output_dir, "score_correlations.csv"))
fwrite(tumour_normal, file.path(output_dir, "tumour_normal_gene_comparisons.csv"))
fwrite(cms_summary, file.path(output_dir, "cms_score_summary.csv"))
fwrite(cms_overall, file.path(output_dir, "cms_kruskal_wallis.csv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
