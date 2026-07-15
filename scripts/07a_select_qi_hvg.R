#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
output_dir <- file.path(task_root, "inputs", "harmony")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

qi_matrix_file <- Sys.getenv("FAP_QI_MATRIX", file.path(task_root, "inputs", "41467_2022_29366_MOESM6_ESM.gz"))
stopifnot(file.exists(qi_matrix_file))

message("Loading Qi et al. sparse count matrix")
qi_counts <- readRDS(gzfile(qi_matrix_file))
stopifnot(inherits(qi_counts, "dgCMatrix"), ncol(qi_counts) == 54103)

sum_counts <- numeric(nrow(qi_counts))
sum_squares <- numeric(nrow(qi_counts))
aggregated_counts <- rowsum(qi_counts@x, group = qi_counts@i + 1L, reorder = FALSE)
aggregated_squares <- rowsum(qi_counts@x^2, group = qi_counts@i + 1L, reorder = FALSE)
sum_counts[as.integer(rownames(aggregated_counts))] <- aggregated_counts[, 1]
sum_squares[as.integer(rownames(aggregated_squares))] <- aggregated_squares[, 1]
gene_mean <- sum_counts / ncol(qi_counts)
gene_variance <- pmax((sum_squares - ncol(qi_counts) * gene_mean^2) / (ncol(qi_counts) - 1), 0)
dispersion <- gene_variance / pmax(gene_mean, 1e-8)

stats <- data.table(
  gene = rownames(qi_counts),
  mean_count = gene_mean,
  variance_count = gene_variance,
  dispersion = dispersion
)
stats[, excluded := grepl("^MT-|^RPL|^RPS", gene) | mean_count < 0.01]
setorder(stats, excluded, -dispersion)
selected <- stats[excluded == FALSE][1:3000, gene]

fwrite(stats, file.path(output_dir, "qi_gene_variability.csv"))
writeLines(selected, file.path(output_dir, "qi_top3000_variable_genes.txt"))
message("Selected ", length(selected), " Qi variable genes")
