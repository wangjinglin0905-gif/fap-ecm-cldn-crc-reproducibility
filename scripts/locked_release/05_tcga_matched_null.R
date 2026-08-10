options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260805)

root <- normalizePath(getwd(), winslash = "/")
expression_file <- file.path(root, "data", "public",
                             "TCGA_COADREAD_expression.txt.gz")
manifest_file <- file.path(root, "results", "analysis",
                           "TCGA_full_primary_manifest.csv")
out_dir <- file.path(root, "results", "analysis")
stopifnot(file.exists(expression_file), file.exists(manifest_file))

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
target_genes <- c(FAP13, matrix4)

spearman <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}

expression <- read.delim(gzfile(expression_file), row.names = 1,
                         check.names = FALSE)
manifest <- read.csv(manifest_file, check.names = FALSE)
stopifnot(all(manifest$sample_id %in% colnames(expression)))
expression <- as.matrix(expression[, manifest$sample_id, drop = FALSE])

gene_mean <- rowMeans(expression)
gene_sd <- apply(expression, 1, sd)
gene_detection <- rowMeans(expression > 0)
eligible <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0 &
  is.finite(gene_detection) & !duplicated(rownames(expression))
expression <- expression[eligible, , drop = FALSE]
gene_mean <- gene_mean[eligible]
gene_sd <- gene_sd[eligible]
gene_detection <- gene_detection[eligible]
stopifnot(all(target_genes %in% rownames(expression)))

z_expression <- t(scale(t(expression)))
z_expression[!is.finite(z_expression)] <- NA_real_
feature_matrix <- cbind(
  mean = as.numeric(scale(gene_mean)),
  log_sd = as.numeric(scale(log(gene_sd))),
  detection = as.numeric(scale(gene_detection))
)
rownames(feature_matrix) <- rownames(expression)

candidate_universe <- setdiff(rownames(expression), target_genes)
candidate_k <- min(500L, length(candidate_universe))
candidate_list <- lapply(target_genes, function(gene) {
  delta <- sweep(feature_matrix[candidate_universe, , drop = FALSE], 2,
                 feature_matrix[gene, ], FUN = "-")
  distance <- rowSums(delta^2)
  candidate_universe[order(distance)[seq_len(candidate_k)]]
})
names(candidate_list) <- target_genes

candidate_diagnostics <- do.call(rbind, lapply(target_genes, function(gene) {
  candidates <- candidate_list[[gene]]
  data.frame(
    target_gene = gene,
    candidate_count = length(candidates),
    target_mean = gene_mean[gene],
    target_sd = gene_sd[gene],
    target_detection = gene_detection[gene],
    candidate_mean_min = min(gene_mean[candidates]),
    candidate_mean_max = max(gene_mean[candidates]),
    candidate_sd_min = min(gene_sd[candidates]),
    candidate_sd_max = max(gene_sd[candidates]),
    candidate_detection_min = min(gene_detection[candidates]),
    candidate_detection_max = max(gene_detection[candidates])
  )
}))
write.csv(candidate_diagnostics,
          file.path(out_dir,
                    "TCGA_matched_null_extended_candidate_diagnostics.csv"),
          row.names = FALSE)

draw_unique <- function(targets, used = character()) {
  selected <- used
  result <- character(length(targets))
  for (index in seq_along(targets)) {
    pool <- setdiff(candidate_list[[targets[index]]], selected)
    if (!length(pool)) stop("Matched candidate pool exhausted")
    result[index] <- sample(pool, 1L)
    selected <- c(selected, result[index])
  }
  result
}

fixed_fap <- colMeans(z_expression[FAP13, , drop = FALSE])
fixed_matrix <- colMeans(z_expression[matrix4, , drop = FALSE])
observed <- spearman(fixed_fap, fixed_matrix)
n_draws <- 10000L
null <- matrix(NA_real_, nrow = n_draws, ncol = 3L,
               dimnames = list(NULL, c("random_vs_random",
                                       "fixed_FAP13_vs_random_matrix4",
                                       "random_FAP13_vs_fixed_matrix4")))
for (iteration in seq_len(n_draws)) {
  random_fap_genes <- draw_unique(FAP13)
  random_matrix_genes <- draw_unique(matrix4, used = random_fap_genes)
  random_fap <- colMeans(
    z_expression[random_fap_genes, , drop = FALSE])
  random_matrix <- colMeans(
    z_expression[random_matrix_genes, , drop = FALSE])
  null[iteration, ] <- c(
    spearman(random_fap, random_matrix),
    spearman(fixed_fap, random_matrix),
    spearman(random_fap, fixed_matrix)
  )
}

null_long <- data.frame(
  iteration = rep(seq_len(n_draws), times = ncol(null)),
  null_type = rep(colnames(null), each = n_draws),
  rho = as.numeric(null)
)
write.csv(null_long,
          file.path(out_dir, "TCGA_matched_null_extended_distribution.csv"),
          row.names = FALSE)

summary_rows <- do.call(rbind, lapply(colnames(null), function(null_type) {
  values <- null[, null_type]
  data.frame(
    null_type = null_type,
    observed_rho = observed,
    random_draws = n_draws,
    null_mean = mean(values),
    null_sd = sd(values),
    null_q025 = unname(quantile(values, 0.025)),
    null_median = median(values),
    null_q975 = unname(quantile(values, 0.975)),
    empirical_p_greater_equal =
      (1 + sum(values >= observed)) / (n_draws + 1),
    empirical_p_two_sided =
      (1 + sum(abs(values) >= abs(observed))) / (n_draws + 1),
    matching_features =
      "mean expression, log standard deviation, and detection rate",
    candidate_pool_per_target_gene = candidate_k
  )
}))
write.csv(summary_rows,
          file.path(out_dir, "TCGA_matched_null_extended_summary.csv"),
          row.names = FALSE)

log_file <- file.path(out_dir, "TCGA_matched_null_extended.log")
sink(log_file)
cat("TCGA MATCHED RANDOM-GENE-SET EMPIRICAL NULL\n")
cat("Run UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
cat("Seed: 20260805\n\n")
print(summary_rows, row.names = FALSE)
cat("\nSession information:\n")
print(sessionInfo())
sink()
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "TCGA_matched_null_extended_sessionInfo.txt"))

