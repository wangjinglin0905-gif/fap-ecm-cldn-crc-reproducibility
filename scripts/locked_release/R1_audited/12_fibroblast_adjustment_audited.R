#!/usr/bin/env Rscript

# Audited fibroblast-content sensitivity analysis.
#
# The primary adjusted statistic is the conventional partial Spearman
# correlation: rank-transform x, y, and covariates; regress ranked x and y on
# ranked covariates; calculate Pearson's correlation between both residuals.
# MCP-counter and EPIC covariates are selected by exact row names, never by
# row position. Percentile confidence intervals use 5,000 row bootstraps.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript fibroblast_adjustment_audited.R INPUT_DIR OUTPUT_DIR")
}
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260810L
n_boot <- 5000L
set.seed(seed)

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
           "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
MATRIX4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
RECEPTOR2 <- c("SDC4", "CD44")
FIB5 <- c("PDGFRA", "PDGFRB", "LUM", "DCN", "COL14A1")

sample_barcode <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(parts) {
    if (length(parts) >= 4L) paste(parts[1:4], collapse = "-") else paste(parts, collapse = "-")
  }, character(1))
}

patient_barcode <- function(x) {
  vapply(strsplit(as.character(x), "-", fixed = TRUE), function(parts) {
    if (length(parts) >= 3L) paste(parts[1:3], collapse = "-") else paste(parts, collapse = "-")
  }, character(1))
}

read_xena <- function(path, project) {
  x <- read.delim(gzfile(path), check.names = FALSE, row.names = 1,
                  stringsAsFactors = FALSE)
  barcodes <- sample_barcode(colnames(x))
  code <- vapply(strsplit(barcodes, "-", fixed = TRUE), function(parts) {
    if (length(parts) >= 4L) substr(parts[[4]], 1, 2) else ""
  }, character(1))
  x <- x[, code == "01", drop = FALSE]
  colnames(x) <- sample_barcode(colnames(x))
  ord <- order(colnames(x))
  x <- x[, ord, drop = FALSE]
  patients <- patient_barcode(colnames(x))
  x <- x[, !duplicated(patients), drop = FALSE]
  manifest <- data.frame(
    sample = colnames(x), patient = patient_barcode(colnames(x)), project = project,
    stringsAsFactors = FALSE
  )
  list(expression = x, manifest = manifest)
}

zscore_mean <- function(x, genes) {
  missing <- setdiff(genes, rownames(x))
  if (length(missing)) stop("Required genes absent: ", paste(missing, collapse = ", "))
  sub <- as.matrix(x[genes, , drop = FALSE])
  z <- t(scale(t(sub), center = TRUE, scale = TRUE))
  colMeans(z)
}

read_named_covariate <- function(path, row_name) {
  x <- read.csv(path, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
  matches <- which(tolower(trimws(rownames(x))) == tolower(row_name))
  if (length(matches) != 1L) {
    stop("Expected exactly one row named '", row_name, "' in ", path)
  }
  v <- as.numeric(x[matches, ])
  names(v) <- sample_barcode(colnames(x))
  v <- v[!duplicated(names(v))]
  v
}

partial_spearman <- function(x, y, z) {
  if (is.null(dim(z))) z <- matrix(z, ncol = 1L)
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  rz <- apply(z, 2, rank, ties.method = "average")
  if (is.null(dim(rz))) rz <- matrix(rz, ncol = 1L)
  res_x <- residuals(lm.fit(cbind(1, rz), rx))
  res_y <- residuals(lm.fit(cbind(1, rz), ry))
  r <- cor(res_x, res_y, method = "pearson")
  df <- length(x) - ncol(rz) - 2L
  t_value <- r * sqrt(df / max(.Machine$double.eps, 1 - r^2))
  p <- 2 * pt(abs(t_value), df = df, lower.tail = FALSE)
  c(rho = unname(r), p_value = unname(p), df = df)
}

bootstrap_ci <- function(x, y, z, local_seed) {
  if (is.null(dim(z))) z <- matrix(z, ncol = 1L)
  set.seed(local_seed)
  n <- length(x)
  vals <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    vals[[i]] <- partial_spearman(x[idx], y[idx], z[idx, , drop = FALSE])[["rho"]]
  }
  unname(quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, type = 7))
}

marginal_row <- function(pair, x, y) {
  keep <- complete.cases(x, y)
  result <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(
    pair = pair, adjustment = "none", method = "Spearman", n = sum(keep),
    rho = unname(result$estimate), p_value = result$p.value, df = sum(keep) - 2L,
    ci_low = NA_real_, ci_high = NA_real_, stringsAsFactors = FALSE
  )
}

adjusted_row <- function(pair, x, y, covars, label, local_seed) {
  common <- Reduce(intersect, c(list(names(x), names(y)), lapply(covars, names)))
  x2 <- x[common]
  y2 <- y[common]
  z <- do.call(cbind, lapply(covars, function(v) v[common]))
  keep <- complete.cases(x2, y2, z)
  x2 <- x2[keep]
  y2 <- y2[keep]
  z <- z[keep, , drop = FALSE]
  est <- partial_spearman(x2, y2, z)
  ci <- bootstrap_ci(x2, y2, z, local_seed)
  data.frame(
    pair = pair, adjustment = label,
    method = "partial Spearman (Pearson correlation of rank residuals)",
    n = length(x2), rho = est[["rho"]], p_value = est[["p_value"]],
    df = est[["df"]], ci_low = ci[[1]], ci_high = ci[[2]],
    stringsAsFactors = FALSE
  )
}

coad <- read_xena(file.path(input_dir, "TCGA.COAD.HiSeqV2.gz"), "TCGA-COAD")
read <- read_xena(file.path(input_dir, "TCGA.READ.HiSeqV2.gz"), "TCGA-READ")
common_genes <- intersect(rownames(coad$expression), rownames(read$expression))
expression <- cbind(coad$expression[common_genes, , drop = FALSE],
                    read$expression[common_genes, , drop = FALSE])
if (anyDuplicated(colnames(expression))) stop("Duplicate sample barcodes remain")
manifest <- rbind(coad$manifest, read$manifest)
write.csv(manifest, file.path(output_dir, "primary_sample_manifest.csv"), row.names = FALSE)

senmayo_source <- readLines(file.path(input_dir, "senmayo_genes.txt"), warn = FALSE)
senmayo_source <- trimws(senmayo_source[nzchar(trimws(senmayo_source))])
senmayo_nonoverlap <- setdiff(senmayo_source, FAP13)
senmayo_alias <- ifelse(senmayo_nonoverlap == "CXCL8", "IL8", senmayo_nonoverlap)
senmayo_present <- senmayo_alias[senmayo_alias %in% rownames(expression)]

scores <- data.frame(
  FAP13 = zscore_mean(expression, FAP13),
  matrix4 = zscore_mean(expression, MATRIX4),
  receptor2 = zscore_mean(expression, RECEPTOR2),
  fib5 = zscore_mean(expression, FIB5),
  SenMayo = zscore_mean(expression, senmayo_present),
  check.names = FALSE
)
rownames(scores) <- colnames(expression)
write.csv(cbind(sample = rownames(scores), scores),
          file.path(output_dir, "tcga_primary_scores.csv"), row.names = FALSE)

mcp_fib <- read_named_covariate(file.path(input_dir, "mcpcounter_tumor_scores.csv"), "Fibroblasts")
epic_caf <- read_named_covariate(file.path(input_dir, "epic_tumor_cell_fractions.csv"), "CAFs")
cov_samples <- union(names(mcp_fib), names(epic_caf))
cov_frame <- data.frame(
  sample = cov_samples,
  MCPcounter_Fibroblasts = mcp_fib[cov_samples],
  EPIC_CAFs = epic_caf[cov_samples],
  stringsAsFactors = FALSE
)
write.csv(cov_frame, file.path(output_dir, "named_fibroblast_covariates.csv"), row.names = FALSE)

series <- lapply(scores, function(v) setNames(v, rownames(scores)))
pairs <- list(
  list("FAP13-matrix4", series$FAP13, series$matrix4),
  list("FAP13-receptor2", series$FAP13, series$receptor2),
  list("SenMayo-FAP13", series$SenMayo, series$FAP13),
  list("SenMayo-matrix4", series$SenMayo, series$matrix4)
)

results <- list()
counter <- 1L
for (pair_index in seq_along(pairs)) {
  item <- pairs[[pair_index]]
  pair <- item[[1]]
  x <- item[[2]]
  y <- item[[3]]
  results[[counter]] <- marginal_row(pair, x, y); counter <- counter + 1L
  results[[counter]] <- adjusted_row(pair, x, y, list(series$fib5),
                                     "fib5 transcript score", seed + pair_index * 10L + 1L); counter <- counter + 1L
  results[[counter]] <- adjusted_row(pair, x, y, list(mcp_fib),
                                     "MCP-counter Fibroblasts", seed + pair_index * 10L + 2L); counter <- counter + 1L
  results[[counter]] <- adjusted_row(pair, x, y, list(epic_caf),
                                     "EPIC CAFs", seed + pair_index * 10L + 3L); counter <- counter + 1L
  results[[counter]] <- adjusted_row(pair, x, y, list(mcp_fib, epic_caf),
                                     "MCP-counter Fibroblasts + EPIC CAFs", seed + pair_index * 10L + 4L); counter <- counter + 1L
}
results <- do.call(rbind, results)
write.csv(results, file.path(output_dir, "fibroblast_adjusted_correlations.csv"), row.names = FALSE)

hashes <- tools::md5sum(list.files(input_dir, full.names = TRUE))
provenance <- c(
  paste0("R_version=", R.version.string),
  paste0("seed=", seed),
  paste0("bootstrap_resamples=", n_boot),
  paste0("primary_samples=", ncol(expression)),
  paste0("TCGA_COAD_samples=", sum(manifest$project == "TCGA-COAD")),
  paste0("TCGA_READ_samples=", sum(manifest$project == "TCGA-READ")),
  paste0("SenMayo_source_genes=", length(senmayo_source)),
  paste0("SenMayo_after_FAP13_overlap_removal=", length(senmayo_nonoverlap)),
  paste0("SenMayo_represented_after_aliasing=", length(senmayo_present)),
  paste0("SenMayo_absent_after_aliasing=", paste(setdiff(senmayo_alias, senmayo_present), collapse = ";")),
  paste0("input_md5_", names(hashes), "=", unname(hashes))
)
writeLines(provenance, file.path(output_dir, "analysis_provenance.txt"), useBytes = TRUE)

print(provenance)
print(results, digits = 6, row.names = FALSE)
