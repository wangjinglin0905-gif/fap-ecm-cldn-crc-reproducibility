#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "Usage: Rscript 02_bulk_target_purged_composition.R",
    "<TCGA/deconvolution input_dir> <pinned MCPcounter genes.txt> <output_dir>"
  ))
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
mcp_gene_path <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(MCPcounter)
  library(EPIC)
})

seed <- 2026090310L
n_boot <- 5000L

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
           "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
MATRIX4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
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
  x <- x[, order(colnames(x)), drop = FALSE]
  patients <- patient_barcode(colnames(x))
  x <- x[, !duplicated(patients), drop = FALSE]
  manifest <- data.frame(
    sample = colnames(x), patient = patient_barcode(colnames(x)), project = project,
    stringsAsFactors = FALSE
  )
  list(expression = x, manifest = manifest)
}

zscore_mean <- function(x, genes) {
  genes <- intersect(genes, rownames(x))
  if (!length(genes)) stop("No score genes represented")
  sub <- as.matrix(x[genes, , drop = FALSE])
  z <- t(scale(t(sub), center = TRUE, scale = TRUE))
  z <- z[apply(z, 1, function(v) all(is.finite(v))), , drop = FALSE]
  colMeans(z)
}

read_named_covariate <- function(path, row_name) {
  x <- read.csv(path, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
  hit <- which(tolower(trimws(rownames(x))) == tolower(row_name))
  if (length(hit) != 1L) stop("Expected one row named ", row_name, " in ", path)
  v <- as.numeric(x[hit, ])
  names(v) <- sample_barcode(colnames(x))
  v[!duplicated(names(v))]
}

partial_spearman <- function(x, y, z) {
  if (is.null(dim(z))) z <- matrix(z, ncol = 1L)
  keep <- complete.cases(x, y, z)
  x <- x[keep]; y <- y[keep]; z <- z[keep, , drop = FALSE]
  rx <- rank(x, ties.method = "average")
  ry <- rank(y, ties.method = "average")
  rz <- apply(z, 2, rank, ties.method = "average")
  if (is.null(dim(rz))) rz <- matrix(rz, ncol = 1L)
  design <- cbind(1, rz)
  ex <- qr.resid(qr(design), rx)
  ey <- qr.resid(qr(design), ry)
  r <- suppressWarnings(cor(ex, ey, method = "pearson"))
  df <- length(x) - ncol(rz) - 2L
  t_value <- r * sqrt(df / max(.Machine$double.eps, 1 - r^2))
  c(rho = unname(r), p_value = 2 * pt(abs(t_value), df = df, lower.tail = FALSE), df = df)
}

bootstrap_partial <- function(x, y, z, local_seed) {
  if (is.null(dim(z))) z <- matrix(z, ncol = 1L)
  keep <- complete.cases(x, y, z)
  x <- x[keep]; y <- y[keep]; z <- z[keep, , drop = FALSE]
  set.seed(local_seed)
  vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(length(x), length(x), replace = TRUE)
    vals[[b]] <- partial_spearman(x[idx], y[idx], z[idx, , drop = FALSE])[["rho"]]
  }
  list(
    ci = unname(quantile(vals, c(0.025, 0.975), type = 6, na.rm = TRUE)),
    draws = vals
  )
}

marginal_row <- function(pair, x, y) {
  keep <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.frame(
    pair = pair, adjustment = "none", proxy_definition = "not applicable",
    n = sum(keep), rho = unname(ct$estimate), p_value = ct$p.value,
    ci_low = NA_real_, ci_high = NA_real_, stringsAsFactors = FALSE
  )
}

adjusted_row <- function(pair, x, y, covariates, label, definition, local_seed) {
  common <- Reduce(intersect, c(list(names(x), names(y)), lapply(covariates, names)))
  x2 <- x[common]; y2 <- y[common]
  z <- do.call(cbind, lapply(covariates, function(v) v[common]))
  keep <- complete.cases(x2, y2, z)
  x2 <- x2[keep]; y2 <- y2[keep]; z <- z[keep, , drop = FALSE]
  est <- partial_spearman(x2, y2, z)
  boot <- bootstrap_partial(x2, y2, z, local_seed)
  data.frame(
    pair = pair, adjustment = label, proxy_definition = definition,
    n = length(x2), rho = est[["rho"]], p_value = est[["p_value"]],
    ci_low = boot$ci[[1]], ci_high = boot$ci[[2]], stringsAsFactors = FALSE
  )
}

coad <- read_xena(file.path(input_dir, "TCGA.COAD.HiSeqV2.gz"), "TCGA-COAD")
read <- read_xena(file.path(input_dir, "TCGA.READ.HiSeqV2.gz"), "TCGA-READ")
common_genes <- intersect(rownames(coad$expression), rownames(read$expression))
expression <- cbind(coad$expression[common_genes, , drop = FALSE],
                    read$expression[common_genes, , drop = FALSE])
if (anyDuplicated(colnames(expression))) stop("Duplicate primary patient barcodes remain")
manifest <- rbind(coad$manifest, read$manifest)
write.csv(manifest, file.path(out_dir, "tcga_primary_sample_manifest.csv"), row.names = FALSE)

sen_source <- trimws(readLines(file.path(input_dir, "senmayo_genes.txt"), warn = FALSE))
sen_source <- sen_source[nzchar(sen_source)]
sen_nonoverlap <- setdiff(sen_source, FAP13)
sen_alias <- ifelse(sen_nonoverlap == "CXCL8", "IL8", sen_nonoverlap)
SENMAYO <- unique(intersect(sen_alias, rownames(expression)))

scores <- data.frame(
  FAP13 = zscore_mean(expression, FAP13),
  matrix4 = zscore_mean(expression, MATRIX4),
  fib5 = zscore_mean(expression, FIB5),
  SenMayo = zscore_mean(expression, SENMAYO),
  check.names = FALSE
)
rownames(scores) <- colnames(expression)
write.csv(cbind(sample = rownames(scores), scores),
          file.path(out_dir, "tcga_primary_target_scores.csv"), row.names = FALSE)

mcp_markers <- read.delim(mcp_gene_path, check.names = FALSE, quote = "\"",
                          stringsAsFactors = FALSE, colClasses = "character")
mcp_fib_markers <- unique(mcp_markers[mcp_markers[["Cell population"]] == "Fibroblasts", "HUGO symbols"])
epic_sig <- unique(EPIC::TRef$sigGenes)

target_sets <- list(
  FAP13_matrix4 = unique(c(FAP13, MATRIX4)),
  SenMayo_FAP13 = unique(c(SENMAYO, FAP13)),
  SenMayo_matrix4 = unique(c(SENMAYO, MATRIX4)),
  global_disjoint = unique(c(SENMAYO, FAP13, MATRIX4))
)

overlap_rows <- list()
k <- 1L
for (nm in names(target_sets)) {
  target <- target_sets[[nm]]
  for (proxy in c("MCP-counter Fibroblasts", "EPIC TRef signature genes")) {
    universe <- if (proxy == "MCP-counter Fibroblasts") mcp_fib_markers else epic_sig
    ov <- intersect(universe, target)
    overlap_rows[[k]] <- data.frame(
      proxy = proxy, target_set = nm, proxy_feature_count = length(universe),
      target_feature_count = length(target), overlap_count = length(ov),
      overlap_genes = paste(ov, collapse = ";"),
      retained_after_purge = length(setdiff(universe, target)),
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
overlap_audit <- do.call(rbind, overlap_rows)
write.csv(overlap_audit, file.path(out_dir, "composition_proxy_feature_overlap_audit.csv"), row.names = FALSE)

run_mcp <- function(purge = character()) {
  markers <- mcp_markers[!(mcp_markers[["HUGO symbols"]] %in% purge), , drop = FALSE]
  estimate <- MCPcounter::MCPcounter.estimate(
    expression, featuresType = "HUGO_symbols", genes = markers
  )
  if (!"Fibroblasts" %in% rownames(estimate)) stop("MCP-counter fibroblast row absent after purge")
  setNames(as.numeric(estimate["Fibroblasts", ]), colnames(estimate))
}

message("Recomputing full and target-purged MCP-counter proxies")
mcp_variants <- list(full = run_mcp())
for (nm in names(target_sets)) mcp_variants[[nm]] <- run_mcp(target_sets[[nm]])

archived_mcp <- read_named_covariate(file.path(input_dir, "mcpcounter_tumor_scores.csv"), "Fibroblasts")
archived_epic <- read_named_covariate(file.path(input_dir, "epic_tumor_cell_fractions.csv"), "CAFs")

compare_vector <- function(candidate, archived, label) {
  common <- intersect(names(candidate), names(archived))
  d <- candidate[common] - archived[common]
  data.frame(
    candidate = label, n_common = length(common),
    pearson = suppressWarnings(cor(candidate[common], archived[common], method = "pearson")),
    spearman = suppressWarnings(cor(candidate[common], archived[common], method = "spearman")),
    rmse = sqrt(mean(d ^ 2)), max_abs_difference = max(abs(d)),
    stringsAsFactors = FALSE
  )
}

mcp_calibration <- compare_vector(mcp_variants$full, archived_mcp, "MCPcounter logged Xena expression")

run_epic_core <- function(bulk, scale_exprs, purge = character()) {
  ref <- EPIC::TRef
  keep_ref <- !(rownames(ref$refProfiles) %in% purge)
  ref$refProfiles <- ref$refProfiles[keep_ref, , drop = FALSE]
  ref$refProfiles.var <- ref$refProfiles.var[keep_ref, , drop = FALSE]
  ref$sigGenes <- setdiff(unique(ref$sigGenes), purge)
  bulk <- bulk[!(rownames(bulk) %in% purge), , drop = FALSE]
  fit <- EPIC::EPIC(bulk = as.matrix(bulk), reference = ref,
                    scaleExprs = scale_exprs, withOtherCells = TRUE)
  cf <- fit$cellFractions
  if (!"CAFs" %in% colnames(cf)) stop("EPIC CAF column absent")
  setNames(as.numeric(cf[, "CAFs"]), rownames(cf))
}

message("Calibrating EPIC preprocessing against the archived covariate")
linearized_expression <- 2 ^ as.matrix(expression) - 1
linearized_expression[linearized_expression < 0] <- 0
epic_inputs <- list(
  logged_scale_TRUE = list(bulk = expression, scale = TRUE),
  logged_scale_FALSE = list(bulk = expression, scale = FALSE),
  linearized_scale_TRUE = list(bulk = linearized_expression, scale = TRUE),
  linearized_scale_FALSE = list(bulk = linearized_expression, scale = FALSE)
)
epic_candidates <- list()
epic_cal_rows <- list()
for (nm in names(epic_inputs)) {
  message("  EPIC candidate: ", nm)
  candidate <- run_epic_core(epic_inputs[[nm]]$bulk, epic_inputs[[nm]]$scale)
  epic_candidates[[nm]] <- candidate
  epic_cal_rows[[nm]] <- compare_vector(candidate, archived_epic, nm)
}
epic_calibration <- do.call(rbind, epic_cal_rows)
epic_calibration <- epic_calibration[order(-epic_calibration$pearson, epic_calibration$rmse), , drop = FALSE]
chosen_epic <- epic_calibration$candidate[[1]]
if (!is.finite(epic_calibration$pearson[[1]]) || epic_calibration$pearson[[1]] < 0.95) {
  stop("No EPIC preprocessing candidate reproduced the archived CAF ranking adequately")
}
write.csv(rbind(mcp_calibration, epic_calibration),
          file.path(out_dir, "deconvolution_preprocessing_calibration.csv"), row.names = FALSE)

message("Chosen EPIC preprocessing: ", chosen_epic)
chosen_input <- epic_inputs[[chosen_epic]]
epic_variants <- list(full = epic_candidates[[chosen_epic]])
for (nm in names(target_sets)) {
  message("  EPIC purged variant: ", nm)
  epic_variants[[nm]] <- run_epic_core(
    chosen_input$bulk, chosen_input$scale, target_sets[[nm]]
  )
}

all_samples <- rownames(scores)
covariates <- data.frame(
  sample = all_samples,
  MCP_full = mcp_variants$full[all_samples],
  MCP_FAP13_matrix4_purged = mcp_variants$FAP13_matrix4[all_samples],
  MCP_SenMayo_FAP13_purged = mcp_variants$SenMayo_FAP13[all_samples],
  MCP_SenMayo_matrix4_purged = mcp_variants$SenMayo_matrix4[all_samples],
  MCP_global_disjoint = mcp_variants$global_disjoint[all_samples],
  EPIC_full = epic_variants$full[all_samples],
  EPIC_FAP13_matrix4_purged = epic_variants$FAP13_matrix4[all_samples],
  EPIC_SenMayo_FAP13_purged = epic_variants$SenMayo_FAP13[all_samples],
  EPIC_SenMayo_matrix4_purged = epic_variants$SenMayo_matrix4[all_samples],
  EPIC_global_disjoint = epic_variants$global_disjoint[all_samples],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write.csv(covariates, file.path(out_dir, "target_purged_composition_covariates.csv"), row.names = FALSE)

series <- lapply(scores, function(v) setNames(v, rownames(scores)))
pairs <- list(
  FAP13_matrix4 = list(x = series$FAP13, y = series$matrix4),
  SenMayo_FAP13 = list(x = series$SenMayo, y = series$FAP13),
  SenMayo_matrix4 = list(x = series$SenMayo, y = series$matrix4)
)

results <- list()
counter <- 1L
for (i in seq_along(pairs)) {
  nm <- names(pairs)[[i]]
  x <- pairs[[i]]$x
  y <- pairs[[i]]$y
  pair_key <- nm
  results[[counter]] <- marginal_row(pair_key, x, y); counter <- counter + 1L
  results[[counter]] <- adjusted_row(
    pair_key, x, y, list(series$fib5), "fib5 transcript score",
    "five genes disjoint from FAP13, matrix4 and SenMayo", seed + i * 100L + 1L
  ); counter <- counter + 1L

  purged_mcp <- mcp_variants[[nm]]
  purged_epic <- epic_variants[[nm]]
  adjustment_specs <- list(
    list("MCP recomputed full", list(mcp_variants$full), "official fibroblast marker score; target overlap retained"),
    list("MCP pair-purged", list(purged_mcp), "official fibroblast marker score after removal of both target sets"),
    list("MCP global-disjoint", list(mcp_variants$global_disjoint), "official fibroblast marker score after removal of all three manuscript target sets"),
    list("EPIC recomputed full", list(epic_variants$full), "TRef CAF fraction; target overlap retained"),
    list("EPIC pair-purged", list(purged_epic), "TRef CAF fraction after removal of both target sets from signatures, reference scaling and bulk scaling"),
    list("EPIC global-disjoint", list(epic_variants$global_disjoint), "TRef CAF fraction after removal of all three manuscript target sets"),
    list("MCP+EPIC recomputed full", list(mcp_variants$full, epic_variants$full), "joint full official proxies"),
    list("MCP+EPIC pair-purged", list(purged_mcp, purged_epic), "joint proxies purged for this target pair"),
    list("MCP+EPIC global-disjoint", list(mcp_variants$global_disjoint, epic_variants$global_disjoint), "joint proxies disjoint from all manuscript target sets"),
    list("MCP archived legacy", list(archived_mcp), "archived covariate; exact primary-sample intersection"),
    list("EPIC archived legacy", list(archived_epic), "archived covariate; exact primary-sample intersection"),
    list("MCP+EPIC archived legacy", list(archived_mcp, archived_epic), "archived covariates; exact primary-sample intersection")
  )
  for (j in seq_along(adjustment_specs)) {
    spec <- adjustment_specs[[j]]
    results[[counter]] <- adjusted_row(
      pair_key, x, y, spec[[2]], spec[[1]], spec[[3]], seed + i * 100L + j + 1L
    )
    counter <- counter + 1L
  }
}

results <- do.call(rbind, results)
write.csv(results, file.path(out_dir, "target_purged_partial_correlations.csv"), row.names = FALSE)

provenance <- c(
  paste0("R_version=", R.version.string),
  paste0("libPaths=", paste(.libPaths(), collapse = " | ")),
  paste0("MCPcounter_version=", as.character(packageVersion("MCPcounter"))),
  paste0("EPIC_version=", as.character(packageVersion("EPIC"))),
  paste0("MCPcounter_source_commit=b6eac73e91c246fcff0bb1a5c68a816cd588fc48"),
  paste0("MCPcounter_genes_sha256=3BF042D9C1637CA2C1A6B8D1D1B2E4A14B97A6A7553B19A94FE47E4F37F3516A"),
  paste0("TCGA_primary_samples=", nrow(scores)),
  paste0("SenMayo_represented_nonoverlap=", length(SENMAYO)),
  paste0("bootstrap_resamples=", n_boot),
  paste0("seed=", seed),
  paste0("chosen_EPIC_preprocessing=", chosen_epic)
)
writeLines(provenance, file.path(out_dir, "analysis_provenance.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_bulk_composition.txt"))

cat("Completed target-purged composition rerun.\n")
print(overlap_audit, row.names = FALSE)
print(rbind(mcp_calibration, epic_calibration), digits = 5, row.names = FALSE)
print(results, digits = 5, row.names = FALSE)
