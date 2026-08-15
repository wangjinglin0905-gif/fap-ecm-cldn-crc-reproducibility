options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript test_signature_specificity.R <GSE132465_seurat.rds> <senmayo_genes.txt> <output_directory>")
}
input_rds <- normalizePath(args[[1]], mustWork = TRUE)
senmayo_file <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fap13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
sasp25 <- c("IL6", "CXCL8", "IL1A", "IL1B", "CCL2", "CCL5", "CXCL1", "CXCL2", "CXCL3", "CXCL10", "MMP1", "MMP3", "MMP9", "MMP10", "MMP13", "SERPINE1", "PLAU", "TIMP2", "VEGFA", "GDF15", "IGFBP3", "TNF", "CSF2", "HGF", "FAS")
senmayo_source <- unique(trimws(readLines(senmayo_file, warn = FALSE)))
senmayo_source <- senmayo_source[nzchar(senmayo_source)]
senmayo_nonoverlap <- setdiff(senmayo_source, fap13)

obj <- readRDS(input_rds)
meta <- obj@meta.data
rna_data <- LayerData(obj[["RNA"]], layer = "data")
fib_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")
is_fib <- meta$Cell_subtype %in% fib_subtypes
is_epi <- meta$Cell_type == "Epithelial cells"
is_tumor <- meta$Class == "Tumor"
cells <- colnames(obj)[is_tumor & (is_fib | is_epi)]
x <- rna_data[, cells, drop = FALSE]

gene_mean <- Matrix::rowMeans(x)
gene_second <- Matrix::rowMeans(x ^ 2)
gene_sd <- sqrt(pmax(gene_second - gene_mean ^ 2, 0))
gene_detect <- Matrix::rowMeans(x > 0)

patients <- sort(unique(meta[cells, "Patient"]))
diff_z <- matrix(NA_real_, nrow = nrow(x), ncol = length(patients), dimnames = list(rownames(x), patients))
for (j in seq_along(patients)) {
  p <- patients[j]
  pcells <- cells[meta[cells, "Patient"] == p]
  fcells <- pcells[meta[pcells, "Cell_subtype"] %in% fib_subtypes]
  ecells <- pcells[meta[pcells, "Cell_type"] == "Epithelial cells"]
  stopifnot(length(fcells) > 0L, length(ecells) > 0L)
  raw_diff <- Matrix::rowMeans(rna_data[, fcells, drop = FALSE]) - Matrix::rowMeans(rna_data[, ecells, drop = FALSE])
  diff_z[, j] <- raw_diff / gene_sd
}

eligible <- is.finite(gene_mean) & is.finite(gene_sd) & gene_sd > 0 & gene_detect >= 0.001
all_signature_genes <- unique(c(senmayo_source, sasp25, fap13, matrix4, "MKI67"))
eligible[names(eligible) %in% all_signature_genes] <- FALSE
candidate_genes <- names(eligible)[eligible]

feature_matrix <- cbind(
  log_mean = log1p(gene_mean),
  log_sd = log1p(gene_sd),
  detection_logit = qlogis(pmin(pmax(gene_detect, 1e-5), 1 - 1e-5))
)
feature_matrix <- scale(feature_matrix)

build_neighbors <- function(target_genes, k = 300L) {
  target_genes <- intersect(target_genes, rownames(feature_matrix))
  target_genes <- target_genes[is.finite(gene_sd[target_genes]) & gene_sd[target_genes] > 0]
  candidate_features <- feature_matrix[candidate_genes, , drop = FALSE]
  neighbors <- lapply(target_genes, function(g) {
    delta <- sweep(candidate_features, 2, feature_matrix[g, ], "-")
    d2 <- rowSums(delta ^ 2)
    candidate_genes[order(d2)[seq_len(min(k, length(d2)))]]
  })
  names(neighbors) <- target_genes
  neighbors
}

draw_unique_set <- function(neighbors) {
  order_idx <- sample(seq_along(neighbors))
  selected <- character(length(neighbors))
  used <- character(0)
  for (idx in order_idx) {
    choices <- setdiff(neighbors[[idx]], used)
    if (!length(choices)) choices <- neighbors[[idx]]
    selected[idx] <- sample(choices, 1L)
    used <- c(used, selected[idx])
  }
  selected
}

test_signature <- function(signature_name, genes, n_draws = 5000L, seed = 20260814L) {
  genes <- intersect(genes, rownames(diff_z))
  genes <- genes[apply(diff_z[genes, , drop = FALSE], 1, function(v) all(is.finite(v)))]
  neighbors <- build_neighbors(genes)
  observed_patient <- colMeans(diff_z[genes, , drop = FALSE])
  observed <- c(
    mean_patient_difference = mean(observed_patient),
    median_patient_difference = median(observed_patient),
    concordant_patients = sum(observed_patient > 0)
  )
  set.seed(seed)
  null <- matrix(NA_real_, nrow = n_draws, ncol = 3L)
  colnames(null) <- names(observed)
  for (b in seq_len(n_draws)) {
    selected <- draw_unique_set(neighbors)
    patient_diff <- colMeans(diff_z[selected, , drop = FALSE])
    null[b, ] <- c(mean(patient_diff), median(patient_diff), sum(patient_diff > 0))
  }
  summary <- data.frame(
    signature = signature_name,
    represented_genes = length(genes),
    patients = length(observed_patient),
    statistic = names(observed),
    observed = as.numeric(observed),
    null_median = apply(null, 2, median),
    null_q025 = apply(null, 2, quantile, probs = 0.025),
    null_q975 = apply(null, 2, quantile, probs = 0.975),
    empirical_p_upper = vapply(seq_along(observed), function(j) (1 + sum(null[, j] >= observed[j])) / (n_draws + 1), numeric(1)),
    percentile = vapply(seq_along(observed), function(j) mean(null[, j] <= observed[j]), numeric(1))
  )
  null_df <- data.frame(signature = signature_name, draw = seq_len(n_draws), null)
  patient_df <- data.frame(signature = signature_name, Patient = names(observed_patient), difference = as.numeric(observed_patient))
  list(summary = summary, null = null_df, patient = patient_df, genes = genes)
}

sen <- test_signature("SenMayo_nonoverlap", senmayo_nonoverlap, seed = 20260814L)
sasp <- test_signature("SASP25", sasp25, seed = 20260815L)
write.csv(rbind(sen$summary, sasp$summary), file.path(out_dir, "signature_matched_null_summary.csv"), row.names = FALSE)
write.csv(rbind(sen$null, sasp$null), file.path(out_dir, "signature_matched_null_draws.csv"), row.names = FALSE)
write.csv(rbind(sen$patient, sasp$patient), file.path(out_dir, "signature_matched_null_patient_differences.csv"), row.names = FALSE)
writeLines(c(
  paste0("SenMayo represented genes: ", length(sen$genes)),
  paste0("SASP represented genes: ", length(sasp$genes)),
  paste0("Candidate pool: ", length(candidate_genes)),
  paste0("Draws per signature: ", nrow(sen$null))
), file.path(out_dir, "signature_matched_null_run_info.txt"))

cat("Completed expression-matched signature specificity tests.\n")
