options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
})

root <- normalizePath(getwd(), winslash = "/")
default_source <- file.path(root, "data", "raw", "GSE132465")
source_root <- Sys.getenv("GSE132465_SOURCE_ROOT", unset = default_source)
annotation_file <- file.path(
  source_root, "data", "GSE132465",
  "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
)
raw_umi_file <- file.path(
  source_root, "data", "GSE132465",
  "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
)
seurat_file <- file.path(source_root, "results", "GSE132465_seurat.rds")
out_dir <- file.path(root, "results", "analysis", "GSE132465")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.exists(c(annotation_file, raw_umi_file, seurat_file))))

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
receptor2 <- c("SDC4", "CD44")
target_genes <- unique(c(FAP13, matrix4, receptor2))
fibroblast_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2",
                         "Stromal 3")

bootstrap_spearman <- function(x, y, seed, reps = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  set.seed(seed)
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  estimates <- replicate(reps, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  estimates <- estimates[is.finite(estimates)]
  test <- suppressWarnings(cor.test(x, y, method = "spearman",
                                    exact = FALSE))
  c(n = n, rho = observed,
    ci_low = unname(quantile(estimates, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(estimates, 0.975, type = 6, na.rm = TRUE)),
    p_value = test$p.value)
}

score_zmean <- function(expression, genes, patients) {
  stopifnot(all(genes %in% rownames(expression)))
  z <- t(scale(t(expression[genes, patients, drop = FALSE])))
  colMeans(z, na.rm = FALSE)
}

make_group_indicator <- function(group, levels) {
  keep <- !is.na(group)
  sparseMatrix(
    i = which(keep),
    j = match(group[keep], levels),
    x = 1,
    dims = c(length(group), length(levels)),
    dimnames = list(NULL, levels)
  )
}

normalise_pseudobulk <- function(counts) {
  library_size <- Matrix::colSums(counts)
  stopifnot(all(library_size > 0))
  cpm <- t(t(as.matrix(counts)) / library_size) * 1e6
  log2(cpm + 0.25)
}

# The raw UMI header is checked independently against the GEO annotation.
annotation <- fread(annotation_file)
raw_connection <- gzfile(raw_umi_file, open = "rt")
raw_header <- readLines(raw_connection, n = 1L, warn = FALSE)
close(raw_connection)
raw_cell_ids <- strsplit(raw_header, "\t", fixed = TRUE)[[1]][-1]
stopifnot(nrow(annotation) == 63689L)
stopifnot(identical(raw_cell_ids, annotation$Index))

object <- readRDS(seurat_file)
counts <- LayerData(object[["RNA"]], layer = "counts")
metadata <- object[[]]
stopifnot(inherits(counts, "sparseMatrix"))
stopifnot(identical(colnames(counts), annotation$Index))
stopifnot(identical(colnames(counts), rownames(metadata)))
stopifnot(all(target_genes %in% rownames(counts)))

# Cross-check annotations embedded in the frozen object against the GEO file.
annotation_columns <- c("Patient", "Class", "Sample", "Cell_type",
                        "Cell_subtype")
annotation_equal <- vapply(annotation_columns, function(variable) {
  identical(as.character(metadata[[variable]]),
            as.character(annotation[[variable]]))
}, logical(1))
stopifnot(all(annotation_equal))

is_fibroblast <- metadata$Cell_type == "Stromal cells" &
  metadata$Cell_subtype %in% fibroblast_subtypes
is_epithelial <- metadata$Cell_type == "Epithelial cells"
compartment <- ifelse(is_fibroblast, "fibroblast_lineage",
                      ifelse(is_epithelial, "epithelial", NA_character_))
group <- ifelse(
  is.na(compartment), NA_character_,
  paste(metadata$Class, compartment, metadata$Patient, sep = "::")
)
group_levels <- unique(group[!is.na(group)])
indicator <- make_group_indicator(group, group_levels)
pseudobulk_counts <- counts %*% indicator

group_parts <- tstrsplit(group_levels, "::", fixed = TRUE)
group_manifest <- data.frame(
  group = group_levels,
  class = group_parts[[1]],
  compartment = group_parts[[2]],
  patient = group_parts[[3]],
  cells = as.integer(Matrix::colSums(indicator)),
  library_size = as.numeric(Matrix::colSums(pseudobulk_counts)),
  stringsAsFactors = FALSE
)
write.csv(group_manifest,
          file.path(out_dir, "GSE132465_pseudobulk_group_manifest.csv"),
          row.names = FALSE)

normalised <- matrix(NA_real_, nrow = nrow(pseudobulk_counts),
                     ncol = ncol(pseudobulk_counts),
                     dimnames = dimnames(pseudobulk_counts))
for (this_class in unique(group_manifest$class)) {
  for (this_compartment in unique(group_manifest$compartment)) {
    selected <- group_manifest$class == this_class &
      group_manifest$compartment == this_compartment
    if (sum(selected) >= 2L) {
      normalised[, selected] <- normalise_pseudobulk(
        pseudobulk_counts[, selected, drop = FALSE]
      )
    }
  }
}

tumour_manifest <- group_manifest[group_manifest$class == "Tumor", ]
cell_counts <- reshape(
  tumour_manifest[, c("patient", "compartment", "cells")],
  idvar = "patient", timevar = "compartment", direction = "wide"
)
names(cell_counts) <- sub("^cells\\.", "n_", names(cell_counts))
eligible <- cell_counts$patient[
  cell_counts$n_fibroblast_lineage >= 20L &
    cell_counts$n_epithelial >= 20L
]
stopifnot(length(eligible) == 15L)

group_name <- function(class, compartment, patient) {
  paste(class, compartment, patient, sep = "::")
}
fib_groups <- group_name("Tumor", "fibroblast_lineage", eligible)
epi_groups <- group_name("Tumor", "epithelial", eligible)
fib_expression <- normalised[, fib_groups, drop = FALSE]
epi_expression <- normalised[, epi_groups, drop = FALSE]
colnames(fib_expression) <- eligible
colnames(epi_expression) <- eligible

patient_scores <- data.frame(
  patient = eligible,
  n_fib = cell_counts$n_fibroblast_lineage[match(eligible,
                                                  cell_counts$patient)],
  n_epi = cell_counts$n_epithelial[match(eligible, cell_counts$patient)],
  fib_FAP_logCPM = as.numeric(fib_expression["FAP", eligible]),
  fib_matrix4_zmean = score_zmean(fib_expression, matrix4, eligible),
  fib_receptor2_zmean = score_zmean(fib_expression, receptor2, eligible),
  epi_receptor2_zmean = score_zmean(epi_expression, receptor2, eligible),
  stringsAsFactors = FALSE
)
write.csv(patient_scores,
          file.path(out_dir, "GSE132465_raw_UMI_patient_scores.csv"),
          row.names = FALSE)

comparisons <- list(
  "Fibroblast-lineage FAP vs fibroblast-lineage matrix4" =
    c("fib_FAP_logCPM", "fib_matrix4_zmean"),
  "Fibroblast-lineage FAP vs epithelial receptor2" =
    c("fib_FAP_logCPM", "epi_receptor2_zmean"),
  "Fibroblast-lineage matrix4 vs epithelial receptor2" =
    c("fib_matrix4_zmean", "epi_receptor2_zmean"),
  "Fibroblast-lineage FAP vs fibroblast-lineage receptor2" =
    c("fib_FAP_logCPM", "fib_receptor2_zmean")
)
correlations <- do.call(rbind, lapply(seq_along(comparisons), function(index) {
  label <- names(comparisons)[index]
  variables <- comparisons[[label]]
  bootstrap_seed <- 20260805L + index - 1L
  data.frame(
    comparison = label,
    bootstrap_seed = bootstrap_seed,
    t(bootstrap_spearman(patient_scores[[variables[1]]],
                         patient_scores[[variables[2]]],
                         seed = bootstrap_seed))
  )
}))
correlations$fdr_bh_four <- p.adjust(correlations$p_value, method = "BH")
rownames(correlations) <- NULL
write.csv(correlations,
          file.path(out_dir, "GSE132465_raw_UMI_correlations.csv"),
          row.names = FALSE)

# Paired tumour-normal analysis is exploratory and uses only patients with
# both classes and at least 20 fibroblast-lineage cells in each class.
paired_counts <- reshape(
  group_manifest[group_manifest$compartment == "fibroblast_lineage",
                 c("patient", "class", "cells")],
  idvar = "patient", timevar = "class", direction = "wide"
)
paired_patients <- paired_counts$patient[
  !is.na(paired_counts$cells.Tumor) & !is.na(paired_counts$cells.Normal) &
    paired_counts$cells.Tumor >= 20L & paired_counts$cells.Normal >= 20L
]
paired_rows <- do.call(rbind, lapply(paired_patients, function(patient) {
  tumour_group <- group_name("Tumor", "fibroblast_lineage", patient)
  normal_group <- group_name("Normal", "fibroblast_lineage", patient)
  data.frame(
    patient = patient,
    tumour_cells = group_manifest$cells[match(tumour_group,
                                               group_manifest$group)],
    normal_cells = group_manifest$cells[match(normal_group,
                                               group_manifest$group)],
    tumour_FAP_logCPM = normalised["FAP", tumour_group],
    normal_FAP_logCPM = normalised["FAP", normal_group],
    tumour_matrix4_mean_logCPM = mean(
      normalised[matrix4, tumour_group], na.rm = FALSE),
    normal_matrix4_mean_logCPM = mean(
      normalised[matrix4, normal_group], na.rm = FALSE)
  )
}))
write.csv(paired_rows,
          file.path(out_dir, "GSE132465_paired_tumour_normal_scores.csv"),
          row.names = FALSE)

paired_tests <- rbind(
  data.frame(
    feature = "FAP logCPM", n_pairs = nrow(paired_rows),
    median_tumour_minus_normal = median(
      paired_rows$tumour_FAP_logCPM - paired_rows$normal_FAP_logCPM),
    p_value = wilcox.test(paired_rows$tumour_FAP_logCPM,
                          paired_rows$normal_FAP_logCPM,
                          paired = TRUE, exact = FALSE)$p.value
  ),
  data.frame(
    feature = "matrix4 mean logCPM", n_pairs = nrow(paired_rows),
    median_tumour_minus_normal = median(
      paired_rows$tumour_matrix4_mean_logCPM -
        paired_rows$normal_matrix4_mean_logCPM),
    p_value = wilcox.test(paired_rows$tumour_matrix4_mean_logCPM,
                          paired_rows$normal_matrix4_mean_logCPM,
                          paired = TRUE, exact = FALSE)$p.value
  )
)
paired_tests$fdr_bh_two <- p.adjust(paired_tests$p_value, method = "BH")
write.csv(paired_tests,
          file.path(out_dir, "GSE132465_paired_tumour_normal_tests.csv"),
          row.names = FALSE)

source_manifest <- data.frame(
  file = c(annotation_file, raw_umi_file, seurat_file),
  bytes = as.numeric(file.info(c(annotation_file, raw_umi_file,
                                 seurat_file))$size),
  md5 = unname(tools::md5sum(c(annotation_file, raw_umi_file, seurat_file))),
  role = c("GEO cell annotation", "GEO raw UMI matrix",
           "Frozen sparse raw-count object used for aggregation"),
  stringsAsFactors = FALSE
)
write.csv(source_manifest,
          file.path(out_dir, "GSE132465_source_manifest.csv"),
          row.names = FALSE)

structure_audit <- data.frame(
  item = c("raw_UMI_header_cells", "annotation_cells", "Seurat_cells",
           "Seurat_genes", "patients", "tumour_cells", "normal_cells",
           "eligible_tumour_patients", "paired_fibroblast_patients",
           paste0("annotation_equal_", annotation_columns)),
  value = c(length(raw_cell_ids), nrow(annotation), ncol(counts), nrow(counts),
            uniqueN(metadata$Patient), sum(metadata$Class == "Tumor"),
            sum(metadata$Class == "Normal"), length(eligible),
            length(paired_patients), as.character(annotation_equal)),
  stringsAsFactors = FALSE
)
write.csv(structure_audit,
          file.path(out_dir, "GSE132465_structure_audit.csv"),
          row.names = FALSE)

log_file <- file.path(out_dir, "GSE132465_raw_pseudobulk.log")
sink(log_file)
cat("GSE132465 RAW-COUNT PSEUDOBULK AUDIT\n")
cat("Run UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n\n")
print(structure_audit, row.names = FALSE)
cat("\nPatient-level correlations:\n")
print(correlations, row.names = FALSE)
cat("\nPaired tumour-normal sensitivity:\n")
print(paired_tests, row.names = FALSE)
cat("\nSession information:\n")
print(sessionInfo())
sink()
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "GSE132465_raw_pseudobulk_sessionInfo.txt"))
