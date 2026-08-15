options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(lme4)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript recompute_single_cell_senescence.R <GSE132465_seurat.rds> <senmayo_genes.txt> <output_directory>")
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
senmayo_overlap <- intersect(senmayo_source, fap13)
senmayo_nonoverlap <- setdiff(senmayo_source, fap13)
markers <- c("CDKN2A", "CDKN2B", "CDKN1A", "LMNB1", "MKI67", "IL6", "CXCL8")
ligands <- c("GDF15", "CCL2", "TGFB1", "MMP3", "MMP1", "CXCL8", "IL6")
receptors <- c("IL6R", "IL6ST", "CXCR1", "CXCR2", "TGFBR1", "TGFBR2", "CCR2", "EGFR")

obj <- readRDS(input_rds)
meta <- obj@meta.data
stopifnot(nrow(meta) == ncol(obj), identical(rownames(meta), colnames(obj)))
rna_counts <- LayerData(obj[["RNA"]], layer = "counts")
rna_data <- LayerData(obj[["RNA"]], layer = "data")

fib_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")
is_fib <- meta$Cell_subtype %in% fib_subtypes
is_epi <- meta$Cell_type == "Epithelial cells"
is_t <- meta$Cell_type == "T cells"
is_myeloid <- meta$Cell_type == "Myeloids"
is_tumor <- meta$Class == "Tumor"

cell_counts <- data.frame(
  set = c("all_cells", "all_fibroblast_lineage", "tumor_fibroblast_lineage", "normal_fibroblast_lineage", "all_epithelial", "tumor_epithelial", "normal_epithelial", "all_FAP_detected", "all_fib_FAP_detected", "tumor_fib_FAP_detected", "normal_fib_FAP_detected"),
  n = c(
    nrow(meta), sum(is_fib), sum(is_fib & is_tumor), sum(is_fib & !is_tumor),
    sum(is_epi), sum(is_epi & is_tumor), sum(is_epi & !is_tumor),
    sum(rna_counts["FAP", ] > 0), sum(rna_counts["FAP", is_fib] > 0),
    sum(rna_counts["FAP", is_fib & is_tumor] > 0), sum(rna_counts["FAP", is_fib & !is_tumor] > 0)
  )
)
write.csv(cell_counts, file.path(out_dir, "cell_counts.csv"), row.names = FALSE)

gene_availability <- data.frame(
  set = c("SenMayo_source", "SenMayo_FAP13_overlap", "SenMayo_nonoverlap_source", "SenMayo_nonoverlap_represented", "SASP25_source", "SASP25_represented", "markers_represented"),
  n = c(length(senmayo_source), length(senmayo_overlap), length(senmayo_nonoverlap), sum(senmayo_nonoverlap %in% rownames(rna_data)), length(sasp25), sum(sasp25 %in% rownames(rna_data)), sum(markers %in% rownames(rna_data))),
  genes = c(
    paste(senmayo_source, collapse = ";"), paste(senmayo_overlap, collapse = ";"), paste(senmayo_nonoverlap, collapse = ";"),
    paste(intersect(senmayo_nonoverlap, rownames(rna_data)), collapse = ";"), paste(sasp25, collapse = ";"),
    paste(intersect(sasp25, rownames(rna_data)), collapse = ";"), paste(intersect(markers, rownames(rna_data)), collapse = ";")
  )
)
write.csv(gene_availability, file.path(out_dir, "gene_availability.csv"), row.names = FALSE)

score_mean <- function(mat, genes, cells) {
  g <- intersect(genes, rownames(mat))
  if (length(g) == 0L) stop("No represented genes")
  as.numeric(Matrix::colMeans(mat[g, cells, drop = FALSE]))
}

score_gene_z <- function(mat, genes, cells) {
  g <- intersect(genes, rownames(mat))
  x <- as.matrix(mat[g, cells, drop = FALSE])
  mu <- rowMeans(x)
  sig <- apply(x, 1, sd)
  keep <- is.finite(sig) & sig > 0
  z <- sweep(x[keep, , drop = FALSE], 1, mu[keep], "-")
  z <- sweep(z, 1, sig[keep], "/")
  colMeans(z)
}

safe_z <- function(x) {
  sx <- sd(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / sx)
}

make_df <- function(cells, compartment = NULL, score_scope = c("current", "global")) {
  score_scope <- match.arg(score_scope)
  d <- data.frame(
    cell = cells,
    Patient = meta[cells, "Patient"],
    Class = meta[cells, "Class"],
    Sample = meta[cells, "Sample"],
    Cell_type = meta[cells, "Cell_type"],
    Cell_subtype = meta[cells, "Cell_subtype"],
    nCount_RNA = meta[cells, "nCount_RNA"],
    stringsAsFactors = FALSE
  )
  d$log_nCount_RNA <- log1p(d$nCount_RNA)
  if (!is.null(compartment)) d$compartment <- compartment
  d$FAP_detected <- as.numeric(rna_counts["FAP", cells] > 0)
  d$FAP_status <- factor(ifelse(d$FAP_detected == 1, "FAP+", "FAP-"), levels = c("FAP-", "FAP+"))
  d$FAP_expr <- as.numeric(rna_data["FAP", cells])
  d$MKI67 <- as.numeric(rna_data["MKI67", cells])
  d$MKI67_z <- safe_z(d$MKI67)
  d$SenMayo_mean <- score_mean(rna_data, senmayo_nonoverlap, cells)
  d$SASP_mean <- score_mean(rna_data, sasp25, cells)
  d$matrix4_mean <- score_mean(rna_data, matrix4, cells)
  d$SenMayo_zmean <- score_gene_z(rna_data, senmayo_nonoverlap, cells)
  d$SASP_zmean <- score_gene_z(rna_data, sasp25, cells)
  d$matrix4_zmean <- score_gene_z(rna_data, matrix4, cells)
  d$SenMayo_mean_z <- safe_z(d$SenMayo_mean)
  d$SASP_mean_z <- safe_z(d$SASP_mean)
  d$matrix4_mean_z <- safe_z(d$matrix4_mean)
  d
}

model_row <- function(fit, model_id, subset_id, outcome) {
  ct <- as.data.frame(summary(fit)$coefficients)
  ct$term <- rownames(ct)
  rownames(ct) <- NULL
  names(ct)[1:5] <- c("estimate", "std_error", "df", "statistic", "p_value")
  ct$ci_low <- ct$estimate - 1.96 * ct$std_error
  ct$ci_high <- ct$estimate + 1.96 * ct$std_error
  vc <- as.data.frame(VarCorr(fit))
  random_var <- vc$vcov[vc$grp == "Patient"][1]
  residual_var <- vc$vcov[vc$grp == "Residual"][1]
  ct$icc <- random_var / (random_var + residual_var)
  ct$model_id <- model_id
  ct$subset <- subset_id
  ct$outcome <- outcome
  ct[, c("model_id", "subset", "outcome", "term", "estimate", "std_error", "df", "statistic", "p_value", "ci_low", "ci_high", "icc")]
}

fit_lmer <- function(formula, df, model_id, subset_id, outcome) {
  fit <- lmerTest::lmer(formula, data = df, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
  model_row(fit, model_id, subset_id, outcome)
}

paired_summary <- function(df, value, group, id = "Patient", group_a, group_b, analysis_id) {
  ag <- aggregate(df[[value]], list(id = df[[id]], group = df[[group]]), mean, na.rm = TRUE)
  names(ag)[3] <- "value"
  wa <- ag[ag$group == group_a, c("id", "value")]
  wb <- ag[ag$group == group_b, c("id", "value")]
  names(wa)[2] <- "a"
  names(wb)[2] <- "b"
  w <- merge(wa, wb, by = "id")
  w$diff_a_minus_b <- w$a - w$b
  wt <- suppressWarnings(wilcox.test(w$a, w$b, paired = TRUE, exact = FALSE, conf.int = FALSE))
  s <- data.frame(
    analysis_id = analysis_id,
    value = value,
    group_a = group_a,
    group_b = group_b,
    n_pairs = nrow(w),
    mean_a = mean(w$a),
    mean_b = mean(w$b),
    difference_of_means = mean(w$a) - mean(w$b),
    median_a = median(w$a),
    median_b = median(w$b),
    difference_of_medians = median(w$a) - median(w$b),
    median_paired_difference = median(w$diff_a_minus_b),
    positive = sum(w$diff_a_minus_b > 0),
    negative = sum(w$diff_a_minus_b < 0),
    zero = sum(w$diff_a_minus_b == 0),
    wilcoxon_W = unname(wt$statistic),
    p_value = wt$p.value
  )
  list(summary = s, patient = w)
}

all_comp_cells <- colnames(obj)[is_fib | is_epi]
all_comp <- make_df(
  all_comp_cells,
  ifelse(meta[all_comp_cells, "Cell_subtype"] %in% fib_subtypes, "Fibroblast", "Epithelial")
)
all_comp$compartment <- factor(all_comp$compartment, levels = c("Epithelial", "Fibroblast"))
all_fib <- all_comp[all_comp$compartment == "Fibroblast", ]

tumor_comp_cells <- colnames(obj)[is_tumor & (is_fib | is_epi)]
tumor_comp <- make_df(
  tumor_comp_cells,
  ifelse(meta[tumor_comp_cells, "Cell_subtype"] %in% fib_subtypes, "Fibroblast", "Epithelial")
)
tumor_comp$compartment <- factor(tumor_comp$compartment, levels = c("Epithelial", "Fibroblast"))
tumor_fib <- tumor_comp[tumor_comp$compartment == "Fibroblast", ]

normal_comp_cells <- colnames(obj)[!is_tumor & (is_fib | is_epi)]
normal_comp <- make_df(
  normal_comp_cells,
  ifelse(meta[normal_comp_cells, "Cell_subtype"] %in% fib_subtypes, "Fibroblast", "Epithelial")
)
normal_comp$compartment <- factor(normal_comp$compartment, levels = c("Epithelial", "Fibroblast"))
normal_fib <- normal_comp[normal_comp$compartment == "Fibroblast", ]

write.csv(all_comp, file.path(out_dir, "cell_scores_all_fibroblast_epithelial.csv"), row.names = FALSE)
write.csv(tumor_comp, file.path(out_dir, "cell_scores_tumor_fibroblast_epithelial.csv"), row.names = FALSE)

model_results <- rbind(
  fit_lmer(SenMayo_mean ~ FAP_status + MKI67 + (1 | Patient), all_fib, "B1_candidate_SenMayo", "all_fibroblast_including_normal", "SenMayo_mean"),
  fit_lmer(SASP_mean ~ FAP_status + MKI67 + (1 | Patient), all_fib, "B1_candidate_SASP", "all_fibroblast_including_normal", "SASP_mean"),
  fit_lmer(SenMayo_mean ~ FAP_status + MKI67 + Class + (1 | Patient), all_fib, "B1_class_adjusted_SenMayo", "all_fibroblast_including_normal", "SenMayo_mean"),
  fit_lmer(SenMayo_mean ~ FAP_status + MKI67 + (1 | Patient), tumor_fib, "B1_tumor_SenMayo", "tumor_fibroblast", "SenMayo_mean"),
  fit_lmer(SASP_mean ~ FAP_status + MKI67 + (1 | Patient), tumor_fib, "B1_tumor_SASP", "tumor_fibroblast", "SASP_mean"),
  fit_lmer(SenMayo_zmean ~ FAP_status + MKI67_z + (1 | Patient), tumor_fib, "B1_tumor_SenMayo_gene_z", "tumor_fibroblast", "SenMayo_zmean"),
  fit_lmer(SenMayo_zmean ~ FAP_status + MKI67_z + log_nCount_RNA + Cell_subtype + (1 | Patient), tumor_fib, "B1_tumor_SenMayo_depth_subtype", "tumor_fibroblast", "SenMayo_zmean"),
  fit_lmer(SASP_zmean ~ FAP_status + MKI67_z + log_nCount_RNA + Cell_subtype + (1 | Patient), tumor_fib, "B1_tumor_SASP_depth_subtype", "tumor_fibroblast", "SASP_zmean"),
  fit_lmer(SenMayo_mean ~ compartment + (1 | Patient), all_comp, "B2_candidate_unadjusted", "all_cells_including_normal", "SenMayo_mean"),
  fit_lmer(SenMayo_mean ~ compartment + MKI67 + (1 | Patient), all_comp, "B2_candidate_MKI67_adjusted", "all_cells_including_normal", "SenMayo_mean"),
  fit_lmer(SASP_mean ~ compartment + MKI67 + (1 | Patient), all_comp, "B2_candidate_SASP_adjusted", "all_cells_including_normal", "SASP_mean"),
  fit_lmer(MKI67 ~ compartment + (1 | Patient), all_comp, "B2_candidate_MKI67", "all_cells_including_normal", "MKI67"),
  fit_lmer(SenMayo_mean ~ compartment + MKI67 + Class + (1 | Patient), all_comp, "B2_class_adjusted", "all_cells_including_normal", "SenMayo_mean"),
  fit_lmer(SenMayo_mean ~ compartment + (1 | Patient), tumor_comp, "B2_tumor_unadjusted", "tumor_cells", "SenMayo_mean"),
  fit_lmer(SenMayo_mean ~ compartment + MKI67 + (1 | Patient), tumor_comp, "B2_tumor_MKI67_adjusted", "tumor_cells", "SenMayo_mean"),
  fit_lmer(SenMayo_zmean ~ compartment + MKI67_z + log_nCount_RNA + (1 | Patient), tumor_comp, "B2_tumor_depth_adjusted", "tumor_cells", "SenMayo_zmean"),
  fit_lmer(SASP_mean ~ compartment + MKI67 + (1 | Patient), tumor_comp, "B2_tumor_SASP_adjusted", "tumor_cells", "SASP_mean"),
  fit_lmer(MKI67 ~ compartment + (1 | Patient), tumor_comp, "B2_tumor_MKI67", "tumor_cells", "MKI67"),
  fit_lmer(matrix4_mean ~ SenMayo_mean + MKI67 + (1 | Patient), tumor_fib, "matrix_senescence_unadjusted_for_FAP", "tumor_fibroblast", "matrix4_mean"),
  fit_lmer(matrix4_mean ~ SenMayo_mean + FAP_expr + MKI67 + (1 | Patient), tumor_fib, "matrix_senescence_direct", "tumor_fibroblast", "matrix4_mean"),
  fit_lmer(matrix4_mean ~ SenMayo_mean + FAP_status + MKI67 + (1 | Patient), tumor_fib, "matrix_senescence_FAPstatus_adjusted", "tumor_fibroblast", "matrix4_mean"),
  fit_lmer(SenMayo_mean ~ matrix4_mean + FAP_status + MKI67 + (1 | Patient), tumor_fib, "senescence_matrix_FAPstatus_adjusted", "tumor_fibroblast", "SenMayo_mean"),
  fit_lmer(matrix4_mean ~ FAP_status + MKI67 + (1 | Patient), tumor_fib, "matrix_FAP_status", "tumor_fibroblast", "matrix4_mean")
)
write.csv(model_results, file.path(out_dir, "mixed_model_results.csv"), row.names = FALSE)

paired_objects <- list(
  paired_summary(all_comp, "SenMayo_mean", "compartment", group_a = "Fibroblast", group_b = "Epithelial", analysis_id = "all_compartment_SenMayo"),
  paired_summary(all_comp, "SASP_mean", "compartment", group_a = "Fibroblast", group_b = "Epithelial", analysis_id = "all_compartment_SASP"),
  paired_summary(all_comp, "MKI67", "compartment", group_a = "Fibroblast", group_b = "Epithelial", analysis_id = "all_compartment_MKI67"),
  paired_summary(tumor_comp, "SenMayo_mean", "compartment", group_a = "Fibroblast", group_b = "Epithelial", analysis_id = "tumor_compartment_SenMayo"),
  paired_summary(tumor_comp, "SASP_mean", "compartment", group_a = "Fibroblast", group_b = "Epithelial", analysis_id = "tumor_compartment_SASP"),
  paired_summary(tumor_comp, "MKI67", "compartment", group_a = "Fibroblast", group_b = "Epithelial", analysis_id = "tumor_compartment_MKI67"),
  paired_summary(tumor_fib, "SenMayo_mean", "FAP_status", group_a = "FAP+", group_b = "FAP-", analysis_id = "tumor_fib_FAP_SenMayo"),
  paired_summary(tumor_fib, "SASP_mean", "FAP_status", group_a = "FAP+", group_b = "FAP-", analysis_id = "tumor_fib_FAP_SASP"),
  paired_summary(tumor_fib, "matrix4_mean", "FAP_status", group_a = "FAP+", group_b = "FAP-", analysis_id = "tumor_fib_FAP_matrix4")
)
paired_summaries <- do.call(rbind, lapply(paired_objects, `[[`, "summary"))
write.csv(paired_summaries, file.path(out_dir, "patient_paired_summaries.csv"), row.names = FALSE)
for (x in paired_objects) write.csv(x$patient, file.path(out_dir, paste0(x$summary$analysis_id, "_patient_values.csv")), row.names = FALSE)

marker_detection <- function(df, subset_id) {
  cells <- df$cell
  out <- list()
  k <- 0L
  for (g in intersect(markers, rownames(rna_counts))) {
    det <- as.numeric(rna_counts[g, cells] > 0)
    expr <- as.numeric(rna_data[g, cells])
    for (st in c("FAP-", "FAP+")) {
      ii <- df$FAP_status == st
      k <- k + 1L
      out[[k]] <- data.frame(subset = subset_id, gene = g, FAP_status = st, n_cells = sum(ii), detection_rate = mean(det[ii]), mean_log_expression = mean(expr[ii]), median_log_expression = median(expr[ii]))
    }
  }
  do.call(rbind, out)
}
marker_results <- rbind(marker_detection(all_fib, "all_fibroblast_including_normal"), marker_detection(tumor_fib, "tumor_fibroblast"), marker_detection(normal_fib, "normal_fibroblast"))
write.csv(marker_results, file.path(out_dir, "marker_detection_by_FAP_status.csv"), row.names = FALSE)

patient_marker <- list()
kk <- 0L
for (g in intersect(markers, rownames(rna_counts))) {
  det <- as.numeric(rna_counts[g, tumor_fib$cell] > 0)
  tmp <- data.frame(Patient = tumor_fib$Patient, FAP_status = tumor_fib$FAP_status, detection = det)
  ag <- aggregate(detection ~ Patient + FAP_status, tmp, mean)
  wa <- ag[ag$FAP_status == "FAP+", c("Patient", "detection")]
  wb <- ag[ag$FAP_status == "FAP-", c("Patient", "detection")]
  names(wa)[2] <- "FAP_plus"
  names(wb)[2] <- "FAP_minus"
  w <- merge(wa, wb, by = "Patient")
  w$difference <- w$FAP_plus - w$FAP_minus
  wt <- suppressWarnings(wilcox.test(w$FAP_plus, w$FAP_minus, paired = TRUE, exact = FALSE))
  kk <- kk + 1L
  patient_marker[[kk]] <- data.frame(gene = g, n_pairs = nrow(w), median_difference = median(w$difference), positive = sum(w$difference > 0), negative = sum(w$difference < 0), wilcoxon_W = unname(wt$statistic), p_value = wt$p.value)
}
patient_marker <- do.call(rbind, patient_marker)
patient_marker$q_value_BH <- p.adjust(patient_marker$p_value, method = "BH")
write.csv(patient_marker, file.path(out_dir, "patient_paired_marker_detection.csv"), row.names = FALSE)

patient_cor <- function(df, subset_id) {
  ag <- aggregate(cbind(FAP_expr, SenMayo_mean, SASP_mean, matrix4_mean, MKI67) ~ Patient, df, mean)
  pairs <- list(c("FAP_expr", "SenMayo_mean"), c("FAP_expr", "SASP_mean"), c("FAP_expr", "matrix4_mean"), c("SenMayo_mean", "matrix4_mean"))
  out <- do.call(rbind, lapply(pairs, function(p) {
    ct <- suppressWarnings(cor.test(ag[[p[1]]], ag[[p[2]]], method = "spearman", exact = FALSE))
    data.frame(subset = subset_id, x = p[1], y = p[2], n = nrow(ag), rho = unname(ct$estimate), p_value = ct$p.value)
  }))
  list(summary = out, patient = ag)
}
pc_all <- patient_cor(all_fib, "all_fibroblast_including_normal")
pc_tumor <- patient_cor(tumor_fib, "tumor_fibroblast")
write.csv(rbind(pc_all$summary, pc_tumor$summary), file.path(out_dir, "patient_level_correlations.csv"), row.names = FALSE)
write.csv(pc_tumor$patient, file.path(out_dir, "tumor_fibroblast_patient_scores.csv"), row.names = FALSE)

# Exploratory tumour-only FAP+ fibroblast ligand x immune receptor screen.
tumor_fap_cells <- tumor_fib$cell[tumor_fib$FAP_status == "FAP+"]
ligand_mat <- lapply(intersect(ligands, rownames(rna_data)), function(g) {
  aggregate(as.numeric(rna_data[g, tumor_fap_cells]), list(Patient = meta[tumor_fap_cells, "Patient"]), mean)
})
names(ligand_mat) <- intersect(ligands, rownames(rna_data))
immune_groups <- list(T_cells = colnames(obj)[is_tumor & is_t], Myeloid = colnames(obj)[is_tumor & is_myeloid])

screen <- list()
k <- 0L
for (lg in names(ligand_mat)) {
  ldat <- ligand_mat[[lg]]
  names(ldat)[2] <- "ligand_value"
  for (grp in names(immune_groups)) {
    cells <- immune_groups[[grp]]
    for (rc in intersect(receptors, rownames(rna_data))) {
      rdat <- aggregate(as.numeric(rna_data[rc, cells]), list(Patient = meta[cells, "Patient"]), mean)
      names(rdat)[2] <- "receptor_value"
      m <- merge(ldat, rdat, by = "Patient")
      ct <- suppressWarnings(cor.test(m$ligand_value, m$receptor_value, method = "spearman", exact = FALSE))
      k <- k + 1L
      screen[[k]] <- data.frame(ligand = lg, receptor = rc, immune_compartment = grp, n = nrow(m), rho = unname(ct$estimate), p_value = ct$p.value)
    }
  }
}
screen <- do.call(rbind, screen)
screen$q_value_BH <- p.adjust(screen$p_value, method = "BH")
write.csv(screen, file.path(out_dir, "tumor_only_SASP_immune_screen.csv"), row.names = FALSE)

# Restrict interpretation to biologically cognate ligand-receptor pairs. The
# all-by-all matrix above is retained only as a transparent exploratory ledger.
cognate_pairs <- data.frame(
  ligand = c("IL6", "IL6", "CXCL8", "CXCL8", "CCL2", "TGFB1", "TGFB1"),
  receptor = c("IL6R", "IL6ST", "CXCR1", "CXCR2", "CCR2", "TGFBR1", "TGFBR2")
)
cognate_screen <- merge(screen[, c("ligand", "receptor", "immune_compartment", "n", "rho", "p_value")], cognate_pairs, by = c("ligand", "receptor"))
cognate_screen$q_value_BH <- p.adjust(cognate_screen$p_value, method = "BH")
write.csv(cognate_screen, file.path(out_dir, "tumor_only_cognate_ligand_receptor_screen.csv"), row.names = FALSE)

sink(file.path(out_dir, "R_sessionInfo.txt"))
cat("Input RDS:", input_rds, "\n")
cat("RDS dimensions:", paste(dim(obj), collapse = " x "), "\n")
cat("SenMayo source genes:", length(senmayo_source), "\n")
cat("SenMayo overlaps removed:", paste(senmayo_overlap, collapse = ", "), "\n")
cat("SenMayo non-overlap represented:", sum(senmayo_nonoverlap %in% rownames(rna_data)), "\n")
cat("SASP25 represented:", sum(sasp25 %in% rownames(rna_data)), "\n")
print(sessionInfo())
sink()

cat("Completed single-cell senescence reanalysis. Outputs:", normalizePath(out_dir), "\n")
