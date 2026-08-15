suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(lmerTest)
})

options(width = 180, digits = 8)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript independent_recompute_single_cell.R <GSE132465_seurat.rds> <senmayo_genes.txt>")
}
obj_path <- normalizePath(args[[1]], mustWork = TRUE)
sen_path <- normalizePath(args[[2]], mustWork = TRUE)

obj <- readRDS(obj_path)
meta <- obj[[]]
counts <- LayerData(obj, assay = "RNA", layer = "counts")
data <- LayerData(obj, assay = "RNA", layer = "data")

sen_source <- unique(trimws(readLines(sen_path)))
sen_source <- sen_source[nzchar(sen_source)]
fap_overlap <- c("CXCL12", "WNT2", "MMP2", "MMP9")
sen_no_overlap <- setdiff(sen_source, fap_overlap)
sen_genes <- intersect(sen_no_overlap, rownames(data))
sen_full_genes <- intersect(sen_source, rownames(data))
sasp_source <- c(
  "IL6", "CXCL8", "IL1A", "IL1B", "CCL2", "CCL5", "CXCL1", "CXCL2",
  "CXCL3", "CXCL10", "MMP1", "MMP3", "MMP9", "MMP10", "MMP13",
  "SERPINE1", "PLAU", "TIMP2", "VEGFA", "GDF15", "IGFBP3", "TNF",
  "CSF2", "HGF", "FAS"
)
sasp_genes <- intersect(sasp_source, rownames(data))

cat("OBJECT AND GENE COVERAGE\n")
cat("object dim:", nrow(obj), "genes x", ncol(obj), "cells\n")
cat("SenMayo source:", length(sen_source), "minus overlaps:", length(sen_no_overlap), "represented:", length(sen_genes), "\n")
cat("SenMayo full represented (including 4 FAP13 overlaps):", length(sen_full_genes), "\n")
cat("missing after overlap removal:", paste(setdiff(sen_no_overlap, rownames(data)), collapse = ","), "\n")
cat("SASP source:", length(sasp_source), "represented:", length(sasp_genes), "missing:", paste(setdiff(sasp_source, rownames(data)), collapse = ","), "\n\n")

fib_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")
is_fib <- meta$Cell_subtype %in% fib_subtypes
is_epi <- meta$Cell_type == "Epithelial cells"
all_target <- is_fib | is_epi

cat("CELL SET CROSS-TABULATIONS\n")
print(with(meta, table(Class, compartment = ifelse(is_fib, "Fibroblast", ifelse(is_epi, "Epithelial", "Other"))))[, c("Fibroblast", "Epithelial")])
cat("fibroblast subtype by class:\n")
print(with(meta[is_fib, ], table(Class, Cell_subtype)))

score_raw <- function(genes, cells) {
  Matrix::colMeans(data[genes, cells, drop = FALSE])
}

score_zmean <- function(genes, cells) {
  m <- as.matrix(data[genes, cells, drop = FALSE])
  sds <- apply(m, 1, sd)
  keep <- is.finite(sds) & sds > 0
  z <- t(scale(t(m[keep, , drop = FALSE]), center = TRUE, scale = TRUE))
  colMeans(z)
}

target_cells <- colnames(data)[all_target]
fib_cells <- colnames(data)[is_fib]
epi_cells <- colnames(data)[is_epi]

cell_df <- meta[all_target, c("Patient", "Class", "Cell_type", "Cell_subtype", "nCount_RNA", "nFeature_RNA"), drop = FALSE]
cell_df$cell <- rownames(cell_df)
cell_df$compartment <- ifelse(cell_df$Cell_subtype %in% fib_subtypes, "Fibroblast", "Epithelial")
cell_df$compartment <- factor(cell_df$compartment, levels = c("Epithelial", "Fibroblast"))
cell_df$SenMayo_rawmean <- score_raw(sen_genes, rownames(cell_df))
cell_df$SASP_rawmean <- score_raw(sasp_genes, rownames(cell_df))
cell_df$SenMayo_zmean <- score_zmean(sen_genes, rownames(cell_df))
cell_df$SenMayo_full_zmean <- score_zmean(sen_full_genes, rownames(cell_df))
cell_df$SASP_zmean <- score_zmean(sasp_genes, rownames(cell_df))
cell_df$FAP <- as.numeric(data["FAP", rownames(cell_df)])
cell_df$FAPpos <- cell_df$FAP > 0
cell_df$MKI67 <- as.numeric(data["MKI67", rownames(cell_df)])
cell_df$logUMI <- log1p(cell_df$nCount_RNA)
cell_df$logFeature <- log1p(cell_df$nFeature_RNA)

fib <- cell_df[cell_df$compartment == "Fibroblast", , drop = FALSE]
fib$SenMayo_fib_zmean <- score_zmean(sen_genes, fib$cell)
fib$SenMayo_fib_full_zmean <- score_zmean(sen_full_genes, fib$cell)
fib$SASP_fib_zmean <- score_zmean(sasp_genes, fib$cell)
fib_t_cells <- fib$cell[fib$Class == "Tumor"]
fib$SenMayo_tumor_zmean <- NA_real_
fib$SASP_tumor_zmean <- NA_real_
fib$SenMayo_tumor_zmean[fib$Class == "Tumor"] <- score_zmean(sen_genes, fib_t_cells)
fib$SASP_tumor_zmean[fib$Class == "Tumor"] <- score_zmean(sasp_genes, fib_t_cells)
cat("\nFAP STATUS CONFOUNDING\n")
print(with(fib, table(Class, FAPpos)))
cat("P(Tumor|FAP+) =", mean(fib$Class[fib$FAPpos] == "Tumor"), "; P(Tumor|FAP-) =", mean(fib$Class[!fib$FAPpos] == "Tumor"), "\n")
cat("depth summaries by FAP status (median):\n")
print(aggregate(cbind(nCount_RNA, nFeature_RNA, logUMI) ~ FAPpos, fib, median))
cat("Spearman score-depth correlations in all fibroblasts:\n")
for (v in c("SenMayo_rawmean", "SASP_rawmean", "SenMayo_zmean", "SASP_zmean")) {
  cat(v, "vs nCount", cor(fib[[v]], fib$nCount_RNA, method = "spearman"), "vs nFeature", cor(fib[[v]], fib$nFeature_RNA, method = "spearman"), "\n")
}

print_lmer <- function(label, fit, term) {
  co <- summary(fit)$coefficients
  cat(label, ": n=", nobs(fit), "; estimate=", co[term, "Estimate"], "; SE=", co[term, "Std. Error"], "; df=", co[term, "df"], "; t=", co[term, "t value"], "; P=", co[term, "Pr(>|t|)"], "\n", sep = "")
}

cat("\nB1 REPORTED PIPELINE REPRODUCTION (raw log-normalized gene means)\n")
cat("raw score medians by FAP status:\n")
print(aggregate(cbind(SenMayo_rawmean, SASP_rawmean) ~ FAPpos, fib, median))
cat("z-mean score medians by FAP status (Methods-stated alternative):\n")
print(aggregate(cbind(SenMayo_zmean, SASP_zmean) ~ FAPpos, fib, median))

cat("fibroblast-only z-mean score medians by FAP status (actual reported scale):\n")
print(aggregate(cbind(SenMayo_fib_zmean, SASP_fib_zmean) ~ FAPpos, fib, median))
cat("fibroblast-only FULL SenMayo (with 4 overlaps) medians:\n")
print(aggregate(SenMayo_fib_full_zmean ~ FAPpos, fib, median))

m_b1_sen <- lmer(SenMayo_rawmean ~ FAPpos + MKI67 + (1 | Patient), data = fib)
m_b1_sasp <- lmer(SASP_rawmean ~ FAPpos + MKI67 + (1 | Patient), data = fib)
print_lmer("SenMayo ~ FAPpos + MKI67 + (1|Patient)", m_b1_sen, "FAPposTRUE")
print_lmer("SASP ~ FAPpos + MKI67 + (1|Patient)", m_b1_sasp, "FAPposTRUE")
m_b1_sen_z <- lmer(SenMayo_fib_zmean ~ FAPpos + MKI67 + (1 | Patient), data = fib)
m_b1_sasp_z <- lmer(SASP_fib_zmean ~ FAPpos + MKI67 + (1 | Patient), data = fib)
print_lmer("REPORTED-SCALE SenMayo fib-zmean ~ FAPpos + MKI67", m_b1_sen_z, "FAPposTRUE")
print_lmer("REPORTED-SCALE SASP fib-zmean ~ FAPpos + MKI67", m_b1_sasp_z, "FAPposTRUE")
m_b1_sen_full_z <- lmer(SenMayo_fib_full_zmean ~ FAPpos + MKI67 + (1 | Patient), data = fib)
print_lmer("FULL SenMayo (4 overlaps retained) fib-zmean ~ FAPpos + MKI67", m_b1_sen_full_z, "FAPposTRUE")

cat("\nB1 SENSITIVITY MODELS\n")
m_b1_sen_class <- lmer(SenMayo_rawmean ~ FAPpos + MKI67 + Class + logUMI + (1 | Patient), data = fib)
m_b1_sasp_class <- lmer(SASP_rawmean ~ FAPpos + MKI67 + Class + logUMI + (1 | Patient), data = fib)
print_lmer("SenMayo adjusted for Class + logUMI", m_b1_sen_class, "FAPposTRUE")
print_lmer("SASP adjusted for Class + logUMI", m_b1_sasp_class, "FAPposTRUE")
m_b1_sen_z_class <- lmer(SenMayo_fib_zmean ~ FAPpos + MKI67 + Class + logUMI + (1 | Patient), data = fib)
m_b1_sasp_z_class <- lmer(SASP_fib_zmean ~ FAPpos + MKI67 + Class + logUMI + (1 | Patient), data = fib)
print_lmer("REPORTED-SCALE SenMayo adjusted Class + logUMI", m_b1_sen_z_class, "FAPposTRUE")
print_lmer("REPORTED-SCALE SASP adjusted Class + logUMI", m_b1_sasp_z_class, "FAPposTRUE")

fib_t <- fib[fib$Class == "Tumor", , drop = FALSE]
m_b1_sen_t <- lmer(SenMayo_rawmean ~ FAPpos + MKI67 + logUMI + (1 | Patient), data = fib_t)
m_b1_sasp_t <- lmer(SASP_rawmean ~ FAPpos + MKI67 + logUMI + (1 | Patient), data = fib_t)
print_lmer("Tumor-only SenMayo adjusted logUMI", m_b1_sen_t, "FAPposTRUE")
print_lmer("Tumor-only SASP adjusted logUMI", m_b1_sasp_t, "FAPposTRUE")
m_b1_sen_t_z <- lmer(SenMayo_tumor_zmean ~ FAPpos + MKI67 + logUMI + (1 | Patient), data = fib_t)
m_b1_sasp_t_z <- lmer(SASP_tumor_zmean ~ FAPpos + MKI67 + logUMI + (1 | Patient), data = fib_t)
print_lmer("Tumor-only RECOMPUTED-z SenMayo adjusted logUMI", m_b1_sen_t_z, "FAPposTRUE")
print_lmer("Tumor-only RECOMPUTED-z SASP adjusted logUMI", m_b1_sasp_t_z, "FAPposTRUE")
m_b1_sen_t_nf <- lmer(SenMayo_tumor_zmean ~ FAPpos + MKI67 + logFeature + (1 | Patient), data = fib_t)
m_b1_sasp_t_nf <- lmer(SASP_tumor_zmean ~ FAPpos + MKI67 + logFeature + (1 | Patient), data = fib_t)
print_lmer("Tumor-only RECOMPUTED-z SenMayo adjusted logFeature", m_b1_sen_t_nf, "FAPposTRUE")
print_lmer("Tumor-only RECOMPUTED-z SASP adjusted logFeature", m_b1_sasp_t_nf, "FAPposTRUE")

cat("\nTABLE 3 FULL-LIBRARY NORMALIZATION REPRODUCTION\n")
fap_t_z <- as.numeric(scale(fib_t$FAP))
matrix4_t_zmean <- score_zmean(c("COL1A1", "COL1A2", "COL3A1", "FN1"), fib_t$cell)
receptor2_t_zmean <- score_zmean(c("SDC4", "CD44"), fib_t$cell)
table3_df <- data.frame(Patient = fib_t$Patient, FAP_z = fap_t_z, matrix4_zmean = matrix4_t_zmean, receptor2_zmean = receptor2_t_zmean)
m_table3_matrix <- lmer(matrix4_zmean ~ FAP_z + (1 | Patient), data = table3_df)
m_table3_receptor <- lmer(receptor2_zmean ~ FAP_z + (1 | Patient), data = table3_df)
print_lmer("full-library matrix4 ~ FAP_z", m_table3_matrix, "FAP_z")
print_lmer("full-library receptor2 ~ FAP_z", m_table3_receptor, "FAP_z")

cat("random-slope sensitivity (may be singular):\n")
m_b1_sen_rs <- suppressWarnings(lmer(SenMayo_rawmean ~ FAPpos + MKI67 + Class + logUMI + (1 + FAPpos | Patient), data = fib, control = lmerControl(check.conv.singular = "ignore")))
m_b1_sasp_rs <- suppressWarnings(lmer(SASP_rawmean ~ FAPpos + MKI67 + Class + logUMI + (1 + FAPpos | Patient), data = fib, control = lmerControl(check.conv.singular = "ignore")))
print_lmer("SenMayo random-slope adjusted", m_b1_sen_rs, "FAPposTRUE")
print_lmer("SASP random-slope adjusted", m_b1_sasp_rs, "FAPposTRUE")
cat("singular SenMayo:", isSingular(m_b1_sen_rs), "SASP:", isSingular(m_b1_sasp_rs), "\n")

patient_status <- aggregate(cbind(SenMayo_rawmean, SASP_rawmean) ~ Patient + FAPpos, fib, mean)
patient_counts <- aggregate(cell ~ Patient + FAPpos, fib, length)
names(patient_counts)[3] <- "n_cells"
patient_status <- merge(patient_status, patient_counts, by = c("Patient", "FAPpos"))
sen_w <- reshape(patient_status[, c("Patient", "FAPpos", "SenMayo_rawmean", "n_cells")], idvar = "Patient", timevar = "FAPpos", direction = "wide")
sasp_w <- reshape(patient_status[, c("Patient", "FAPpos", "SASP_rawmean", "n_cells")], idvar = "Patient", timevar = "FAPpos", direction = "wide")
eligible_sen <- complete.cases(sen_w) & sen_w$n_cells.FALSE >= 5 & sen_w$n_cells.TRUE >= 5
eligible_sasp <- complete.cases(sasp_w) & sasp_w$n_cells.FALSE >= 5 & sasp_w$n_cells.TRUE >= 5
sen_diff <- sen_w$SenMayo_rawmean.TRUE[eligible_sen] - sen_w$SenMayo_rawmean.FALSE[eligible_sen]
sasp_diff <- sasp_w$SASP_rawmean.TRUE[eligible_sasp] - sasp_w$SASP_rawmean.FALSE[eligible_sasp]
cat("patient-paired FAP+/FAP- (>=5 cells/group) SenMayo n=", length(sen_diff), " median diff=", median(sen_diff), " positive=", sum(sen_diff > 0), " Wilcoxon P=", wilcox.test(sen_diff, exact = FALSE)$p.value, "\n", sep = "")
cat("patient-paired FAP+/FAP- (>=5 cells/group) SASP n=", length(sasp_diff), " median diff=", median(sasp_diff), " positive=", sum(sasp_diff > 0), " Wilcoxon P=", wilcox.test(sasp_diff, exact = FALSE)$p.value, "\n", sep = "")

patient_status_t <- aggregate(cbind(SenMayo_tumor_zmean, SASP_tumor_zmean) ~ Patient + FAPpos, fib_t, mean)
patient_counts_t <- aggregate(cell ~ Patient + FAPpos, fib_t, length)
names(patient_counts_t)[3] <- "n_cells"
patient_status_t <- merge(patient_status_t, patient_counts_t, by = c("Patient", "FAPpos"))
sen_w_t <- reshape(patient_status_t[, c("Patient", "FAPpos", "SenMayo_tumor_zmean", "n_cells")], idvar = "Patient", timevar = "FAPpos", direction = "wide")
sasp_w_t <- reshape(patient_status_t[, c("Patient", "FAPpos", "SASP_tumor_zmean", "n_cells")], idvar = "Patient", timevar = "FAPpos", direction = "wide")
eligible_sen_t <- complete.cases(sen_w_t) & sen_w_t$n_cells.FALSE >= 5 & sen_w_t$n_cells.TRUE >= 5
eligible_sasp_t <- complete.cases(sasp_w_t) & sasp_w_t$n_cells.FALSE >= 5 & sasp_w_t$n_cells.TRUE >= 5
sen_diff_t <- sen_w_t$SenMayo_tumor_zmean.TRUE[eligible_sen_t] - sen_w_t$SenMayo_tumor_zmean.FALSE[eligible_sen_t]
sasp_diff_t <- sasp_w_t$SASP_tumor_zmean.TRUE[eligible_sasp_t] - sasp_w_t$SASP_tumor_zmean.FALSE[eligible_sasp_t]
cat("TUMOR-ONLY patient-paired FAP+/FAP- (>=5/group) SenMayo n=", length(sen_diff_t), " median diff=", median(sen_diff_t), " positive=", sum(sen_diff_t > 0), " Wilcoxon P=", wilcox.test(sen_diff_t, exact = FALSE)$p.value, "\n", sep = "")
cat("TUMOR-ONLY patient-paired FAP+/FAP- (>=5/group) SASP n=", length(sasp_diff_t), " median diff=", median(sasp_diff_t), " positive=", sum(sasp_diff_t > 0), " Wilcoxon P=", wilcox.test(sasp_diff_t, exact = FALSE)$p.value, "\n", sep = "")

cat("\nMARKER DETECTION RATES\n")
for (g in c("CDKN2A", "CDKN2B", "LMNB1")) {
  fib[[paste0(g, "_det")]] <- as.numeric(counts[g, fib$cell]) > 0
  rates <- aggregate(fib[[paste0(g, "_det")]], list(FAPpos = fib$FAPpos), mean)
  cat(g, " rates:", paste(rates$FAPpos, round(rates$x, 5), collapse = "; "), "\n")
  gm <- glmer(as.formula(paste0(g, "_det ~ FAPpos + Class + logUMI + (1 | Patient)")), data = fib, family = binomial, control = glmerControl(optimizer = "bobyqa"))
  co <- summary(gm)$coefficients["FAPposTRUE", ]
  cat(g, " adjusted OR=", exp(co["Estimate"]), " z=", co["z value"], " P=", co["Pr(>|z|)"], "\n", sep = "")
}

summarize_compartment <- function(df, min_cells = 1) {
  means <- aggregate(cbind(SenMayo_rawmean, SASP_rawmean, SenMayo_zmean, SenMayo_full_zmean, SASP_zmean, MKI67) ~ Patient + compartment, df, mean)
  ns <- aggregate(cell ~ Patient + compartment, df, length)
  names(ns)[3] <- "n_cells"
  means <- merge(means, ns, by = c("Patient", "compartment"))
  wide <- reshape(means, idvar = "Patient", timevar = "compartment", direction = "wide")
  eligible <- complete.cases(wide) & wide$n_cells.Epithelial >= min_cells & wide$n_cells.Fibroblast >= min_cells
  out <- list(wide = wide[eligible, , drop = FALSE])
  for (v in c("SenMayo_rawmean", "SASP_rawmean", "SenMayo_zmean", "SenMayo_full_zmean", "SASP_zmean", "MKI67")) {
    out[[v]] <- out$wide[[paste0(v, ".Fibroblast")]] - out$wide[[paste0(v, ".Epithelial")]]
  }
  out
}

print_compartment <- function(label, obj) {
  cat(label, " n patients=", nrow(obj$wide), "\n", sep = "")
  for (v in c("SenMayo_rawmean", "SASP_rawmean", "SenMayo_zmean", "SenMayo_full_zmean", "SASP_zmean", "MKI67")) {
    d <- obj[[v]]
    wt <- wilcox.test(d, exact = TRUE)
    cat("  ", v, ": median paired diff=", median(d), "; positive=", sum(d > 0), "/", length(d), "; P=", wt$p.value, "\n", sep = "")
  }
}

boot_median_ci <- function(x, n_boot = 10000, seed = 20260814) {
  set.seed(seed)
  vals <- replicate(n_boot, median(sample(x, length(x), replace = TRUE)))
  unname(quantile(vals, c(0.025, 0.975)))
}

cat("\nB2 PATIENT-LEVEL COMPARTMENT CONTRASTS\n")
b2_all <- summarize_compartment(cell_df, min_cells = 1)
b2_tumor_all <- summarize_compartment(cell_df[cell_df$Class == "Tumor", ], min_cells = 1)
b2_tumor_20 <- summarize_compartment(cell_df[cell_df$Class == "Tumor", ], min_cells = 20)
print_compartment("All tumor+normal cells (reported population)", b2_all)
print_compartment("Tumor-only, no minimum", b2_tumor_all)
print_compartment("Tumor-only, >=20 cells/compartment", b2_tumor_20)
cat("B2 bootstrap median CIs (fixed combined-cell no-overlap score scale):\n")
for (nm in c("all23", "tumor23", "tumor15")) {
  ob <- switch(nm, all23 = b2_all, tumor23 = b2_tumor_all, tumor15 = b2_tumor_20)
  for (v in c("SenMayo_zmean", "SASP_zmean", "MKI67")) {
    ci <- boot_median_ci(ob[[v]], seed = 20260814 + match(nm, c("all23", "tumor23", "tumor15")) + match(v, c("SenMayo_zmean", "SASP_zmean", "MKI67")))
    cat("  ", nm, " ", v, " CI=", ci[1], " to ", ci[2], "\n", sep = "")
  }
}

cat("\nB2 CELL-LEVEL MIXED MODELS\n")
m_b2_sen <- lmer(SenMayo_rawmean ~ compartment + (1 | Patient), data = cell_df)
m_b2_sasp <- lmer(SASP_rawmean ~ compartment + (1 | Patient), data = cell_df)
m_b2_mki <- lmer(MKI67 ~ compartment + (1 | Patient), data = cell_df)
m_b2_sen_adj <- lmer(SenMayo_rawmean ~ compartment + MKI67 + (1 | Patient), data = cell_df)
print_lmer("All cells SenMayo compartment", m_b2_sen, "compartmentFibroblast")
print_lmer("All cells SASP compartment", m_b2_sasp, "compartmentFibroblast")
print_lmer("All cells MKI67 compartment", m_b2_mki, "compartmentFibroblast")
print_lmer("All cells SenMayo compartment + MKI67", m_b2_sen_adj, "compartmentFibroblast")
m_b2_sen_z <- lmer(SenMayo_zmean ~ compartment + (1 | Patient), data = cell_df)
m_b2_sasp_z <- lmer(SASP_zmean ~ compartment + (1 | Patient), data = cell_df)
m_b2_sen_z_adj <- lmer(SenMayo_zmean ~ compartment + MKI67 + (1 | Patient), data = cell_df)
print_lmer("REPORTED-SCALE All cells SenMayo combined-z compartment", m_b2_sen_z, "compartmentFibroblast")
print_lmer("REPORTED-SCALE All cells SASP combined-z compartment", m_b2_sasp_z, "compartmentFibroblast")
print_lmer("REPORTED-SCALE All cells SenMayo combined-z + MKI67", m_b2_sen_z_adj, "compartmentFibroblast")
m_b2_sen_full_z <- lmer(SenMayo_full_zmean ~ compartment + (1 | Patient), data = cell_df)
m_b2_sen_full_z_adj <- lmer(SenMayo_full_zmean ~ compartment + MKI67 + (1 | Patient), data = cell_df)
print_lmer("FULL SenMayo (4 overlaps retained) combined-z compartment", m_b2_sen_full_z, "compartmentFibroblast")
print_lmer("FULL SenMayo (4 overlaps retained) combined-z + MKI67", m_b2_sen_full_z_adj, "compartmentFibroblast")

cell_t <- cell_df[cell_df$Class == "Tumor", , drop = FALSE]
cell_t$SenMayo_tumorcombined_zmean <- score_zmean(sen_genes, cell_t$cell)
cell_t$SASP_tumorcombined_zmean <- score_zmean(sasp_genes, cell_t$cell)
m_b2_sen_t <- lmer(SenMayo_rawmean ~ compartment + MKI67 + logUMI + (1 | Patient), data = cell_t)
m_b2_sasp_t <- lmer(SASP_rawmean ~ compartment + MKI67 + logUMI + (1 | Patient), data = cell_t)
print_lmer("Tumor-only SenMayo compartment + MKI67 + logUMI", m_b2_sen_t, "compartmentFibroblast")
print_lmer("Tumor-only SASP compartment + MKI67 + logUMI", m_b2_sasp_t, "compartmentFibroblast")
m_b2_sen_t_z <- lmer(SenMayo_tumorcombined_zmean ~ compartment + MKI67 + logUMI + (1 | Patient), data = cell_t)
m_b2_sasp_t_z <- lmer(SASP_tumorcombined_zmean ~ compartment + MKI67 + logUMI + (1 | Patient), data = cell_t)
print_lmer("Tumor-only RECOMPUTED-z SenMayo compartment + MKI67 + logUMI", m_b2_sen_t_z, "compartmentFibroblast")
print_lmer("Tumor-only RECOMPUTED-z SASP compartment + MKI67 + logUMI", m_b2_sasp_t_z, "compartmentFibroblast")

tumor_z_means <- aggregate(cbind(SenMayo_tumorcombined_zmean, SASP_tumorcombined_zmean, MKI67) ~ Patient + compartment, cell_t, mean)
tumor_z_ns <- aggregate(cell ~ Patient + compartment, cell_t, length)
names(tumor_z_ns)[3] <- "n_cells"
tumor_z_wide <- reshape(merge(tumor_z_means, tumor_z_ns, by = c("Patient", "compartment")), idvar = "Patient", timevar = "compartment", direction = "wide")
for (minimum in c(1, 20)) {
  elig <- complete.cases(tumor_z_wide) & tumor_z_wide$n_cells.Epithelial >= minimum & tumor_z_wide$n_cells.Fibroblast >= minimum
  cat("Tumor-only recomputed-z patient contrasts, min cells=", minimum, " n=", sum(elig), "\n", sep = "")
  for (v in c("SenMayo_tumorcombined_zmean", "SASP_tumorcombined_zmean")) {
    d <- tumor_z_wide[[paste0(v, ".Fibroblast")]][elig] - tumor_z_wide[[paste0(v, ".Epithelial")]][elig]
    ci <- boot_median_ci(d, seed = 20260830 + minimum + match(v, c("SenMayo_tumorcombined_zmean", "SASP_tumorcombined_zmean")))
    cat("  ", v, ": median=", median(d), " CI=", ci[1], " to ", ci[2], " positive=", sum(d > 0), "/", length(d), " P=", wilcox.test(d, exact = TRUE)$p.value, "\n", sep = "")
  }
}

cat("\nPATIENT-LEVEL FAP-SENMAYO CORRELATION\n")
fib_patient_all <- aggregate(cbind(FAP, SenMayo_rawmean, SenMayo_fib_zmean, SenMayo_fib_full_zmean) ~ Patient, fib, mean)
fib_patient_t <- aggregate(cbind(FAP, SenMayo_rawmean, SenMayo_tumor_zmean) ~ Patient, fib_t, mean)
ct_all <- cor.test(fib_patient_all$FAP, fib_patient_all$SenMayo_fib_zmean, method = "spearman", exact = FALSE)
ct_all_full <- cor.test(fib_patient_all$FAP, fib_patient_all$SenMayo_fib_full_zmean, method = "spearman", exact = FALSE)
ct_t <- cor.test(fib_patient_t$FAP, fib_patient_t$SenMayo_tumor_zmean, method = "spearman", exact = FALSE)
cat("all tumor+normal fibroblasts (fib-z score): rho=", unname(ct_all$estimate), " P=", ct_all$p.value, " n=", nrow(fib_patient_all), "\n", sep = "")
cat("all tumor+normal fibroblasts (FULL overlapping score): rho=", unname(ct_all_full$estimate), " P=", ct_all_full$p.value, " n=", nrow(fib_patient_all), "\n", sep = "")
cat("tumor-only fibroblasts (tumor-z score): rho=", unname(ct_t$estimate), " P=", ct_t$p.value, " n=", nrow(fib_patient_t), "\n", sep = "")

cat("\nDONE\n")
