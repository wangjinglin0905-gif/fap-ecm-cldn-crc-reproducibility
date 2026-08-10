# =============================================================================
# REVIEW_R6: (7) Nomogram for LN metastasis (T stage + StromalScore)
#            (8) Mixed-effects model for single-cell (cells nested in patients)
#            (5) Left/right colon stratification sensitivity
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(rms); library(pROC); library(lme4)
})

PROJ <- "."
OUT  <- file.path(PROJ, "output", "review_r6")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TCGA_EXPR <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"

# ---- ESTIMATE-style stromal signature (same as REVIEW_R5) ----
ESTIMATE_STROMAL_CORE <- unique(c(
  "COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL11A1","COL12A1",
  "COL15A1","COL16A1","COL18A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2",
  "FN1","SPARC","DCN","LUM","BGN","POSTN","THY1","PDPN","TAGLN","ACTA2",
  "FAP","MMP2","MMP11","MMP14","TIMP1","TIMP2","TIMP3","LOXL1","LOXL2",
  "LOX","PLOD1","PLOD2","PLOD3","VCAN","FBN1","FBLN1","FBLN2","ELN","MFAP4",
  "IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP7","TGFBI","CTGF","CYR61",
  "PDGFRB","PDGFRA","CXCL12","SDC1","SDC2","SDC4","CD44","VIM","S100A4",
  "ITGA1","ITGA2","ITGA5","ITGB1","ITGB3","ITGB5","EMP3","MGP","CNN1",
  "CALD1","LMOD1","MYLK","TPM2","FLNA","FLNB","ACTN1","ACTG2","DES","VWF"
))
FAP13 <- c("FAP","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
MATRIX4 <- c("COL1A1","COL1A2","COL3A1","FN1")

zmean <- function(df, genes) {
  avail <- genes[genes %in% colnames(df)]
  if (length(avail) < 5) return(rep(NA, nrow(df)))
  sub <- as.data.frame(lapply(avail, function(g) {
    v <- df[[g]]
    if (all(is.na(v))) return(v)
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }))
  names(sub) <- avail
  rowMeans(sub, na.rm = TRUE)
}

# =============================================================================
# PART A: TCGA data prep (shared by Nomogram + left/right)
# =============================================================================
cat("=== PART A: TCGA prep ===\n")
expr_raw <- read.table(gzfile(TCGA_EXPR), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
tcga <- read.csv(file.path(PROJ, "data", "A1_tcga_coad_merged.csv"))
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) vals[i] <- expr_mat[gene, sid]
    else {
      pid <- substr(sid, 1, 12)
      m <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(m) > 0) vals[i] <- expr_mat[gene, m[1]]
    }
  }
  vals
}
for (g in unique(c(FAP13, MATRIX4, ESTIMATE_STROMAL_CORE))) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga <- tcga[tcga$sample_scope == "Tumor", ]
tcga$FAP13 <- zmean(tcga, FAP13)
tcga$matrix4 <- zmean(tcga, MATRIX4)
tcga$StromalScore <- zmean(tcga, ESTIMATE_STROMAL_CORE)
tcga$N_pos <- ifelse(grepl("N[12]", tcga$ajcc_pathologic_n), 1,
              ifelse(grepl("N0", tcga$ajcc_pathologic_n), 0, NA))
tcga$T34 <- ifelse(grepl("^T[34]", tcga$ajcc_pathologic_t), 1,
            ifelse(grepl("^T[12]$", tcga$ajcc_pathologic_t), 0, NA))
tcga$Tnum <- ifelse(grepl("^T1$", tcga$ajcc_pathologic_t), 1,
             ifelse(grepl("^T2$", tcga$ajcc_pathologic_t), 2,
             ifelse(grepl("^T3$", tcga$ajcc_pathologic_t), 3,
             ifelse(grepl("^T4", tcga$ajcc_pathologic_t), 4, NA))))
# left/right colon
site <- tcga$site_of_resection_or_biopsy
tcga$side <- ifelse(site %in% c("Cecum","Ascending colon","Hepatic flexure of colon","Transverse colon"), "Right",
             ifelse(site %in% c("Splenic flexure of colon","Descending colon","Sigmoid colon"), "Left", NA))
cat("side dist:\n"); print(table(tcga$side, useNA = "ifany"))

# =============================================================================
# PART B (7): Nomogram - LN metastasis (T stage + StromalScore)
# =============================================================================
cat("\n=== PART B (7): Nomogram ===\n")
nomo_data <- tcga[!is.na(tcga$N_pos) & !is.na(tcga$StromalScore) & !is.na(tcga$Tnum), ]
cat("Nomogram n =", nrow(nomo_data), "(N0:", sum(nomo_data$N_pos == 0), "N+:", sum(nomo_data$N_pos == 1), ")\n")

# Keep only modeling variables for datadist (avoid non-numeric columns error)
nomo_model <- nomo_data[, c("N_pos", "Tnum", "StromalScore")]
dd <- datadist(nomo_model)
options(datadist = "dd")
fit <- lrm(N_pos ~ Tnum + StromalScore, data = nomo_model, x = TRUE, y = TRUE)
print(fit)

# calibration
cal <- calibrate(fit, B = 200)
png(file.path(OUT, "R6_nomogram_calibration.png"), width = 7, height = 7, units = "in", res = 300)
plot(cal, main = "Calibration curve (200 bootstrap resamples)")
dev.off()

# Nomogram plot
png(file.path(OUT, "R6_nomogram.png"), width = 10, height = 7, units = "in", res = 300)
plot(nomogram(fit, fun = plogis,
              funlabel = "Probability of LN metastasis",
              lp = FALSE),
     main = "Nomogram: LN metastasis (T stage + StromalScore)")
dev.off()

# ROC comparison: T alone vs T+Stromal
roc_t <- roc(nomo_data$N_pos, nomo_data$Tnum, quiet = TRUE)
roc_ts <- roc(nomo_data$N_pos, predict(fit), quiet = TRUE)
cat(sprintf("AUC(T alone) = %.3f | AUC(T + StromalScore) = %.3f\n", auc(roc_t), auc(roc_ts)))
p_auc <- roc.test(roc_t, roc_ts, method = "delong")
cat(sprintf("DeLong P = %.4f\n", p_auc$p.value))

# =============================================================================
# PART C (5): Left/right colon stratification
# =============================================================================
cat("\n=== PART C (5): Left/right stratification (TCGA) ===\n")
for (sd in c("Right", "Left")) {
  sub <- tcga[tcga$side == sd & !is.na(tcga$FAP13) & !is.na(tcga$matrix4), ]
  if (nrow(sub) >= 8) {
    cc <- cor.test(sub$FAP13, sub$matrix4, method = "spearman")
    cat(sprintf("%s colon: n=%d | FAP13-matrix4 rho=%.3f P=%.4f\n",
        sd, nrow(sub), cc$estimate, cc$p.value))
  }
}
# also FAP13-receptor2 by side (if receptor genes available - not extracted here,
# keep matrix4 as the key coupling test)

# =============================================================================
# PART D (8): Mixed-effects model - single cells nested in patients
# =============================================================================
cat("\n=== PART D (8): Mixed-effects model ===\n")
suppressPackageStartupMessages({library(monocle3); library(SummarizedExperiment)})
cds <- readRDS(file.path(PROJ, "output", "A5_monocle3_cds.rds"))
md  <- as.data.frame(colData(cds))
expr <- assay(cds)

FIB_LINEAGE <- c("Myofibroblasts","Stromal 1","Stromal 2","Stromal 3")
RECEPTOR2 <- c("SDC4","CD44")
MATRIX4 <- c("COL1A1","COL1A2","COL3A1","FN1")

keep <- md$Class == "Tumor" & md$Cell_subtype %in% FIB_LINEAGE
mdk <- md[keep, ]
present <- unique(c("FAP", MATRIX4, RECEPTOR2))
present <- present[present %in% rownames(expr)]
sub <- as.matrix(expr[present, keep, drop = FALSE])
tot <- colSums(sub); tot[tot == 0] <- NA
norm <- log1p(t(sub) / tot * 10000); norm <- t(norm)
keep2 <- !is.na(colSums(norm))
norm <- norm[, keep2]; mdk <- mdk[keep2, ]
cat("fib-lineage cells with scores:", ncol(norm), "| patients:", length(unique(mdk$Patient)), "\n")

# cell-level scores
cell_z <- function(gene) {
  v <- norm[gene, ]
  (v - mean(v)) / sd(v)
}
dat <- data.frame(
  Patient = mdk$Patient,
  FAP = cell_z("FAP"),
  matrix4 = rowMeans(sapply(MATRIX4[MATRIX4 %in% rownames(norm)], cell_z)),
  receptor2 = rowMeans(sapply(RECEPTOR2[RECEPTOR2 %in% rownames(norm)], cell_z)),
  stringsAsFactors = FALSE
)
cat("cells:", nrow(dat), "| patients:", length(unique(dat$Patient)), "\n")

# Mixed-effects model: receptor2 ~ FAP + (1|Patient)
m_fap <- lmer(receptor2 ~ FAP + (1 | Patient), data = dat)
cat("\nMixed model: receptor2 ~ FAP + (1|Patient)\n")
s <- summary(m_fap)
print(s$coefficients)
cat("Random effect SD (Patient):", round(as.data.frame(VarCorr(m_fap))$sdcor[1], 3), "\n")
cat("Residual SD:", round(s$sigma, 3), "\n")

m_mat <- lmer(matrix4 ~ FAP + (1 | Patient), data = dat)
cat("\nMixed model: matrix4 ~ FAP + (1|Patient)\n")
print(summary(m_mat)$coefficients)

# Patient-level (compare with Spearman from REVIEW_R2b for consistency)
pat <- dat %>% group_by(Patient) %>%
  summarise(FAP = mean(FAP), matrix4 = mean(matrix4), receptor2 = mean(receptor2), n = n()) %>%
  filter(n >= 20)
cat("\nPatient-level (>=20 cells), n =", nrow(pat), "\n")
c1 <- cor.test(pat$FAP, pat$matrix4, method = "spearman")
c2 <- cor.test(pat$FAP, pat$receptor2, method = "spearman")
cat(sprintf("Spearman: FAP-matrix4 rho=%.3f P=%.4f | FAP-receptor2 rho=%.3f P=%.4f\n",
    c1$estimate, c1$p.value, c2$estimate, c2$p.value))

# =============================================================================
# Summary table
# =============================================================================
sum_tab <- data.frame(
  analysis = c("Nomogram: T alone AUC", "Nomogram: T + Stromal AUC", "DeLong P",
               "LR chi2 (T+Stromal)", "LR P",
               "Left colon FAP13-matrix4 rho", "Right colon FAP13-matrix4 rho",
               "LMM receptor2 ~ FAP slope", "LMM receptor2 ~ FAP P",
               "LMM matrix4 ~ FAP slope", "LMM matrix4 ~ FAP P",
               "Patient-level FAP-matrix4 rho", "Patient-level FAP-receptor2 rho"),
  value = c(sprintf("%.3f", auc(roc_t)), sprintf("%.3f", auc(roc_ts)),
            sprintf("%.4f", p_auc$p.value),
            sprintf("%.1f", fit$stats["Model L.R."]), sprintf("%.2e", fit$stats["P"]),
            sprintf("rho=%.3f", {
              sL <- tcga[tcga$side == "Left" & !is.na(tcga$FAP13) & !is.na(tcga$matrix4), ]
              if (nrow(sL) >= 8) cor.test(sL$FAP13, sL$matrix4, method="spearman")$estimate else NA
            }),
            sprintf("rho=%.3f", {
              sR <- tcga[tcga$side == "Right" & !is.na(tcga$FAP13) & !is.na(tcga$matrix4), ]
              if (nrow(sR) >= 8) cor.test(sR$FAP13, sR$matrix4, method="spearman")$estimate else NA
            }),
            sprintf("%.3f", fixef(m_fap)["FAP"]),
            sprintf("%.4f", summary(m_fap)$coefficients["FAP", 5]),
            sprintf("%.3f", fixef(m_mat)["FAP"]),
            sprintf("%.4f", summary(m_mat)$coefficients["FAP", 5]),
            sprintf("%.3f", c1$estimate), sprintf("%.3f", c2$estimate)),
  stringsAsFactors = FALSE
)
write.csv(sum_tab, file.path(OUT, "review_r6_summary.csv"), row.names = FALSE)
print(sum_tab)
cat("\n=== [DONE] REVIEW_R6 ===\n")
