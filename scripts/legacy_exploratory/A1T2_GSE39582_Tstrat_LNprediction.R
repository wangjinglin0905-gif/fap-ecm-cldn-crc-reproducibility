# =============================================================================
# A1T2: GSE39582 T-stage stratification (external validation of A1-T) +
#       Lymph node metastasis prediction (TCGA N0/N+ + GSE39582 N0/N+)
# Rationale from thesis: FAP/TSR/TB are LN-metastasis risk factors (FAP OR=7.84);
# T-stage gradient in TCGA needs external confirmation; N prediction tests
# whether FAP-CAF score can stratify LN-positive CRC (translational angle).
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggpubr)
  library(survival); library(pROC)
})

PROJ_ROOT <- "."
FIG_DIR   <- file.path(PROJ_ROOT, "output", "figures")
TAB_DIR   <- file.path(PROJ_ROOT, "output", "tables")
DATA_DIR  <- file.path(PROJ_ROOT, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

ECM_SDC4_CD44 <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")
CLDN_CORE     <- c("CLDN1","CLDN2","CLDN4")

compute_zmean <- function(df, genes) {
  avail <- genes[genes %in% colnames(df)]
  sub <- as.data.frame(lapply(avail, function(g) {
    v <- df[[g]]
    if (all(is.na(v))) return(v)
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }))
  names(sub) <- avail
  rowMeans(sub, na.rm = TRUE)
}

# =============================================================================
# PART 1: GSE39582 T-stage stratification
# =============================================================================
cat("=== PART 1: GSE39582 T-stage stratification ===\n")
# Load GSE39582 expression (already parsed: 585 samples x 70 probes)
expr <- read.csv(file.path(DATA_DIR, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
pheno <- read.csv(file.path(DATA_DIR, "GSE39582", "GSE39582_pheno.csv"))
cat("expr:", dim(expr), "| pheno:", dim(pheno), "\n")

# Probe -> gene mapping
probe_map <- read.csv(file.path(DATA_DIR, "GSE39582", "gene2probe.csv"))
# Collapse probes to genes (mean of probes per gene).
# expr columns: col 1 = sample ID (GSM...), cols 2.. = probe expression.
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) {
    gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
  }
}
cat("Gene-level matrix:", ncol(gene_expr), "cols\n")

# Re-parse tnm.t/tnm.n from series matrix directly
sm_file <- file.path(DATA_DIR, "GSE39582", "GSE39582_series_matrix.txt.gz")
con <- gzfile(sm_file, "rt")
field_rows <- list(); sample_names <- NULL
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_characteristics_ch1")) {
    cells <- strsplit(substr(line, nchar("!Sample_characteristics_ch1") + 1, nchar(line)), "\t")[[1]]
    cells <- gsub('"', '', trimws(cells))
    field_rows[[length(field_rows) + 1]] <- cells
} else if (startsWith(line, "!Sample_geo_accession")) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar("!Sample_geo_accession") + 1, nchar(line)), "\t")[[1]])
    sample_names <- sample_names[sample_names != ""]
}
}
close(con)
# field names from first cells: first cell may be empty (tab prefix),
# field name is in the first cell containing ":" (usually cell #2)
field_names <- sapply(field_rows, function(r) {
  # find first cell containing ":"
  hit <- which(grepl(":", r))[1]
  if (!is.na(hit)) sub(":.*", "", r[hit]) else ""
})
# values: strip field-name prefix from all cells (fixed=TRUE avoids regex metachar issues with "()")
get_field <- function(fname) {
  idx <- which(field_names == fname)
  if (length(idx) == 0) return(NULL)
  vals <- field_rows[[idx[1]]]
  vals <- vals[vals != ""]
  vals <- sub(paste0(fname, ":"), "", vals, fixed = TRUE)
  trimws(vals)  # strip leading/trailing whitespace (" T4" -> "T4")
}
tnm_t <- get_field("tnm.t"); tnm_n <- get_field("tnm.n"); tnm_m <- get_field("tnm.m")
os_ev <- get_field("os.event"); os_delay <- get_field("os.delay (months)")
rfs_ev <- get_field("rfs.event"); rfs_delay <- get_field("rfs.delay")
sex_f <- get_field("Sex"); age_f <- get_field("age.at.diagnosis (year)")
mmr_f <- get_field("mmr.status")

n <- length(sample_names)
gse <- data.frame(
  sample = sample_names,
  tnm_t = if (length(tnm_t) == n) tnm_t else NA,
  tnm_n = if (length(tnm_n) == n) tnm_n else NA,
  tnm_m = if (length(tnm_m) == n) tnm_m else NA,
  os_event = as.numeric(if (length(os_ev) == n) os_ev else NA),
  os_delay = as.numeric(if (length(os_delay) == n) os_delay else NA),
  rfs_event = as.numeric(if (length(rfs_ev) == n) rfs_ev else NA),
  rfs_delay = as.numeric(if (length(rfs_delay) == n) rfs_delay else NA),
  Sex = if (length(sex_f) == n) sex_f else NA,
  age = as.numeric(if (length(age_f) == n) age_f else NA),
  mmr = if (length(mmr_f) == n) mmr_f else NA,
  stringsAsFactors = FALSE
)
cat("GSE39582 parsed:", nrow(gse), "samples\n")
cat("tnm.t dist:\n"); print(table(gse$tnm_t, useNA = "ifany"))
cat("tnm.n dist:\n"); print(table(gse$tnm_n, useNA = "ifany"))

# Merge expression scores
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
# Normalize T/N groups
gse$T_group <- ifelse(grepl("^T[12]$", gse$tnm_t), "T1-2",
               ifelse(grepl("^T[34]", gse$tnm_t), "T3-4", NA))
gse$T_grad <- ifelse(grepl("^T1$", gse$tnm_t), "T1",
              ifelse(grepl("^T2$", gse$tnm_t), "T2",
              ifelse(grepl("^T3$", gse$tnm_t), "T3",
              ifelse(grepl("^T4", gse$tnm_t), "T4", NA))))
gse$N_pos <- ifelse(grepl("N[12]", gse$tnm_n), 1, ifelse(grepl("N0", gse$tnm_n), 0, NA))

# Scores
gse$FAP_CAF <- compute_zmean(gse, c("FAP","COL1A1","COL1A2","COL3A1","FN1","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A"))
gse$ECM_SDC4_CD44 <- compute_zmean(gse, ECM_SDC4_CD44)
gse$CLDN_core <- compute_zmean(gse, CLDN_CORE)

# 1.1 T-stratified correlations
gse_t <- gse[!is.na(gse$T_group), ]
cat("\nGSE39582 T-stratified samples:", nrow(gse_t), "\n")
cat("  T1-2:", sum(gse_t$T_group == "T1-2"), "| T3-4:", sum(gse_t$T_group == "T3-4"), "\n")

cor_t12 <- cor.test(gse_t$FAP_CAF[gse_t$T_group == "T1-2"], gse_t$ECM_SDC4_CD44[gse_t$T_group == "T1-2"], method = "spearman")
cor_t34 <- cor.test(gse_t$FAP_CAF[gse_t$T_group == "T3-4"], gse_t$ECM_SDC4_CD44[gse_t$T_group == "T3-4"], method = "spearman")
cat(sprintf("FAP-ECM: T1-2 rho=%.3f P=%.3f | T3-4 rho=%.3f P=%.3f\n",
    cor_t12$estimate, cor_t12$p.value, cor_t34$estimate, cor_t34$p.value))

corc_t12 <- cor.test(gse_t$FAP_CAF[gse_t$T_group == "T1-2"], gse_t$CLDN_core[gse_t$T_group == "T1-2"], method = "spearman")
corc_t34 <- cor.test(gse_t$FAP_CAF[gse_t$T_group == "T3-4"], gse_t$CLDN_core[gse_t$T_group == "T3-4"], method = "spearman")
cat(sprintf("FAP-CLDN: T1-2 rho=%.3f P=%.3f | T3-4 rho=%.3f P=%.3f\n",
    corc_t12$estimate, corc_t12$p.value, corc_t34$estimate, corc_t34$p.value))

# 1.2 T gradient
cat("\nGSE39582 FAP-CLDN gradient (T1/T2/T3/T4):\n")
grad <- data.frame(t_grad = c("T1","T2","T3","T4"), n = NA, rho = NA, P = NA)
for (i in seq_len(nrow(grad))) {
  g <- grad$t_grad[i]
  idx <- gse_t$T_grad == g
  grad$n[i] <- sum(idx)
  if (sum(idx) >= 5) {
    cc <- cor.test(gse_t$FAP_CAF[idx], gse_t$CLDN_core[idx], method = "spearman")
    grad$rho[i] <- round(cc$estimate, 3); grad$P[i] <- signif(cc$p.value, 3)
  }
}
print(grad)

# 1.3 Survival T1-2 vs T3-4 (OS)
gse_t$os_time <- gse_t$os_delay
gse_s <- gse_t[!is.na(gse_t$os_time) & gse_t$os_time > 0, ]
cat("\nGSE39582 survival:", nrow(gse_s), "patients,", sum(gse_s$os_event), "OS events\n")
cat("  T1-2 events:", sum(gse_s$os_event[gse_s$T_group == "T1-2"]), "/", sum(gse_s$T_group == "T1-2"), "\n")
cat("  T3-4 events:", sum(gse_s$os_event[gse_s$T_group == "T3-4"]), "/", sum(gse_s$T_group == "T3-4"), "\n")
sd_t <- survdiff(Surv(os_time, os_event) ~ T_group, data = gse_s)
p_lr <- pchisq(sd_t$chisq, 1, lower.tail = FALSE)
cox_t <- coxph(Surv(os_time, os_event) ~ T_group, data = gse_s)
cat(sprintf("Log-rank P = %.4f | Cox T3-4 HR = %.2f (%.2f-%.2f) P=%.4f\n",
    p_lr, exp(coef(cox_t)), exp(confint(cox_t))[1], exp(confint(cox_t))[2],
    summary(cox_t)$coefficients[5]))

# KM plot
km_fit <- survfit(Surv(os_time, os_event) ~ T_group, data = gse_s)
png(file.path(FIG_DIR, "A1T2_GSE39582_KM_Tgroups.png"), width = 10, height = 7.5, units = "in", res = 300)
plot(km_fit, col = c("#377EB8", "#E41A1C"), lwd = 2, xlab = "Overall survival (months)",
     ylab = "Survival probability", main = sprintf("GSE39582: OS by T Group (n=%d, log-rank P=%.4f)", nrow(gse_s), p_lr))
legend("topright", c("T1-2", "T3-4"), col = c("#377EB8", "#E41A1C"), lwd = 2)
dev.off()

# Gradient plot
grad_plot <- grad[!is.na(grad$rho), ]
p1 <- ggplot(grad_plot, aes(x = factor(t_grad, levels = c("T1","T2","T3","T4")), y = rho)) +
  geom_line(aes(group = 1), color = "grey40", linetype = 2) +
  geom_point(aes(size = n), color = "#E64B35") +
  geom_text(aes(label = sprintf("rho=%.2f", rho)), vjust = -1.2, size = 3.6) +
  scale_size_continuous(range = c(3, 9)) +
  labs(x = "T stage (GSE39582)", y = "Spearman rho (FAP-CAF vs CLDN core)",
       title = "FAP-CLDN Correlation Gradient Across T Stages (GSE39582 external)") +
  theme_classic(base_size = 13) + ylim(-0.6, 0.6)
ggsave(file.path(FIG_DIR, "A1T2_GSE39582_FAP_CLDN_gradient.png"), p1, width = 7, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A1T2_GSE39582_FAP_CLDN_gradient.tiff"), p1, width = 7, height = 6, dpi = 600, compression = "lzw")

# =============================================================================
# PART 2: Lymph node metastasis prediction
# =============================================================================
cat("\n=== PART 2: Lymph node metastasis prediction ===\n")

# 2.1 TCGA: N0 vs N+
tcga <- read.csv(file.path(DATA_DIR, "A1_tcga_coad_merged.csv"))
tcga$N_pos <- ifelse(grepl("N[12]", tcga$ajcc_pathologic_n), 1,
              ifelse(grepl("N0", tcga$ajcc_pathologic_n), 0, NA))
expr_file <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"
expr_raw <- read.table(gzfile(expr_file), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) vals[i] <- expr_mat[gene, sid]
    else {
      pid <- substr(sid, 1, 12)
      matches <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(matches) > 0) vals[i] <- expr_mat[gene, matches[1]]
    }
  }
  vals
}
for (g in c(ECM_SDC4_CD44, CLDN_CORE)) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga$ECM_SDC4_CD44_score <- compute_zmean(tcga, ECM_SDC4_CD44)
tcga$CLDN_core_score <- compute_zmean(tcga, CLDN_CORE)
tcga <- tcga[!is.na(tcga$N_pos) & tcga$sample_scope == "Tumor", ]
cat("TCGA N analysis:", nrow(tcga), "tumors (N0:", sum(tcga$N_pos == 0), "N+:", sum(tcga$N_pos == 1), ")\n")

wt_fap <- wilcox.test(tcga$FAP_CAF[tcga$N_pos == 1], tcga$FAP_CAF[tcga$N_pos == 0])
wt_ecm <- wilcox.test(tcga$ECM_SDC4_CD44_score[tcga$N_pos == 1], tcga$ECM_SDC4_CD44_score[tcga$N_pos == 0])
cat(sprintf("TCGA FAP_CAF N+ vs N0: P=%.4f (median %.3f vs %.3f)\n", wt_fap$p.value,
    median(tcga$FAP_CAF[tcga$N_pos == 1]), median(tcga$FAP_CAF[tcga$N_pos == 0])))
cat(sprintf("TCGA ECM N+ vs N0: P=%.4f (median %.3f vs %.3f)\n", wt_ecm$p.value,
    median(tcga$ECM_SDC4_CD44_score[tcga$N_pos == 1]), median(tcga$ECM_SDC4_CD44_score[tcga$N_pos == 0])))

# AUC
roc_fap <- roc(tcga$N_pos, tcga$FAP_CAF, quiet = TRUE)
roc_ecm <- roc(tcga$N_pos, tcga$ECM_SDC4_CD44_score, quiet = TRUE)
cat(sprintf("AUC FAP_CAF = %.3f | AUC ECM-SDC4/CD44 = %.3f\n", auc(roc_fap), auc(roc_ecm)))

# Logistic N prediction
log_n <- glm(N_pos ~ FAP_CAF + ECM_SDC4_CD44_score, data = tcga, family = binomial)
cat("Logistic (FAP+ECM) summary:\n")
print(summary(log_n)$coefficients)

# 2.2 GSE39582: N0 vs N+
gse_n <- gse[!is.na(gse$N_pos), ]
cat("\nGSE39582 N analysis:", nrow(gse_n), "samples (N0:", sum(gse_n$N_pos == 0), "N+:", sum(gse_n$N_pos == 1), ")\n")
wt_fap_g <- wilcox.test(gse_n$FAP_CAF[gse_n$N_pos == 1], gse_n$FAP_CAF[gse_n$N_pos == 0])
wt_ecm_g <- wilcox.test(gse_n$ECM_SDC4_CD44[gse_n$N_pos == 1], gse_n$ECM_SDC4_CD44[gse_n$N_pos == 0])
cat(sprintf("GSE39582 FAP_CAF N+ vs N0: P=%.4f (median %.3f vs %.3f)\n", wt_fap_g$p.value,
    median(gse_n$FAP_CAF[gse_n$N_pos == 1]), median(gse_n$FAP_CAF[gse_n$N_pos == 0])))
cat(sprintf("GSE39582 ECM N+ vs N0: P=%.4f (median %.3f vs %.3f)\n", wt_ecm_g$p.value,
    median(gse_n$ECM_SDC4_CD44[gse_n$N_pos == 1]), median(gse_n$ECM_SDC4_CD44[gse_n$N_pos == 0])))
roc_fap_g <- roc(gse_n$N_pos, gse_n$FAP_CAF, quiet = TRUE)
roc_ecm_g <- roc(gse_n$N_pos, gse_n$ECM_SDC4_CD44, quiet = TRUE)
cat(sprintf("AUC FAP_CAF = %.3f | AUC ECM-SDC4/CD44 = %.3f\n", auc(roc_fap_g), auc(roc_ecm_g)))

# N+ prediction plot (TCGA + GSE39582)
n_plot <- data.frame(
  cohort = c(rep("TCGA-COAD", nrow(tcga)), rep("GSE39582", nrow(gse_n))),
  N_status = c(tcga$N_pos, gse_n$N_pos),
  FAP_CAF = c(tcga$FAP_CAF, gse_n$FAP_CAF),
  stringsAsFactors = FALSE
)
n_plot$N_status <- factor(ifelse(n_plot$N_status == 1, "LN+", "LN0"), levels = c("LN0", "LN+"))
p2 <- ggplot(n_plot, aes(x = cohort, y = FAP_CAF, fill = N_status)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.75) +
  stat_compare_means(aes(group = N_status), method = "wilcox.test", label = "p.format", size = 3.2) +
  scale_fill_manual(values = c("LN0" = "#377EB8", "LN+" = "#E41A1C")) +
  labs(x = "", y = "FAP-CAF score",
       title = "FAP-CAF Score by Lymph Node Status (TCGA + GSE39582)") +
  theme_classic(base_size = 13)
ggsave(file.path(FIG_DIR, "A1T2_LN_prediction_boxplot.png"), p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A1T2_LN_prediction_boxplot.tiff"), p2, width = 8, height = 6, dpi = 600, compression = "lzw")

# Save summary
sum_tab <- data.frame(
  analysis = c("GSE39582 FAP-ECM T1-2", "GSE39582 FAP-ECM T3-4",
               "GSE39582 FAP-CLDN T1-2", "GSE39582 FAP-CLDN T3-4",
               "GSE39582 KM logrank T1-2 vs T3-4",
               "TCGA FAP N+ vs N0 P", "TCGA AUC FAP", "TCGA AUC ECM",
               "GSE39582 FAP N+ vs N0 P", "GSE39582 AUC FAP", "GSE39582 AUC ECM"),
  value = c(sprintf("rho=%.3f P=%.3f", cor_t12$estimate, cor_t12$p.value),
            sprintf("rho=%.3f P=%.3f", cor_t34$estimate, cor_t34$p.value),
            sprintf("rho=%.3f P=%.3f", corc_t12$estimate, corc_t12$p.value),
            sprintf("rho=%.3f P=%.3f", corc_t34$estimate, corc_t34$p.value),
            sprintf("P=%.4f HR=%.2f", p_lr, exp(coef(cox_t))),
            sprintf("P=%.4f", wt_fap$p.value), sprintf("%.3f", auc(roc_fap)), sprintf("%.3f", auc(roc_ecm)),
            sprintf("P=%.4f", wt_fap_g$p.value), sprintf("%.3f", auc(roc_fap_g)), sprintf("%.3f", auc(roc_ecm_g))),
  stringsAsFactors = FALSE)
write.csv(sum_tab, file.path(TAB_DIR, "A1T2_summary.csv"), row.names = FALSE)
cat("\n=== [DONE] A1T2 complete ===\n")
print(sum_tab)
