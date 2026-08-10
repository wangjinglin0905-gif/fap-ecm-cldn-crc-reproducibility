# =============================================================================
# A1 deepening + A4: Cox + KM + Claudin individual analysis
# Data: TCGA-COAD expression + clinical (112 patients with complete data)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggpubr); library(ggtext)
  library(survival); library(survminer)
  library(pheatmap); library(corrplot)
})

PROJ_ROOT <- "."
FIG_DIR   <- file.path(PROJ_ROOT, "output", "figures")
TAB_DIR   <- file.path(PROJ_ROOT, "output", "tables")
DATA_DIR  <- file.path(PROJ_ROOT, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

CLAUDIN_GENES <- c("CLDN1","CLDN2","CLDN3","CLDN4","CLDN7",
                   "CLDN8","CLDN12","CLDN15","CLDN18","CLDN23")
ECM_SDC4_CD44 <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")

message("=== Loading data ===")

# 1) Load expression matrix (gzipped)
expr_file <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"
cat("Loading expression from:", expr_file, "\n")
expr_raw <- read.table(gzfile(expr_file), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
cat("Expression matrix:", nrow(expr_raw), "genes x", ncol(expr_raw), "samples\n")

# 2) Load merged clinical+scores
merged <- read.csv(file.path(DATA_DIR, "A1_tcga_coad_merged.csv"))
cat("Clinical+Scores:", nrow(merged), "patients\n")

# 3) Match expression samples to clinical
# Expression colnames: TCGA-XX-XXXX-XX format; clinical has submitter_id
avail_genes <- CLAUDIN_GENES[CLAUDIN_GENES %in% rownames(expr_raw)]
avail_ecm   <- ECM_SDC4_CD44[ECM_SDC4_CD44 %in% rownames(expr_raw)]
cat("Available Claudin genes:", paste(avail_genes, collapse=", "), "\n")
cat("Available ECM/SDC4/CD44:", paste(avail_ecm, collapse=", "), "\n")

# Extract per-sample expression for Claudins & ECM genes
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  # match by sample_id (TCGA-XX-XXXX-01)
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    # Expression colnames may have different formats; try exact match then fuzzy
    if (sid %in% colnames(expr_mat)) {
      vals[i] <- expr_mat[gene, sid]
    } else {
      # Try matching by patient portion
      pid <- substr(sid, 1, 12)
      matches <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(matches) > 0) vals[i] <- expr_mat[gene, matches[1]]
    }
  }
  vals
}

# Add gene expression columns
for (gene in c(avail_genes, avail_ecm)) {
  merged[[gene]] <- extract_gene(gene, expr_raw, merged)
}
cat("Gene expression columns added.\n")

# Compute stage_group from ajcc_pathologic_stage
merged$stage_group <- ifelse(merged$ajcc_pathologic_stage == "Stage I", "Early (I)",
  ifelse(merged$ajcc_pathologic_stage %in% c("Stage II","Stage III","Stage IV"), "Late (II-IV)", NA))
# Filter to tumor samples with complete data
merged <- merged[!is.na(merged$stage_group) & merged$sample_scope == "Tumor", ]
cat("Analyzing", nrow(merged), "tumor samples: Early(I) =",
    sum(merged$stage_group=="Early (I)"), "Late(II-IV) =", sum(merged$stage_group=="Late (II-IV)"), "\n")

# =============================================================================
# A1 DEEPENING
# =============================================================================
message("\n=== A1: Cox Regression ===")

# Prepare survival
merged$os_time <- as.numeric(merged$days_to_death)
merged$os_time[is.na(merged$os_time)] <- as.numeric(merged$days_to_last_follow_up)[is.na(merged$os_time)]
merged$os_event <- ifelse(merged$vital_status == "Dead", 1, 0)
merged$age <- as.numeric(merged$age_at_diagnosis) / 365.25  # days -> years

# Exclude patients with missing critical data
cox_data <- merged[!is.na(merged$os_time) & merged$os_time > 0 & !is.na(merged$FAP_CAF), ]
cat("Cox analysis on", nrow(cox_data), "patients with", sum(cox_data$os_event), "events\n")

# Univariate KM: FAP_High vs Low in early stage
cox_data$FAP_binary <- ifelse(cox_data$FAP_CAF > median(cox_data$FAP_CAF, na.rm=TRUE), "High", "Low")

# All stages KM
fit_all <- survfit(Surv(os_time, os_event) ~ FAP_binary, data = cox_data)
p_km_all <- ggsurvplot(fit_all, data = cox_data, pval = TRUE, risk.table = TRUE,
  legend.title = "FAP-CAF Score", legend.labs = c("High","Low"),
  title = "All Stages: OS by FAP-CAF Score",
  palette = c("#E41A1C","#377EB8"))
png(file.path(FIG_DIR, "Fig1A_KM_all_stages.png"), width = 10, height = 8, units = "in", res = 300)
print(p_km_all)
dev.off()

# Early stage KM
cox_early <- cox_data[cox_data$stage_group == "Early (I)", ]
cox_early$FAP_binary <- ifelse(cox_early$FAP_CAF > median(cox_early$FAP_CAF, na.rm=TRUE), "High", "Low")
cat("Early stage KM:", nrow(cox_early), "patients,", sum(cox_early$os_event), "events\n")

if (sum(cox_early$os_event) >= 3 && length(unique(cox_early$FAP_binary)) == 2) {
  fit_early <- survfit(Surv(os_time, os_event) ~ FAP_binary, data = cox_early)
  p_km_early <- ggsurvplot(fit_early, data = cox_early, pval = TRUE, risk.table = TRUE,
    legend.title = "FAP-CAF Score", legend.labs = c("High","Low"),
    title = "Stage I: OS by FAP-CAF Score",
    palette = c("#E41A1C","#377EB8"))
  png(file.path(FIG_DIR, "Fig1B_KM_stage_I.png"), width = 10, height = 8, units = "in", res = 300)
  print(p_km_early)
  dev.off()
} else {
  cat("  Insufficient events for stage I KM (need >=3 events)\n")
}

# Multivariate Cox
cox_multi <- coxph(Surv(os_time, os_event) ~ FAP_CAF + age + stage_group + sex_at_birth, data = cox_data)
cox_summary <- summary(cox_multi)
cat("\nMultivariate Cox:\n")
print(cox_summary$coefficients)
cat("PH test P:", min(cox.zph(cox_multi)$table[,"p"], na.rm=TRUE), "\n")

# Forest plot
p_forest <- ggforest(cox_multi, data = cox_data, main = "Multivariate Cox: Overall Survival")
png(file.path(FIG_DIR, "Fig1C_Cox_forest.png"), width = 10, height = 6, units = "in", res = 300)
print(p_forest)
dev.off()

# Save Cox table
cox_tab <- as.data.frame(cox_summary$coefficients)
cox_tab$Variable <- rownames(cox_tab)
write.csv(cox_tab, file.path(TAB_DIR, "A1_multivariate_cox.csv"), row.names = FALSE)

# =============================================================================
# ECM-SDC4/CD44 score and correlation
# =============================================================================
message("\n=== A1: ECM-SDC4/CD44 correlation ===")

# Compute ECM-SDC4/CD44 score as mean z-score of available genes
compute_score <- function(data, genes) {
  available <- genes[genes %in% colnames(data)]
  if (length(available) < 2) return(rep(NA, nrow(data)))
  sub <- data[, available, drop = FALSE]
  sub_z <- scale(sub)
  rowMeans(sub_z, na.rm = TRUE)
}

merged$ECM_SDC4_CD44_score <- compute_score(merged, avail_ecm)

# Stage-stratified correlations
cor_early <- cor.test(
  merged$FAP_CAF[merged$stage_group == "Early (I)"],
  merged$ECM_SDC4_CD44_score[merged$stage_group == "Early (I)"],
  method = "spearman")
cor_late <- cor.test(
  merged$FAP_CAF[merged$stage_group == "Late (II-IV)"],
  merged$ECM_SDC4_CD44_score[merged$stage_group == "Late (II-IV)"],
  method = "spearman")
cat(sprintf("FAP-CAF vs ECM-SDC4/CD44: Early(I) rho=%.3f P=%.4f, Late(II-IV) rho=%.3f P=%.4f\n",
            cor_early$estimate, cor_early$p.value, cor_late$estimate, cor_late$p.value))

# Correlation scatter plot (stage-stratified)
p_scatter <- ggplot(merged, aes(x = FAP_CAF, y = ECM_SDC4_CD44_score, color = stage_group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  stat_cor(method = "spearman", label.x.npc = "left", size = 3.5) +
  facet_wrap(~stage_group) +
  scale_color_manual(values = c("Early (I)" = "#377EB8", "Late (II-IV)" = "#E41A1C")) +
  labs(x = "FAP-CAF Score", y = "ECM-SDC4/CD44 Score",
       title = "FAP Stromal Program vs ECM-SDC4/CD44 Axis (Stage-Stratified)",
       subtitle = paste0("Spearman: Early(I) \u03c1=", round(cor_early$estimate,3), 
                         ", Late(II-IV) \u03c1=", round(cor_late$estimate,3))) +
  theme_minimal(base_size = 12)
png(file.path(FIG_DIR, "Fig1D_correlation_stage_stratified.png"), width = 12, height = 6, units = "in", res = 300)
print(p_scatter)
dev.off()

# =============================================================================
# A4: CLAUDIN INDIVIDUAL ANALYSIS
# =============================================================================
message("\n=== A4: Claudin Individual Analysis ===")

# Compute claudin correlations with FAP_CAF, stratified
claudin_results <- data.frame()
for (gene in avail_genes) {
  if (sum(!is.na(merged[[gene]])) < 10) next
  
  # Early stage
  e_rows <- merged$stage_group == "Early (I)" & !is.na(merged[[gene]])
  if (sum(e_rows) >= 10) {
    ce <- cor.test(merged$FAP_CAF[e_rows], merged[[gene]][e_rows], method = "spearman")
  } else { ce <- list(estimate = NA, p.value = NA) }
  
  # Late stage
  l_rows <- merged$stage_group == "Late (II-IV)" & !is.na(merged[[gene]])
  if (sum(l_rows) >= 10) {
    cl <- cor.test(merged$FAP_CAF[l_rows], merged[[gene]][l_rows], method = "spearman")
  } else { cl <- list(estimate = NA, p.value = NA) }
  
  # Overall
  co <- cor.test(merged$FAP_CAF[!is.na(merged[[gene]])],
                 merged[[gene]][!is.na(merged[[gene]])], method = "spearman")
  
  claudin_results <- rbind(claudin_results, data.frame(
    Gene = gene,
    rho_overall = co$estimate,  p_overall = co$p.value,
    rho_early = ce$estimate,    p_early = ce$p.value,
    rho_late = cl$estimate,     p_late = cl$p.value,
    rho_delta = (ce$estimate %||% 0) - (cl$estimate %||% 0),
    n_early = sum(e_rows), n_late = sum(l_rows),
    early_mean = mean(merged[[gene]][e_rows], na.rm = TRUE),
    late_mean  = mean(merged[[gene]][l_rows], na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

# Sort by absolute delta
claudin_results <- claudin_results[order(abs(claudin_results$rho_delta), decreasing = TRUE), ]
write.csv(claudin_results, file.path(TAB_DIR, "A4_claudin_fap_correlation.csv"), row.names = FALSE)
cat("Claudin results:\n")
print(claudin_results[, c("Gene","rho_early","rho_late","rho_delta","p_early","p_late")])

# Waterfall plot: delta rho
p_waterfall <- ggplot(claudin_results, aes(x = reorder(Gene, rho_delta), y = rho_delta, 
                                           fill = rho_delta > 0)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", rho_delta), 
                hjust = ifelse(rho_delta > 0, -0.2, 1.2)), size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"), guide = "none") +
  coord_flip() +
  labs(x = "", y = expression(Delta*rho ~ "(Early - Late)"),
       title = "Claudin-FAP Correlation: Stage-Specific Divergence",
       subtitle = expression("Positive " * Delta * rho * " = stronger FAP-CAF correlation in early-stage CRC")) +
  theme_minimal(base_size = 13)
png(file.path(FIG_DIR, "Fig4A_claudin_waterfall.png"), width = 8, height = 6, units = "in", res = 300)
print(p_waterfall)
dev.off()

# Heatmap: rho_early vs rho_late
heat_data <- as.matrix(claudin_results[, c("rho_early","rho_late","rho_delta")])
rownames(heat_data) <- claudin_results$Gene
colnames(heat_data) <- c("Early (I)", "Late (II-IV)", "Delta")
p_heat <- pheatmap(heat_data, cluster_rows = TRUE, cluster_cols = FALSE,
  display_numbers = TRUE, number_format = "%.3f",
  color = colorRampPalette(c("#377EB8","white","#E41A1C"))(100),
  main = "Claudin-FAP Stromal Correlation (Stage-Stratified Spearman \u03c1)",
  fontsize_number = 9, fontsize = 11)
png(file.path(FIG_DIR, "Fig4B_claudin_heatmap.png"), width = 7, height = 5, units = "in", res = 300)
print(p_heat)
dev.off()

# Save merged data with all gene columns
write.csv(merged, file.path(DATA_DIR, "A1_A4_merged_full.csv"), row.names = FALSE)

# =============================================================================
# OUTPUT SUMMARY
# =============================================================================
cat("\n========================================\n")
cat("A1 + A4 Analysis Complete!\n")
cat("========================================\n")
cat("\n--- KEY FINDINGS ---\n")
cat(sprintf("1. FAP-CAF vs ECM-SDC4/CD44 (stage-stratified):\n"))
cat(sprintf("   Early(I) rho=%.3f P=%.4f\n", cor_early$estimate, cor_early$p.value))
cat(sprintf("   Late(II-IV) rho=%.3f P=%.4f\n", cor_late$estimate, cor_late$p.value))
cat(sprintf("2. Multivariate Cox: FAP_CAF HR=%.2f (%.2f-%.2f), P=%.4f\n",
            exp(coef(cox_multi)["FAP_CAF"]),
            exp(confint(cox_multi)["FAP_CAF",1]),
            exp(confint(cox_multi)["FAP_CAF",2]),
            summary(cox_multi)$coefficients["FAP_CAF","Pr(>|z|)"]))
cat(sprintf("3. Claudin with largest delta: %s (delta=%.3f)\n",
            claudin_results$Gene[1], claudin_results$rho_delta[1]))
cat(sprintf("4. Early events: %d/%d patients\n",
            sum(cox_early$os_event), nrow(cox_early)))
cat("\nFigures saved to:", FIG_DIR, "\n")
cat("Tables saved to:", TAB_DIR, "\n")
cat("========================================\n")
