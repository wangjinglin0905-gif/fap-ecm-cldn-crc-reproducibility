# =============================================================================
# A1-T: T-stage stratified analysis (T1-2 vs T3-4, plus T1/T2/T3/T4 gradient)
# Rationale: pT1a/pT2 concept from endoscopy-pathology --- the stroma transition
# threshold localizes to the submucosal 1000um boundary, NOT the AJCC stage
# boundary. TCGA Stage I is mostly pT2 (45/48), so AJCC stratification cannot
# capture the "zero metastasis potential" window. T-stage re-stratification is
# the best available transcriptomic proxy.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(ggpubr); library(ggtext)
  library(survival); library(survminer)
})

PROJ_ROOT <- "."
FIG_DIR   <- file.path(PROJ_ROOT, "output", "figures")
TAB_DIR   <- file.path(PROJ_ROOT, "output", "tables")
DATA_DIR  <- file.path(PROJ_ROOT, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

ECM_SDC4_CD44 <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")
CLDN_CORE     <- c("CLDN1","CLDN2","CLDN4")

message("=== Loading data ===")
expr_file <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"
expr_raw <- read.table(gzfile(expr_file), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
cat("Expression matrix:", nrow(expr_raw), "genes x", ncol(expr_raw), "samples\n")

merged <- read.csv(file.path(DATA_DIR, "A1_tcga_coad_merged.csv"))
cat("Clinical+Scores:", nrow(merged), "patients\n")

# --- Recompute scores with same definitions as AJCR/A1_A4 --------------------
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) {
      vals[i] <- expr_mat[gene, sid]
    } else {
      pid <- substr(sid, 1, 12)
      matches <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(matches) > 0) vals[i] <- expr_mat[gene, matches[1]]
    }
  }
  vals
}

for (gene in c(ECM_SDC4_CD44, CLDN_CORE)) {
  merged[[gene]] <- extract_gene(gene, expr_raw, merged)
}

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

merged$ECM_SDC4_CD44_score <- compute_zmean(merged, ECM_SDC4_CD44)
merged$CLDN_core_score      <- compute_zmean(merged, CLDN_CORE)

# --- T-stage grouping --------------------------------------------------------
# Normalize T stage: T1/T2 -> 'T1-2', T3/T4/T4a/T4b -> 'T3-4'
merged$t_raw <- merged$ajcc_pathologic_t
merged$t_group <- ifelse(grepl("^T[12]$", merged$t_raw), "T1-2",
                  ifelse(grepl("^T[34]", merged$t_raw), "T3-4", NA))
# Gradient grouping T1/T2/T3/T4 (T4a/T4b -> T4)
merged$t_grad <- ifelse(grepl("^T1$", merged$t_raw), "T1",
                  ifelse(grepl("^T2$", merged$t_raw), "T2",
                  ifelse(grepl("^T3$", merged$t_raw), "T3",
                  ifelse(grepl("^T4", merged$t_raw), "T4", NA))))

# Filter tumor samples with T stage
tdata <- merged[!is.na(merged$t_group) & merged$sample_scope == "Tumor", ]
cat("T-stratified tumor samples:", nrow(tdata), "\n")
print(table(tdata$t_group, useNA = "ifany"))
cat("\nT gradient:\n")
print(table(tdata$t_grad, useNA = "ifany"))
cat("\nT1-2 vs T3-4 by AJCC stage:\n")
print(table(tdata$t_group, tdata$ajcc_pathologic_stage, useNA = "ifany"))

# --- 1. FAP-CAF vs ECM-SDC4/CD44 correlation by T group -----------------------
cat("\n=== A1-T: FAP_CAF vs ECM-SDC4/CD44 (T-stratified) ===\n")
cor_t12 <- cor.test(tdata$FAP_CAF[tdata$t_group == "T1-2"],
                    tdata$ECM_SDC4_CD44_score[tdata$t_group == "T1-2"],
                    method = "spearman")
cor_t34 <- cor.test(tdata$FAP_CAF[tdata$t_group == "T3-4"],
                    tdata$ECM_SDC4_CD44_score[tdata$t_group == "T3-4"],
                    method = "spearman")
cat(sprintf("  T1-2: rho=%.3f P=%.4f (n=%d)\n", cor_t12$estimate, cor_t12$p.value,
            sum(tdata$t_group == "T1-2")))
cat(sprintf("  T3-4: rho=%.3f P=%.4f (n=%d)\n", cor_t34$estimate, cor_t34$p.value,
            sum(tdata$t_group == "T3-4")))

# --- 2. FAP-CAF vs CLDN core by T group ---------------------------------------
cat("\n=== A1-T: FAP_CAF vs CLDN_core (T-stratified) ===\n")
corc_t12 <- cor.test(tdata$FAP_CAF[tdata$t_group == "T1-2"],
                     tdata$CLDN_core_score[tdata$t_group == "T1-2"],
                     method = "spearman")
corc_t34 <- cor.test(tdata$FAP_CAF[tdata$t_group == "T3-4"],
                     tdata$CLDN_core_score[tdata$t_group == "T3-4"],
                     method = "spearman")
cat(sprintf("  T1-2: rho=%.3f P=%.4f (n=%d)\n", corc_t12$estimate, corc_t12$p.value,
            sum(tdata$t_group == "T1-2")))
cat(sprintf("  T3-4: rho=%.3f P=%.4f (n=%d)\n", corc_t34$estimate, corc_t34$p.value,
            sum(tdata$t_group == "T3-4")))

# --- 3. FAP-CAF vs CLDN core gradient across T1/T2/T3/T4 -----------------------
cat("\n=== A1-T: FAP_CAF vs CLDN_core gradient (T1/T2/T3/T4) ===\n")
grad_tab <- data.frame(t_grad = c("T1","T2","T3","T4"),
                       n = NA, rho = NA, P = NA, stringsAsFactors = FALSE)
for (i in seq_len(nrow(grad_tab))) {
  g <- grad_tab$t_grad[i]
  idx <- tdata$t_grad == g
  grad_tab$n[i] <- sum(idx)
  if (sum(idx) >= 5) {
    cc <- cor.test(tdata$FAP_CAF[idx], tdata$CLDN_core_score[idx], method = "spearman")
    grad_tab$rho[i] <- round(cc$estimate, 3); grad_tab$P[i] <- signif(cc$p.value, 3)
  }
}
print(grad_tab)

# --- 4. Survival: T1-2 vs T3-4 KM + Cox ---------------------------------------
cat("\n=== A1-T: Survival T1-2 vs T3-4 ===\n")
tdata$os_time <- as.numeric(tdata$days_to_death)
tdata$os_time[is.na(tdata$os_time)] <- as.numeric(tdata$days_to_last_follow_up)[is.na(tdata$os_time)]
tdata$os_event <- ifelse(tdata$vital_status == "Dead", 1, 0)
tdata$age <- as.numeric(tdata$age_at_diagnosis) / 365.25

sdat <- tdata[!is.na(tdata$os_time) & tdata$os_time > 0, ]
cat("Survival analysis on", nrow(sdat), "patients,", sum(sdat$os_event), "events\n")
cat("  T1-2 events:", sum(sdat$os_event[sdat$t_group == "T1-2"]), "/", sum(sdat$t_group == "T1-2"), "\n")
cat("  T3-4 events:", sum(sdat$os_event[sdat$t_group == "T3-4"]), "/", sum(sdat$t_group == "T3-4"), "\n")

# KM by T group
fit_t <- survfit(Surv(os_time, os_event) ~ t_group, data = sdat)
sd_t <- survdiff(Surv(os_time, os_event) ~ t_group, data = sdat)
p_lr <- pchisq(sd_t$chisq, 1, lower.tail = FALSE)
cat(sprintf("  Log-rank P = %.4f\n", p_lr))

# Cox univariate T1-2 vs T3-4
cox_t <- coxph(Surv(os_time, os_event) ~ t_group, data = sdat)
cat("  Cox T3-4 vs T1-2: HR =", round(exp(coef(cox_t)), 2),
    "95%CI", round(exp(confint(cox_t))[1], 2), "-", round(exp(confint(cox_t))[2], 2),
    "P =", signif(summary(cox_t)$coefficients[5], 3), "\n")

# Cox multivariable: FAP_CAF + age + t_group
sdat2 <- sdat[complete.cases(sdat$FAP_CAF, sdat$age, sdat$os_time, sdat$os_event), ]
cox_multi <- coxph(Surv(os_time, os_event) ~ FAP_CAF + age + t_group, data = sdat2)
co <- summary(cox_multi)$coefficients; ci <- summary(cox_multi)$conf.int
multi_tab <- data.frame(
  term = rownames(co),
  HR = round(ci[, "exp(coef)"], 2),
  lower = round(ci[, "lower .95"], 2),
  upper = round(ci[, "upper .95"], 2),
  P = signif(co[, "Pr(>|z|)"], 3),
  stringsAsFactors = FALSE)
cat("\nMultivariable Cox (FAP_CAF + age + T group):\n")
print(multi_tab)
write.csv(multi_tab, file.path(TAB_DIR, "A1T_multivariable_Cox.csv"), row.names = FALSE)

# --- 5. Figures ----------------------------------------------------------------
# Fig: scatter FAP-CAF vs ECM by T group
p1 <- ggplot(tdata, aes(x = FAP_CAF, y = ECM_SDC4_CD44_score, color = t_group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15) +
  stat_cor(method = "spearman", label.x.npc = "left", size = 3.5) +
  scale_color_manual(values = c("T1-2" = "#377EB8", "T3-4" = "#E41A1C")) +
  labs(x = "FAP-CAF Score", y = "ECM-SDC4/CD44 Score",
       title = "FAP Stromal Program vs ECM-SDC4/CD44 Axis (T-Stratified, TCGA-COAD)",
       color = "T stage") +
  theme_classic(base_size = 12)
ggsave(file.path(FIG_DIR, "A1T_FAP_ECM_corr_Tstrat.png"), p1, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(FIG_DIR, "A1T_FAP_ECM_corr_Tstrat.tiff"), p1, width = 10, height = 6.5, dpi = 600, compression = "lzw")

# Fig: correlation barplot FAP-CLDN by T group vs gradient
cor_df <- data.frame(
  group = c("T1-2","T3-4"),
  rho = c(corc_t12$estimate, corc_t34$estimate),
  P = c(corc_t12$p.value, corc_t34$p.value),
  stringsAsFactors = FALSE)
p2 <- ggplot(cor_df, aes(x = group, y = rho, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("rho=%.2f\nP=%.3f", rho, P)), vjust = -0.4, size = 3.8) +
  scale_fill_manual(values = c("T1-2" = "#377EB8", "T3-4" = "#E41A1C")) +
  ylim(-0.5, 0.8) +
  labs(x = "T stage group", y = "Spearman rho (FAP-CAF vs CLDN core)",
       title = "FAP-CAF vs CLDN Core: T-Stratified Correlation (TCGA-COAD)") +
  theme_classic(base_size = 13) + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "A1T_FAP_CLDN_corr_Tstrat.png"), p2, width = 7, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A1T_FAP_CLDN_corr_Tstrat.tiff"), p2, width = 7, height = 6, dpi = 600, compression = "lzw")

# Fig: KM by T group
km_plot <- ggsurvplot(fit_t, data = sdat, pval = TRUE, risk.table = TRUE,
  palette = c("#377EB8", "#E41A1C"),
  legend.title = "T stage", legend.labs = c("T1-2", "T3-4"),
  xlab = "Overall survival (days)", ylab = "Survival probability",
  title = sprintf("OS by T Stage Group - TCGA-COAD (n=%d)", nrow(sdat)),
  ggtheme = theme_classic())
png(file.path(FIG_DIR, "A1T_KM_Tgroups.png"), width = 10, height = 7.5, units = "in", res = 300)
print(km_plot)
dev.off()
tiff(file.path(FIG_DIR, "A1T_KM_Tgroups.tiff"), width = 10, height = 7.5, units = "in", res = 600, compression = "lzw")
print(km_plot)
dev.off()

# Fig: FAP-CLDN gradient across T1/T2/T3/T4 (mirror of thesis Table 6)
grad_plot_df <- data.frame(
  t_grad = factor(grad_tab$t_grad, levels = c("T1","T2","T3","T4")),
  rho = grad_tab$rho, n = grad_tab$n, stringsAsFactors = FALSE)
p3 <- ggplot(grad_plot_df, aes(x = t_grad, y = rho)) +
  geom_line(aes(group = 1), color = "grey40", linetype = 2) +
  geom_point(aes(size = n), color = "#E64B35") +
  geom_text(aes(label = sprintf("rho=%.2f", rho)), vjust = -1.2, size = 3.6) +
  scale_size_continuous(range = c(3, 9)) +
  labs(x = "T stage (TCGA-COAD)", y = "Spearman rho (FAP-CAF vs CLDN core)",
       title = "FAP-CLDN Correlation Gradient Across T Stages",
       subtitle = "Thesis Table 6 (IHC): FAP(++) 15%->54%->82%->86%->100% (HGIN+pT1a->pT1b->pT2->pT3->pT4)") +
  theme_classic(base_size = 13) + ylim(-0.6, 0.6)
ggsave(file.path(FIG_DIR, "A1T_FAP_CLDN_gradient.png"), p3, width = 7, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A1T_FAP_CLDN_gradient.tiff"), p3, width = 7, height = 6, dpi = 600, compression = "lzw")

# --- 6. Save summary table -----------------------------------------------------
summary_tab <- data.frame(
  analysis = c("FAP-ECM T1-2", "FAP-ECM T3-4", "FAP-CLDN T1-2", "FAP-CLDN T3-4",
               "KM logrank T1-2 vs T3-4"),
  rho_or_HR = c(round(cor_t12$estimate,3), round(cor_t34$estimate,3),
                round(corc_t12$estimate,3), round(corc_t34$estimate,3),
                round(exp(coef(cox_t)),2)),
  P = c(signif(cor_t12$p.value,3), signif(cor_t34$p.value,3),
        signif(corc_t12$p.value,3), signif(corc_t34$p.value,3),
        signif(p_lr,3)),
  stringsAsFactors = FALSE)
write.csv(summary_tab, file.path(TAB_DIR, "A1T_summary.csv"), row.names = FALSE)
cat("\n=== [DONE] A1-T analysis complete ===\n")
print(summary_tab)
