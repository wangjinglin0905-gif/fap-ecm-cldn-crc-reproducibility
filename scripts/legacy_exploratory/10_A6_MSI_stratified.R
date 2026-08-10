# =============================================================================
# A6: MSI-status stratification in GSE39582 (pMMR vs dMMR)
# Question: does the FAP-ECM / FAP-CLDN relationship differ by MMR status?
# - FAP_CAF / ECM_SDC4_CD44 / CLDN_core scores: pMMR vs dMMR (Wilcoxon)
# - FAP-ECM & FAP-CLDN Spearman correlations within each MMR group
# - OS survival by MMR group (sanity check: dMMR known to have better OS)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggpubr)
  library(survival)
})

PROJ_ROOT <- "."
FIG_DIR   <- file.path(PROJ_ROOT, "output", "figures")
TAB_DIR   <- file.path(PROJ_ROOT, "output", "tables")
DATA_DIR  <- file.path(PROJ_ROOT, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

ECM_SDC4_CD44 <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")
CLDN_CORE     <- c("CLDN1","CLDN2","CLDN4")
FAP_CAF_17    <- c("FAP","COL1A1","COL1A2","COL3A1","FN1","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")

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

cat("=== A6: MSI stratification (GSE39582) ===\n")

# ---- Load expression ----
expr <- read.csv(file.path(DATA_DIR, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(DATA_DIR, "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}

# ---- Re-parse MMR + clinical fields from series matrix ----
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
field_names <- sapply(field_rows, function(r) {
  hit <- which(grepl(":", r))[1]
  if (!is.na(hit)) sub(":.*", "", r[hit]) else ""
})
get_field <- function(fname) {
  idx <- which(field_names == fname)
  if (length(idx) == 0) return(NULL)
  vals <- field_rows[[idx[1]]]
  vals <- vals[vals != ""]
  vals <- sub(paste0(fname, ":"), "", vals, fixed = TRUE)
  trimws(vals)
}
mmr_f    <- get_field("mmr.status")
tnm_t    <- get_field("tnm.t")
tnm_n    <- get_field("tnm.n")
os_ev    <- get_field("os.event")
os_delay <- get_field("os.delay (months)")
rfs_ev   <- get_field("rfs.event")
rfs_delay<- get_field("rfs.delay")

n <- length(sample_names)
gse <- data.frame(
  sample = sample_names,
  mmr = if (length(mmr_f) == n) mmr_f else NA,
  tnm_t = if (length(tnm_t) == n) tnm_t else NA,
  tnm_n = if (length(tnm_n) == n) tnm_n else NA,
  os_event = as.numeric(if (length(os_ev) == n) os_ev else NA),
  os_delay = as.numeric(if (length(os_delay) == n) os_delay else NA),
  rfs_event = as.numeric(if (length(rfs_ev) == n) rfs_ev else NA),
  rfs_delay = as.numeric(if (length(rfs_delay) == n) rfs_delay else NA),
  stringsAsFactors = FALSE
)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
gse$FAP_CAF    <- compute_zmean(gse, FAP_CAF_17)
gse$ECM_SDC4_CD44 <- compute_zmean(gse, ECM_SDC4_CD44)
gse$CLDN_core  <- compute_zmean(gse, CLDN_CORE)
gse$MSI_group  <- ifelse(gse$mmr == "dMMR", "dMMR (MSI)",
                  ifelse(gse$mmr == "pMMR", "pMMR (MSS)", NA))
cat("MMR distribution:\n"); print(table(gse$mmr, useNA = "ifany"))

msi <- gse[!is.na(gse$MSI_group), ]
cat("\nAnalyzable:", nrow(msi), "(dMMR:", sum(msi$MSI_group == "dMMR (MSI)"),
    "| pMMR:", sum(msi$MSI_group == "pMMR (MSS)"), ")\n")

# ===== 1. Score comparison: pMMR vs dMMR =====
cat("\n--- Score comparison (Wilcoxon) ---\n")
score_list <- list(FAP_CAF = "FAP_CAF", ECM_SDC4_CD44 = "ECM_SDC4_CD44", CLDN_core = "CLDN_core")
score_tab <- data.frame(score = names(score_list), p_value = NA, dMMR_median = NA, pMMR_median = NA, stringsAsFactors = FALSE)
for (i in seq_along(score_list)) {
  s <- score_list[[i]]
  dm <- msi[[s]][msi$MSI_group == "dMMR (MSI)"]
  pm <- msi[[s]][msi$MSI_group == "pMMR (MSS)"]
  w <- wilcox.test(dm, pm)
  score_tab$p_value[i] <- signif(w$p.value, 3)
  score_tab$dMMR_median[i] <- round(median(dm, na.rm = TRUE), 3)
  score_tab$pMMR_median[i] <- round(median(pm, na.rm = TRUE), 3)
  cat(sprintf("%s: dMMR median=%.3f vs pMMR median=%.3f | P=%.4f\n",
      s, median(dm, na.rm=TRUE), median(pm, na.rm=TRUE), w$p.value))
}

# ===== 2. Correlations within each MMR group =====
cat("\n--- Correlations within MMR groups (Spearman) ---\n")
cor_tab <- data.frame(pair = c("FAP-ECM", "FAP-CLDN"), group = NA, rho = NA, P = NA, n = NA)
cor_ecm_d <- cor.test(msi$FAP_CAF[msi$MSI_group == "dMMR (MSI)"], msi$ECM_SDC4_CD44[msi$MSI_group == "dMMR (MSI)"], method = "spearman")
cor_ecm_p <- cor.test(msi$FAP_CAF[msi$MSI_group == "pMMR (MSS)"], msi$ECM_SDC4_CD44[msi$MSI_group == "pMMR (MSS)"], method = "spearman")
cor_cld_d <- cor.test(msi$FAP_CAF[msi$MSI_group == "dMMR (MSI)"], msi$CLDN_core[msi$MSI_group == "dMMR (MSI)"], method = "spearman")
cor_cld_p <- cor.test(msi$FAP_CAF[msi$MSI_group == "pMMR (MSS)"], msi$CLDN_core[msi$MSI_group == "pMMR (MSS)"], method = "spearman")
cat(sprintf("FAP-ECM dMMR: rho=%.3f P=%.4f | pMMR: rho=%.3f P=%.4f\n",
    cor_ecm_d$estimate, cor_ecm_d$p.value, cor_ecm_p$estimate, cor_ecm_p$p.value))
cat(sprintf("FAP-CLDN dMMR: rho=%.3f P=%.4f | pMMR: rho=%.3f P=%.4f\n",
    cor_cld_d$estimate, cor_cld_d$p.value, cor_cld_p$estimate, cor_cld_p$p.value))

# Fisher r-to-z test for correlation difference between groups
fz <- function(r1, n1, r2, n2) {
  z1 <- atanh(r1); z2 <- atanh(r2)
  se <- sqrt(1/(n1 - 3) + 1/(n2 - 3))
  z <- (z1 - z2) / se
  2 * pnorm(-abs(z))
}
n_d <- sum(msi$MSI_group == "dMMR (MSI)"); n_p <- sum(msi$MSI_group == "pMMR (MSS)")
p_diff_ecm <- fz(cor_ecm_d$estimate, n_d, cor_ecm_p$estimate, n_p)
p_diff_cld <- fz(cor_cld_d$estimate, n_d, cor_cld_p$estimate, n_p)
cat(sprintf("Correlation diff (Fisher z): FAP-ECM P=%.4f | FAP-CLDN P=%.4f\n", p_diff_ecm, p_diff_cld))

# ===== 3. OS survival by MMR (sanity check) =====
msi_s <- msi[!is.na(msi$os_delay) & msi$os_delay > 0 & !is.na(msi$MSI_group), ]
cat("\n--- OS by MMR ---\n")
cat("n =", nrow(msi_s), "| events:", sum(msi_s$os_event), "\n")
cat("  dMMR events:", sum(msi_s$os_event[msi_s$MSI_group == "dMMR (MSI)"]), "/", sum(msi_s$MSI_group == "dMMR (MSI)"), "\n")
cat("  pMMR events:", sum(msi_s$os_event[msi_s$MSI_group == "pMMR (MSS)"]), "/", sum(msi_s$MSI_group == "pMMR (MSS)"), "\n")
sd_msi <- survdiff(Surv(os_delay, os_event) ~ MSI_group, data = msi_s)
p_msi <- pchisq(sd_msi$chisq, 1, lower.tail = FALSE)
cox_msi <- coxph(Surv(os_delay, os_event) ~ MSI_group, data = msi_s)
cat(sprintf("Log-rank P=%.4f | Cox dMMR vs pMMR HR=%.2f (%.2f-%.2f) P=%.4f\n",
    p_msi, exp(coef(cox_msi)), exp(confint(cox_msi))[1], exp(confint(cox_msi))[2],
    summary(cox_msi)$coefficients[5]))

# ===== 4. MSI x T-stage distribution =====
msi$T_group <- ifelse(grepl("^T[12]$", msi$tnm_t), "T1-2", ifelse(grepl("^T[34]", msi$tnm_t), "T3-4", NA))
cat("\n--- MSI x T-stage ---\n")
print(table(msi$MSI_group, msi$T_group, useNA = "ifany"))
ft <- fisher.test(table(msi$MSI_group, msi$T_group))
cat(sprintf("Fisher P = %.4f\n", ft$p.value))

# ===== 5. MSI x N-stage =====
msi$N_pos <- ifelse(grepl("N[12]", msi$tnm_n), 1, ifelse(grepl("N0", msi$tnm_n), 0, NA))
cat("\n--- MSI x N-stage ---\n")
print(table(msi$MSI_group, msi$N_pos, useNA = "ifany"))
ftn <- fisher.test(table(msi$MSI_group, msi$N_pos))
cat(sprintf("Fisher P = %.4f\n", ftn$p.value))

# ===== 6. Plots =====
# 6.1 Score boxplot by MMR
plot_df <- tidyr::pivot_longer(msi[, c("sample", "MSI_group", "FAP_CAF", "ECM_SDC4_CD44", "CLDN_core")],
                               cols = c("FAP_CAF", "ECM_SDC4_CD44", "CLDN_core"),
                               names_to = "score", values_to = "value")
plot_df$score <- factor(plot_df$score, levels = c("FAP_CAF", "ECM_SDC4_CD44", "CLDN_core"),
                        labels = c("FAP-CAF", "ECM-SDC4/CD44", "CLDN core"))
p1 <- ggplot(plot_df, aes(x = MSI_group, y = value, fill = MSI_group)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.75) +
  stat_compare_means(aes(group = MSI_group), method = "wilcox.test", label = "p.format", size = 3.2) +
  facet_wrap(~ score, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("dMMR (MSI)" = "#E64B35", "pMMR (MSS)" = "#4DBBD5")) +
  labs(x = "", y = "z-mean score",
       title = sprintf("GSE39582: Stromal / tight-junction scores by MMR status (n=%d)", nrow(msi))) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(FIG_DIR, "A6_GSE39582_MSI_score_boxplot.png"), p1, width = 10, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A6_GSE39582_MSI_score_boxplot.tiff"), p1, width = 10, height = 6, dpi = 600, compression = "lzw")

# 6.2 KM by MMR
km_msi <- survfit(Surv(os_delay, os_event) ~ MSI_group, data = msi_s)
png(file.path(FIG_DIR, "A6_GSE39582_KM_MMR.png"), width = 9, height = 7, units = "in", res = 300)
plot(km_msi, col = c("#E64B35", "#4DBBD5"), lwd = 2, xlab = "Overall survival (months)",
     ylab = "Survival probability",
     main = sprintf("GSE39582: OS by MMR status (n=%d, log-rank P=%.4f)", nrow(msi_s), p_msi))
legend("topright", c("dMMR (MSI)", "pMMR (MSS)"), col = c("#E64B35", "#4DBBD5"), lwd = 2)
dev.off()

# 6.3 Correlation comparison plot
cor_plot <- data.frame(
  pair = rep(c("FAP-ECM", "FAP-CLDN"), each = 2),
  group = rep(c("dMMR (MSI)", "pMMR (MSS)"), 2),
  rho = c(cor_ecm_d$estimate, cor_ecm_p$estimate, cor_cld_d$estimate, cor_cld_p$estimate),
  stringsAsFactors = FALSE
)
p3 <- ggplot(cor_plot, aes(x = pair, y = rho, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.2f", rho)), position = position_dodge(0.75), vjust = -0.4, size = 3.6) +
  scale_fill_manual(values = c("dMMR (MSI)" = "#E64B35", "pMMR (MSS)" = "#4DBBD5")) +
  labs(x = "", y = "Spearman rho (FAP-CAF vs score)",
       title = "GSE39582: Correlation of FAP-CAF with ECM / CLDN by MMR status",
       subtitle = sprintf("Fisher z-test: FAP-ECM diff P=%.3f | FAP-CLDN diff P=%.3f", p_diff_ecm, p_diff_cld)) +
  theme_classic(base_size = 13) + ylim(-0.3, 1.0)
ggsave(file.path(FIG_DIR, "A6_GSE39582_MSI_correlation.png"), p3, width = 8, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A6_GSE39582_MSI_correlation.tiff"), p3, width = 8, height = 6, dpi = 600, compression = "lzw")

# ===== 7. Summary table =====
sum_tab <- data.frame(
  analysis = c("n dMMR / pMMR",
               "FAP_CAF dMMR vs pMMR (Wilcoxon)",
               "ECM_SDC4_CD44 dMMR vs pMMR (Wilcoxon)",
               "CLDN_core dMMR vs pMMR (Wilcoxon)",
               "FAP-ECM rho dMMR",
               "FAP-ECM rho pMMR",
               "FAP-ECM correlation diff (Fisher z)",
               "FAP-CLDN rho dMMR",
               "FAP-CLDN rho pMMR",
               "FAP-CLDN correlation diff (Fisher z)",
               "OS log-rank dMMR vs pMMR",
               "MSI x T-stage Fisher",
               "MSI x N-stage Fisher"),
  value = c(sprintf("%d / %d", sum(msi$MSI_group == "dMMR (MSI)"), sum(msi$MSI_group == "pMMR (MSS)")),
            sprintf("P=%.4f (median %.3f vs %.3f)", score_tab$p_value[1], score_tab$dMMR_median[1], score_tab$pMMR_median[1]),
            sprintf("P=%.4f (median %.3f vs %.3f)", score_tab$p_value[2], score_tab$dMMR_median[2], score_tab$pMMR_median[2]),
            sprintf("P=%.4f (median %.3f vs %.3f)", score_tab$p_value[3], score_tab$dMMR_median[3], score_tab$pMMR_median[3]),
            sprintf("rho=%.3f P=%.4f", cor_ecm_d$estimate, cor_ecm_d$p.value),
            sprintf("rho=%.3f P=%.4f", cor_ecm_p$estimate, cor_ecm_p$p.value),
            sprintf("P=%.4f", p_diff_ecm),
            sprintf("rho=%.3f P=%.4f", cor_cld_d$estimate, cor_cld_d$p.value),
            sprintf("rho=%.3f P=%.4f", cor_cld_p$estimate, cor_cld_p$p.value),
            sprintf("P=%.4f", p_diff_cld),
            sprintf("P=%.4f HR=%.2f", p_msi, exp(coef(cox_msi))),
            sprintf("P=%.4f", ft$p.value),
            sprintf("P=%.4f", ftn$p.value)),
  stringsAsFactors = FALSE
)
write.csv(sum_tab, file.path(TAB_DIR, "A6_MSI_summary.csv"), row.names = FALSE)
cat("\n=== [DONE] A6 MSI stratification complete ===\n")
print(sum_tab)
