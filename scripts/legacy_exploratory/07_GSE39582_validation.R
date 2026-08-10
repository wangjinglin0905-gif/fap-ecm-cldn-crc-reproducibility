#!/usr/bin/env Rscript
# =====================================================================
# 07_GSE39582_validation_FAP_ECM.R
# GSE39582 (n=585, Marisa 2013, GPL570) 外部验证
# 复算 AJCR 评分体系：FAP_CAF(17基因) / CLDN_core(CLDN1/2/4) / ECM-SDC4-CD44(6基因)
# 分析：分期分层相关 + FAP梯度 + OS/RFS 生存（KM + 多因素Cox）
# 输出：figures/FigGSE39582_*.png(tiff) + tables/GSE39582_validation_*.csv
# =====================================================================
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})

OUTD <- "./data/GSE39582"
FIGD <- "./output/figures"
TABD <- "./output/tables"
dir.create(FIGD, recursive = TRUE, showWarnings = FALSE)
dir.create(TABD, recursive = TRUE, showWarnings = FALSE)

# ---- 1. load probe-gene map ----
g2p <- read.csv(file.path(OUTD, "gene2probe.csv"), stringsAsFactors = FALSE)

# ---- 2. load expression + pheno ----
expr <- read.csv(file.path(OUTD, "GSE39582_expr.csv"), stringsAsFactors = FALSE,
                 check.names = FALSE)
pheno <- read.csv(file.path(OUTD, "GSE39582_pheno.csv"), stringsAsFactors = FALSE)
cat("expr:", nrow(expr), "x", ncol(expr) - 1, "\n")
cat("pheno:", nrow(pheno), "\n")

# ---- 3. probe -> gene (mean of probes per gene) ----
probe_cols <- colnames(expr)[-1]
collapse_to_gene <- function(expr_df, g2p, gene) {
  probes <- g2p$probe[g2p$gene == gene]
  probes <- intersect(probes, probe_cols)
  if (length(probes) == 0) return(rep(NA_real_, nrow(expr_df)))
  mat <- as.matrix(expr_df[, probes, drop = FALSE])
  storage.mode(mat) <- "numeric"
  if (length(probes) == 1) return(mat[, 1]) else return(rowMeans(mat, na.rm = TRUE))
}

genes_all <- unique(g2p$gene)
gene_mat <- sapply(genes_all, function(g) collapse_to_gene(expr, g2p, g))
rownames(gene_mat) <- expr$sample
cat("gene matrix:", nrow(gene_mat), "x", ncol(gene_mat), "\n")

# ---- 4. scores (same definition as AJCR/TCGA L0) ----
score_zmean <- function(mat, genes) {
  genes <- intersect(genes, colnames(mat))
  if (length(genes) < 3) stop("too few genes for score")
  z <- scale(mat[, genes, drop = FALSE])
  rowMeans(z, na.rm = TRUE)
}

fap_genes <- c("FAP","COL1A1","COL1A2","COL3A1","FN1","POSTN","THY1","PDPN",
               "TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
cldn_core <- c("CLDN1","CLDN2","CLDN4")
ecm_genes <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")

sc <- data.frame(
  sample = expr$sample,
  FAP_CAF   = score_zmean(gene_mat, fap_genes),
  CLDN_core = score_zmean(gene_mat, cldn_core),
  ECM_SDC4_CD44 = score_zmean(gene_mat, ecm_genes)
)

# ---- 5. merge pheno ----
m <- merge(sc, pheno, by = "sample", all.x = TRUE)
m$os_time <- as.numeric(m$os_delay_months)
m$os_ev   <- as.numeric(m$os_event)
m$rfs_time <- as.numeric(m$rfs_delay)
m$rfs_ev   <- as.numeric(m$rfs_event)
m$stage   <- as.numeric(m$tnm_stage)
m$age     <- as.numeric(m$age_at_diagnosis_year)
m$sex     <- m$Sex
m$stage_grp <- ifelse(m$stage == 1, "Stage I", ifelse(m$stage %in% 2:4, "Stage II-IV", NA))
cat("valid stage:", sum(!is.na(m$stage_grp)), " (I:", sum(m$stage_grp == "Stage I", na.rm=TRUE),
    " II-IV:", sum(m$stage_grp == "Stage II-IV", na.rm=TRUE), ")\n")

# ---- 6. stage-stratified correlations ----
spear <- function(x, y) {
  k <- complete.cases(x, y)
  ct <- suppressWarnings(cor.test(x[k], y[k], method = "spearman", exact = FALSE))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = sum(k))
}
cor_tab <- rbind(
  data.frame(comparison = "FAP_CAF_vs_CLDN_core", stage = "All", spear(m$FAP_CAF, m$CLDN_core)),
  data.frame(comparison = "FAP_CAF_vs_CLDN_core", stage = "Stage I", spear(m$FAP_CAF[m$stage_grp=="Stage I"], m$CLDN_core[m$stage_grp=="Stage I"])),
  data.frame(comparison = "FAP_CAF_vs_CLDN_core", stage = "Stage II-IV", spear(m$FAP_CAF[m$stage_grp=="Stage II-IV"], m$CLDN_core[m$stage_grp=="Stage II-IV"])),
  data.frame(comparison = "FAP_CAF_vs_ECM", stage = "All", spear(m$FAP_CAF, m$ECM_SDC4_CD44)),
  data.frame(comparison = "FAP_CAF_vs_ECM", stage = "Stage I", spear(m$FAP_CAF[m$stage_grp=="Stage I"], m$ECM_SDC4_CD44[m$stage_grp=="Stage I"])),
  data.frame(comparison = "FAP_CAF_vs_ECM", stage = "Stage II-IV", spear(m$FAP_CAF[m$stage_grp=="Stage II-IV"], m$ECM_SDC4_CD44[m$stage_grp=="Stage II-IV"]))
)
write.csv(cor_tab, file.path(TABD, "GSE39582_stage_stratified_correlations.csv"), row.names = FALSE)
cat("\n=== Stage-stratified correlations ===\n"); print(cor_tab, row.names = FALSE)

# ---- 7. FAP_CAF high/low KM (OS + RFS, whole cohort + Stage I) ----
km_panel <- function(dat, time, ev, score, label, fname_png, fname_tiff) {
  d <- dat[complete.cases(dat[[time]], dat[[ev]], dat[[score]]), ]
  d$grp <- ifelse(d[[score]] > median(d[[score]]), "FAP-high", "FAP-low")
  s <- Surv(d[[time]], d[[ev]])
  fit <- survfit(s ~ grp, data = d)
  sd <- survdiff(s ~ grp, data = d)
  p_lr <- pchisq(sd$chisq, 1, lower.tail = FALSE)
  cox <- coxph(s ~ d[[score]])
  hr <- exp(coef(cox)); ci <- exp(confint(cox))
  p_cox <- summary(cox)$coefficients[5]
  do_km_plot <- function(file, res) {
    png(file, width = 7, height = 5.5, units = "in", res = res)
    par(mar = c(4.5, 4.5, 3, 1.2))
    plot(fit, col = c("#3C5488", "#E64B35"), lwd = 2.2,
         xlab = paste0(label, " (months)"), ylab = "Survival probability",
         main = paste0(label, " by FAP-CAF score - GSE39582 (n=", nrow(d), ")"),
         mark.time = TRUE)
    legend("bottomleft",
           legend = c(paste0("FAP-low (n=", sum(d$grp == "FAP-low"), ")"),
                      paste0("FAP-high (n=", sum(d$grp == "FAP-high"), ")")),
           col = c("#3C5488", "#E64B35"), lwd = 2.2, bty = "n")
    text(par("usr")[2] * 0.55, 0.15,
         sprintf("log-rank P = %.4f\nHR = %.3f [%.3f-%.3f]",
                 p_lr, hr, ci[1], ci[2]), cex = 0.85, adj = 0)
    dev.off()
  }
  do_km_plot(fname_png, 300)
  do_km_plot(fname_tiff, 600)
  data.frame(endpoint = label, n = nrow(d), events = sum(d[[ev]]),
             logrank_p = p_lr, HR = hr, lower = ci[1], upper = ci[2], cox_p = p_cox)
}

km_res <- rbind(
  km_panel(m, "os_time", "os_ev", "FAP_CAF", "OS", file.path(FIGD, "FigGSE39582_OS_KM.png"), file.path(FIGD, "FigGSE39582_OS_KM.tiff")),
  km_panel(m, "rfs_time", "rfs_ev", "FAP_CAF", "RFS", file.path(FIGD, "FigGSE39582_RFS_KM.png"), file.path(FIGD, "FigGSE39582_RFS_KM.tiff"))
)
# Stage I subset (RFS preferred - more events)
mI <- m[!is.na(m$stage_grp) & m$stage_grp == "Stage I", ]
km_res <- rbind(km_res,
  km_panel(mI, "os_time", "os_ev", "FAP_CAF", "OS_StageI", file.path(FIGD, "FigGSE39582_OS_StageI_KM.png"), file.path(FIGD, "FigGSE39582_OS_StageI_KM.tiff")),
  km_panel(mI, "rfs_time", "rfs_ev", "FAP_CAF", "RFS_StageI", file.path(FIGD, "FigGSE39582_RFS_StageI_KM.png"), file.path(FIGD, "FigGSE39582_RFS_StageI_KM.tiff"))
)
write.csv(km_res, file.path(TABD, "GSE39582_KM_results.csv"), row.names = FALSE)
cat("\n=== KM results ===\n"); print(km_res, row.names = FALSE)

# ---- 8. multivariable Cox (OS: FAP_CAF + age + sex + stage[continuous]) ----
s <- m[complete.cases(m$os_time, m$os_ev, m$FAP_CAF, m$age, m$sex, m$stage), ]
s$stage_cont <- as.numeric(s$stage)
mul <- coxph(Surv(os_time, os_ev) ~ FAP_CAF + age + sex + stage_cont, data = s)
co <- summary(mul)$coefficients; ci <- summary(mul)$conf.int
mul_tab <- data.frame(term = rownames(co),
                      HR = ci[, "exp(coef)"], lower = ci[, "lower .95"],
                      upper = ci[, "upper .95"], p = co[, "Pr(>|z|)"], n = nrow(s))
write.csv(mul_tab, file.path(TABD, "GSE39582_multivariable_Cox_OS.csv"), row.names = FALSE)
cat("\n=== Multivariable Cox (OS, stage continuous) ===\n"); print(mul_tab, row.names = FALSE)

# RFS multivariable
s2 <- m[complete.cases(m$rfs_time, m$rfs_ev, m$FAP_CAF, m$age, m$sex, m$stage), ]
s2$stage_cont <- as.numeric(s2$stage)
mul2 <- coxph(Surv(rfs_time, rfs_ev) ~ FAP_CAF + age + sex + stage_cont, data = s2)
co2 <- summary(mul2)$coefficients; ci2 <- summary(mul2)$conf.int
mul2_tab <- data.frame(term = rownames(co2),
                       HR = ci2[, "exp(coef)"], lower = ci2[, "lower .95"],
                       upper = ci2[, "upper .95"], p = co2[, "Pr(>|z|)"], n = nrow(s2))
write.csv(mul2_tab, file.path(TABD, "GSE39582_multivariable_Cox_RFS.csv"), row.names = FALSE)
cat("\n=== Multivariable Cox (RFS, stage continuous) ===\n"); print(mul2_tab, row.names = FALSE)

# ---- 9. FAP gradient vs ECM (correlation fig) ----
p <- ggplot(m, aes(FAP_CAF, ECM_SDC4_CD44)) +
  geom_point(alpha = 0.4, color = "#3C5488") +
  geom_smooth(method = "lm", se = TRUE, color = "#E64B35") +
  labs(title = "FAP-CAF vs ECM-SDC4/CD44 - GSE39582",
       x = "FAP-CAF score", y = "ECM-SDC4/CD44 score") +
  theme_classic()
ggsave(file.path(FIGD, "FigGSE39582_FAP_ECM_corr.png"), p, width = 6, height = 5, dpi = 300)
ggsave(file.path(FIGD, "FigGSE39582_FAP_ECM_corr.tiff"), p, width = 6, height = 5, dpi = 600, compression = "lzw")

cat("\n[DONE] GSE39582 validation complete\n")
