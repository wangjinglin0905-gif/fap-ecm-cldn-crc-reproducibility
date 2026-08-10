# =============================================================================
# A4: Claudin 家族成员逐个分析
# 扩展 A1 的 Claudin 输出 + CPTAC 蛋白验证 + GSE39582 外部验证
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggpubr)
  library(survival); library(survminer); library(corrplot)
})

PROJ_ROOT <- "."
OUT_DIR   <- file.path(PROJ_ROOT, "output"); FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR   <- file.path(OUT_DIR, "tables"); DATA_DIR <- file.path(PROJ_ROOT, "data")

CLAUDIN_GENES <- c("CLDN1","CLDN2","CLDN3","CLDN4","CLDN7",
                   "CLDN8","CLDN12","CLDN15","CLDN18","CLDN23")

message("=== A4: Claudin Individual Analysis ===")

# ---- Load A1 results ----
clinical <- readRDS(file.path(DATA_DIR, "A1_clinical_scored.rds"))
res_df   <- readRDS(file.path(DATA_DIR, "A1_DESeq2_results.rds"))

# ---- Part 1: Correlation matrix: Claudin x FAP_CAF_score (stage-stratified) ----
message("[A4-1] Claudin-FAP stromal correlation matrix...")

claudin_expr <- clinical[, CLAUDIN_GENES[CLAUDIN_GENES %in% colnames(clinical)]]
claudin_expr <- as.matrix(claudin_expr)

cor_matrix_early <- cor(claudin_expr[clinical$stage_group == "Early (I)", ],
                        clinical$FAP_CAF_score[clinical$stage_group == "Early (I)"],
                        method = "spearman", use = "complete.obs")
cor_matrix_late  <- cor(claudin_expr[clinical$stage_group == "Late (II-IV)", ],
                        clinical$FAP_CAF_score[clinical$stage_group == "Late (II-IV)"],
                        method = "spearman", use = "complete.obs")

claudin_summary <- data.frame(
  Gene = CLAUDIN_GENES[CLAUDIN_GENES %in% rownames(cor_matrix_early)],
  rho_early = cor_matrix_early[,1],
  rho_late  = cor_matrix_late[,1],
  rho_delta = cor_matrix_early[,1] - cor_matrix_late[,1]
)
claudin_summary <- claudin_summary[order(abs(claudin_summary$rho_delta), decreasing = TRUE), ]

# Merge with DESeq2 results
des_in_common <- res_df[res_df$gene %in% claudin_summary$Gene, c("gene","log2FoldChange","padj")]
colnames(des_in_common)[1] <- "Gene"
claudin_summary <- merge(claudin_summary, des_in_common, by = "Gene", all.x = TRUE)

write.csv(claudin_summary, file.path(TAB_DIR, "A4_claudin_fap_correlation_stage_stratified.csv"),
          row.names = FALSE)

# Waterfall plot
p_waterfall <- ggplot(claudin_summary, aes(x = reorder(Gene, rho_delta), y = rho_delta, fill = rho_delta > 0)) +
  geom_col() +
  geom_text(aes(label = sprintf("%.3f", rho_delta)), hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  coord_flip() +
  labs(x = "Claudin Gene", y = expression(Delta*rho ~ "(Early - Late)"),
       title = "Claudin-FAP Correlation: Early vs Late CRC",
       subtitle = expression("Positive " * Delta * rho * " = stronger correlation in early-stage CRC")) +
  theme_minimal()
ggsave(file.path(FIG_DIR, "A4_claudin_waterfall.png"), p_waterfall, width = 8, height = 6, dpi = 300)

# ---- Part 2: Survival analysis per Claudin (early stage) ----
message("[A4-2] Survival analysis per Claudin (early stage only)...")

clinical_early <- clinical[clinical$stage_group == "Early (I)", ]
surv_results <- data.frame()

for (gene in CLAUDIN_GENES) {
  if (!gene %in% colnames(clinical_early)) next
  val <- clinical_early[[gene]]
  med <- median(val, na.rm = TRUE)
  clinical_early$claudin_high <- val > med
  fit <- survfit(Surv(os_time, os_event) ~ claudin_high, data = clinical_early)
  diff <- survdiff(Surv(os_time, os_event) ~ claudin_high, data = clinical_early)
  surv_results <- rbind(surv_results, data.frame(
    Gene = gene, chi2 = diff$chisq, p = 1 - pchisq(diff$chisq, 1)
  ))
}
write.csv(surv_results, file.path(TAB_DIR, "A4_claudin_survival_early.csv"), row.names = FALSE)

# ---- Part 3: CPTAC protein validation (placeholder) ----
message("[A4-3] CPTAC protein validation — deferred to dedicated script (04b_A4_CPTAC.R)")
message("CPTAC requires downloading from https://proteomics.cancer.gov/ or cBioPortal.")

# ---- Part 4: GSE39582 external validation (placeholder) ----
message("[A4-4] GSE39582 external validation — deferred to dedicated script (04c_A4_GSE39582.R)")
message("GSE39582 requires downloading from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39582).")

# ---- Part 5: Heatmap summary ----
message("[A4-5] Generating Claudin-FAP heatmap...")

heat_data <- claudin_summary[, c("rho_early", "rho_late", "rho_delta")]
rownames(heat_data) <- claudin_summary$Gene
p_heat <- pheatmap(heat_data, cluster_rows = TRUE, cluster_cols = FALSE,
                   display_numbers = TRUE, number_format = "%.3f",
                   color = colorRampPalette(c("blue", "white", "red"))(100),
                   main = "Claudin-FAP Stromal Program Correlation\nStage-Stratified",
                   filename = file.path(FIG_DIR, "A4_claudin_heatmap.png"),
                   width = 8, height = 5)

message("A4 complete. Key finding: check A4_claudin_fap_correlation_stage_stratified.csv")
message("Positive delta_rho → Claudin member more strongly associated with FAP in early CRC")
