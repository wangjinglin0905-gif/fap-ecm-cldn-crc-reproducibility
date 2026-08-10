# =============================================================================
# A1: TCGA-COAD 分期分层 bulk 分析
# 检验：FAP 间质-ECM-SDC4/CD44 轴活性在 I 期 vs II-IV 期 CRC 中的差异
# =============================================================================
source(file.path(dirname(sys.frame(1)$ofile) %||% ".", "00_config.R"), 
       local = FALSE, echo = FALSE, chdir = TRUE)

# ... fallback when sys.frame doesn't give the script path
tryCatch({
  source(file.path(PROJ_ROOT, "scripts", "00_config.R"))
}, error = function(e) {
  message("Running A1 in standalone mode; loading packages directly")
  for (pkg in c("TCGAbiolinks","SummarizedExperiment","DESeq2","edgeR",
                "GSVA","survival","survminer","rstatix","ggplot2","ggpubr",
                "pheatmap","dplyr","tidyr","tibble","stringr")) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  PROJ_ROOT <- "."
  OUT_DIR   <- file.path(PROJ_ROOT, "output")
  FIG_DIR   <- file.path(OUT_DIR, "figures")
  TAB_DIR   <- file.path(OUT_DIR, "tables")
  DATA_DIR  <- file.path(PROJ_ROOT, "data")
  FDR_THRESHOLD   <- 0.05
  LOGFC_THRESHOLD <- 0.5
  FAP_CAF_GENES <- c("FAP","COL1A1","COL1A2","FN1","ACTA2")
  ECM_SDC4_CD44_GENES <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")
  PROLIFERATION_CONFOUNDER <- "MKI67"
  CLAUDIN_GENES <- c("CLDN1","CLDN2","CLDN3","CLDN4","CLDN7",
                     "CLDN8","CLDN12","CLDN15","CLDN18","CLDN23")
  STAGE_EARLY <- "I"
  STAGE_LATE  <- c("II","III","IV")
})

message("=== A1: TCGA-COAD Staging-Stratified Bulk Analysis ===")
message("Output directory: ", OUT_DIR)

# ---- Step 0: Download / Load TCGA-COAD data ----
rda_path <- file.path(DATA_DIR, "tcga_coad_query.rda")
count_path <- file.path(DATA_DIR, "tcga_coad_counts.rds")
clin_path  <- file.path(DATA_DIR, "tcga_coad_clinical.rds")

if (!file.exists(count_path)) {
  message("[A1-1] Downloading TCGA-COAD RNA-seq (STAR-counts) via GDC...")
  query <- GDCquery(
    project = "TCGA-COAD",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  saveRDS(query, rda_path)
  GDCdownload(query, method = "api", directory = file.path(DATA_DIR, "GDCdata"))
  
  message("[A1-2] Preparing count matrix and clinical data...")
  tcga <- GDCprepare(query, directory = file.path(DATA_DIR, "GDCdata"))
  
  # Extract counts (unstranded)
  counts <- assay(tcga, "unstranded")
  rownames(counts) <- rowData(tcga)$gene_name
  
  # Extract clinical
  clinical <- as.data.frame(colData(tcga))
  
  saveRDS(counts, count_path)
  saveRDS(clinical, clin_path)
  message("[A1] Data saved.")
} else {
  message("[A1] Loading cached TCGA-COAD data...")
  counts   <- readRDS(count_path)
  clinical <- readRDS(clin_path)
}

# ---- Step 1: Stage-stratify clinical ----
message("[A1-3] Stratifying by AJCC stage (I vs II-IV)...")
clinical$stage_group <- ifelse(
  clinical$ajcc_pathologic_stage %in% STAGE_EARLY, "Early (I)",
  ifelse(clinical$ajcc_pathologic_stage %in% STAGE_LATE, "Late (II-IV)", NA)
)
clinical <- clinical[!is.na(clinical$stage_group), ]

table_stage <- table(clinical$stage_group)
message("  Early (I)  : ", table_stage["Early (I)"])
message("  Late (II-IV): ", table_stage["Late (II-IV)"])

# ---- Step 2: DESeq2 differential expression (I vs II-IV) ----
message("[A1-4] Running DESeq2: I vs II-IV...")
keep_samples <- clinical$barcode[clinical$barcode %in% colnames(counts)]
clinical_filtered <- clinical[clinical$barcode %in% keep_samples, ]
counts_filtered <- counts[, keep_samples, drop = FALSE]

# Remove low-expression genes
keep_genes <- rowSums(counts_filtered >= 10) >= 10
counts_filtered <- counts_filtered[keep_genes, ]
counts_filtered <- round(counts_filtered)  # ensure integer

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData   = clinical_filtered,
  design    = ~ stage_group
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("stage_group", "Late (II-IV)", "Early (I)"), alpha = FDR_THRESHOLD)
res <- lfcShrink(dds, contrast = c("stage_group", "Late (II-IV)", "Early (I)"), 
                 res = res, type = "apeglm")
res_df <- as.data.frame(res[order(res$padj), ])
res_df$gene <- rownames(res_df)

# Save
write.csv(res_df, file.path(TAB_DIR, "A1_DESeq2_I_vs_II-IV.csv"), row.names = FALSE)

# Volcano for key genes
key_genes <- unique(c(FAP_CAF_GENES, ECM_SDC4_CD44_GENES, CLAUDIN_GENES, PROLIFERATION_CONFOUNDER))
res_df$is_key <- res_df$gene %in% key_genes
res_df$neg_log10_padj <- -log10(res_df$padj + 1e-300)

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = is_key, size = is_key), alpha = 0.6) +
  geom_vline(xintercept = c(-LOGFC_THRESHOLD, LOGFC_THRESHOLD), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed", alpha = 0.3) +
  geom_text(data = subset(res_df, is_key), aes(label = gene), size = 3, vjust = -0.5) +
  scale_color_manual(values = c("grey70", "red")) +
  scale_size_manual(values = c(0.5, 2)) +
  labs(x = "log2 Fold Change (Late vs Early)",
       y = expression(-log[10](adjusted~P)),
       title = "TCGA-COAD: I vs II-IV Differential Expression") +
  theme_minimal()
ggsave(file.path(FIG_DIR, "A1_volcano_I_vs_II-IV.png"), p_volcano, width = 10, height = 8, dpi = 300)

message("[A1-4] DESeq2 complete. Significant genes: ", sum(res_df$padj < 0.05, na.rm = TRUE))

# ---- Step 3: ssGSEA scoring ----
message("[A1-5] Calculating ssGSEA scores...")
vst <- vst(dds, blind = FALSE)
expr_vst <- assay(vst)
rownames(expr_vst) <- rownames(counts_filtered)

gene_sets <- list(
  "FAP_CAF_score"  = FAP_CAF_GENES,
  "ECM_SDC4_CD44"  = ECM_SDC4_CD44_GENES
)

ssgsea_res <- gsva(as.matrix(expr_vst), gene_sets, method = "ssgsea",
                   kcdf = "Gaussian", min.sz = 3, max.sz = 500)

# Merge into clinical
clinical_filtered$FAP_CAF_score <- as.numeric(ssgsea_res["FAP_CAF_score", ])
clinical_filtered$ECM_SDC4_CD44_score <- as.numeric(ssgsea_res["ECM_SDC4_CD44", ])

# ESTIMATE stromal score (simple proxy: mean of COL1A1 + COL1A2)
stromal_proxy <- colMeans(expr_vst[c("COL1A1", "COL1A2"), , drop = FALSE])
clinical_filtered$stromal_proxy <- as.numeric(stromal_proxy)

# ---- Step 4: Correlations (stage-stratified) ----
message("[A1-6] Spearman correlations: FAP-CAF vs ECM-SDC4/CD44 (stage-stratified)...")

# Overall
cor_overall <- cor.test(clinical_filtered$FAP_CAF_score,
                        clinical_filtered$ECM_SDC4_CD44_score, method = "spearman")

# Stage-stratified
cor_early <- cor.test(
  clinical_filtered$FAP_CAF_score[clinical_filtered$stage_group == "Early (I)"],
  clinical_filtered$ECM_SDC4_CD44_score[clinical_filtered$stage_group == "Early (I)"],
  method = "spearman")
cor_late <- cor.test(
  clinical_filtered$FAP_CAF_score[clinical_filtered$stage_group == "Late (II-IV)"],
  clinical_filtered$ECM_SDC4_CD44_score[clinical_filtered$stage_group == "Late (II-IV)"],
  method = "spearman")

# Partial correlation correcting for stromal proxy
clinical_filtered$resid_FAP <- resid(lm(FAP_CAF_score ~ stromal_proxy, data = clinical_filtered))
clinical_filtered$resid_ECM  <- resid(lm(ECM_SDC4_CD44_score ~ stromal_proxy, data = clinical_filtered))
cor_partial_early <- cor.test(
  clinical_filtered$resid_FAP[clinical_filtered$stage_group == "Early (I)"],
  clinical_filtered$resid_ECM[clinical_filtered$stage_group == "Early (I)"],
  method = "spearman")

cor_results <- data.frame(
  Comparison = c("Overall", "Early (I)", "Late (II-IV)", "Early (I) - stromal corrected"),
  rho = c(cor_overall$estimate, cor_early$estimate, cor_late$estimate, cor_partial_early$estimate),
  p_value = c(cor_overall$p.value, cor_early$p.value, cor_late$p.value, cor_partial_early$p.value)
)
write.csv(cor_results, file.path(TAB_DIR, "A1_correlations.csv"), row.names = FALSE)

# Scatter plot
p_scatter <- ggplot(clinical_filtered, aes(x = FAP_CAF_score, y = ECM_SDC4_CD44_score, color = stage_group)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~stage_group) +
  stat_cor(method = "spearman", label.x.npc = "left") +
  labs(x = "FAP-CAF Score (ssGSEA)", y = "ECM-SDC4/CD44 Score (ssGSEA)",
       title = "FAP Stromal Program vs ECM-SDC4/CD44 Axis (Stage-Stratified)") +
  theme_minimal()
ggsave(file.path(FIG_DIR, "A1_correlation_stage_stratified.png"), p_scatter, width = 12, height = 6, dpi = 300)

# ---- Step 5: Survival Analysis ----
message("[A1-7] Kaplan-Meier + Cox regression...")

# Add survival data
clinical_filtered$os_time <- as.numeric(clinical_filtered$days_to_death)
clinical_filtered$os_event <- ifelse(clinical_filtered$vital_status == "Dead", 1, 0)
# Use PFI if OS events too few
clinical_filtered$pfi_time <- as.numeric(clinical_filtered$days_to_new_tumor_event_after_initial_treatment)
clinical_filtered$pfi_event <- ifelse(!is.na(clinical_filtered$pfi_time), 1, 0)

# Dichotomize by median FAP-CAF score
clinical_filtered$FAP_high <- clinical_filtered$FAP_CAF_score > median(clinical_filtered$FAP_CAF_score, na.rm = TRUE)

# KM (all stages)
fit_all <- survfit(Surv(os_time, os_event) ~ FAP_high, data = clinical_filtered)
p_km_all <- ggsurvplot(fit_all, data = clinical_filtered, pval = TRUE,
                        title = "All Stages: OS by FAP-CAF Score", risk.table = TRUE)

# KM (early stage only)
clinical_early <- clinical_filtered[clinical_filtered$stage_group == "Early (I)", ]
clinical_early$FAP_high_early <- clinical_early$FAP_CAF_score > median(clinical_early$FAP_CAF_score, na.rm = TRUE)
fit_early <- survfit(Surv(os_time, os_event) ~ FAP_high_early, data = clinical_early)
p_km_early <- ggsurvplot(fit_early, data = clinical_early, pval = TRUE,
                          title = "Early Stage (I): OS by FAP-CAF Score", risk.table = TRUE)

# Multivariate Cox
cox_multi <- coxph(Surv(os_time, os_event) ~ FAP_CAF_score + 
                     age_at_index + gender + ajcc_pathologic_stage + 
                     stromal_proxy,
                   data = clinical_filtered)
cox_ph_test <- cox.zph(cox_multi)
message("  PH assumption test p = ", format(min(cox_ph_test$table[, "p"]), digits = 3))

cox_summary <- summary(cox_multi)
cox_table <- as.data.frame(cox_summary$coefficients)
cox_table$Variable <- rownames(cox_table)
write.csv(cox_table, file.path(TAB_DIR, "A1_multivariate_cox.csv"), row.names = FALSE)

# Forest plot
p_forest <- ggforest(cox_multi, data = clinical_filtered,
                     main = "Multivariate Cox: Overall Survival")
ggsave(file.path(FIG_DIR, "A1_cox_forest.png"), p_forest, width = 10, height = 6, dpi = 300)

# ---- Step 6: Claudin individual analysis ----
message("[A1-8] Claudin family individual expression...")
claudin_in_data <- CLAUDIN_GENES[CLAUDIN_GENES %in% rownames(res_df)]

claudin_table <- data.frame(
  Gene = character(),
  log2FC = numeric(),
  pvalue = numeric(),
  padj = numeric(),
  Early_mean = numeric(),
  Late_mean = numeric(),
  stringsAsFactors = FALSE
)

for (gene in claudin_in_data) {
  claudin_table <- rbind(claudin_table, data.frame(
    Gene   = gene,
    log2FC = res_df[gene, "log2FoldChange"],
    pvalue = res_df[gene, "pvalue"],
    padj   = res_df[gene, "padj"],
    Early_mean = mean(expr_vst[gene, clinical_filtered$stage_group == "Early (I)"]),
    Late_mean  = mean(expr_vst[gene, clinical_filtered$stage_group == "Late (II-IV)"]),
    stringsAsFactors = FALSE
  ))
}
write.csv(claudin_table, file.path(TAB_DIR, "A1_claudin_individual.csv"), row.names = FALSE)

# ---- Step 7: Output summary ----
message("\n========================================")
message("A1 Analysis Complete!")
message("========================================")
message("Figures saved to: ", FIG_DIR)
message("  - A1_volcano_I_vs_II-IV.png")
message("  - A1_correlation_stage_stratified.png")
message("  - A1_cox_forest.png")
message("Tables saved to: ", TAB_DIR)
message("  - A1_DESeq2_I_vs_II-IV.csv")
message("  - A1_correlations.csv")
message("  - A1_multivariate_cox.csv")
message("  - A1_claudin_individual.csv")
message("========================================")

# ---- Key findings summary ----
cat("\n--- A1 Key Findings ---\n")
cat("FAP log2FC (Late vs Early): ", signif(res_df["FAP", "log2FoldChange"], 3), 
    ", padj = ", signif(res_df["FAP", "padj"], 3), "\n")
cat("SDC4 log2FC (Late vs Early): ", 
    signif(res_df["SDC4", "log2FoldChange"], 3), 
    ", padj = ", signif(res_df["SDC4", "padj"], 3), "\n")
cat("Overall correlation rho: ", signif(cor_overall$estimate, 3),
    ", P = ", signif(cor_overall$p.value, 3), "\n")
cat("Early (I) correlation rho: ", signif(cor_early$estimate, 3),
    ", P = ", signif(cor_early$p.value, 3), "\n")
cat("Late (II-IV) correlation rho: ", signif(cor_late$estimate, 3),
    ", P = ", signif(cor_late$p.value, 3), "\n")
cat("Early (I) partial rho (stroma-corrected): ", signif(cor_partial_early$estimate, 3),
    ", P = ", signif(cor_partial_early$p.value, 3), "\n")

saveRDS(clinical_filtered, file.path(DATA_DIR, "A1_clinical_scored.rds"))
saveRDS(res_df, file.path(DATA_DIR, "A1_DESeq2_results.rds"))
message("Data saved. A1 done.")
