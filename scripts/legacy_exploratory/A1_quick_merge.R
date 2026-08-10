# Quick A1 analysis: merge TCGA scores with clinical, stage-stratify, output summary
suppressPackageStartupMessages({
  library(dplyr); library(TCGAbiolinks)
})

# Read pre-computed scores
tcga <- read.csv("legacy_external/tcga_recomputed_scores.csv")
tcga$submitter_id <- substr(tcga$patient_id, 1, 12)
cat("TCGA scores:", nrow(tcga), "rows\n")

# Get clinical
clinical <- GDCquery_clinic(project = "TCGA-COAD", type = "clinical")
cat("Clinical:", nrow(clinical), "rows\n")

# Merge
merged <- merge(tcga, clinical, by = "submitter_id", all.x = TRUE)
cat("Merged:", nrow(merged), "rows\n")

# Stage
cat("\nStage distribution:\n")
print(table(merged$ajcc_pathologic_stage, useNA = "ifany"))

# Stage group
merged$stage_group <- ifelse(merged$ajcc_pathologic_stage == "Stage I", "Early (I)",
  ifelse(merged$ajcc_pathologic_stage %in% c("Stage II", "Stage III", "Stage IV"), "Late (II-IV)", NA))
merged <- merged[!is.na(merged$stage_group), ]
cat("\nN per group:", table(merged$stage_group), "\n")

# Summary stats
merged %>%
  filter(!is.na(ajcc_pathologic_stage)) %>%
  group_by(ajcc_pathologic_stage) %>%
  summarise(
    n = n(),
    FAP_CAF_mean = mean(FAP_CAF, na.rm = TRUE),
    CLDN_core_mean = mean(CLDN_core, na.rm = TRUE),
    .groups = "drop"
  ) %>% print()

# Wilcoxon: FAP_CAF I vs II-IV
early <- merged$FAP_CAF[merged$stage_group == "Early (I)"]
late  <- merged$FAP_CAF[merged$stage_group == "Late (II-IV)"]
wt <- wilcox.test(early, late)
cat(sprintf("\nFAP_CAF: Early mean=%.3f, Late mean=%.3f, Wilcoxon P=%.4f\n",
            mean(early), mean(late), wt$p.value))

# Spearman: FAP vs CLDN, by stage
ce <- cor.test(early, merged$CLDN_core[merged$stage_group == "Early (I)"], method = "spearman")
cl <- cor.test(late,  merged$CLDN_core[merged$stage_group == "Late (II-IV)"], method = "spearman")
cat(sprintf("FAP-CAF vs CLDN_core rho: Early(I)=%.3f P=%.4f, Late(II-IV)=%.3f P=%.4f\n",
            ce$estimate, ce$p.value, cl$estimate, cl$p.value))

# Survival
merged$os_time <- as.numeric(merged$days_to_death)
merged$os_event <- ifelse(merged$vital_status == "Dead", 1, 0)
cat(sprintf("Survival: %d events / %d patients\n",
            sum(merged$os_event, na.rm = TRUE), nrow(merged)))

# Save
out <- "./data/A1_tcga_coad_merged.csv"
write.csv(merged, out, row.names = FALSE)
cat("Saved to", out, "\n")
