options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(lme4)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript domain_marker_check.R <GSE132465_seurat.rds> <output_directory>")
}
input_rds <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_text <- file.path(out_dir, "domain_tumor_marker_check.txt")

obj <- readRDS(input_rds)
meta <- obj@meta.data
counts <- LayerData(obj[["RNA"]], layer = "counts")
fib_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3")
cells <- rownames(meta)[meta$Class == "Tumor" & meta$Cell_subtype %in% fib_subtypes]
markers <- c("CDKN2A", "CDKN2B", "CDKN1A", "LMNB1", "MKI67")

d <- data.frame(
  cell = cells,
  Patient = factor(meta[cells, "Patient"]),
  Cell_subtype = factor(meta[cells, "Cell_subtype"]),
  log_nCount = log1p(meta[cells, "nCount_RNA"]),
  log_nFeature = log1p(meta[cells, "nFeature_RNA"]),
  FAP_status = factor(ifelse(as.numeric(counts["FAP", cells] > 0), "FAP+", "FAP-"), levels = c("FAP-", "FAP+")),
  stringsAsFactors = FALSE
)

for (g in markers) d[[g]] <- as.numeric(counts[g, cells] > 0)

sink(out_text)
cat("Tumor-only fibroblast-lineage marker audit\n")
cat("n cells:", nrow(d), "; patients:", nlevels(d$Patient), "\n")
print(table(d$FAP_status))

for (g in markers) {
  cat("\n###", g, "\n")
  tab <- aggregate(d[[g]], list(FAP_status = d$FAP_status), mean)
  names(tab)[2] <- "detection_rate"
  print(tab)

  per_patient <- aggregate(d[[g]], list(Patient = d$Patient, FAP_status = d$FAP_status), mean)
  names(per_patient)[3] <- "rate"
  plus <- per_patient[per_patient$FAP_status == "FAP+", c("Patient", "rate")]
  minus <- per_patient[per_patient$FAP_status == "FAP-", c("Patient", "rate")]
  names(plus)[2] <- "plus"
  names(minus)[2] <- "minus"
  paired <- merge(plus, minus, by = "Patient")
  wt <- suppressWarnings(wilcox.test(paired$plus, paired$minus, paired = TRUE, exact = FALSE))
  cat("patient paired n:", nrow(paired), "; median difference:", median(paired$plus - paired$minus),
      "; positive/negative/zero:", sum(paired$plus > paired$minus), "/", sum(paired$plus < paired$minus), "/", sum(paired$plus == paired$minus),
      "; P:", wt$p.value, "\n")

  form <- as.formula(paste0(g, " ~ FAP_status + log_nCount + Cell_subtype + (1 | Patient)"))
  fit <- glmer(form, data = d, family = binomial(), control = glmerControl(optimizer = "bobyqa"))
  cf <- coef(summary(fit))["FAP_statusFAP+", ]
  cat("library-depth/subtype-adjusted GLMM: logOR =", cf[["Estimate"]], "; OR =", exp(cf[["Estimate"]]),
      "; 95% CI =", exp(cf[["Estimate"]] - 1.96 * cf[["Std. Error"]]), "to", exp(cf[["Estimate"]] + 1.96 * cf[["Std. Error"]]),
      "; P =", cf[["Pr(>|z|)"]], "; singular =", isSingular(fit), "\n")
}

score_file <- file.path(out_dir, "cell_scores_tumor_fibroblast_epithelial.csv")
scores <- read.csv(score_file, check.names = FALSE)
scores <- scores[scores$compartment == "Fibroblast", ]
score_df <- merge(
  scores,
  d[, c("cell", "log_nCount", "log_nFeature")],
  by = "cell",
  all.x = TRUE,
  sort = FALSE
)
score_df$Patient <- factor(score_df$Patient)
score_df$FAP_status <- factor(score_df$FAP_status, levels = c("FAP-", "FAP+"))
score_df$Cell_subtype <- factor(score_df$Cell_subtype)

cat("\n### SenMayo score, tumour-only fibroblasts with technical/subtype adjustment\n")
score_fit <- lmerTest::lmer(
  SenMayo_zmean ~ FAP_status + MKI67_z + log_nCount + Cell_subtype + (1 | Patient),
  data = score_df,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
print(coef(summary(score_fit)))
cat("singular:", isSingular(score_fit), "\n")

cat("\n### SenMayo score, tumour-only fibroblasts with depth adjustment only\n")
score_fit_depth <- lmerTest::lmer(
  SenMayo_zmean ~ FAP_status + MKI67_z + log_nCount + (1 | Patient),
  data = score_df,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
print(coef(summary(score_fit_depth)))

cat("\n### SenMayo score, tumour-only fibroblasts with subtype adjustment only\n")
score_fit_subtype <- lmerTest::lmer(
  SenMayo_zmean ~ FAP_status + MKI67_z + Cell_subtype + (1 | Patient),
  data = score_df,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
print(coef(summary(score_fit_subtype)))

cat("\n### SASP score, tumour-only fibroblasts with technical/subtype adjustment\n")
sasp_fit <- lmerTest::lmer(
  SASP_zmean ~ FAP_status + MKI67_z + log_nCount + Cell_subtype + (1 | Patient),
  data = score_df,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
print(coef(summary(sasp_fit)))
cat("singular:", isSingular(sasp_fit), "\n")

all_score_df <- read.csv(score_file, check.names = FALSE)
all_score_df$log_nCount <- log1p(meta[all_score_df$cell, "nCount_RNA"])
all_score_df$Patient <- factor(all_score_df$Patient)
all_score_df$compartment <- factor(all_score_df$compartment, levels = c("Epithelial", "Fibroblast"))

cat("\n### SenMayo score, tumour-only compartment contrast with depth adjustment\n")
comp_fit <- lmerTest::lmer(
  SenMayo_zmean ~ compartment + MKI67_z + log_nCount + (1 | Patient),
  data = all_score_df,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
print(coef(summary(comp_fit)))
cat("singular:", isSingular(comp_fit), "\n")
sink()

cat("Wrote", normalizePath(out_text), "\n")
