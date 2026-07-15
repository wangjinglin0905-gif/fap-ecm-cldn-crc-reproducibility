#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 300)
set.seed(20260714)

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
output_dir <- file.path(task_root, "results", "L0_TCGA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

local_lib <- file.path(getwd(), "work", "CMSlib")
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(CMSclassifier)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(jsonlite)
})

expression_file <- Sys.getenv("FAP_TCGA_EXPRESSION", file.path(task_root, "inputs", "TCGA_COADREAD_expression.txt.gz"))
legacy_cms_file <- Sys.getenv("FAP_LEGACY_CMS", file.path(task_root, "inputs", "CMS_stratified_data.csv"))
legacy_msi_file <- Sys.getenv("FAP_LEGACY_MSI", file.path(task_root, "inputs", "msi_status_clean.csv"))

stopifnot(file.exists(expression_file), file.exists(legacy_cms_file), file.exists(legacy_msi_file))

message("Loading TCGA COAD/READ expression matrix")
raw_table <- fread(expression_file)
gene_symbols <- raw_table[[1]]
expression_matrix <- as.matrix(raw_table[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_symbols

sample_ids <- colnames(expression_matrix)
sample_type <- substr(sample_ids, 14, 15)
sample_scope <- ifelse(sample_type %in% c("01", "02", "06"), "Tumor",
                       ifelse(sample_type == "11", "Adjacent_normal", "Other"))
sample_manifest <- data.frame(
  sample_id = sample_ids,
  patient_id = substr(sample_ids, 1, 12),
  sample_type_code = sample_type,
  sample_scope = sample_scope
)
fwrite(sample_manifest, file.path(output_dir, "tcga_sample_manifest.csv"))
fwrite(as.data.frame(table(sample_type_code = sample_type, sample_scope = sample_scope)),
       file.path(output_dir, "tcga_sample_type_counts.csv"))

if (sum(sample_scope == "Tumor") != 383L || sum(sample_scope == "Adjacent_normal") != 51L) {
  stop("Unexpected TCGA sample composition")
}

score_zmean <- function(matrix_input, genes) {
  missing_genes <- setdiff(genes, rownames(matrix_input))
  if (length(missing_genes)) stop("Missing score genes: ", paste(missing_genes, collapse = ", "))
  gene_z <- t(scale(t(matrix_input[genes, , drop = FALSE]), center = TRUE, scale = TRUE))
  colMeans(gene_z, na.rm = TRUE)
}

fap_genes <- c("FAP", "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "THY1", "PDPN",
               "TAGLN", "ACTA2", "MMP2", "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
deligand_genes <- setdiff(fap_genes, c("TGFB1", "INHBA", "WNT2", "WNT5A"))
cldn_core_genes <- c("CLDN1", "CLDN2", "CLDN4")

score_table <- data.frame(
  sample_id = sample_ids,
  patient_id = substr(sample_ids, 1, 12),
  sample_scope = sample_scope,
  FAP_CAF = score_zmean(expression_matrix, fap_genes),
  FAP_CAF_de_ligand = score_zmean(expression_matrix, deligand_genes),
  CLDN_core = score_zmean(expression_matrix, cldn_core_genes)
)
fwrite(score_table, file.path(output_dir, "tcga_recomputed_scores.csv"))

spearman_ci <- function(x, y) {
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(test$estimate)
  n <- sum(complete.cases(x, y))
  z <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1, 1) * qnorm(0.975) * se)
  data.frame(n = n, rho = rho, ci_lower = ci[1], ci_upper = ci[2], p_value = test$p.value)
}

robustness <- rbind(
  cbind(comparison = "FAP_CAF_vs_de_ligand_all_profiles",
        spearman_ci(score_table$FAP_CAF, score_table$FAP_CAF_de_ligand)),
  cbind(comparison = "FAP_CAF_vs_CLDN_core_all_profiles",
        spearman_ci(score_table$FAP_CAF, score_table$CLDN_core)),
  cbind(comparison = "de_ligand_vs_CLDN_core_all_profiles",
        spearman_ci(score_table$FAP_CAF_de_ligand, score_table$CLDN_core))
)
fwrite(robustness, file.path(output_dir, "score_correlations.csv"))

tumor_normal_genes <- c("FAP", "CLDN1", "CLDN2", "CLDN4")
tumor_normal <- rbindlist(lapply(tumor_normal_genes, function(gene) {
  tumor_values <- expression_matrix[gene, sample_scope == "Tumor"]
  normal_values <- expression_matrix[gene, sample_scope == "Adjacent_normal"]
  test <- wilcox.test(tumor_values, normal_values, exact = FALSE)
  data.frame(
    gene = gene,
    n_tumor = length(tumor_values),
    n_adjacent_normal = length(normal_values),
    median_tumor = median(tumor_values),
    median_adjacent_normal = median(normal_values),
    median_difference = median(tumor_values) - median(normal_values),
    p_value = test$p.value
  )
}))
tumor_normal$fdr_bh <- p.adjust(tumor_normal$p_value, method = "BH")
fwrite(tumor_normal, file.path(output_dir, "tumor_normal_gene_comparisons.csv"))

message("Running CMSclassifier on 383 tumor profiles")
tumor_expression <- expression_matrix[, sample_scope == "Tumor", drop = FALSE]
symbol_to_entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(tumor_expression),
  keytype = "SYMBOL",
  column = "ENTREZID",
  multiVals = "first"
)
mapped <- !is.na(symbol_to_entrez)
entrez_expression <- tumor_expression[mapped, , drop = FALSE]
rownames(entrez_expression) <- unname(symbol_to_entrez[mapped])
entrez_expression <- rowsum(entrez_expression, group = rownames(entrez_expression), reorder = FALSE)
entrez_expression <- as.data.frame(entrez_expression, check.names = FALSE)

cms_fit <- classifyCMS(entrez_expression, method = "SSP")
cms_table <- data.frame(
  sample_id = rownames(cms_fit$predictedCMS),
  SSP_predicted = cms_fit$predictedCMS$SSP,
  SSP_nearest = cms_fit$nearestCMS$SSP
)

legacy_cms <- fread(legacy_cms_file)
legacy_cms <- legacy_cms[Sample %in% cms_table$sample_id, .(sample_id = Sample, legacy_CMS = CMS)]
cms_table <- merge(cms_table, legacy_cms, by = "sample_id", all.x = TRUE)
cms_table <- merge(cms_table, score_table[, c("sample_id", "FAP_CAF", "FAP_CAF_de_ligand", "CLDN_core")],
                   by = "sample_id", all.x = TRUE)
fwrite(cms_table, file.path(output_dir, "cms_classifier_results.csv"))

cms_agreement <- rbind(
  data.frame(comparison = "SSP_vs_legacy", n_complete = sum(complete.cases(cms_table$SSP_predicted, cms_table$legacy_CMS)),
             agreement = mean(cms_table$SSP_predicted == cms_table$legacy_CMS, na.rm = TRUE))
)
fwrite(cms_agreement, file.path(output_dir, "cms_label_agreement.csv"))

cms_test_one <- function(labels, method_name) {
  keep <- !is.na(labels)
  dat <- data.frame(score = cms_table$FAP_CAF[keep], CMS = labels[keep])
  kw <- kruskal.test(score ~ CMS, data = dat)
  pairwise <- pairwise.wilcox.test(dat$score, dat$CMS, p.adjust.method = "bonferroni", exact = FALSE)
  summary <- rbindlist(lapply(split(dat$score, dat$CMS), function(values) {
    data.frame(n = length(values), median = median(values), q1 = quantile(values, 0.25), q3 = quantile(values, 0.75))
  }), idcol = "CMS")
  summary$method <- method_name
  list(summary = summary, kw = data.frame(method = method_name, n = nrow(dat), statistic = unname(kw$statistic),
                                           df = unname(kw$parameter), p_value = kw$p.value),
       pairwise = as.data.frame(as.table(pairwise$p.value), stringsAsFactors = FALSE))
}

cms_results <- lapply(list(SSP = cms_table$SSP_predicted, legacy = cms_table$legacy_CMS),
                      function(x) cms_test_one(x, "placeholder"))
names(cms_results) <- c("SSP", "legacy")
for (name in names(cms_results)) {
  cms_results[[name]]$summary$method <- name
  cms_results[[name]]$kw$method <- name
  cms_results[[name]]$pairwise$method <- name
}
fwrite(rbindlist(lapply(cms_results, `[[`, "summary")), file.path(output_dir, "cms_score_summary.csv"))
fwrite(rbindlist(lapply(cms_results, `[[`, "kw")), file.path(output_dir, "cms_kruskal_wallis.csv"))
fwrite(rbindlist(lapply(cms_results, `[[`, "pairwise"), fill = TRUE), file.path(output_dir, "cms_pairwise_bonferroni.csv"))

message("Using MANTIS scores fetched from cBioPortal for the 383 tumor profiles")
mantis_file <- file.path(task_root, "inputs", "cbioportal_mantis_scores.csv")
stopifnot(file.exists(mantis_file))
mantis <- fread(mantis_file)
mantis$mantis_score <- as.numeric(mantis$mantis_score)
mantis$mantis_group_0.4 <- ifelse(mantis$mantis_score > 0.4, "MSI-H", "MSS")

legacy_msi <- fread(legacy_msi_file)
names(legacy_msi) <- c("patient_id", "legacy_MSI_group")
msi_table <- merge(score_table[score_table$sample_scope == "Tumor", ],
                   mantis[, c("sample_id", "mantis_score", "mantis_group_0.4")],
                   by = "sample_id", all.x = TRUE)
msi_table <- merge(msi_table, legacy_msi, by = "patient_id", all.x = TRUE)
fwrite(msi_table, file.path(output_dir, "msi_matched_scores.csv"))

msi_tests <- rbindlist(lapply(c("FAP_CAF", "FAP_CAF_de_ligand", "CLDN_core"), function(score_name) {
  complete <- !is.na(msi_table$mantis_group_0.4) & msi_table$mantis_group_0.4 %in% c("MSS", "MSI-H")
  values <- msi_table[[score_name]][complete]
  groups <- msi_table$mantis_group_0.4[complete]
  test <- wilcox.test(values ~ groups, exact = FALSE)
  data.frame(score = score_name, n = length(values), n_MSI_H = sum(groups == "MSI-H"), n_MSS = sum(groups == "MSS"),
             median_MSI_H = median(values[groups == "MSI-H"]), median_MSS = median(values[groups == "MSS"]),
             p_value = test$p.value)
}))
msi_tests$fdr_bh_three_scores <- p.adjust(msi_tests$p_value, method = "BH")
fwrite(msi_tests, file.path(output_dir, "msi_score_comparisons.csv"))

overlap <- complete.cases(msi_table$legacy_MSI_group, msi_table$mantis_group_0.4)
msi_agreement <- data.frame(
  n_overlap = sum(overlap),
  agreement = mean(msi_table$legacy_MSI_group[overlap] == msi_table$mantis_group_0.4[overlap])
)
fwrite(msi_agreement, file.path(output_dir, "msi_label_agreement.csv"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("L0 TCGA analysis complete: ", output_dir)
