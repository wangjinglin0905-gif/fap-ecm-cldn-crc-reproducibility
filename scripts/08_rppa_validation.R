#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(data.table))

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
output_dir <- file.path(task_root, "results", "RPPA")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- Sys.getenv("FAP_TCGA_EXPRESSION", file.path(task_root, "inputs", "TCGA_COADREAD_expression.txt.gz"))
rppa_file <- file.path(task_root, "inputs", "cbioportal_coadread_tcga_rppa_smad3_smad4_cdh1.csv")
legacy_merged_file <- Sys.getenv("FAP_LEGACY_RPPA", file.path(task_root, "inputs", "tcga_rna_rppa_merged.csv"))
stopifnot(file.exists(expression_file), file.exists(rppa_file))

expression_table <- fread(expression_file)
gene_symbols <- expression_table[[1]]
expression_matrix <- as.matrix(expression_table[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_symbols

tgf_beta_genes <- c("TGFB1", "TGFB2", "TGFB3", "SMAD2", "SMAD3", "SMAD4", "SMAD7")
stopifnot(all(tgf_beta_genes %in% rownames(expression_matrix)))
gene_z <- t(scale(t(expression_matrix[tgf_beta_genes, , drop = FALSE])))
tgf_beta_score <- colMeans(gene_z, na.rm = TRUE)
rna_scores <- data.table(sample_id = names(tgf_beta_score), TGF_beta_RNA = as.numeric(tgf_beta_score))

rppa_long <- fread(rppa_file)
rppa_wide <- dcast(rppa_long, sample_id + patient_id ~ gene, value.var = "value")
setnames(rppa_wide, c("SMAD3", "SMAD4", "CDH1"), c("SMAD3_RPPA", "SMAD4_RPPA", "CDH1_RPPA"))
merged <- merge(rna_scores, rppa_wide, by = "sample_id", all = FALSE)
fwrite(merged, file.path(output_dir, "tcga_rna_rppa_merged_from_raw_api.csv"))

correlation_results <- rbindlist(lapply(c("SMAD3_RPPA", "SMAD4_RPPA", "CDH1_RPPA"), function(protein) {
  complete <- complete.cases(merged$TGF_beta_RNA, merged[[protein]])
  test <- suppressWarnings(cor.test(merged$TGF_beta_RNA[complete], merged[[protein]][complete],
                                    method = "spearman", exact = FALSE))
  data.table(
    protein = protein,
    n = sum(complete),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}))
correlation_results[, fdr_bh_three_proteins := p.adjust(p_value, method = "BH")]
fwrite(correlation_results, file.path(output_dir, "rppa_tgf_beta_correlations_from_raw_api.csv"))

if (file.exists(legacy_merged_file)) {
  legacy <- fread(legacy_merged_file)
  audit <- merge(
    merged,
    legacy,
    by.x = "sample_id",
    by.y = "sampleId",
    suffixes = c("_raw_api", "_legacy")
  )
  comparison <- rbindlist(lapply(c("TGF_beta_RNA", "SMAD3_RPPA", "SMAD4_RPPA", "CDH1_RPPA"), function(variable) {
    raw_column <- paste0(variable, "_raw_api")
    legacy_column <- paste0(variable, "_legacy")
    data.table(
      variable = variable,
      n_overlap = sum(complete.cases(audit[[raw_column]], audit[[legacy_column]])),
      pearson_agreement = cor(audit[[raw_column]], audit[[legacy_column]], use = "complete.obs"),
      max_absolute_difference = max(abs(audit[[raw_column]] - audit[[legacy_column]]), na.rm = TRUE)
    )
  }))
  fwrite(comparison, file.path(output_dir, "legacy_processed_table_agreement.csv"))
}

manifest <- data.table(
  item = c("study", "molecular_profile", "sample_list", "api_retrieval_date", "rna_score_genes"),
  value = c("coadread_tcga", "coadread_tcga_rppa", "coadread_tcga_rppa",
            as.character(Sys.Date()), paste(tgf_beta_genes, collapse = ";"))
)
fwrite(manifest, file.path(output_dir, "rppa_download_and_score_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("RPPA raw-API validation complete: ", output_dir)
