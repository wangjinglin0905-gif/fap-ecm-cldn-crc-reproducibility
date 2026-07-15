#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
input_file <- file.path(task_root, "inputs", "cbioportal_coad_cptac_2019_selected_proteins.csv")
output_dir <- file.path(task_root, "results", "CPTAC_protein")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))

protein_long <- fread(input_file)
protein_wide <- dcast(protein_long, sample_id + patient_id ~ gene, value.var = "value")
fwrite(protein_wide, file.path(output_dir, "cptac_selected_proteins_wide.csv"))

availability <- data.table(
  gene = c("FAP", "COL1A1", "COL1A2", "FN1", "CLDN1", "CLDN2", "CLDN4"),
  n_samples = sapply(c("FAP", "COL1A1", "COL1A2", "FN1", "CLDN1", "CLDN2", "CLDN4"),
                     function(gene) if (gene %in% names(protein_wide)) sum(!is.na(protein_wide[[gene]])) else 0)
)
fwrite(availability, file.path(output_dir, "cptac_protein_availability.csv"))

ecm_genes <- intersect(c("FAP", "COL1A1", "COL1A2", "FN1"), names(protein_wide))
ecm_z <- scale(as.matrix(protein_wide[, ..ecm_genes]))
protein_wide[, ECM_protein_score := rowMeans(ecm_z, na.rm = TRUE)]

comparisons <- list(
  c("FAP", "COL1A1"),
  c("FAP", "COL1A2"),
  c("FAP", "FN1"),
  c("FAP", "CLDN4"),
  c("ECM_protein_score", "CLDN4"),
  c("ECM_protein_score", "CLDN1")
)
correlations <- rbindlist(lapply(comparisons, function(pair) {
  if (!all(pair %in% names(protein_wide))) {
    return(data.table(x = pair[1], y = pair[2], n = 0, rho = NA_real_, p_value = NA_real_))
  }
  complete <- complete.cases(protein_wide[[pair[1]]], protein_wide[[pair[2]]])
  if (sum(complete) < 5) {
    return(data.table(x = pair[1], y = pair[2], n = sum(complete), rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(protein_wide[[pair[1]]][complete], protein_wide[[pair[2]]][complete],
                                    method = "spearman", exact = FALSE))
  data.table(x = pair[1], y = pair[2], n = sum(complete), rho = unname(test$estimate), p_value = test$p.value)
}))
correlations[, fdr_bh_six_tests := p.adjust(p_value, method = "BH")]
fwrite(correlations, file.path(output_dir, "cptac_prespecified_correlations.csv"))

protein_wide[, FAP_group := factor(ifelse(FAP > median(FAP, na.rm = TRUE), "High", "Low"),
                                   levels = c("Low", "High"))]
group_results <- rbindlist(lapply(intersect(c("CLDN1", "CLDN4"), names(protein_wide)), function(gene) {
  complete <- complete.cases(protein_wide[[gene]], protein_wide$FAP_group)
  high <- protein_wide[[gene]][complete & protein_wide$FAP_group == "High"]
  low <- protein_wide[[gene]][complete & protein_wide$FAP_group == "Low"]
  test <- wilcox.test(high, low, exact = FALSE)
  data.table(
    protein = gene,
    n_high = length(high),
    n_low = length(low),
    median_high = median(high),
    median_low = median(low),
    median_difference = median(high) - median(low),
    p_value = test$p.value
  )
}))
group_results[, fdr_bh_two_tests := p.adjust(p_value, method = "BH")]
fwrite(group_results, file.path(output_dir, "cptac_cldn_by_fap_group.csv"))

manifest <- data.table(
  item = c("study", "molecular_profile", "sample_list", "scope", "retrieval_date"),
  value = c("coad_cptac_2019", "coad_cptac_2019_protein_quantification",
            "coad_cptac_2019_protein_quantification", "tumor-only proteomics", as.character(Sys.Date()))
)
fwrite(manifest, file.path(output_dir, "cptac_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("CPTAC protein validation complete: ", output_dir)
