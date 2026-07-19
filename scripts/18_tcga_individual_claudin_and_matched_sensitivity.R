#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 300)
set.seed(20260714)

suppressPackageStartupMessages(library(data.table))

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = TRUE)
output_dir <- file.path(task_root, "results", "L0_TCGA_claudin_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- Sys.getenv(
  "FAP_TCGA_EXPRESSION",
  file.path(task_root, "inputs", "TCGA_COADREAD_expression.txt.gz")
)
stopifnot(file.exists(expression_file))

raw_table <- fread(expression_file)
gene_symbols <- raw_table[[1]]
expression_matrix <- as.matrix(raw_table[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_symbols

sample_ids <- colnames(expression_matrix)
patient_ids <- substr(sample_ids, 1, 12)
sample_type <- substr(sample_ids, 14, 15)

score_zmean <- function(matrix_input, genes) {
  missing <- setdiff(genes, rownames(matrix_input))
  if (length(missing)) stop("Missing genes: ", paste(missing, collapse = ", "))
  gene_z <- t(scale(t(matrix_input[genes, , drop = FALSE]), center = TRUE, scale = TRUE))
  colMeans(gene_z, na.rm = TRUE)
}

fap_genes <- c(
  "FAP", "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "THY1", "PDPN",
  "TAGLN", "ACTA2", "MMP2", "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A"
)
reduced_fap_genes <- setdiff(fap_genes, c("TGFB1", "INHBA", "WNT2", "WNT5A"))
claudin_genes <- c("CLDN1", "CLDN2", "CLDN4")

scores <- data.table(
  sample_id = sample_ids,
  patient_id = patient_ids,
  sample_type = sample_type,
  reduced_FAP_CAF = score_zmean(expression_matrix, reduced_fap_genes),
  CLDN_core = score_zmean(expression_matrix, claudin_genes),
  CLDN1 = as.numeric(expression_matrix["CLDN1", ]),
  CLDN2 = as.numeric(expression_matrix["CLDN2", ]),
  CLDN4 = as.numeric(expression_matrix["CLDN4", ])
)

correlation_row <- function(scope, predictor, outcome, x, y) {
  complete <- complete.cases(x, y)
  x <- x[complete]; y <- y[complete]
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(test$estimate)
  n <- length(x)
  z <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  ci <- tanh(z + c(-1, 1) * qnorm(0.975) * se)
  data.table(
    scope = scope,
    predictor = predictor,
    outcome = outcome,
    n = n,
    spearman_rho = rho,
    ci_low = ci[1],
    ci_high = ci[2],
    p_value = test$p.value
  )
}

scope_masks <- list(
  all_profiles = scores$sample_type %in% c("01", "02", "06", "11"),
  all_tumours = scores$sample_type %in% c("01", "02", "06"),
  primary_tumours = scores$sample_type == "01"
)

correlations <- rbindlist(lapply(names(scope_masks), function(scope) {
  subset <- scores[scope_masks[[scope]]]
  outcomes <- c("CLDN1", "CLDN2", "CLDN4", "CLDN_core")
  result <- rbindlist(lapply(outcomes, function(outcome) {
    correlation_row(
      scope,
      "reduced_FAP_CAF",
      outcome,
      subset$reduced_FAP_CAF,
      subset[[outcome]]
    )
  }))
  result[, fdr_bh_four_outcomes := p.adjust(p_value, method = "BH")]
  result
}))

normal_indices <- which(sample_type == "11")
primary_indices <- which(sample_type == "01")
primary_by_patient <- split(primary_indices, patient_ids[primary_indices])
matched <- rbindlist(lapply(normal_indices, function(normal_index) {
  patient <- patient_ids[normal_index]
  candidates <- primary_by_patient[[patient]]
  if (is.null(candidates) || !length(candidates)) return(NULL)
  tumour_index <- candidates[1]
  data.table(
    patient_id = patient,
    tumour_sample = sample_ids[tumour_index],
    normal_sample = sample_ids[normal_index],
    FAP_difference = expression_matrix["FAP", tumour_index] - expression_matrix["FAP", normal_index],
    CLDN1_difference = expression_matrix["CLDN1", tumour_index] - expression_matrix["CLDN1", normal_index],
    CLDN2_difference = expression_matrix["CLDN2", tumour_index] - expression_matrix["CLDN2", normal_index],
    CLDN4_difference = expression_matrix["CLDN4", tumour_index] - expression_matrix["CLDN4", normal_index]
  )
}))
if (!nrow(matched)) stop("No matched primary tumour-adjacent-normal pairs found")

paired_results <- rbindlist(lapply(c("FAP", claudin_genes), function(gene) {
  differences <- matched[[paste0(gene, "_difference")]]
  test <- suppressWarnings(wilcox.test(
    differences,
    mu = 0,
    alternative = "two.sided",
    exact = FALSE,
    conf.int = TRUE,
    conf.level = 0.95
  ))
  data.table(
    gene = gene,
    paired_patients = length(differences),
    median_tumour_minus_normal = median(differences),
    hodges_lehmann_estimate = unname(test$estimate),
    ci_low = unname(test$conf.int[1]),
    ci_high = unname(test$conf.int[2]),
    statistic = unname(test$statistic),
    p_value = test$p.value
  )
}))
paired_results[, fdr_bh_four_genes := p.adjust(p_value, method = "BH")]

fwrite(scores, file.path(output_dir, "tcga_scores_and_individual_claudins.csv"))
fwrite(correlations, file.path(output_dir, "reduced_fap_caf_individual_claudin_correlations.csv"))
fwrite(matched, file.path(output_dir, "matched_primary_normal_differences.csv"))
fwrite(paired_results, file.path(output_dir, "matched_primary_normal_tests.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

print(correlations)
print(paired_results)
