options(stringsAsFactors = FALSE)

output_dir <- file.path("work", "reproducibility", "results", "L1_TISCH_GSE166555")
scores <- read.csv(file.path(output_dir, "patient_cross_compartment_scores.csv"), check.names = FALSE)
pseudobulk <- read.csv(file.path(output_dir, "patient_group_pseudobulk.csv"), check.names = FALSE)

bootstrap_spearman <- function(predictor, outcome, iterations = 10000L, seed = 20260714L) {
  set.seed(seed)
  sample_size <- length(predictor)
  estimates <- replicate(iterations, {
    indices <- sample.int(sample_size, sample_size, replace = TRUE)
    sampled_predictor <- predictor[indices]
    sampled_outcome <- outcome[indices]
    if (length(unique(sampled_predictor)) < 2L || length(unique(sampled_outcome)) < 2L) {
      return(NA_real_)
    }
    suppressWarnings(cor(sampled_predictor, sampled_outcome, method = "spearman"))
  })
  quantile(estimates, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}

association_tests <- function(minimum_caf_cells) {
  eligible <- scores[
    scores$caf_cells >= minimum_caf_cells & scores$epithelial_cells >= 100,
    ,
    drop = FALSE
  ]
  predictors <- c("FAP_log_cpm", "FAP_positive_fraction", "FAP_CAF", "FAP_CAF_de_ligand")
  outcomes <- c(
    "epithelial_CLDN1_log_cpm",
    "epithelial_CLDN2_log_cpm",
    "epithelial_CLDN4_log_cpm",
    "epithelial_CLDN_core"
  )
  records <- list()
  record_index <- 1L
  for (predictor in predictors) {
    for (outcome in outcomes) {
      complete <- eligible[complete.cases(eligible[, c(predictor, outcome)]), , drop = FALSE]
      test <- suppressWarnings(cor.test(
        complete[[predictor]],
        complete[[outcome]],
        method = "spearman",
        exact = FALSE
      ))
      confidence_interval <- bootstrap_spearman(
        complete[[predictor]],
        complete[[outcome]]
      )
      records[[record_index]] <- data.frame(
        minimum_caf_cells = minimum_caf_cells,
        predictor = predictor,
        outcome = outcome,
        patients = nrow(complete),
        spearman_rho = unname(test$estimate),
        p_value = test$p.value,
        bootstrap_ci_low = confidence_interval[1],
        bootstrap_ci_high = confidence_interval[2],
        prespecified_primary = predictor == "FAP_CAF_de_ligand" && outcome == "epithelial_CLDN_core"
      )
      record_index <- record_index + 1L
    }
  }
  result <- do.call(rbind, records)
  result$fdr_bh <- p.adjust(result$p_value, method = "BH")
  result
}

associations <- rbind(association_tests(20L), association_tests(5L))
write.csv(
  associations,
  file.path(output_dir, "patient_level_associations.csv"),
  row.names = FALSE
)

eligible <- scores[scores$caf_cells >= 20 & scores$epithelial_cells >= 100, , drop = FALSE]
leave_one_out_records <- list()
record_index <- 1L
for (omitted_patient in eligible$patient) {
  retained <- eligible[eligible$patient != omitted_patient, , drop = FALSE]
  for (predictor in c("FAP_log_cpm", "FAP_CAF", "FAP_CAF_de_ligand")) {
    test <- suppressWarnings(cor.test(
      retained[[predictor]],
      retained$epithelial_CLDN_core,
      method = "spearman",
      exact = FALSE
    ))
    leave_one_out_records[[record_index]] <- data.frame(
      omitted_patient = omitted_patient,
      predictor = predictor,
      patients = nrow(retained),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value
    )
    record_index <- record_index + 1L
  }
}
write.csv(
  do.call(rbind, leave_one_out_records),
  file.path(output_dir, "leave_one_patient_out_CLDN_core.csv"),
  row.names = FALSE
)

paired_specifications <- data.frame(
  analysis_group = c("CAF_like", "Epithelial", "Epithelial", "Epithelial"),
  gene = c("FAP", "CLDN1", "CLDN2", "CLDN4")
)
paired_records <- list()
for (row_index in seq_len(nrow(paired_specifications))) {
  analysis_group <- paired_specifications$analysis_group[row_index]
  gene <- paired_specifications$gene[row_index]
  subset <- pseudobulk[
    pseudobulk$analysis_group == analysis_group & pseudobulk$gene == gene,
    c("patient", "sample_origin", "log_cpm"),
    drop = FALSE
  ]
  wide <- reshape(subset, idvar = "patient", timevar = "sample_origin", direction = "wide")
  complete <- wide[complete.cases(wide[, c("log_cpm.Normal", "log_cpm.Tumor")]), , drop = FALSE]
  test <- if (nrow(complete) > 0L) {
    suppressWarnings(wilcox.test(
      complete$log_cpm.Tumor,
      complete$log_cpm.Normal,
      paired = TRUE,
      exact = FALSE
    ))
  } else {
    NULL
  }
  paired_records[[row_index]] <- data.frame(
    analysis_group = analysis_group,
    gene = gene,
    paired_patients = nrow(complete),
    median_tumor_minus_normal = if (nrow(complete) > 0L) {
      median(complete$log_cpm.Tumor - complete$log_cpm.Normal)
    } else {
      NA_real_
    },
    wilcoxon_statistic = if (is.null(test)) NA_real_ else unname(test$statistic),
    p_value = if (is.null(test)) NA_real_ else test$p.value
  )
}
paired_tests <- do.call(rbind, paired_records)
paired_tests$fdr_bh <- p.adjust(paired_tests$p_value, method = "BH")
write.csv(
  paired_tests,
  file.path(output_dir, "paired_tumor_normal_tests.csv"),
  row.names = FALSE
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
