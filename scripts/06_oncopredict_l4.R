#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260714)

suppressPackageStartupMessages({
  library(data.table)
  library(oncoPredict)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
output_dir <- file.path(task_root, "results", "L4_oncoPredict")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- Sys.getenv("FAP_TCGA_EXPRESSION", file.path(task_root, "inputs", "TCGA_COADREAD_expression.txt.gz"))
group_file <- file.path(task_root, "results", "L4_immune", "tumor_fap_caf_groups.csv")
stopifnot(file.exists(expression_file), file.exists(group_file))

raw_table <- fread(expression_file)
gene_symbols <- raw_table[[1]]
expression_matrix <- as.matrix(raw_table[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_symbols

groups <- fread(group_file)
tumor_expression <- expression_matrix[, groups$sample_id, drop = FALSE]

training_expression <- readRDS(system.file("extdata", "GDSC2_Expr_short.rds", package = "oncoPredict"))
training_response <- readRDS(system.file("extdata", "GDSC2_Res.rds", package = "oncoPredict"))

message("Running oncoPredict with GDSC2 training data on 383 tumor profiles")
prediction_warnings <- character()
predictions <- withCallingHandlers(
  calcPhenotype(
    trainingExprData = training_expression,
    trainingPtype = training_response,
    testExprData = tumor_expression,
    batchCorrect = "eb",
    powerTransformPhenotype = TRUE,
    removeLowVaryingGenes = 0.2,
    minNumSamples = 10,
    selection = 1,
    printOutput = TRUE,
    pcr = FALSE,
    removeLowVaringGenesFrom = "homogenizeData",
    report_pc = FALSE,
    cc = FALSE,
    parallel = FALSE
  ),
  warning = function(warning_condition) {
    prediction_warnings <<- c(prediction_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)

predictions <- as.matrix(predictions)
if (!all(groups$sample_id %in% rownames(predictions)) && all(groups$sample_id %in% colnames(predictions))) {
  predictions <- t(predictions)
}
stopifnot(all(groups$sample_id %in% rownames(predictions)))
predictions <- predictions[groups$sample_id, , drop = FALSE]
fwrite(as.data.table(predictions, keep.rownames = "sample_id"), file.path(output_dir, "gdsc2_predicted_response_tumors.csv"))

results <- rbindlist(lapply(colnames(predictions), function(drug) {
  values <- predictions[, drug]
  high <- values[groups$FAP_CAF_group == "High"]
  low <- values[groups$FAP_CAF_group == "Low"]
  test <- wilcox.test(high, low, exact = FALSE)
  data.table(
    drug = drug,
    n_high = length(high),
    n_low = length(low),
    median_high = median(high, na.rm = TRUE),
    median_low = median(low, na.rm = TRUE),
    median_difference = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    mean_high = mean(high, na.rm = TRUE),
    mean_low = mean(low, na.rm = TRUE),
    mean_difference = mean(high, na.rm = TRUE) - mean(low, na.rm = TRUE),
    p_value = test$p.value
  )
}))
results[, fdr_bh := p.adjust(p_value, method = "BH")]
setorder(results, mean_difference)
results[, sensitivity_rank := seq_len(.N)]
fwrite(results, file.path(output_dir, "gdsc2_high_vs_low_all_drugs.csv"))
fwrite(head(results, 25), file.path(output_dir, "gdsc2_top25_lower_predicted_response.csv"))

quality_audit <- data.table(
  drug = colnames(predictions),
  max_absolute_prediction = apply(abs(predictions), 2, max, na.rm = TRUE)
)
quality_audit[, numerically_implausible := !is.finite(max_absolute_prediction) | max_absolute_prediction > 100]
fwrite(quality_audit, file.path(output_dir, "gdsc2_prediction_quality_audit.csv"))

message("Running no-power-transform sensitivity analysis")
no_power_warnings <- character()
predictions_no_power <- withCallingHandlers(
  calcPhenotype(
    trainingExprData = training_expression,
    trainingPtype = training_response,
    testExprData = tumor_expression,
    batchCorrect = "eb",
    powerTransformPhenotype = FALSE,
    removeLowVaryingGenes = 0.2,
    minNumSamples = 10,
    selection = 1,
    printOutput = FALSE,
    pcr = FALSE,
    removeLowVaringGenesFrom = "homogenizeData",
    report_pc = FALSE,
    cc = FALSE,
    parallel = FALSE
  ),
  warning = function(warning_condition) {
    no_power_warnings <<- c(no_power_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)
predictions_no_power <- as.matrix(predictions_no_power)
if (!all(groups$sample_id %in% rownames(predictions_no_power)) &&
    all(groups$sample_id %in% colnames(predictions_no_power))) {
  predictions_no_power <- t(predictions_no_power)
}
predictions_no_power <- predictions_no_power[groups$sample_id, , drop = FALSE]
fwrite(as.data.table(predictions_no_power, keep.rownames = "sample_id"),
       file.path(output_dir, "gdsc2_predicted_response_no_power_transform.csv"))

no_power_results <- rbindlist(lapply(colnames(predictions_no_power), function(drug) {
  values <- predictions_no_power[, drug]
  high <- values[groups$FAP_CAF_group == "High"]
  low <- values[groups$FAP_CAF_group == "Low"]
  test <- wilcox.test(high, low, exact = FALSE)
  data.table(
    drug = drug,
    mean_difference = mean(high, na.rm = TRUE) - mean(low, na.rm = TRUE),
    median_difference = median(high, na.rm = TRUE) - median(low, na.rm = TRUE),
    p_value = test$p.value
  )
}))
no_power_results[, fdr_bh := p.adjust(p_value, method = "BH")]
setorder(no_power_results, mean_difference)
no_power_results[, sensitivity_rank := seq_len(.N)]
fwrite(no_power_results, file.path(output_dir, "gdsc2_no_power_transform_high_vs_low.csv"))
writeLines(unique(no_power_warnings), file.path(output_dir, "oncoPredict_no_power_warnings.txt"))

writeLines(unique(prediction_warnings), file.path(output_dir, "oncoPredict_warnings.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Tumor-only oncoPredict analysis complete: ", output_dir)
