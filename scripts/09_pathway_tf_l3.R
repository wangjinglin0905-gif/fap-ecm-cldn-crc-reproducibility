#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(progeny)
  library(decoupleR)
  library(dorothea)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
input_dir <- file.path(task_root, "results", "L1_single_cell")
output_dir <- file.path(task_root, "results", "L3_pathway_TF")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(input_dir, "myofibroblast_pseudobulk_counts.tsv.gz")
burden_file <- file.path(input_dir, "patient_sender_burden.csv")
stopifnot(file.exists(count_file), file.exists(burden_file))

count_table <- fread(count_file)
genes <- count_table[[1]]
counts <- as.matrix(count_table[, -1])
storage.mode(counts) <- "integer"
rownames(counts) <- genes

burden <- fread(burden_file)
eligible <- burden[tumor_myofibroblasts >= 20]
eligible[, group := factor(ifelse(strict_sender_fraction > median(strict_sender_fraction), "High", "Low"),
                           levels = c("Low", "High"))]
counts <- counts[, eligible$patient, drop = FALSE]

dge <- DGEList(counts = counts)
keep <- filterByExpr(dge, group = eligible$group, min.count = 10)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)

compare_scores <- function(score_matrix, method_name) {
  results <- rbindlist(lapply(colnames(score_matrix), function(feature) {
    values <- score_matrix[, feature]
    high <- values[eligible$group == "High"]
    low <- values[eligible$group == "Low"]
    group_test <- wilcox.test(high, low, exact = FALSE)
    correlation_test <- suppressWarnings(cor.test(values, eligible$strict_sender_fraction,
                                                   method = "spearman", exact = FALSE))
    data.table(
      method = method_name,
      feature = feature,
      n_patients = length(values),
      high_mean = mean(high),
      low_mean = mean(low),
      mean_difference = mean(high) - mean(low),
      wilcoxon_p = group_test$p.value,
      burden_rho = unname(correlation_test$estimate),
      burden_p = correlation_test$p.value
    )
  }))
  results[, wilcoxon_fdr_bh := p.adjust(wilcoxon_p, method = "BH")]
  results[, burden_fdr_bh := p.adjust(burden_p, method = "BH")]
  results[]
}

message("Running PROGENy on myofibroblast pseudobulk expression")
progeny_scores <- progeny(log_cpm, scale = TRUE, organism = "Human", top = 500, perm = 1, verbose = TRUE)
fwrite(as.data.table(progeny_scores, keep.rownames = "patient"), file.path(output_dir, "progeny_patient_scores.csv"))
progeny_results <- compare_scores(progeny_scores, "PROGENy")
fwrite(progeny_results, file.path(output_dir, "progeny_high_vs_low_and_burden.csv"))

message("Running DoRothEA A+B regulons with ULM")
network <- as.data.table(dorothea_hs)[confidence %in% c("A", "B"), .(source = tf, target, mor)]
tf_long <- as.data.table(run_ulm(log_cpm, network, .source = source, .target = target, .mor = mor, minsize = 5))
score_column <- if ("score" %in% names(tf_long)) "score" else "estimate"
tf_wide <- dcast(tf_long, condition ~ source, value.var = score_column)
setnames(tf_wide, "condition", "patient")
tf_matrix <- as.matrix(tf_wide[, -1])
rownames(tf_matrix) <- tf_wide$patient
tf_matrix <- tf_matrix[eligible$patient, , drop = FALSE]
fwrite(as.data.table(tf_matrix, keep.rownames = "patient"), file.path(output_dir, "dorothea_patient_scores.csv"))
tf_results <- compare_scores(tf_matrix, "DoRothEA_ULM_AB")
fwrite(tf_results, file.path(output_dir, "dorothea_high_vs_low_and_burden.csv"))

emt_tfs <- c("SMAD3", "SMAD4", "TWIST1", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TCF7L2", "LEF1")
fwrite(tf_results[feature %in% emt_tfs], file.path(output_dir, "dorothea_selected_emt_tfs.csv"))
fwrite(eligible, file.path(output_dir, "eligible_patient_groups.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("L3 pathway and TF analysis complete: ", output_dir)
