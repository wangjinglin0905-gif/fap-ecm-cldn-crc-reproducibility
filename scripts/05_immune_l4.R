#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260714)

suppressPackageStartupMessages({
  library(data.table)
  library(MCPcounter)
  library(EPIC)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
output_dir <- file.path(task_root, "results", "L4_immune")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- Sys.getenv("FAP_TCGA_EXPRESSION", file.path(task_root, "inputs", "TCGA_COADREAD_expression.txt.gz"))
score_file <- file.path(task_root, "results", "L0_TCGA", "tcga_recomputed_scores.csv")
stopifnot(file.exists(expression_file), file.exists(score_file))

raw_table <- fread(expression_file)
gene_symbols <- raw_table[[1]]
expression_matrix <- as.matrix(raw_table[, -1])
storage.mode(expression_matrix) <- "numeric"
rownames(expression_matrix) <- gene_symbols

scores <- fread(score_file)
scores <- scores[sample_scope == "Tumor"]
tumor_expression <- expression_matrix[, scores$sample_id, drop = FALSE]
group_median <- median(scores$FAP_CAF)
scores[, FAP_CAF_group := factor(ifelse(FAP_CAF > group_median, "High", "Low"), levels = c("Low", "High"))]
fwrite(scores[, .(sample_id, FAP_CAF, FAP_CAF_group)], file.path(output_dir, "tumor_fap_caf_groups.csv"))

compare_groups <- function(score_matrix, group_table, method_name) {
  results <- rbindlist(lapply(rownames(score_matrix), function(cell_type) {
    values <- as.numeric(score_matrix[cell_type, group_table$sample_id])
    high <- values[group_table$FAP_CAF_group == "High"]
    low <- values[group_table$FAP_CAF_group == "Low"]
    test <- wilcox.test(high, low, exact = FALSE)
    data.table(
      method = method_name,
      cell_type = cell_type,
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
  results[]
}

message("Running MCP-counter on 383 tumor profiles")
mcp_genes_file <- system.file("extdata", "genes.txt", package = "IOBR")
if (!nzchar(mcp_genes_file)) stop("IOBR MCP-counter signature file not found")
mcp_genes <- fread(mcp_genes_file, data.table = FALSE, check.names = FALSE)
mcp_scores <- MCPcounter.estimate(
  tumor_expression,
  featuresType = "HUGO_symbols",
  genes = mcp_genes
)
fwrite(as.data.table(mcp_scores, keep.rownames = "cell_type"), file.path(output_dir, "mcpcounter_tumor_scores.csv"))
mcp_results <- compare_groups(mcp_scores, scores, "MCPcounter")
fwrite(mcp_results, file.path(output_dir, "mcpcounter_high_vs_low.csv"))

message("Running EPIC on back-transformed tumor expression")
epic_input <- pmax(2^tumor_expression - 1, 0)
epic_warnings <- character()
epic_fit <- withCallingHandlers(
  EPIC(bulk = epic_input, reference = "TRef", scaleExprs = TRUE, withOtherCells = TRUE),
  warning = function(warning_condition) {
    epic_warnings <<- c(epic_warnings, conditionMessage(warning_condition))
    invokeRestart("muffleWarning")
  }
)
epic_fractions <- t(as.matrix(epic_fit$cellFractions))
colnames(epic_fractions) <- rownames(epic_fit$cellFractions)
fwrite(as.data.table(epic_fractions, keep.rownames = "cell_type"), file.path(output_dir, "epic_tumor_cell_fractions.csv"))
epic_results <- compare_groups(epic_fractions, scores, "EPIC")
fwrite(epic_results, file.path(output_dir, "epic_high_vs_low.csv"))
writeLines(unique(epic_warnings), file.path(output_dir, "epic_warnings.txt"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("L4 immune deconvolution complete: ", output_dir)
