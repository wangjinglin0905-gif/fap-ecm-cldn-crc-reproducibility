#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(paste(
    "Usage: Rscript 03_correct_spatial_matched_null.R",
    "<Route-B results_full directory> <output directory>"
  ))
}

source_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- c("GSE280315", "Valdeolivas_Visium")

correct_one <- function(cohort) {
  cohort_dir <- file.path(source_dir, cohort)
  summary_path <- file.path(cohort_dir, "expression_matched_null_summary.csv")
  draws_path <- file.path(cohort_dir, "expression_matched_null_draws.csv")
  old <- read.csv(summary_path, check.names = FALSE)
  draws_table <- read.csv(draws_path, check.names = FALSE)
  if (nrow(old) != 1L || ncol(draws_table) < 2L) {
    stop("Unexpected matched-null schema for ", cohort)
  }

  observed <- as.numeric(old$observed_median_patient_mean_difference[[1]])
  draws <- as.numeric(draws_table[[2]])
  draws <- draws[is.finite(draws)]
  b <- length(draws)
  center <- median(draws)

  upper <- (1 + sum(draws >= observed)) / (b + 1)
  lower <- (1 + sum(draws <= observed)) / (b + 1)
  doubled_tail <- min(1, 2 * min(upper, lower))
  centered_abs <- (1 + sum(abs(draws - center) >= abs(observed - center))) / (b + 1)

  data.frame(
    cohort = cohort,
    observed_median_patient_mean_difference = observed,
    null_center_median = center,
    null_q025 = unname(quantile(draws, 0.025, type = 6)),
    null_q975 = unname(quantile(draws, 0.975, type = 6)),
    observed_minus_null_center = observed - center,
    empirical_upper_tail_p = upper,
    empirical_doubled_tail_p = doubled_tail,
    empirical_centered_absolute_p = centered_abs,
    original_reported_two_sided_p = as.numeric(old$empirical_two_sided_p[[1]]),
    draws = b,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(cohorts, correct_one))
write.csv(
  results,
  file.path(out_dir, "spatial_matched_null_corrected.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

notes <- c(
  "Correction rule",
  "The original two-sided implementation compared abs(null draw) with abs(observed), implicitly centering the null at zero.",
  "The corrected primary empirical P value compares absolute deviations around the empirical null median:",
  "(1 + sum(abs(null - median(null)) >= abs(observed - median(null)))) / (B + 1).",
  "The one-sided upper-tail and doubled-tail values are retained as sensitivity quantities.",
  "No Route-B source output was overwritten."
)
writeLines(notes, file.path(out_dir, "spatial_matched_null_correction_notes.txt"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_spatial_null.txt"))

print(results, digits = 6, row.names = FALSE)
