#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 13a_lock_scRNA_endpoint_estimands.R <aligned scores.csv> <output.csv>")
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
output_path <- args[[2]]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_path, check.names = FALSE)

definition <- c(
  `GSE132465|SenMayo_score` = "mean cell-level gene-wise z score; 119 source-available overlap-purged genes",
  `GSE132465|SASP_score` = "mean cell-level gene-wise z score; SASP25",
  `GSE132465|MKI67_score` = "mean cell-level log1p(counts per 10,000) MKI67 expression",
  `GSE166555|SenMayo_score` = "patient-compartment pseudobulk gene-wise z score; 114 source-available overlap-purged genes",
  `GSE166555|SASP_score` = "patient-compartment pseudobulk gene-wise z score; SASP25",
  `GSE166555|MKI67_score` = "patient-compartment pseudobulk log1p counts-per-million MKI67 expression"
)

endpoint_label <- c(SenMayo_score = "SenMayo", SASP_score = "SASP25", MKI67_score = "MKI67")

bootstrap_ci <- function(x, statistic = c("median", "mean"), b = 10000L, seed = 1L) {
  statistic <- match.arg(statistic)
  set.seed(seed)
  n <- length(x)
  draws <- replicate(b, {
    z <- x[sample.int(n, n, replace = TRUE)]
    if (statistic == "median") median(z) else mean(z)
  })
  unname(quantile(draws, c(0.025, 0.975), names = FALSE, type = 7))
}

rows <- list()
k <- 1L
for (cohort in c("GSE132465", "GSE166555")) {
  for (endpoint in c("SenMayo_score", "SASP_score", "MKI67_score")) {
    x <- d[d$cohort == cohort, c("patient", "compartment", endpoint)]
    w <- reshape(x, idvar = "patient", timevar = "compartment", direction = "wide")
    diff <- w[[paste0(endpoint, ".Fibroblast")]] - w[[paste0(endpoint, ".Epithelial")]]
    seed_base <- if (cohort == "GSE132465") {
      c(SenMayo_score = 20260817L, SASP_score = 20260818L, MKI67_score = 20260819L)[[endpoint]]
    } else {
      c(SenMayo_score = 20260815L, SASP_score = 20260816L, MKI67_score = 20260817L)[[endpoint]]
    }
    ci_median <- bootstrap_ci(diff, "median", seed = seed_base)
    ci_mean <- bootstrap_ci(diff, "mean", seed = seed_base + 100L)
    wt <- suppressWarnings(wilcox.test(diff, mu = 0, paired = FALSE, exact = TRUE,
                                       alternative = "two.sided", conf.int = FALSE))
    key <- paste(cohort, endpoint, sep = "|")
    rows[[k]] <- data.frame(
      cohort = cohort,
      endpoint = endpoint_label[[endpoint]],
      score_definition = definition[[key]],
      n_patients = length(diff),
      primary_summary = "median fibroblast-minus-epithelial patient difference",
      median_difference = median(diff),
      median_ci_low = ci_median[[1]],
      median_ci_high = ci_median[[2]],
      mean_difference = mean(diff),
      mean_ci_low = ci_mean[[1]],
      mean_ci_high = ci_mean[[2]],
      fibroblast_higher = sum(diff > 0),
      epithelial_higher = sum(diff < 0),
      exact_wilcoxon_p = wt$p.value,
      bootstrap_resamples = 10000L,
      bootstrap_seed = seed_base,
      source_artifact = input_path,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

out <- do.call(rbind, rows)
write.csv(out, output_path, row.names = FALSE)
print(out)
