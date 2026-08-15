options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript compute_minor_revision_ci.R <GSE132465_result_directory> <output_directory>")
}
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_mean <- function(file_name, outcome, seed, n_boot = 10000L) {
  path <- file.path(input_dir, file_name)
  dat <- read.csv(path, check.names = FALSE)
  delta <- dat[["diff_a_minus_b"]]
  stopifnot(length(delta) == 23L, all(is.finite(delta)))

  set.seed(seed)
  boot <- replicate(n_boot, mean(sample(delta, length(delta), replace = TRUE)))
  ci <- unname(quantile(boot, c(0.025, 0.975), names = FALSE))
  exact_test <- suppressWarnings(wilcox.test(delta, mu = 0, exact = TRUE,
                                              alternative = "two.sided",
                                              correct = FALSE))

  data.frame(
    outcome = outcome,
    estimand = "mean paired fibroblast-minus-epithelial difference",
    n_patients = length(delta),
    estimate = mean(delta),
    ci_level = 0.95,
    ci_method = paste0(n_boot, " patient-level percentile bootstrap resamples"),
    ci_low = ci[1],
    ci_high = ci[2],
    seed = seed,
    exact_wilcoxon_p = unname(exact_test$p.value),
    input_file = file.path("results", "GSE132465", file_name),
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  bootstrap_mean("tumor_compartment_SASP_patient_values.csv", "SASP25", 20260818L),
  bootstrap_mean("tumor_compartment_MKI67_patient_values.csv", "MKI67", 20260819L)
)

write.csv(results, file.path(output_dir, "minor_revision_paired_effect_cis.csv"),
          row.names = FALSE, quote = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "R_sessionInfo_minor_revision_ci.txt"))

print(results)
