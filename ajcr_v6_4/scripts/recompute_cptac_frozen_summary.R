options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: Rscript recompute_cptac_frozen_summary.R ",
    "<cptac_protein_sdc4_cd44.csv> <output_directory>"
  )
}

input_file <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(input_file, check.names = FALSE)
required <- c("sample_id", "FAP", "COL1A1", "COL1A2", "FN1", "SDC4", "CD44")
stopifnot(all(required %in% names(dat)), nrow(dat) == 97L)

for (gene in setdiff(required, "sample_id")) {
  dat[[paste0(gene, "_z")]] <- as.numeric(scale(dat[[gene]]))
}
dat$matrix3 <- rowMeans(dat[, c("COL1A1_z", "COL1A2_z", "FN1_z")])
dat$receptor2 <- rowMeans(dat[, c("SDC4_z", "CD44_z")])

bootstrap_spearman <- function(x, y, seed, reps = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  set.seed(seed)
  boot <- replicate(reps, {
    idx <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    n = length(x),
    rho = unname(test$estimate),
    ci_low = unname(quantile(boot, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(boot, 0.975, type = 6, na.rm = TRUE)),
    p_value = test$p.value,
    bootstrap_resamples = reps,
    bootstrap_seed = seed
  )
}

result <- rbind(
  cbind(
    comparison = "FAP protein vs FAP-excluded matrix3 protein score",
    bootstrap_spearman(dat$FAP_z, dat$matrix3, 2026081505L)
  ),
  cbind(
    comparison = "matrix3 protein score vs receptor2 protein score",
    bootstrap_spearman(dat$matrix3, dat$receptor2, 2026081506L)
  )
)
result$fdr_bh_two <- p.adjust(result$p_value, method = "BH")
result$interpretation <- c(
  "FAP-ECM protein covariation; not a SASP or senescence assay",
  "No coordinated receptor co-induction"
)

write.csv(
  result,
  file.path(out_dir, "cptac_fap_ecm_and_receptor_frozen.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "R_sessionInfo_cptac_frozen.txt")
)

message("Wrote frozen CPTAC summary to: ", normalizePath(out_dir))
