options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026080603)

root <- normalizePath(getwd(), winslash = "/")
input_file <- file.path(root, "data", "public",
                        "cptac_protein_sdc4_cd44.csv")
out_dir <- file.path(root, "results", "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))

dat <- read.csv(input_file, check.names = FALSE)
required <- c("sample_id", "FAP", "COL1A1", "COL1A2", "FN1",
              "SDC4", "CD44")
stopifnot(all(required %in% names(dat)), nrow(dat) == 97L)
for (gene in setdiff(required, "sample_id")) dat[[paste0(gene, "_z")]] <- as.numeric(scale(dat[[gene]]))
dat$matrix3 <- rowMeans(dat[, c("COL1A1_z", "COL1A2_z", "FN1_z")])
dat$receptor2 <- rowMeans(dat[, c("SDC4_z", "CD44_z")])

bootstrap_cor <- function(x, y, seed, reps = 5000L) {
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  set.seed(seed)
  boot <- replicate(reps, { i <- sample.int(length(x), replace = TRUE); suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(n = length(x), rho = unname(test$estimate),
             ci_low = unname(quantile(boot, 0.025, type = 6, na.rm = TRUE)),
             ci_high = unname(quantile(boot, 0.975, type = 6, na.rm = TRUE)),
             p_value = test$p.value)
}
matrix_result <- bootstrap_cor(dat$FAP_z, dat$matrix3, 2026080611L)
receptor_result <- bootstrap_cor(dat$matrix3, dat$receptor2, 2026080612L)
result <- rbind(
  cbind(comparison = "FAP vs matrix3", matrix_result),
  cbind(comparison = "matrix3 vs receptor2", receptor_result)
)
result$fdr_bh_two <- p.adjust(result$p_value, method = "BH")
write.csv(result, file.path(out_dir, "CPTAC_primary_correlations.csv"), row.names = FALSE)
write.csv(dat[, c("sample_id", "FAP_z", "matrix3", "receptor2")],
          file.path(out_dir, "Figure1_CPTAC_source_data.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "CPTAC_sessionInfo.txt"))
