options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript domain_tumor_only_check.R <single_cell_result_directory>")
}
in_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_file <- file.path(in_dir, "domain_tumor_only_check.txt")

all_cells <- read.csv(file.path(in_dir, "cell_scores_all_fibroblast_epithelial.csv"), check.names = FALSE)
tumor_cells <- read.csv(file.path(in_dir, "cell_scores_tumor_fibroblast_epithelial.csv"), check.names = FALSE)

prep <- function(d) {
  d$Patient <- factor(d$Patient)
  d$Class <- factor(d$Class)
  d$Cell_subtype <- factor(d$Cell_subtype)
  d$compartment <- factor(d$compartment, levels = c("Epithelial", "Fibroblast"))
  d$FAP_status <- factor(d$FAP_status, levels = c("FAP-", "FAP+"))
  d
}

all_cells <- prep(all_cells)
tumor_cells <- prep(tumor_cells)
all_fib <- droplevels(all_cells[all_cells$compartment == "Fibroblast", ])
tumor_fib <- droplevels(tumor_cells[tumor_cells$compartment == "Fibroblast", ])

fit_and_print <- function(label, formula, data) {
  cat("\n###", label, "\n")
  cat("n =", nrow(data), "; patients =", nlevels(droplevels(data$Patient)), "\n")
  fit <- lmerTest::lmer(
    formula,
    data = data,
    REML = FALSE,
    control = lme4::lmerControl(optimizer = "bobyqa")
  )
  print(coef(summary(fit)))
  cat("singular:", lme4::isSingular(fit), "\n")
  invisible(fit)
}

paired_print <- function(label, data, value) {
  ag <- aggregate(data[[value]], list(Patient = data$Patient, group = data$compartment), mean)
  names(ag)[3] <- "value"
  a <- ag[ag$group == "Fibroblast", c("Patient", "value")]
  b <- ag[ag$group == "Epithelial", c("Patient", "value")]
  names(a)[2] <- "fibroblast"
  names(b)[2] <- "epithelial"
  w <- merge(a, b, by = "Patient")
  d <- w$fibroblast - w$epithelial
  wt <- suppressWarnings(wilcox.test(w$fibroblast, w$epithelial, paired = TRUE, exact = FALSE))
  cat("\n###", label, "\n")
  cat("n pairs =", nrow(w), "; mean diff =", mean(d), "; median diff =", median(d),
      "; positive/negative/zero =", sum(d > 0), "/", sum(d < 0), "/", sum(d == 0),
      "; W =", unname(wt$statistic), "; P =", wt$p.value, "\n")
}

sink(out_file)
cat("Cell-composition audit\n")
print(with(all_fib, table(Class, FAP_status)))

fit_and_print(
  "B1 all fibroblasts, manuscript-like",
  SenMayo_mean ~ FAP_status + MKI67 + (1 | Patient),
  all_fib
)
fit_and_print(
  "B1 all fibroblasts, Class-adjusted",
  SenMayo_mean ~ FAP_status + MKI67 + Class + (1 | Patient),
  all_fib
)
fit_and_print(
  "B1 tumor-only fibroblasts",
  SenMayo_mean ~ FAP_status + MKI67 + (1 | Patient),
  tumor_fib
)
fit_and_print(
  "B1 tumor-only fibroblasts, subtype-adjusted",
  SenMayo_mean ~ FAP_status + MKI67 + Cell_subtype + (1 | Patient),
  tumor_fib
)
fit_and_print(
  "B1 tumor-only fibroblasts, gene-z score",
  SenMayo_zmean ~ FAP_status + MKI67_z + (1 | Patient),
  tumor_fib
)
fit_and_print(
  "B1 tumor-only SASP",
  SASP_mean ~ FAP_status + MKI67 + (1 | Patient),
  tumor_fib
)

fit_and_print(
  "B2 all cells, manuscript-like",
  SenMayo_mean ~ compartment + MKI67 + (1 | Patient),
  all_cells
)
fit_and_print(
  "B2 all cells, Class-adjusted",
  SenMayo_mean ~ compartment + MKI67 + Class + (1 | Patient),
  all_cells
)
fit_and_print(
  "B2 tumor-only cells",
  SenMayo_mean ~ compartment + MKI67 + (1 | Patient),
  tumor_cells
)
fit_and_print(
  "B2 all cells, gene-z score",
  SenMayo_zmean ~ compartment + MKI67_z + (1 | Patient),
  all_cells
)
fit_and_print(
  "B2 tumor-only cells, gene-z score",
  SenMayo_zmean ~ compartment + MKI67_z + (1 | Patient),
  tumor_cells
)
fit_and_print(
  "B2 tumor-only SASP",
  SASP_mean ~ compartment + MKI67 + (1 | Patient),
  tumor_cells
)
fit_and_print(
  "B2 tumor-only MKI67",
  MKI67 ~ compartment + (1 | Patient),
  tumor_cells
)

paired_print("B2 all-compartment SenMayo paired", all_cells, "SenMayo_mean")
paired_print("B2 tumor-only SenMayo paired", tumor_cells, "SenMayo_mean")
paired_print("B2 all-compartment SenMayo gene-z paired", all_cells, "SenMayo_zmean")
paired_print("B2 tumor-only SenMayo gene-z paired", tumor_cells, "SenMayo_zmean")
paired_print("B2 all-compartment SASP paired", all_cells, "SASP_mean")
paired_print("B2 tumor-only SASP paired", tumor_cells, "SASP_mean")
paired_print("B2 all-compartment MKI67 paired", all_cells, "MKI67")
paired_print("B2 tumor-only MKI67 paired", tumor_cells, "MKI67")
sink()

cat("Wrote", normalizePath(out_file), "\n")
