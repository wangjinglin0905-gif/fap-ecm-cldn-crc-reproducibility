#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: Rscript 13b_recover_gse39582_bulk_scores.R",
    "<GSE39582_expr.csv> <gene2probe.csv> <gse39582_sample_dataset.csv>",
    "<tcga_primary_target_scores.csv> <output_dir>"
  ))
}

expr_path <- normalizePath(args[[1]], mustWork = TRUE)
map_path <- normalizePath(args[[2]], mustWork = TRUE)
dataset_path <- normalizePath(args[[3]], mustWork = TRUE)
tcga_path <- normalizePath(args[[4]], mustWork = TRUE)
out_dir <- args[[5]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expr <- read.csv(expr_path, check.names = FALSE)
probe_map <- read.csv(map_path, check.names = FALSE)
dataset_map <- read.csv(dataset_path, check.names = FALSE)
tcga <- read.csv(tcga_path, check.names = FALSE)

stopifnot(names(expr)[1] == "sample")
stopifnot(all(c("gene", "probe") %in% names(probe_map)))
stopifnot(all(c("sample", "dataset") %in% names(dataset_map)))
stopifnot(all(c("sample", "FAP13", "matrix4") %in% names(tcga)))

FAP13 <- c(
  "FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
  "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A"
)
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
receptor2 <- c("SDC4", "CD44")
required_genes <- unique(c(FAP13, matrix4, receptor2))

gene_expr <- data.frame(sample = as.character(expr[[1]]), check.names = FALSE)
for (gene in required_genes) {
  probes <- intersect(probe_map$probe[probe_map$gene == gene], names(expr))
  if (!length(probes)) stop("No mapped probe was available for ", gene)
  gene_expr[[gene]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}

gse <- merge(dataset_map, gene_expr, by = "sample", all.x = TRUE, sort = FALSE)
gse <- gse[!is.na(gse$dataset) & gse$dataset != "Non Tumoral", ]
stopifnot(nrow(gse) == 566L)

zmean <- function(data, genes) {
  rowMeans(scale(data[, genes, drop = FALSE]), na.rm = TRUE)
}
gse$FAP13 <- zmean(gse, FAP13)
gse$matrix4 <- zmean(gse, matrix4)
gse$receptor2 <- zmean(gse, receptor2)

bootstrap_spearman <- function(data, cohort, seed, reps = 5000L) {
  keep <- is.finite(data$FAP13) & is.finite(data$matrix4)
  x <- data$FAP13[keep]
  y <- data$matrix4[keep]
  n <- length(x)
  estimate <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  boot <- replicate(reps, {
    index <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    cohort = cohort,
    n = n,
    rho = estimate,
    ci_low = unname(quantile(boot, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(boot, 0.975, type = 6, na.rm = TRUE)),
    p_value = test$p.value,
    bootstrap_replicates = reps,
    stringsAsFactors = FALSE
  )
}

marginal <- rbind(
  bootstrap_spearman(tcga, "TCGA-COAD/READ", 2026081501L),
  bootstrap_spearman(gse, "GSE39582", 2026081503L)
)

stopifnot(abs(marginal$rho[marginal$cohort == "TCGA-COAD/READ"] - 0.929585676530053) < 1e-10)
stopifnot(abs(marginal$rho[marginal$cohort == "GSE39582"] - 0.9134) < 1e-3)

write.csv(
  gse[, c("sample", "dataset", "FAP13", "matrix4", "receptor2")],
  file.path(out_dir, "gse39582_tumour_sample_scores.csv"),
  row.names = FALSE
)
write.csv(
  marginal,
  file.path(out_dir, "bulk_marginal_correlations_bootstrap.csv"),
  row.names = FALSE
)

provenance <- data.frame(
  role = c("expression", "probe map", "tumour/non-tumour map", "TCGA scores"),
  source_file = basename(c(expr_path, map_path, dataset_path, tcga_path)),
  input_rows = c(nrow(expr), nrow(probe_map), nrow(dataset_map), nrow(tcga)),
  note = c(
    "Public GSE39582 probe-level expression matrix; locally archived input",
    "Frozen gene-to-probe mapping used by the locked release pipeline",
    "Frozen sample classification; 19 non-tumour samples excluded",
    "Frozen TCGA primary tumour signature scores"
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(out_dir, "bulk_marginal_source_provenance.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_gse39582_recovery.txt"))

cat("Recovered GSE39582 sample-level scores and bootstrap estimates.\n")
print(marginal)
