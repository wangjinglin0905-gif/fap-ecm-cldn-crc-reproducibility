options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260805)

root <- normalizePath(getwd(), winslash = "/")
public_dir <- file.path(root, "data", "public")
derived_dir <- file.path(root, "results", "analysis")
runtime_dir <- file.path(root, "results", "analysis")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)

tcga_expression_file <- file.path(public_dir,
                                   "TCGA_COADREAD_expression.txt.gz")
tcga_manifest_file <- file.path(derived_dir,
                                "TCGA_full_primary_manifest.csv")
absolute_file <- file.path(public_dir,
                           "TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
gse_expression_file <- file.path(public_dir, "GSE39582",
                                 "GSE39582_expr.csv")
gse_probe_file <- file.path(public_dir, "GSE39582", "gene2probe.csv")
gse_meta_file <- file.path(derived_dir, "GSE39582_parsed_metadata.csv")
cptac_file <- file.path(public_dir, "cptac_protein_sdc4_cd44.csv")
patient_file <- file.path(public_dir, "GSE132465",
                          "GSE132465_raw_UMI_patient_scores.csv")

inputs <- c(tcga_expression_file, tcga_manifest_file, absolute_file,
            gse_expression_file, gse_probe_file, gse_meta_file,
            cptac_file, patient_file)
stopifnot(all(file.exists(inputs)))

absolute_expected_md5 <- "8ea2ca92c8ae58350538999dfa1174da"
absolute_observed_md5 <- unname(tools::md5sum(absolute_file))
if (!identical(tolower(absolute_observed_md5), absolute_expected_md5)) {
  stop("The frozen GDC ABSOLUTE file failed its MD5 integrity check.")
}

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
all_score_genes <- c(FAP13, matrix4)

score_zmean <- function(data, genes) {
  stopifnot(all(genes %in% names(data)))
  z <- scale(as.matrix(data[, genes, drop = FALSE]))
  rowMeans(z, na.rm = FALSE)
}

spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

bootstrap_spearman <- function(x, y, reps = 5000L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  estimates <- replicate(reps, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  estimates <- estimates[is.finite(estimates)]
  c(n = n, rho = spearman(x, y),
    ci_low = unname(quantile(estimates, 0.025, na.rm = TRUE, type = 6)),
    ci_high = unname(quantile(estimates, 0.975, na.rm = TRUE, type = 6)))
}

partial_rank_cor <- function(data, adjust_purity = TRUE) {
  keep <- complete.cases(data[, c("FAP13", "matrix4", "project_id",
                                  "purity")])
  d <- data[keep, , drop = FALSE]
  rx <- rank(d$FAP13, ties.method = "average")
  ry <- rank(d$matrix4, ties.method = "average")
  rp <- rank(d$purity, ties.method = "average")
  project <- factor(d$project_id)
  if (adjust_purity) {
    ex <- residuals(lm(rx ~ rp + project))
    ey <- residuals(lm(ry ~ rp + project))
  } else {
    ex <- residuals(lm(rx ~ project))
    ey <- residuals(lm(ry ~ project))
  }
  suppressWarnings(cor(ex, ey, method = "pearson"))
}

bootstrap_partial_rank <- function(data, adjust_purity = TRUE,
                                   reps = 5000L) {
  keep <- complete.cases(data[, c("FAP13", "matrix4", "project_id",
                                  "purity")])
  d <- data[keep, , drop = FALSE]
  observed <- partial_rank_cor(d, adjust_purity = adjust_purity)
  estimates <- replicate(reps, {
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    partial_rank_cor(d[idx, , drop = FALSE],
                     adjust_purity = adjust_purity)
  })
  c(n = nrow(d), rho = observed,
    ci_low = unname(quantile(estimates, 0.025, na.rm = TRUE, type = 6)),
    ci_high = unname(quantile(estimates, 0.975, na.rm = TRUE, type = 6)))
}

leave_one_sample_out <- function(dataset, ids, x, y) {
  ok <- is.finite(x) & is.finite(y) & !is.na(ids)
  ids <- as.character(ids[ok])
  x <- x[ok]
  y <- y[ok]
  observed <- spearman(x, y)
  estimates <- vapply(seq_along(x), function(i) {
    spearman(x[-i], y[-i])
  }, numeric(1))
  data.frame(
    dataset = dataset,
    n = length(x),
    observed_rho = observed,
    loo_min = min(estimates, na.rm = TRUE),
    loo_max = max(estimates, na.rm = TRUE),
    max_absolute_delta = max(abs(estimates - observed), na.rm = TRUE),
    most_influential_id = ids[which.max(abs(estimates - observed))],
    most_influential_loo_rho = estimates[which.max(abs(estimates - observed))]
  )
}

gene_loo_rows <- function(dataset, data) {
  full_fap <- score_zmean(data, FAP13)
  full_matrix <- score_zmean(data, matrix4)
  full_rho <- spearman(full_fap, full_matrix)
  fap_rows <- do.call(rbind, lapply(FAP13, function(gene) {
    estimate <- spearman(score_zmean(data, setdiff(FAP13, gene)),
                         full_matrix)
    data.frame(dataset = dataset, perturbed_score = "FAP13",
               omitted_gene = gene, remaining_genes = length(FAP13) - 1L,
               n = nrow(data), full_rho = full_rho, loo_rho = estimate,
               delta_from_full = estimate - full_rho)
  }))
  matrix_rows <- do.call(rbind, lapply(matrix4, function(gene) {
    estimate <- spearman(full_fap,
                         score_zmean(data, setdiff(matrix4, gene)))
    data.frame(dataset = dataset, perturbed_score = "matrix4",
               omitted_gene = gene, remaining_genes = length(matrix4) - 1L,
               n = nrow(data), full_rho = full_rho, loo_rho = estimate,
               delta_from_full = estimate - full_rho)
  }))
  rbind(fap_rows, matrix_rows)
}

# Reconstruct cohort-level gene expression used by the frozen scores.
tcga_manifest <- read.csv(tcga_manifest_file, check.names = FALSE)
tcga_full <- read.delim(gzfile(tcga_expression_file), row.names = 1,
                        check.names = FALSE)
stopifnot(all(all_score_genes %in% rownames(tcga_full)))
stopifnot(all(tcga_manifest$sample_id %in% colnames(tcga_full)))
tcga_gene <- as.data.frame(t(as.matrix(
  tcga_full[all_score_genes, tcga_manifest$sample_id, drop = FALSE]
)))

gse_probe <- read.csv(gse_expression_file, check.names = FALSE)
gse_map <- read.csv(gse_probe_file, check.names = FALSE)
gse_gene <- data.frame(sample = as.character(gse_probe[[1]]))
for (gene in unique(gse_map$gene)) {
  probes <- intersect(gse_map$probe[gse_map$gene == gene], names(gse_probe))
  if (length(probes)) {
    gse_gene[[gene]] <- rowMeans(gse_probe[, probes, drop = FALSE],
                                 na.rm = TRUE)
  }
}
gse_meta <- read.csv(gse_meta_file, check.names = FALSE)
gse_gene <- merge(gse_meta[, c("sample", "dataset")], gse_gene,
                  by = "sample", all.x = TRUE, sort = FALSE)
gse_gene <- gse_gene[!is.na(gse_gene$dataset) &
                       gse_gene$dataset != "Non Tumoral", , drop = FALSE]
stopifnot(nrow(gse_gene) == 566L)
stopifnot(all(all_score_genes %in% names(gse_gene)))

# Gene-level leave-one-out sensitivity in both bulk cohorts.
gene_loo <- rbind(gene_loo_rows("TCGA-COAD/READ", tcga_gene),
                  gene_loo_rows("GSE39582", gse_gene))
write.csv(gene_loo,
          file.path(derived_dir, "score_gene_leave_one_out.csv"),
          row.names = FALSE)

gene_loo_summary <- do.call(rbind, lapply(
  split(gene_loo, interaction(gene_loo$dataset,
                              gene_loo$perturbed_score, drop = TRUE)),
  function(d) data.frame(
    dataset = d$dataset[1],
    perturbed_score = d$perturbed_score[1],
    n_omissions = nrow(d),
    full_rho = d$full_rho[1],
    minimum_loo_rho = min(d$loo_rho),
    maximum_loo_rho = max(d$loo_rho),
    largest_absolute_delta = max(abs(d$delta_from_full)),
    most_influential_gene = d$omitted_gene[which.max(
      abs(d$delta_from_full))]
  )
))
rownames(gene_loo_summary) <- NULL
write.csv(gene_loo_summary,
          file.path(derived_dir, "score_gene_leave_one_out_summary.csv"),
          row.names = FALSE)

# Patient/sample influence in each retained orthogonal data layer.
tcga_fap13 <- score_zmean(tcga_gene, FAP13)
tcga_matrix4 <- score_zmean(tcga_gene, matrix4)
gse_fap13 <- score_zmean(gse_gene, FAP13)
gse_matrix4 <- score_zmean(gse_gene, matrix4)

cptac <- read.csv(cptac_file, check.names = FALSE)
cptac_matrix3 <- score_zmean(cptac, c("COL1A1", "COL1A2", "FN1"))
patient <- read.csv(patient_file, check.names = FALSE)
patient <- patient[patient$n_fib >= 20 & patient$n_epi >= 20 &
                     complete.cases(patient[, c("fib_FAP_logCPM",
                                                 "fib_matrix4_zmean")]),
                   , drop = FALSE]

sample_loo <- rbind(
  leave_one_sample_out("TCGA-COAD/READ", tcga_manifest$sample_id,
                       tcga_fap13, tcga_matrix4),
  leave_one_sample_out("GSE39582", gse_gene$sample,
                       gse_fap13, gse_matrix4),
  leave_one_sample_out("CPTAC-COAD", cptac$sample_id,
                       cptac$FAP, cptac_matrix3),
  leave_one_sample_out("GSE132465 supportive", patient$patient,
                       patient$fib_FAP_logCPM, patient$fib_matrix4_zmean)
)
write.csv(sample_loo,
          file.path(derived_dir, "score_sample_leave_one_out_summary.csv"),
          row.names = FALSE)

# FAP-only sensitivity is independent of the other 12 FAP13 members.
fap_only <- rbind(
  data.frame(dataset = "TCGA-COAD/READ",
             t(bootstrap_spearman(tcga_gene$FAP, tcga_matrix4))),
  data.frame(dataset = "GSE39582",
             t(bootstrap_spearman(gse_gene$FAP, gse_matrix4)))
)
write.csv(fap_only, file.path(derived_dir, "FAP_only_matrix4_sensitivity.csv"),
          row.names = FALSE)

# Formal COAD-versus-READ heterogeneity test by permuting project labels while
# preserving project sample sizes. The estimand is the difference in Spearman
# correlations, not a causal interaction coefficient.
project_ok <- !is.na(tcga_manifest$project_id)
project <- tcga_manifest$project_id[project_ok]
hx <- tcga_fap13[project_ok]
hy <- tcga_matrix4[project_ok]
rho_coad <- spearman(hx[project == "TCGA-COAD"],
                     hy[project == "TCGA-COAD"])
rho_read <- spearman(hx[project == "TCGA-READ"],
                     hy[project == "TCGA-READ"])
observed_delta <- rho_coad - rho_read
set.seed(2026080501)
n_perm <- 100000L
permuted_delta <- replicate(n_perm, {
  shuffled <- sample(project, replace = FALSE)
  spearman(hx[shuffled == "TCGA-COAD"], hy[shuffled == "TCGA-COAD"]) -
    spearman(hx[shuffled == "TCGA-READ"], hy[shuffled == "TCGA-READ"])
})
set.seed(2026080502)
bootstrap_delta <- replicate(5000L, {
  coad_idx <- sample(which(project == "TCGA-COAD"),
                     sum(project == "TCGA-COAD"), replace = TRUE)
  read_idx <- sample(which(project == "TCGA-READ"),
                     sum(project == "TCGA-READ"), replace = TRUE)
  spearman(hx[coad_idx], hy[coad_idx]) -
    spearman(hx[read_idx], hy[read_idx])
})
heterogeneity <- data.frame(
  n_coad = sum(project == "TCGA-COAD"),
  n_read = sum(project == "TCGA-READ"),
  rho_coad = rho_coad,
  rho_read = rho_read,
  delta_coad_minus_read = observed_delta,
  delta_ci_low = unname(quantile(bootstrap_delta, 0.025, type = 6)),
  delta_ci_high = unname(quantile(bootstrap_delta, 0.975, type = 6)),
  permutations = n_perm,
  permutation_p_two_sided =
    (1 + sum(abs(permuted_delta) >= abs(observed_delta))) / (n_perm + 1)
)
write.csv(heterogeneity,
          file.path(derived_dir, "TCGA_COAD_READ_heterogeneity.csv"),
          row.names = FALSE)

# Orthogonal TCGA ABSOLUTE purity sensitivity. The high-confidence analysis
# uses rows marked "called". An all-finite estimate including legacy calls is
# retained as a secondary sensitivity result.
absolute <- read.delim(absolute_file, check.names = FALSE)
purity <- merge(tcga_manifest,
                absolute[, c("array", "sample", "call status", "purity")],
                by.x = "sample_id", by.y = "array", all.x = TRUE,
                sort = FALSE)
purity$source_uuid <- "4f277128-f793-4354-a13d-30cc7fe9f6b5"
write.csv(purity[, c("sample_id", "patient_id", "project_id", "FAP13",
                     "matrix4", "call status", "purity", "source_uuid")],
          file.path(derived_dir, "TCGA_ABSOLUTE_purity_manifest.csv"),
          row.names = FALSE)

make_purity_rows <- function(d, label) {
  unadjusted <- bootstrap_spearman(d$FAP13, d$matrix4)
  project_adjusted <- bootstrap_partial_rank(d, adjust_purity = FALSE)
  purity_adjusted <- bootstrap_partial_rank(d, adjust_purity = TRUE)
  rbind(
    data.frame(population = label, model = "Unadjusted Spearman",
               t(unadjusted)),
    data.frame(population = label,
               model = "Rank residuals adjusted for project",
               t(project_adjusted)),
    data.frame(population = label,
               model = "Rank residuals adjusted for ABSOLUTE purity and project",
               t(purity_adjusted))
  )
}

called <- purity[purity[["call status"]] == "called" &
                   is.finite(purity$purity) & !is.na(purity$project_id),
                 , drop = FALSE]
all_finite <- purity[is.finite(purity$purity) & !is.na(purity$project_id),
                     , drop = FALSE]
stopifnot(nrow(called) == 351L, nrow(all_finite) == 371L)
purity_results <- rbind(
  make_purity_rows(called, "ABSOLUTE called"),
  make_purity_rows(all_finite, "All finite ABSOLUTE estimates")
)
write.csv(purity_results,
          file.path(derived_dir, "TCGA_ABSOLUTE_purity_adjustment.csv"),
          row.names = FALSE)

purity_score_correlations <- rbind(
  data.frame(population = "ABSOLUTE called", score = "FAP13",
             t(bootstrap_spearman(called$FAP13, called$purity))),
  data.frame(population = "ABSOLUTE called", score = "matrix4",
             t(bootstrap_spearman(called$matrix4, called$purity)))
)
write.csv(purity_score_correlations,
          file.path(derived_dir, "TCGA_ABSOLUTE_score_purity_correlations.csv"),
          row.names = FALSE)

# Machine-readable interpretation gates used for manuscript wording.
gates <- data.frame(
  analysis = c("Gene leave-one-out", "Sample leave-one-out",
               "ABSOLUTE purity adjustment", "COAD/READ heterogeneity"),
  prespecified_gate = c(
    "All leave-one-gene-out rho estimates remain >=0.70 with no sign reversal",
    "All leave-one-sample-out rho estimates remain >=0.70; report max delta",
    "Adjusted estimate remains positive with 95% bootstrap CI excluding zero; do not infer cell-intrinsic coupling",
    "Two-sided 100000-permutation P <0.05 indicates detectable heterogeneity; borderline results require restrained wording"
  ),
  passed = c(
    all(gene_loo$loo_rho >= 0.70),
    all(sample_loo$loo_min >= 0.70),
    with(subset(purity_results,
                population == "ABSOLUTE called" &
                  grepl("purity and project", model)),
         rho > 0 & ci_low > 0),
    heterogeneity$permutation_p_two_sided >= 0.05
  ),
  interpretation_if_passed = c(
    "The bulk association is not dependent on any single score gene",
    "No single patient or sample explains the retained cross-cohort association",
    "Covariation persists after orthogonal global-purity adjustment, but composition confounding is not eliminated",
    "No detectable COAD-versus-READ difference in the correlation"
  )
)
write.csv(gates, file.path(derived_dir, "robustness_interpretation_gates.csv"),
          row.names = FALSE)

provenance <- data.frame(
  item = c("source_page", "download_url", "GDC_UUID", "file_name",
           "file_bytes", "MD5", "download_date", "analysis_language",
           "random_seed"),
  value = c(
    "https://gdc.cancer.gov/about-data/publications/pancanatlas",
    "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
    "4f277128-f793-4354-a13d-30cc7fe9f6b5",
    basename(absolute_file), as.character(file.info(absolute_file)$size),
    absolute_observed_md5, "2026-08-05", R.version.string, "20260805"
  )
)
write.csv(provenance,
          file.path(derived_dir, "TCGA_ABSOLUTE_purity_provenance.csv"),
          row.names = FALSE)

log_file <- file.path(runtime_dir, "additional_robustness_analysis.log")
sink(log_file)
cat("ADDITIONAL ROBUSTNESS ANALYSIS\n")
cat("Run UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n\n")
cat("Gene leave-one-out summary:\n")
print(gene_loo_summary, row.names = FALSE)
cat("\nSample leave-one-out summary:\n")
print(sample_loo, row.names = FALSE)
cat("\nFAP-only sensitivity:\n")
print(fap_only, row.names = FALSE)
cat("\nCOAD/READ heterogeneity:\n")
print(heterogeneity, row.names = FALSE)
cat("\nABSOLUTE purity adjustment:\n")
print(purity_results, row.names = FALSE)
cat("\nScore-purity correlations:\n")
print(purity_score_correlations, row.names = FALSE)
cat("\nInterpretation gates:\n")
print(gates, row.names = FALSE)
cat("\nSession information:\n")
print(sessionInfo())
sink()

writeLines(capture.output(sessionInfo()),
           file.path(runtime_dir,
                     "additional_robustness_analysis_sessionInfo.txt"))
