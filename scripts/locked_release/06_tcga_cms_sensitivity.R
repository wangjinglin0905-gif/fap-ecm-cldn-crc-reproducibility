options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/")
suppressPackageStartupMessages({
  library(CMSclassifier)
})

set.seed(2026080601)

input_dir <- file.path(root, "data", "public")
derived_dir <- file.path(root, "results", "analysis")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

expression_file <- file.path(input_dir, "TCGA_COADREAD_expression.txt.gz")
manifest_file <- file.path(root, "results", "analysis",
                           "TCGA_full_primary_manifest.csv")
gene_info_file <- file.path(input_dir, "Homo_sapiens.gene_info.gz")
stopifnot(all(file.exists(c(expression_file, manifest_file, gene_info_file))))

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")

score_zmean <- function(expression, genes) {
  stopifnot(all(genes %in% rownames(expression)))
  colMeans(t(scale(t(expression[genes, , drop = FALSE]))), na.rm = FALSE)
}

bootstrap_spearman <- function(x, y, seed, reps = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  estimates <- replicate(reps, {
    index <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(n = n, rho = observed,
    ci_low = unname(quantile(estimates, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(estimates, 0.975, type = 6, na.rm = TRUE)),
    p_value = test$p.value)
}

interaction_test <- function(data, label, permutations = 20000L,
                             seed = 2026080602L) {
  data <- data[complete.cases(data[, c("FAP13", "matrix4", "CMS",
                                       "project_id")]), , drop = FALSE]
  data$CMS <- droplevels(factor(data$CMS))
  data$project_id <- factor(data$project_id)
  data$x_rank <- rank(data$FAP13, ties.method = "average")
  data$y_rank <- rank(data$matrix4, ties.method = "average")
  reduced <- lm(y_rank ~ x_rank + CMS + project_id, data = data)
  full <- lm(y_rank ~ x_rank * CMS + project_id, data = data)
  observed <- anova(reduced, full)$F[2]
  parametric_p <- anova(reduced, full)$`Pr(>F)`[2]

  residuals_reduced <- residuals(reduced)
  fitted_reduced <- fitted(reduced)
  project_blocks <- split(seq_len(nrow(data)), data$project_id)
  set.seed(seed)
  permuted_f <- replicate(permutations, {
    permuted_residuals <- residuals_reduced
    for (indices in project_blocks) {
      permuted_residuals[indices] <- sample(residuals_reduced[indices],
                                             replace = FALSE)
    }
    permuted_y <- fitted_reduced + permuted_residuals
    reduced_p <- lm(permuted_y ~ x_rank + CMS + project_id, data = data)
    full_p <- lm(permuted_y ~ x_rank * CMS + project_id, data = data)
    anova(reduced_p, full_p)$F[2]
  })
  empirical_p <- (1 + sum(permuted_f >= observed, na.rm = TRUE)) /
    (permutations + 1)
  data.frame(
    classifier = label,
    n = nrow(data),
    cms_levels = nlevels(data$CMS),
    observed_f = observed,
    numerator_df = anova(reduced, full)$Df[2],
    parametric_p = parametric_p,
    permutations = permutations,
    freedman_lane_p = empirical_p,
    permutation_seed = seed
  )
}

manifest <- read.csv(manifest_file, check.names = FALSE)
expression <- read.delim(gzfile(expression_file), row.names = 1,
                         check.names = FALSE)
stopifnot(nrow(manifest) == 380L)
stopifnot(all(manifest$sample_id %in% colnames(expression)))
expression <- as.matrix(expression[, manifest$sample_id, drop = FALSE])

manifest$project_id_source <- ifelse(is.na(manifest$project_id),
                                     "TSS majority mapping", "GDC clinical")
missing_project <- which(is.na(manifest$project_id))
for (i in missing_project) {
  tss <- sub("^TCGA-([^-]+)-.*$", "\\1", manifest$sample_id[i])
  same_tss <- grepl(paste0("^TCGA-", tss, "-"), manifest$sample_id)
  known_projects <- unique(na.omit(manifest$project_id[same_tss]))
  stopifnot(length(known_projects) == 1L)
  manifest$project_id[i] <- known_projects
}

manifest$FAP13 <- score_zmean(expression, FAP13)
manifest$matrix4 <- score_zmean(expression, matrix4)

gene_info <- read.delim(gzfile(gene_info_file), quote = "", comment.char = "",
                        check.names = FALSE)
gene_info <- gene_info[gene_info$`#tax_id` == 9606,
                       c("GeneID", "Symbol", "type_of_gene"), drop = FALSE]
gene_info$GeneID <- as.character(gene_info$GeneID)
gene_info <- gene_info[!duplicated(gene_info$GeneID), , drop = FALSE]

model_genes <- unique(c(listModelGenes("RF"), listModelGenes("SSP")))
model_map <- gene_info[gene_info$GeneID %in% model_genes, , drop = FALSE]
model_map$in_expression <- model_map$Symbol %in% rownames(expression)
mapped <- model_map[model_map$in_expression, , drop = FALSE]
mapped <- mapped[!duplicated(mapped$GeneID), , drop = FALSE]
classifier_expression <- expression[mapped$Symbol, , drop = FALSE]
rownames(classifier_expression) <- mapped$GeneID
# CMSclassifier 1.0.0's naImpute() expects named data-frame columns when it
# appends model genes that are absent from the input expression matrix.
classifier_expression <- as.data.frame(classifier_expression,
                                       check.names = FALSE)
stopifnot(identical(colnames(classifier_expression), manifest$sample_id))

score_gene_ids <- unique(na.omit(
  gene_info$GeneID[match(c(FAP13, matrix4), gene_info$Symbol)]
))
classifier_expression_no_overlap <- classifier_expression[
  !rownames(classifier_expression) %in% score_gene_ids, , drop = FALSE
]

rf <- classifyCMS.RF(classifier_expression, center = TRUE, minPosterior = 0.5)
ssp <- classifyCMS.SSP(classifier_expression, minCor = 0.15, minDelta = 0.06)
rf_no_overlap <- classifyCMS.RF(classifier_expression_no_overlap, center = TRUE,
                                minPosterior = 0.5)
ssp_no_overlap <- classifyCMS.SSP(classifier_expression_no_overlap,
                                  minCor = 0.15, minDelta = 0.06)
stopifnot(all(vapply(list(rf, ssp, rf_no_overlap, ssp_no_overlap),
                     function(x) identical(rownames(x), manifest$sample_id),
                     logical(1))))

labels <- cbind(
  manifest[, c("sample_id", "patient_id", "project_id", "project_id_source",
               "FAP13", "matrix4")],
  rf, ssp,
  setNames(rf_no_overlap, paste0(names(rf_no_overlap), ".overlap_omitted")),
  setNames(ssp_no_overlap, paste0(names(ssp_no_overlap), ".overlap_omitted"))
)
write.csv(labels, file.path(derived_dir, "TCGA_CMS_classifier_labels.csv"),
          row.names = FALSE)

coverage <- data.frame(
  method = rep(c("RF", "SSP"), 2),
  input_variant = rep(c("full model input", "score-gene overlap omitted"),
                      each = 2),
  model_genes = rep(c(length(listModelGenes("RF")),
                      length(listModelGenes("SSP"))), 2),
  genes_mapped_in_TCGA = c(
    sum(listModelGenes("RF") %in% rownames(classifier_expression)),
    sum(listModelGenes("SSP") %in% rownames(classifier_expression)),
    sum(listModelGenes("RF") %in% rownames(classifier_expression_no_overlap)),
    sum(listModelGenes("SSP") %in% rownames(classifier_expression_no_overlap))
  ),
  direct_overlap_with_FAP13_matrix4 = c(
    sum(listModelGenes("RF") %in% score_gene_ids),
    sum(listModelGenes("SSP") %in% score_gene_ids), 0, 0
  ),
  high_confidence_labels = c(sum(!is.na(rf$RF.predictedCMS)),
                             sum(!is.na(ssp$SSP.predictedCMS)),
                             sum(!is.na(rf_no_overlap$RF.predictedCMS)),
                             sum(!is.na(ssp_no_overlap$SSP.predictedCMS))),
  nearest_labels = 380L,
  stringsAsFactors = FALSE
)
write.csv(coverage, file.path(derived_dir, "TCGA_CMS_classifier_coverage.csv"),
          row.names = FALSE)

make_stratified <- function(cms, classifier, seed_base) {
  data <- manifest
  data$CMS <- cms
  data <- data[!is.na(data$CMS), , drop = FALSE]
  groups <- sort(unique(data$CMS))
  rows <- do.call(rbind, lapply(seq_along(groups), function(i) {
    group <- groups[i]
    selected <- data$CMS == group
    data.frame(classifier = classifier, CMS = group,
               t(bootstrap_spearman(data$FAP13[selected],
                                    data$matrix4[selected],
                                    seed = seed_base + i)))
  }))
  rows$fdr_bh_within_classifier <- p.adjust(rows$p_value, method = "BH")
  rows
}

stratified <- rbind(
  make_stratified(rf$RF.predictedCMS, "RF high-confidence", 2026080610L),
  make_stratified(ssp$SSP.predictedCMS, "SSP high-confidence", 2026080620L),
  make_stratified(rf_no_overlap$RF.predictedCMS,
                  "RF overlap-omitted high-confidence", 2026080630L),
  make_stratified(ssp_no_overlap$SSP.predictedCMS,
                  "SSP overlap-omitted high-confidence", 2026080640L)
)
write.csv(stratified,
          file.path(derived_dir, "TCGA_CMS_stratified_correlations.csv"),
          row.names = FALSE)

rf_data <- transform(manifest, CMS = rf$RF.predictedCMS)
ssp_data <- transform(manifest, CMS = ssp$SSP.predictedCMS)
rf_no_overlap_data <- transform(manifest, CMS = rf_no_overlap$RF.predictedCMS)
ssp_no_overlap_data <- transform(manifest,
                                 CMS = ssp_no_overlap$SSP.predictedCMS)
interactions <- rbind(
  interaction_test(rf_data, "RF high-confidence", seed = 2026080631L),
  interaction_test(ssp_data, "SSP high-confidence", seed = 2026080632L),
  interaction_test(rf_no_overlap_data,
                   "RF overlap-omitted high-confidence", seed = 2026080633L),
  interaction_test(ssp_no_overlap_data,
                   "SSP overlap-omitted high-confidence", seed = 2026080634L)
)
interactions$fdr_bh_four <- p.adjust(interactions$freedman_lane_p,
                                     method = "BH")
write.csv(interactions,
          file.path(derived_dir, "TCGA_CMS_global_interaction.csv"),
          row.names = FALSE)

counts <- rbind(
  data.frame(classifier = "RF high-confidence",
             CMS = names(table(rf$RF.predictedCMS, useNA = "ifany")),
             n = as.integer(table(rf$RF.predictedCMS, useNA = "ifany"))),
  data.frame(classifier = "SSP high-confidence",
             CMS = names(table(ssp$SSP.predictedCMS, useNA = "ifany")),
             n = as.integer(table(ssp$SSP.predictedCMS, useNA = "ifany"))),
  data.frame(classifier = "RF overlap-omitted high-confidence",
             CMS = names(table(rf_no_overlap$RF.predictedCMS,
                               useNA = "ifany")),
             n = as.integer(table(rf_no_overlap$RF.predictedCMS,
                                  useNA = "ifany"))),
  data.frame(classifier = "SSP overlap-omitted high-confidence",
             CMS = names(table(ssp_no_overlap$SSP.predictedCMS,
                               useNA = "ifany")),
             n = as.integer(table(ssp_no_overlap$SSP.predictedCMS,
                                  useNA = "ifany")))
)
write.csv(counts, file.path(derived_dir, "TCGA_CMS_counts.csv"),
          row.names = FALSE)

agreement_row <- function(full, omitted, method) {
  common <- !is.na(full) & !is.na(omitted)
  data.frame(method = method, labels_compared = sum(common),
             exact_agreement_n = sum(full[common] == omitted[common]),
             exact_agreement_fraction = mean(full[common] == omitted[common]))
}
agreement <- rbind(
  agreement_row(rf$RF.predictedCMS, rf_no_overlap$RF.predictedCMS,
                "RF high-confidence"),
  agreement_row(ssp$SSP.predictedCMS, ssp_no_overlap$SSP.predictedCMS,
                "SSP high-confidence"),
  agreement_row(rf$RF.nearestCMS, rf_no_overlap$RF.nearestCMS, "RF nearest"),
  agreement_row(ssp$SSP.nearestCMS, ssp_no_overlap$SSP.nearestCMS,
                "SSP nearest")
)
write.csv(agreement, file.path(derived_dir, "TCGA_CMS_overlap_omission_agreement.csv"),
          row.names = FALSE)

provenance <- data.frame(
  item = c("CMSclassifier_package", "CMSclassifier_version",
           "CMSclassifier_source", "CMSclassifier_mirror_commit",
           "NCBI_gene_info_file", "NCBI_gene_info_MD5",
           "missing_project_resolution", "analysis_date", "R_version"),
  value = c("CMSclassifier", as.character(packageVersion("CMSclassifier")),
            "Original package metadata: Aurelien de Reynies and Justin Guinney",
            "https://github.com/sgosline/CMSclassifier/tree/9639060",
            basename(gene_info_file), unname(tools::md5sum(gene_info_file)),
            paste0(length(missing_project),
                   " sample assigned from unanimous known project within TSS"),
            "2026-08-06", R.version.string)
)
write.csv(provenance,
          file.path(derived_dir, "TCGA_CMS_provenance.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(derived_dir, "TCGA_CMS_sessionInfo.txt"))

cat("CMS sensitivity analysis complete.\n")
print(coverage, row.names = FALSE)
print(counts, row.names = FALSE)
print(stratified, row.names = FALSE)
print(interactions, row.names = FALSE)
print(agreement, row.names = FALSE)
