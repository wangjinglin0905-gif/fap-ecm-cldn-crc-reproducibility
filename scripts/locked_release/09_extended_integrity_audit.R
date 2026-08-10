options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/")
input_dir <- file.path(root, "results", "analysis")
out_dir <- input_dir

read_output <- function(name) {
  path <- file.path(input_dir, name)
  stopifnot(file.exists(path))
  read.csv(path, check.names = FALSE)
}

cms_coverage <- read_output("TCGA_CMS_classifier_coverage.csv")
cms_counts <- read_output("TCGA_CMS_counts.csv")
cms_cor <- read_output("TCGA_CMS_stratified_correlations.csv")
cms_interaction <- read_output("TCGA_CMS_global_interaction.csv")
cms_agreement <- read_output("TCGA_CMS_overlap_omission_agreement.csv")

stopifnot(all(tapply(cms_counts$n, cms_counts$classifier, sum) == 380L))
recomputed_cms_fdr <- ave(cms_cor$p_value, cms_cor$classifier,
                          FUN = function(x) p.adjust(x, method = "BH"))
stopifnot(isTRUE(all.equal(recomputed_cms_fdr,
                           cms_cor$fdr_bh_within_classifier,
                           tolerance = 1e-12)))
stopifnot(all(cms_interaction$n == tapply(
  cms_counts$n[cms_counts$CMS != "<NA>"],
  cms_counts$classifier[cms_counts$CMS != "<NA>"], sum
)[cms_interaction$classifier]))
stopifnot(isTRUE(all.equal(
  p.adjust(cms_interaction$freedman_lane_p, method = "BH"),
  cms_interaction$fdr_bh_four, tolerance = 1e-12
)))

cms_summary <- data.frame(
  metric = c("within_CMS_rho_min", "within_CMS_rho_max",
             "overlap_omitted_rho_min", "overlap_omitted_rho_max",
             "full_RF_interaction_P", "full_SSP_interaction_P",
             "omitted_RF_interaction_P", "omitted_SSP_interaction_P",
             "RF_high_confidence_agreement", "SSP_high_confidence_agreement"),
  value = c(
    min(cms_cor$rho[!grepl("overlap-omitted", cms_cor$classifier)]),
    max(cms_cor$rho[!grepl("overlap-omitted", cms_cor$classifier)]),
    min(cms_cor$rho[grepl("overlap-omitted", cms_cor$classifier)]),
    max(cms_cor$rho[grepl("overlap-omitted", cms_cor$classifier)]),
    cms_interaction$freedman_lane_p[cms_interaction$classifier ==
                                      "RF high-confidence"],
    cms_interaction$freedman_lane_p[cms_interaction$classifier ==
                                      "SSP high-confidence"],
    cms_interaction$freedman_lane_p[cms_interaction$classifier ==
                                      "RF overlap-omitted high-confidence"],
    cms_interaction$freedman_lane_p[cms_interaction$classifier ==
                                      "SSP overlap-omitted high-confidence"],
    cms_agreement$exact_agreement_fraction[cms_agreement$method ==
                                             "RF high-confidence"],
    cms_agreement$exact_agreement_fraction[cms_agreement$method ==
                                             "SSP high-confidence"]
  )
)
write.csv(cms_summary, file.path(out_dir, "extended_audit_CMS_summary.csv"),
          row.names = FALSE)

patients <- read_output("GSE132465_cross_compartment_patient_scores.csv")
pathway_cor <- read_output(
  "GSE132465_fibFAP_epithelial_pathway_correlations.csv"
)
fgsea <- read_output("GSE132465_FAP_continuous_hallmark_fgsea.csv")
fgsea_adjusted <- read_output(
  "GSE132465_FAP_cellcount_adjusted_hallmark_fgsea.csv"
)
partial_cor <- read_output(
  "GSE132465_fibFAP_epithelial_pathway_partial_correlations.csv"
)
candidate <- read_output("GSE132465_candidate_ligand_receptor_screen.csv")
gene_models <- read_output("GSE132465_FAP_continuous_gene_models.csv")
gene_models_adjusted <- read_output(
  "GSE132465_FAP_cellcount_adjusted_gene_models.csv"
)

gene_universe_matches <- all(vapply(
  intersect(unique(gene_models$compartment),
            unique(gene_models_adjusted$compartment)),
  function(compartment) {
    setequal(gene_models$gene[gene_models$compartment == compartment],
             gene_models_adjusted$gene[
               gene_models_adjusted$compartment == compartment])
  }, logical(1)
))

stopifnot(nrow(patients) == 15L, !anyDuplicated(patients$patient))
stopifnot(all(pathway_cor$n == 15L), all(candidate$n == 15L))

primary <- pathway_cor[pathway_cor$sensitivity == "primary", , drop = FALSE]
hallmark_primary <- primary[primary$method == "Epithelial Hallmark mean-z", ]
progeny_primary <- primary[primary$method == "PROGENy top500", ]
stopifnot(isTRUE(all.equal(
  p.adjust(hallmark_primary$p_value, method = "BH"),
  hallmark_primary$fdr_bh_family, tolerance = 1e-12
)))
partial_split <- split(partial_cor, partial_cor$method)
stopifnot(all(vapply(partial_split, function(x) isTRUE(all.equal(
  p.adjust(x$partial_p_value, method = "BH"), x$fdr_bh_within_method,
  tolerance = 1e-12
)), logical(1))))
stopifnot(isTRUE(all.equal(
  p.adjust(progeny_primary$p_value, method = "BH"),
  progeny_primary$fdr_bh_family, tolerance = 1e-12
)))

cell_count_cor <- do.call(rbind, lapply(c("n_fib", "n_epi"), function(v) {
  test <- suppressWarnings(cor.test(patients$fib_FAP_logCPM, patients[[v]],
                                    method = "spearman", exact = FALSE))
  data.frame(variable = v, n = nrow(patients),
             rho = unname(test$estimate), p_value = test$p.value)
}))
cell_count_cor$fdr_bh <- p.adjust(cell_count_cor$p_value, method = "BH")
write.csv(cell_count_cor,
          file.path(out_dir, "GSE132465_FAP_cell_count_sensitivity.csv"),
          row.names = FALSE)

fgsea_key <- merge(
  fgsea[!fgsea$score_genes_omitted,
        c("pathway", "compartment", "NES", "padj")],
  fgsea[fgsea$score_genes_omitted,
        c("pathway", "compartment", "NES", "padj")],
  by = c("pathway", "compartment"), suffixes = c("_full", "_omitted")
)
fgsea_key$direction_concordant <- sign(fgsea_key$NES_full) ==
  sign(fgsea_key$NES_omitted)
fgsea_key$delta_NES_omitted_minus_full <- fgsea_key$NES_omitted -
  fgsea_key$NES_full
write.csv(fgsea_key,
          file.path(out_dir, "GSE132465_fgsea_score_gene_omission_audit.csv"),
          row.names = FALSE)

fgsea_adjusted_key <- merge(
  fgsea_adjusted[!fgsea_adjusted$score_genes_omitted,
                 c("pathway", "compartment", "NES", "padj")],
  fgsea_adjusted[fgsea_adjusted$score_genes_omitted,
                 c("pathway", "compartment", "NES", "padj")],
  by = c("pathway", "compartment"), suffixes = c("_full", "_omitted")
)
fgsea_adjusted_key$direction_concordant <-
  sign(fgsea_adjusted_key$NES_full) == sign(fgsea_adjusted_key$NES_omitted)
fgsea_adjusted_key$delta_NES_omitted_minus_full <-
  fgsea_adjusted_key$NES_omitted - fgsea_adjusted_key$NES_full
write.csv(
  fgsea_adjusted_key,
  file.path(out_dir,
            "GSE132465_adjusted_fgsea_score_gene_omission_audit.csv"),
  row.names = FALSE
)

progeny_sensitivity <- pathway_cor[
  grepl("^PROGENy", pathway_cor$method),
  c("method", "pathway", "n", "rho", "ci_low", "ci_high", "p_value",
    "fdr_bh_family", "sensitivity")
]
write.csv(progeny_sensitivity,
          file.path(out_dir, "GSE132465_PROGENy_footprint_sensitivity.csv"),
          row.names = FALSE)

summarise_gene_models <- function(data, adjustment) {
  do.call(rbind, lapply(split(data, data$compartment), function(x) {
    data.frame(adjustment = adjustment, compartment = unique(x$compartment),
               genes_tested = nrow(x),
               fdr_lt_0_05 = sum(x$adj.P.Val < 0.05, na.rm = TRUE),
               positive_fdr_lt_0_05 = sum(x$adj.P.Val < 0.05 & x$logFC > 0,
                                          na.rm = TRUE),
               negative_fdr_lt_0_05 = sum(x$adj.P.Val < 0.05 & x$logFC < 0,
                                          na.rm = TRUE))
  }))
}
gene_model_summary <- rbind(
  summarise_gene_models(gene_models, "unadjusted"),
  summarise_gene_models(gene_models_adjusted,
                        "compartment_specific_cell_counts")
)
write.csv(gene_model_summary,
          file.path(out_dir, "GSE132465_gene_model_audit_summary.csv"),
          row.names = FALSE)

candidate_summary <- aggregate(
  cbind(candidate_FDR = candidate$fdr_bh_candidate_family < 0.05,
        partial_candidate_FDR =
          candidate$partial_fdr_bh_candidate_family < 0.05,
        genomewide_FDR = candidate$model_fdr_genomewide < 0.05,
        adjusted_genomewide_FDR =
          candidate$adjusted_model_fdr_genomewide < 0.05),
  by = list(family = candidate$family), FUN = sum, na.rm = TRUE
)
write.csv(candidate_summary,
          file.path(out_dir, "GSE132465_candidate_screen_audit_summary.csv"),
          row.names = FALSE)

checks <- data.frame(
  check = c("CMS classifier totals equal 380",
            "CMS BH values reproduce",
            "CMS interaction denominators match high-confidence counts",
            "CMS interaction BH values reproduce",
            "GSE132465 has 15 unique patient replication units",
            "Adjusted and unadjusted gene universes match",
            "Pathway-correlation BH values reproduce",
            "Partial pathway-correlation BH values reproduce",
            "Full versus score-gene-omitted fgsea directions concordant",
            "Adjusted full versus omitted fgsea directions concordant"),
  passed = c(TRUE, TRUE, TRUE, TRUE, TRUE, gene_universe_matches, TRUE,
             TRUE, all(fgsea_key$direction_concordant),
             all(fgsea_adjusted_key$direction_concordant)),
  stringsAsFactors = FALSE
)
write.csv(checks, file.path(out_dir, "extended_analysis_integrity_checks.csv"),
          row.names = FALSE)
stopifnot(all(checks$passed))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "extended_analysis_audit_sessionInfo.txt"))

cat("Extended analysis audit complete.\n")
print(cms_summary, row.names = FALSE)
print(cell_count_cor, row.names = FALSE)
print(gene_model_summary, row.names = FALSE)
print(candidate_summary, row.names = FALSE)
