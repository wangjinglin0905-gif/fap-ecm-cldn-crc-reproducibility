#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
stopifnot(requireNamespace("data.table", quietly = TRUE))

result_root <- "work/reproducibility/results"
cellchat_root <- "work/cellchat_reanalysis/results"
out_dir <- "work/revision_round3/reports"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(folder, file) {
  data.table::fread(file.path(result_root, folder, file), data.table = FALSE)
}

fmt <- function(x, digits = 3L) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) {
  ifelse(x < 0.001, formatC(x, format = "e", digits = 2), fmt(x, 3L))
}

ledger <- list()
add <- function(domain, finding, estimate, uncertainty, multiplicity, interpretation,
                source_file, claim_level) {
  ledger[[length(ledger) + 1L]] <<- data.frame(
    domain = domain,
    finding = finding,
    estimate = estimate,
    uncertainty = uncertainty,
    multiplicity = multiplicity,
    interpretation = interpretation,
    source_file = source_file,
    claim_level = claim_level
  )
}

tcga_cor <- read_result("L0_TCGA", "score_correlations.csv")
row <- tcga_cor[tcga_cor$comparison == "FAP_CAF_vs_de_ligand_all_profiles", ]
add(
  "TCGA bulk RNA", "Original versus ligand-excluded FAP-CAF score",
  sprintf("rho=%s; n=%d", fmt(row$rho), row$n),
  sprintf("95%% CI %s to %s; P=%s", fmt(row$ci_lower), fmt(row$ci_upper), fmt_p(row$p_value)),
  "Not part of a confirmatory test family; score robustness check",
  "Removing four candidate ligands does not materially alter the stromal/ECM score.",
  "L0_TCGA/score_correlations.csv", "validated"
)
row <- tcga_cor[tcga_cor$comparison == "de_ligand_vs_CLDN_core_all_profiles", ]
add(
  "TCGA bulk RNA", "Ligand-excluded FAP-CAF score versus CLDN core across all profiles",
  sprintf("rho=%s; n=%d", fmt(row$rho), row$n),
  sprintf("95%% CI %s to %s; P=%s", fmt(row$ci_lower), fmt(row$ci_upper), fmt_p(row$p_value)),
  "Unadjusted prespecified correlation; primary-tumour sensitivity reported separately",
  "No positive bulk association was detected; the interval excludes a moderate positive correlation.",
  "L0_TCGA/score_correlations.csv", "observed"
)

primary_cor <- read_result("L0_TCGA_primary_sensitivity", "score_correlations.csv")
row <- primary_cor[primary_cor$comparison == "de_ligand_vs_CLDN_core_primary_only", ]
add(
  "TCGA bulk RNA", "Primary-tumour-only ligand-excluded FAP-CAF versus CLDN core",
  sprintf("rho=%s; n=%d", fmt(row$rho), row$n),
  sprintf("P=%s; BH-FDR=%s", fmt_p(row$p_value), fmt_p(row$fdr_bh_three_tests)),
  "BH correction across three prespecified sensitivity correlations",
  "A weak inverse, not positive, association was detected in primary tumours.",
  "L0_TCGA_primary_sensitivity/score_correlations.csv", "observed"
)

tn <- read_result("L0_TCGA", "tumor_normal_gene_comparisons.csv")
add(
  "TCGA bulk RNA", "Tumour versus adjacent-normal expression",
  paste(sprintf("%s: median difference %s", tn$gene, fmt(tn$median_difference, 2L)), collapse = "; "),
  paste(sprintf("%s FDR=%s", tn$gene, fmt_p(tn$fdr_bh)), collapse = "; "),
  "BH correction across FAP, CLDN1, CLDN2 and CLDN4",
  "FAP, CLDN1 and CLDN2 were higher in tumour; CLDN4 was modestly lower.",
  "L0_TCGA/tumor_normal_gene_comparisons.csv", "observed"
)

cms <- read_result("L0_TCGA", "cms_kruskal_wallis.csv")
cms_summary <- read_result("L0_TCGA", "cms_score_summary.csv")
cms4 <- cms_summary[cms_summary$method == "SSP" & cms_summary$CMS == "CMS4", ]
ssp <- cms[cms$method == "SSP", ]
add(
  "TCGA CMS", "FAP-CAF score across single-sample predictor CMS groups",
  sprintf("n=%d classified tumours; CMS4 median=%s", ssp$n, fmt(cms4$median)),
  sprintf("Kruskal-Wallis P=%s", fmt_p(ssp$p_value)),
  "Overall test followed by Bonferroni-adjusted pairwise tests",
  "The FAP-CAF programme is enriched in CMS4; this is contextual rather than novel subtype discovery.",
  "L0_TCGA/cms_kruskal_wallis.csv; cms_score_summary.csv", "validated"
)

msi <- read_result("L0_TCGA", "msi_score_comparisons.csv")
msi_fap <- msi[msi$score == "FAP_CAF_de_ligand", ]
msi_cldn <- msi[msi$score == "CLDN_core", ]
add(
  "TCGA MSI", "MANTIS-defined MSI context",
  sprintf("n=%d; ligand-excluded FAP-CAF FDR=%s; CLDN core FDR=%s",
          msi_fap$n, fmt_p(msi_fap$fdr_bh_three_scores), fmt_p(msi_cldn$fdr_bh_three_scores)),
  sprintf("MSI-H=%d; MSS=%d", msi_fap$n_MSI_H, msi_fap$n_MSS),
  "BH correction across three score comparisons",
  "The stromal score was not associated with MSI after correction, whereas the CLDN core was.",
  "L0_TCGA/msi_score_comparisons.csv", "observed"
)

cptac <- read_result("CPTAC_protein", "cptac_prespecified_correlations.csv")
ecm <- cptac[cptac$x == "FAP" & cptac$y %in% c("COL1A1", "COL1A2", "FN1"), ]
add(
  "CPTAC proteomics", "FAP protein versus ECM proteins",
  paste(sprintf("%s rho=%s (n=%d)", ecm$y, fmt(ecm$rho), ecm$n), collapse = "; "),
  paste(sprintf("%s FDR=%s", ecm$y, fmt_p(ecm$fdr_bh_six_tests)), collapse = "; "),
  "BH correction across six prespecified protein correlations",
  "FAP protein covaries strongly with collagen/fibronectin proteins.",
  "CPTAC_protein/cptac_prespecified_correlations.csv", "validated"
)
cldn_protein <- cptac[grepl("CLDN", cptac$y), ]
add(
  "CPTAC proteomics", "FAP/ECM protein measures versus available CLDN proteins",
  paste(sprintf("%s vs %s: rho=%s (n=%d)", cldn_protein$x, cldn_protein$y,
                fmt(cldn_protein$rho), cldn_protein$n), collapse = "; "),
  paste(sprintf("FDR=%s", fmt_p(cldn_protein$fdr_bh_six_tests)), collapse = "; "),
  "BH correction across six prespecified protein correlations; pairwise complete cases",
  "No positive association was detected, but CLDN1 coverage was sparse and CLDN2 was unavailable.",
  "CPTAC_protein/cptac_prespecified_correlations.csv", "observed"
)

paired <- read_result("L1_TISCH_GSE166555", "paired_tumor_normal_tests.csv")
add(
  "GSE166555 single-cell pseudobulk", "Paired tumour-normal compartment-specific expression",
  paste(sprintf("%s %s: median difference %s (n=%d)", paired$analysis_group, paired$gene,
                fmt(paired$median_tumor_minus_normal), paired$paired_patients), collapse = "; "),
  paste(sprintf("%s FDR=%s", paired$gene, fmt_p(paired$fdr_bh)), collapse = "; "),
  "BH correction across four paired tests",
  "CAF-like FAP, epithelial CLDN1 and CLDN2 increased; epithelial CLDN4 decreased.",
  "L1_TISCH_GSE166555/paired_tumor_normal_tests.csv", "validated"
)

assoc <- read_result("L1_TISCH_GSE166555", "patient_level_associations.csv")
primary <- assoc[assoc$minimum_caf_cells == 20 & assoc$prespecified_primary, ]
add(
  "GSE166555 single-cell pseudobulk", "Primary cross-compartment patient-level association",
  sprintf("rho=%s; n=%d", fmt(primary$spearman_rho), primary$patients),
  sprintf("bootstrap 95%% CI %s to %s; P=%s; BH-FDR=%s",
          fmt(primary$bootstrap_ci_low), fmt(primary$bootstrap_ci_high),
          fmt_p(primary$p_value), fmt_p(primary$fdr_bh)),
  "BH correction across 16 predictor-outcome combinations at the 20-cell threshold",
  "No positive patient-level coupling was detected; the small cohort leaves a wide interval.",
  "L1_TISCH_GSE166555/patient_level_associations.csv", "observed"
)
sensitivity <- assoc[assoc$minimum_caf_cells == 5 & assoc$prespecified_primary, ]
add(
  "GSE166555 single-cell pseudobulk", "Lower CAF-cell-threshold sensitivity association",
  sprintf("rho=%s; n=%d", fmt(sensitivity$spearman_rho), sensitivity$patients),
  sprintf("bootstrap 95%% CI %s to %s; P=%s; BH-FDR=%s",
          fmt(sensitivity$bootstrap_ci_low), fmt(sensitivity$bootstrap_ci_high),
          fmt_p(sensitivity$p_value), fmt_p(sensitivity$fdr_bh)),
  "BH correction across 16 predictor-outcome combinations at the 5-cell threshold",
  "The direction remained inverse and non-significant after including one additional patient.",
  "L1_TISCH_GSE166555/patient_level_associations.csv", "sensitivity"
)

qi_moran <- read_result("L2_spatial", "spatial_morans_i.csv")
qi_co <- read_result("L2_spatial", "spatial_fap_epithelial_cooccurrence.csv")
add(
  "Qi spatial cohort", "FAP spatial autocorrelation and FAP-high/epithelial-high co-occurrence",
  sprintf("Moran I range %s-%s; co-occurrence %s%%, %s%% and %s%%",
          fmt(min(qi_moran$moran_i)), fmt(max(qi_moran$moran_i)),
          fmt(qi_co$percent_all_spots[1], 2L), fmt(qi_co$percent_all_spots[2], 2L),
          fmt(qi_co$percent_all_spots[3], 2L)),
  "Each Moran permutation P=0.001 with 999 permutations",
  "Per-patient spatial tests; no pooled spot-level inference",
  "FAP-rich regions were organised, but upper-quartile overlap with epithelial-high spots was limited.",
  "L2_spatial/spatial_morans_i.csv; spatial_fap_epithelial_cooccurrence.csv", "validated"
)

spatial <- read_result("L2_Valdeolivas_spatial", "spatial_patient_summary.csv")
for (metric in spatial$metric) {
  row <- spatial[spatial$metric == metric, ]
  add(
    "Valdeolivas raw spatial", metric,
    sprintf("median=%s; n=%d; positive=%d; negative=%d", fmt(row$median, 4L), row$patients,
            row$positive_patients, row$negative_patients),
    sprintf("two-sided exact Wilcoxon P=%s; BH-FDR=%s", fmt_p(row$wilcoxon_p), fmt_p(row$fdr_bh)),
    "BH correction across three patient-level spatial metrics",
    if (metric == "raw_bivariate_moran") {
      "Raw scores did not show a consistent direction."
    } else {
      "No positive spatial coupling was detected; negative estimates are suggestive, not confirmatory, of compartmental separation."
    },
    "L2_Valdeolivas_spatial/spatial_patient_summary.csv", "observed"
  )
}

full_chat <- read.csv(file.path(cellchat_root, "cellchat_reanalysis_summary.csv"))
full_chat <- full_chat[full_chat$analysis == "unweighted", ]
add(
  "GSE132465 CellChat", "Full pooled FAP-high myofibroblast-to-epithelial analysis",
  sprintf("%d significant pairs; %d collagen; %d FN1",
          full_chat$pair_count, full_chat$collagen_pair_count, full_chat$fibronectin_pair_count),
  "100 bootstrap iterations; CellChat threshold P<0.05",
  "Exploratory pooled cell-level inference; no patient-level replication",
  "Co-expression-compatible communication was dominated by ECM-receptor pairs.",
  "work/cellchat_reanalysis/results/cellchat_reanalysis_summary.csv", "exploratory"
)

excluded <- read.csv(file.path(cellchat_root, "cellchat_FAPhigh_to_epithelial_significant_excluding_SMC20.csv"))
add(
  "GSE132465 CellChat", "SMC20-excluded sensitivity analysis",
  sprintf("%d significant pairs; %d COLLAGEN-pathway; %d FN1-pathway",
          nrow(excluded), sum(excluded$pathway_name == "COLLAGEN"), sum(excluded$pathway_name == "FN1")),
  "162 sender cells and 16,945 epithelial cells; 100 bootstrap iterations",
  "Sensitivity exclusion because SMC20 contributed 122/284 candidate senders",
  "All six prespecified collagen/fibronectin-SDC4/CD44 pairs remained among the strongest signals; no TGFB-family, INHBA or CCL2 ligand was detected in this direction.",
  "work/cellchat_reanalysis/results/cellchat_FAPhigh_to_epithelial_significant_excluding_SMC20.csv", "sensitivity"
)

prevalence <- read_result("L3_spatial_CellChat", "spatial_cellchat_pair_prevalence.csv")
add(
  "Valdeolivas spatial CellChat", "Prespecified pair reproducibility across paired sections",
  paste(sprintf("%s: both %d/%d; any %d/%d", prevalence$interaction_name_2,
                prevalence$both_replicates_patients, prevalence$eligible_patients,
                prevalence$any_replicate_patients, prevalence$eligible_patients), collapse = "; "),
  "12 sections from six patients; 100 bootstrap iterations per section",
  "BH correction across all six prespecified pairs within each section",
  "SDC4 pairs reproduced in both sections for 6/6 patients; CD44 pairs were less consistent.",
  "L3_spatial_CellChat/spatial_cellchat_pair_prevalence.csv", "computationally_reproduced"
)

nichenet <- read_result("L3_NicheNet", "epithelial_DE_summary.csv")
add(
  "NicheNet", "Receiver differential-expression gate",
  sprintf("%s genes tested; %d positive genes at FDR<0.10", format(nichenet$tested_genes, big.mark = ","),
          nichenet$fdr_lt_0.10_positive),
  sprintf("eligible patients=%d (%d high, %d low)", nichenet$eligible_patients,
          nichenet$high_patients, nichenet$low_patients),
  "edgeR BH-FDR; ligand ranking used top 200 unadjusted positive genes only",
  "Ligand ranking is hypothesis-generating and cannot support a direct mechanism.",
  "L3_NicheNet/epithelial_DE_summary.csv", "negative"
)

progeny <- read_result("L3_pathway_TF", "progeny_high_vs_low_and_burden.csv")
dorothea <- read_result("L3_pathway_TF", "dorothea_selected_emt_tfs.csv")
add(
  "PROGENy/DoRothEA", "Patient-level pathway and EMT-related TF contrasts",
  sprintf("PROGENy minimum group FDR=%s; selected DoRothEA minimum group FDR=%s",
          fmt_p(min(progeny$wilcoxon_fdr_bh)), fmt_p(min(dorothea$wilcoxon_fdr_bh))),
  "n=13 patients",
  "BH correction performed separately for group contrasts and sender-burden correlations",
  "No tested pathway or selected EMT-related TF survived correction.",
  "L3_pathway_TF/progeny_high_vs_low_and_burden.csv; dorothea_selected_emt_tfs.csv", "negative"
)

rppa <- read_result("RPPA", "rppa_tgf_beta_correlations_from_raw_api.csv")
add(
  "TCGA RPPA", "TGF-beta RNA score versus SMAD3, SMAD4 and E-cadherin protein",
  paste(sprintf("%s rho=%s", rppa$protein, fmt(rppa$rho)), collapse = "; "),
  paste(sprintf("FDR=%s", fmt_p(rppa$fdr_bh_three_proteins)), collapse = "; "),
  "BH correction across three proteins; n=316",
  "The expected TGF-beta/SMAD context was internally consistent, but direction from FAP-high CAFs was not tested.",
  "RPPA/rppa_tgf_beta_correlations_from_raw_api.csv", "contextual"
)

ualcan <- read_result("UALCAN_CPTAC", "ualcan_cptac_primary_vs_normal.csv")
add(
  "UALCAN CPTAC interface", "Primary tumour versus normal display summaries",
  paste(sprintf("%s: %s (median difference %s)", ualcan$gene, ualcan$direction,
                fmt(ualcan$median_difference)), collapse = "; "),
  paste(sprintf("P=%s", fmt_p(ualcan$p_value)), collapse = "; "),
  "Website-reported unpaired t tests; same CPTAC source as the primary protein analysis",
  "A reproducibility check of the public interface, not an independent cohort.",
  "UALCAN_CPTAC/ualcan_cptac_primary_vs_normal.csv", "contextual"
)

ledger_table <- data.table::rbindlist(ledger, fill = TRUE)
data.table::fwrite(ledger_table, file.path(out_dir, "verified_evidence_ledger_r461.csv"), na = "")

md <- c(
  "# Verified evidence ledger (R 4.6.1)",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "All values below were read from the completed R 4.6.1 result tables. Claim levels distinguish observed, validated/reproduced, sensitivity, contextual, exploratory and negative findings.",
  ""
)
for (i in seq_len(nrow(ledger_table))) {
  md <- c(
    md,
    paste0("## ", ledger_table$domain[i], ": ", ledger_table$finding[i]),
    "",
    paste0("- Estimate: ", ledger_table$estimate[i]),
    paste0("- Uncertainty: ", ledger_table$uncertainty[i]),
    paste0("- Multiplicity/unit: ", ledger_table$multiplicity[i]),
    paste0("- Interpretation: ", ledger_table$interpretation[i]),
    paste0("- Claim level: `", ledger_table$claim_level[i], "`"),
    paste0("- Source: `", ledger_table$source_file[i], "`"),
    ""
  )
}
writeLines(md, file.path(out_dir, "verified_evidence_ledger_r461.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "verified_evidence_ledger_sessionInfo.txt"))
cat("Evidence ledger rows:", nrow(ledger_table), "\n")
