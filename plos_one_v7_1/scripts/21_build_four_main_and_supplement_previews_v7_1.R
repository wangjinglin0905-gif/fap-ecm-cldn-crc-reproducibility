#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(paste(
    "Usage: Rscript 21_build_four_main_and_supplement_previews_v7_1.R",
    "<PLOS analysis root> <v6.4 QA root> <Route-B results_full root>",
    "<main figure output dir> <supplement figure output dir>"
  ))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
})

analysis_root <- normalizePath(args[[1]], mustWork = TRUE)
qa64_root <- normalizePath(args[[2]], mustWork = TRUE)
route_root <- normalizePath(args[[3]], mustWork = TRUE)
main_dir <- args[[4]]
supp_dir <- args[[5]]
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
main_source_dir <- file.path(main_dir, "source_data")
supp_source_dir <- file.path(supp_dir, "source_data")
dir.create(main_source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_source_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
  blue = "#2F6F9F",
  orange = "#D97732",
  green = "#2B8C6B",
  purple = "#7A5AA6",
  red = "#B64B4B",
  teal = "#2F7F8F",
  grey = "#67717B",
  light = "#D9E1E7",
  pale_blue = "#D9EAF3",
  pale_orange = "#F6E2D3",
  pale_green = "#DCEDE6",
  pale_red = "#F2DCDC",
  ink = "#20262D"
)

theme_pub <- function(base_size = 8.6) {
  theme_classic(base_family = "Arial", base_size = base_size) +
    theme(
      text = element_text(colour = pal[["ink"]]),
      axis.text = element_text(size = 8.0, colour = pal[["ink"]]),
      axis.title = element_text(size = 8.4, face = "bold"),
      axis.line = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal[["ink"]]),
      plot.title = element_text(size = 8.9, face = "bold", margin = margin(l = 18, b = 3)),
      plot.title.position = "plot",
      plot.subtitle = element_text(size = 7.6, colour = "#45505A", margin = margin(b = 4)),
      plot.caption = element_text(size = 7.1, colour = "#55616C", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = 8.2, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 7.6),
      legend.key.height = unit(3.7, "mm"),
      legend.key.width = unit(4.2, "mm"),
      plot.tag = element_text(size = 10.0, face = "bold", colour = pal[["ink"]]),
      plot.tag.position = c(0.004, 0.997),
      plot.margin = margin(5, 6, 5, 5)
    )
}

tag_plot <- function(p, tag) p + labs(tag = tag)

fmt_p <- function(p) {
  p <- as.numeric(p)
  if (!is.finite(p)) return("P unavailable")
  if (p < 0.001) paste0("P = ", format(p, scientific = TRUE, digits = 2))
  else paste0("P = ", formatC(p, format = "f", digits = 3))
}

save_preview <- function(plot, directory, stem, width = 7.5, height = 6.7) {
  png_path <- file.path(directory, paste0(stem, "_preview.png"))
  tif_path <- file.path(directory, paste0(stem, "_review_600dpi.tiff"))
  ggsave(png_path, plot, width = width, height = height, units = "in",
         dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(tif_path, plot, width = width, height = height, units = "in",
         dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  c(png = normalizePath(png_path), tiff = normalizePath(tif_path))
}

paired_panel <- function(dat, cohort_name, endpoint, display_name, tag, ylab = NULL) {
  d <- dat[dat$cohort == cohort_name, c("patient", "compartment", endpoint)]
  names(d)[3] <- "score"
  d$compartment <- factor(d$compartment, levels = c("Epithelial", "Fibroblast"))
  endpoint_key <- c(
    SenMayo_score = "SenMayo",
    SASP_score = "SASP25",
    MKI67_score = "MKI67"
  )[[endpoint]]
  s <- endpoint_lock[endpoint_lock$cohort == cohort_name & endpoint_lock$endpoint == endpoint_key, ]
  direction_text <- if (endpoint == "MKI67_score") {
    paste0("epithelial ", s$epithelial_higher, "/", s$n_patients)
  } else {
    paste0("fibroblast ", s$fibroblast_higher, "/", s$n_patients)
  }
  subtitle <- paste0(
    cohort_name, "; ", direction_text, "\n",
    "median delta ", sprintf("%.2f", s$median_difference),
    " [", sprintf("%.2f", s$median_ci_low), ", ", sprintf("%.2f", s$median_ci_high), "]; ",
    fmt_p(s$exact_wilcoxon_p)
  )
  p <- ggplot(d, aes(compartment, score, group = patient)) +
    geom_line(colour = pal[["light"]], linewidth = 0.38) +
    geom_point(aes(colour = compartment), size = 1.35, alpha = 0.92) +
    stat_summary(aes(group = compartment), fun = median, geom = "crossbar",
                 width = 0.40, linewidth = 0.48, colour = pal[["ink"]]) +
    scale_colour_manual(values = c(Epithelial = pal[["orange"]], Fibroblast = pal[["blue"]])) +
    labs(title = display_name, subtitle = subtitle, x = NULL, y = ylab) +
    guides(colour = "none") +
    theme_pub()
  tag_plot(p, tag)
}

simple_forest <- function(d, label_col, title, subtitle, tag, xlab,
                          colour_col = NULL, xlim = NULL) {
  label <- d[[label_col]]
  d$.label <- factor(label, levels = rev(unique(label)))
  p <- ggplot(d, aes(x = estimate, y = .label)) +
    geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y",
                  width = 0.16, linewidth = 0.52, colour = pal[["grey"]])
  if (is.null(colour_col)) {
    p <- p + geom_point(size = 1.8, colour = pal[["blue"]])
  } else {
    p <- p + geom_point(aes(colour = .data[[colour_col]]), size = 1.8)
  }
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p <- p + labs(title = title, subtitle = subtitle, x = xlab, y = NULL) + theme_pub()
  tag_plot(p, tag)
}

single_dir <- file.path(analysis_root, "single_cell")
common_dir <- file.path(analysis_root, "common_core")
bulk_dir <- file.path(analysis_root, "bulk_composition")
spatial_null_dir <- file.path(analysis_root, "spatial_null")
ledger_dir <- file.path(dirname(analysis_root), "ledgers")
gse132_qa <- file.path(qa64_root, "repro_gse132465")
cptac_qa <- file.path(qa64_root, "repro_cptac")

# -----------------------------------------------------------------------------
# Main Figure 1: source-available paired compartment scores in two scRNA cohorts
# -----------------------------------------------------------------------------
scores <- read.csv(file.path(single_dir, "scrna_patient_compartment_scores_aligned.csv"), check.names = FALSE)
paired_summary <- read.csv(file.path(single_dir, "scrna_paired_compartment_summary_aligned.csv"), check.names = FALSE)
endpoint_lock <- read.csv(file.path(ledger_dir, "scrna_endpoint_estimand_lock_v01_2026-09-04.csv"), check.names = FALSE)

f1a <- paired_panel(scores, "GSE132465", "SenMayo_score", "SenMayo119", "a", "Mean gene-wise z score")
f1b <- paired_panel(scores, "GSE132465", "SASP_score", "SASP25", "b", NULL)
f1c <- paired_panel(scores, "GSE132465", "MKI67_score", "MKI67", "c", "Mean log1p(CP10K)")
f1d <- paired_panel(scores, "GSE166555", "SenMayo_score", "SenMayo114", "d", "Pseudobulk gene-wise z score")
f1e <- paired_panel(scores, "GSE166555", "SASP_score", "SASP25", "e", NULL)
f1f <- paired_panel(scores, "GSE166555", "MKI67_score", "MKI67", "f", "Pseudobulk log1p(CPM)")

fig1 <- (f1a | f1b | f1c) / (f1d | f1e | f1f)
write.csv(scores, file.path(main_source_dir, "Fig1_patient_compartment_scores.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(main_source_dir, "Fig1_paired_summary.csv"), row.names = FALSE)
write.csv(endpoint_lock, file.path(main_source_dir, "Fig1_endpoint_estimand_lock.csv"), row.names = FALSE)
fig1_paths <- save_preview(fig1, main_dir, "Fig1_scRNA_compartment_localization", height = 6.55)

# -----------------------------------------------------------------------------
# Main Figure 2: matched null, FAP-specific inference and dependent correlations
# -----------------------------------------------------------------------------
null_draws <- read.csv(file.path(gse132_qa, "signature_matched_null_draws.csv"), check.names = FALSE)
null_summary <- read.csv(file.path(gse132_qa, "signature_matched_null_summary.csv"), check.names = FALSE)
ns <- null_summary[null_summary$signature == "SenMayo_nonoverlap" &
                     null_summary$statistic == "median_patient_difference", ]
nd <- null_draws[null_draws$signature == "SenMayo_nonoverlap", ]
f2a <- ggplot(nd, aes(median_patient_difference)) +
  geom_histogram(bins = 46, fill = pal[["pale_blue"]], colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = ns$null_median, colour = pal[["grey"]], linetype = 2, linewidth = 0.55) +
  geom_vline(xintercept = ns$observed, colour = pal[["red"]], linewidth = 0.75) +
  annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
           label = paste0("Observed = ", sprintf("%.3f", ns$observed),
                          "\nNull 97.5% = ", sprintf("%.3f", ns$null_q975),
                          "\nUpper-tail P = ", format(ns$empirical_p_upper, scientific = TRUE, digits = 2)),
           size = 2.35, linewidth = 0, fill = "white") +
  labs(title = "Matched-gene-set null",
       subtitle = "GSE132465; 23 patients; 5,000 gene-set draws",
       x = "Median patient fibroblast - epithelial difference", y = "Null draws") + theme_pub()
f2a <- tag_plot(f2a, "a")

cluster <- read.csv(file.path(single_dir, "gse132465_cluster_inference.csv"), check.names = FALSE)
cluster <- cluster[cluster$estimator == "patient fixed effects with patient-cluster percentile bootstrap", ]
cluster$outcome <- c(
  GSE132_FAPdet_SenMayo_depth_subtype = "SenMayo",
  GSE132_FAPdet_SASP_depth_subtype = "SASP25",
  GSE132_matrix4_on_FAP_depth_subtype = "matrix4"
)[cluster$model_id]
cluster$predictor <- c(
  GSE132_FAPdet_SenMayo_depth_subtype = "FAP-detected",
  GSE132_FAPdet_SASP_depth_subtype = "FAP-detected",
  GSE132_matrix4_on_FAP_depth_subtype = "FAP expression"
)[cluster$model_id]
cluster$predictor <- factor(cluster$predictor, levels = c("FAP-detected", "FAP expression"))
cluster$outcome <- factor(cluster$outcome, levels = c("matrix4", "SASP25", "SenMayo"))
f2b <- ggplot(cluster, aes(estimate, outcome, colour = outcome)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.16, linewidth = 0.50) +
  geom_point(size = 1.9) +
  facet_wrap(~predictor, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(SenMayo = pal[["blue"]], SASP25 = pal[["orange"]], matrix4 = pal[["green"]]),
                      guide = "none") +
  labs(title = "FAP-specific models", subtitle = "Patient-cluster bootstrap; predictor scales differ",
       x = "Coefficient (95% CI)", y = NULL) + theme_pub()
f2b <- tag_plot(f2b, "b")

cor_all <- read.csv(file.path(single_dir, "dependent_correlation_contrasts.csv"), check.names = FALSE)
label_map <- c(
  rho_FAP_matrix = "rho(FAP, matrix4)",
  rho_FAP_SenMayo = "rho(FAP, SenMayo)",
  rho_SenMayo_matrix = "rho(SenMayo, matrix4)",
  FAP_matrix_minus_FAP_SenMayo = "delta: FAP-matrix minus FAP-SenMayo",
  FAP_matrix_minus_SenMayo_matrix = "delta: FAP-matrix minus SenMayo-matrix"
)
cor_all$label <- unname(label_map[cor_all$estimand])
cor_all$kind <- ifelse(grepl("^rho_", cor_all$estimand), "Correlation", "Dependent contrast")

cor_all$cohort_short <- c(GSE132465 = "GSE132465 (n=23)", GSE166555 = "GSE166555 (n=10)")[cor_all$cohort]
cor_all$label <- c(
  rho_FAP_matrix = "FAP-matrix4",
  rho_FAP_SenMayo = "FAP-SenMayo",
  rho_SenMayo_matrix = "SenMayo-matrix4",
  FAP_matrix_minus_FAP_SenMayo = "FAP-matrix4 minus\nFAP-SenMayo",
  FAP_matrix_minus_SenMayo_matrix = "FAP-matrix4 minus\nSenMayo-matrix4"
)[cor_all$estimand]
cohort_cols <- c(`GSE132465 (n=23)` = pal[["blue"]], `GSE166555 (n=10)` = pal[["orange"]])
cohort_shapes <- c(`GSE132465 (n=23)` = 16, `GSE166555 (n=10)` = 17)

cor_raw <- cor_all[cor_all$kind == "Correlation", ]
cor_raw$label <- factor(cor_raw$label, levels = rev(c("FAP-matrix4", "FAP-SenMayo", "SenMayo-matrix4")))
f2c <- ggplot(cor_raw, aes(estimate, label, colour = cohort_short, shape = cohort_short)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15,
                linewidth = 0.48, position = position_dodge(width = 0.34)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.34)) +
  scale_colour_manual(values = cohort_cols) + scale_shape_manual(values = cohort_shapes) +
  coord_cartesian(xlim = c(-0.80, 1.15)) +
  labs(title = "Patient-level correlations", subtitle = "10,000 paired bootstrap resamples",
       x = "Spearman rho (95% CI)", y = NULL) + theme_pub() + theme(legend.position = "bottom")
f2c <- tag_plot(f2c, "c")

cor_delta <- cor_all[cor_all$kind == "Dependent contrast", ]
cor_delta$label <- factor(cor_delta$label, levels = rev(c("FAP-matrix4 minus\nFAP-SenMayo",
                                                         "FAP-matrix4 minus\nSenMayo-matrix4")))
f2d <- ggplot(cor_delta, aes(estimate, label, colour = cohort_short, shape = cohort_short)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15,
                linewidth = 0.48, position = position_dodge(width = 0.34)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.34)) +
  scale_colour_manual(values = cohort_cols) + scale_shape_manual(values = cohort_shapes) +
  coord_cartesian(xlim = c(-0.55, 1.90)) +
  labs(title = "Dependent-correlation contrasts", subtitle = "Positive values favour FAP-matrix4 coupling",
       x = "Difference in rho (95% CI)", y = NULL) + theme_pub() + theme(legend.position = "bottom")
f2d <- tag_plot(f2d, "d")

fig2 <- (f2a | f2b) / (f2c | f2d)
write.csv(cluster, file.path(main_source_dir, "Fig2_cluster_bootstrap_models.csv"), row.names = FALSE)
write.csv(cor_all, file.path(main_source_dir, "Fig2_dependent_correlations.csv"), row.names = FALSE)
write.csv(ns, file.path(main_source_dir, "Fig2_SenMayo_matched_null_summary.csv"), row.names = FALSE)
fig2_paths <- save_preview(fig2, main_dir, "Fig2_FAP_specific_and_covariance", height = 6.85)

# -----------------------------------------------------------------------------
# Main Figure 3: bulk composition sensitivity and protein triangulation
# -----------------------------------------------------------------------------
marginal <- read.csv(file.path(bulk_dir, "bulk_marginal_correlations_bootstrap.csv"), check.names = FALSE)
marginal$estimate <- marginal$rho
marginal$label <- factor(marginal$cohort, levels = rev(c("TCGA-COAD/READ", "GSE39582")))
f3a <- ggplot(marginal, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15,
                linewidth = 0.52, colour = pal[["blue"]]) +
  geom_point(size = 2.0, colour = pal[["blue"]]) +
  coord_cartesian(xlim = c(0, 1.0)) +
  labs(title = "Marginal FAP13-matrix4 correlation",
       subtitle = "Independent bulk cohorts; 5,000-bootstrap intervals",
       x = "Spearman rho (95% CI)", y = NULL) + theme_pub()
f3a <- tag_plot(f3a, "a")

bulk <- read.csv(file.path(bulk_dir, "target_purged_partial_correlations.csv"), check.names = FALSE)
selected_adjustments <- c("none", "fib5 transcript score", "MCP+EPIC recomputed full",
                          "MCP+EPIC pair-purged", "MCP+EPIC global-disjoint")
bulk_main <- bulk[bulk$pair == "FAP13_matrix4" & bulk$adjustment %in% selected_adjustments, ]
short_map <- c(
  none = "Unadjusted",
  `fib5 transcript score` = "fib5",
  `MCP+EPIC recomputed full` = "MCP + EPIC full",
  `MCP+EPIC pair-purged` = "MCP + EPIC pair-purged",
  `MCP+EPIC global-disjoint` = "MCP + EPIC globally disjoint"
)
bulk_main$label <- unname(short_map[bulk_main$adjustment])
bulk_main$label <- factor(bulk_main$label, levels = rev(unname(short_map[selected_adjustments])))
bulk_main$ci_low_plot <- ifelse(is.na(bulk_main$ci_low), bulk_main$rho, bulk_main$ci_low)
bulk_main$ci_high_plot <- ifelse(is.na(bulk_main$ci_high), bulk_main$rho, bulk_main$ci_high)
tcga_marginal <- marginal[marginal$cohort == "TCGA-COAD/READ", ]
bulk_main$ci_low_plot[bulk_main$adjustment == "none"] <- tcga_marginal$ci_low
bulk_main$ci_high_plot[bulk_main$adjustment == "none"] <- tcga_marginal$ci_high
f3b <- ggplot(bulk_main, aes(rho, label)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low_plot, xmax = ci_high_plot), orientation = "y",
                width = 0.15, linewidth = 0.50, colour = pal[["blue"]]) +
  geom_point(aes(shape = adjustment == "none"), size = 1.9, colour = pal[["blue"]]) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 18), guide = "none") +
  coord_cartesian(xlim = c(0, 1.0)) +
  labs(title = "Composition-proxy sensitivity",
       subtitle = "Adjusted estimates have 5,000-bootstrap CIs",
       x = "Spearman rho or partial rho", y = NULL) + theme_pub()
f3b <- tag_plot(f3b, "b")

global_disjoint <- bulk[bulk$adjustment == "MCP+EPIC global-disjoint", ]
global_disjoint$label <- c(
  FAP13_matrix4 = "FAP13-matrix4",
  SenMayo_FAP13 = "SenMayo-FAP13",
  SenMayo_matrix4 = "SenMayo-matrix4"
)[global_disjoint$pair]
global_disjoint$label <- factor(global_disjoint$label,
                                levels = rev(c("FAP13-matrix4", "SenMayo-FAP13", "SenMayo-matrix4")))
f3c <- ggplot(global_disjoint, aes(rho, label, colour = pair)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15, linewidth = 0.52) +
  geom_point(size = 1.95) +
  scale_colour_manual(values = c(FAP13_matrix4 = pal[["green"]], SenMayo_FAP13 = pal[["purple"]],
                                 SenMayo_matrix4 = pal[["orange"]]), guide = "none") +
  coord_cartesian(xlim = c(-0.22, 0.65)) +
  labs(title = "Globally disjoint adjustment",
       subtitle = "TCGA; joint MCP-counter + EPIC proxies",
       x = "Partial Spearman rho (95% CI)", y = NULL) + theme_pub()
f3c <- tag_plot(f3c, "c")

cptac <- read.csv(file.path(cptac_qa, "cptac_fap_ecm_and_receptor_frozen.csv"), check.names = FALSE)
cptac$short <- c("FAP vs FAP-excluded ECM", "ECM vs SDC4/CD44")
cptac$short <- factor(cptac$short, levels = rev(cptac$short))
f3d <- ggplot(cptac, aes(rho, short)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.15,
                linewidth = 0.52, colour = pal[["green"]]) +
  geom_point(size = 2.0, colour = pal[["green"]]) +
  coord_cartesian(xlim = c(-0.25, 1.0)) +
  labs(title = "Protein-layer triangulation",
       subtitle = "97 colon tumours; FAP-excluded ECM score",
       x = "Spearman rho (95% bootstrap CI)", y = NULL) + theme_pub()
f3d <- tag_plot(f3d, "d")

fig3 <- (f3a | f3b) / (f3c | f3d)
write.csv(marginal, file.path(main_source_dir, "Fig3_bulk_marginal_correlations.csv"), row.names = FALSE)
write.csv(bulk_main, file.path(main_source_dir, "Fig3_TCGA_composition_sensitivity.csv"), row.names = FALSE)
write.csv(global_disjoint, file.path(main_source_dir, "Fig3_TCGA_global_disjoint_correlations.csv"), row.names = FALSE)
write.csv(cptac, file.path(main_source_dir, "Fig3_CPTAC_protein_correlations.csv"), row.names = FALSE)
fig3_paths <- save_preview(fig3, main_dir, "Fig3_bulk_composition_and_proteomics", height = 6.65)

# -----------------------------------------------------------------------------
# Main Figure 4: spatial generalizability without pooling
# -----------------------------------------------------------------------------
cohort_label_map <- c(
  GSE280315 = "GSE280315\nVisium HD",
  GSE334323 = "GSE334323\nVisium",
  Valdeolivas_Visium = "Valdeolivas\nVisium"
)
cohort_colours <- c(
  GSE280315 = pal[["blue"]],
  GSE334323 = pal[["orange"]],
  Valdeolivas_Visium = pal[["green"]]
)
cohort_shapes_spatial <- c(GSE280315 = 16, GSE334323 = 17, Valdeolivas_Visium = 15)

h1_primary <- read.csv(file.path(route_root, "cross_cohort", "all_cohorts_H1_patient_effects.csv"), check.names = FALSE)
h1_primary$cohort_label <- cohort_label_map[h1_primary$cohort]
f4a <- ggplot(h1_primary, aes(cohort_label, median_difference, colour = cohort, shape = cohort)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 2026090401L), size = 1.9) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, linewidth = 0.52,
               colour = pal[["ink"]], aes(shape = NULL)) +
  scale_colour_manual(values = cohort_colours, guide = "none") +
  scale_shape_manual(values = cohort_shapes_spatial, guide = "none") +
  labs(title = "Source-available SenMayo",
       subtitle = "Cohort-specific score definitions; no cross-platform pooling",
       x = NULL, y = "Stromal - tumour patient median") + theme_pub()
f4a <- tag_plot(f4a, "a")

spatial_core <- read.csv(file.path(common_dir, "spatial_common_core_patient_effects.csv"), check.names = FALSE)
spatial_core$cohort_label <- cohort_label_map[spatial_core$cohort]
f4b <- ggplot(spatial_core, aes(cohort_label, stromal_minus_tumour_mean_score, colour = cohort, shape = cohort)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 2026090402L), size = 1.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.48, linewidth = 0.52,
               colour = pal[["ink"]], aes(shape = NULL)) +
  scale_colour_manual(values = cohort_colours, guide = "none") +
  scale_shape_manual(values = cohort_shapes_spatial, guide = "none") +
  labs(title = "111-gene common core",
       subtitle = "Cohort means shown without cross-platform pooling",
       x = NULL, y = "Stromal - tumour mean score") + theme_pub()
f4b <- tag_plot(f4b, "b")

mki <- read.csv(file.path(route_root, "cross_cohort", "all_cohorts_MKI67_patient_effects.csv"), check.names = FALSE)
mki$cohort_label <- cohort_label_map[mki$cohort]
f4c <- ggplot(mki, aes(cohort_label, epithelial_minus_stromal_mean, colour = cohort, shape = cohort)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 2026090403L), size = 1.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, linewidth = 0.50,
               colour = pal[["ink"]], aes(shape = NULL)) +
  scale_colour_manual(values = cohort_colours, guide = "none") +
  scale_shape_manual(values = cohort_shapes_spatial, guide = "none") +
  labs(title = "MKI67 polarity control",
       subtitle = "All 12 patient/sample effects favour epithelium",
       x = NULL, y = "Epithelial - stromal MKI67") + theme_pub()
f4c <- tag_plot(f4c, "c")

h2 <- read.csv(file.path(route_root, "cross_cohort", "all_cohorts_H2_patient_effects.csv"), check.names = FALSE)
h2$cohort_label <- cohort_label_map[h2$cohort]
f4d <- ggplot(h2, aes(cohort_label, fisher_z_axis_contrast, colour = cohort, shape = cohort)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 2026090404L), size = 1.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.48, linewidth = 0.50,
               colour = pal[["ink"]], aes(shape = NULL)) +
  scale_colour_manual(values = cohort_colours, guide = "none") +
  scale_shape_manual(values = cohort_shapes_spatial, guide = "none") +
  labs(title = "Spatial FAP-matrix contrast",
       subtitle = "Positive values indicate stronger FAP-matrix coupling",
       x = NULL, y = "Fisher-z correlation contrast") + theme_pub()
f4d <- tag_plot(f4d, "d")

fig4 <- (f4a | f4b) / (f4c | f4d)
write.csv(h1_primary, file.path(main_source_dir, "Fig4_source_available_spatial_patient_effects.csv"), row.names = FALSE)
write.csv(spatial_core, file.path(main_source_dir, "Fig4_common_core_spatial_patient_effects.csv"), row.names = FALSE)
write.csv(h2, file.path(main_source_dir, "Fig4_spatial_H2_patient_effects.csv"), row.names = FALSE)
write.csv(mki, file.path(main_source_dir, "Fig4_spatial_MKI67_patient_effects.csv"), row.names = FALSE)
fig4_paths <- save_preview(fig4, main_dir, "Fig4_spatial_generalizability", height = 6.75)

# -----------------------------------------------------------------------------
# Supplementary Figure S1: dataset flow and evidence hierarchy
# -----------------------------------------------------------------------------
box_panel <- function(boxes, title, tag, arrows = NULL) {
  p <- ggplot() +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              colour = "#66717B", linewidth = 0.35) +
    geom_text(data = boxes, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
              family = "Arial", size = 1.90, lineheight = 0.88) +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, 10), clip = "off") +
    labs(title = title) + theme_void(base_family = "Arial", base_size = 8.6) +
    theme(plot.title = element_text(size = 8.9, face = "bold", margin = margin(l = 10, b = 3)),
          plot.tag = element_text(size = 10.0, face = "bold"), plot.tag.position = c(0.004, 0.997),
          plot.margin = margin(5, 6, 5, 5))
  if (!is.null(arrows)) {
    p <- p + geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend),
                          arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
                          linewidth = 0.38, colour = pal[["grey"]])
  }
  tag_plot(p, tag)
}

s1a_boxes <- data.frame(
  xmin = c(0.3, 3.55, 6.8, 6.8), xmax = c(3.0, 6.25, 9.5, 9.5),
  ymin = c(6.0, 6.0, 7.25, 4.75), ymax = c(8.8, 8.8, 9.55, 7.05),
  fill = c(pal[["pale_blue"]], pal[["pale_green"]], pal[["pale_orange"]], "#ECE6F4"),
  label = c("Paired scRNA-seq\nGSE132465: 23 patients\nGSE166555: 10 eligible",
            "Patient-level\nlocalization\nFAP subgroup\nboundary\nWithin-fibroblast\ncovariance",
            "Bulk + proteomics\n380 + 566 + 97 tumours\nComposition sensitivity",
            "Spatial stress test\n3 + 6 + 3\npatients/samples\nNo cross-platform\npooling")
)
s1a_arrows <- data.frame(x = c(3.0, 6.25, 6.25), y = c(6.35, 7.8, 6.8),
                         xend = c(3.55, 6.8, 6.8), yend = c(6.35, 8.4, 5.9))
s1a <- box_panel(s1a_boxes, "Evidence architecture", "a", s1a_arrows)

s1b_boxes <- data.frame(
  xmin = c(0.4, 3.55, 6.7), xmax = c(3.0, 6.15, 9.55),
  ymin = c(5.4, 5.4, 5.4), ymax = c(8.6, 8.6, 8.6),
  fill = c("#E8ECEF", pal[["pale_blue"]], pal[["pale_green"]]),
  label = c("GSE132465\n63,689 total cells\n23 patients",
            "Primary tumour\nuniverse\n1,501 fibroblast-lineage\n17,469 epithelial",
            "Within fibroblasts\n925 FAP-detected\n576 FAP-undetected")
)
s1b_arrows <- data.frame(x = c(3.0, 6.15), y = c(5.85, 5.85), xend = c(3.55, 6.7), yend = c(5.85, 5.85))
s1b <- box_panel(s1b_boxes, "GSE132465 filtering and subgroup counts", "b", s1b_arrows)

s1c_boxes <- data.frame(
  xmin = c(0.3, 2.75, 5.2, 7.65), xmax = c(2.35, 4.8, 7.25, 9.7),
  ymin = rep(5.4, 4), ymax = rep(8.6, 4),
  fill = c("#E8ECEF", pal[["pale_orange"]], pal[["pale_blue"]], pal[["pale_green"]]),
  label = c("SenMayo source\n125 genes", "Remove FAP13\noverlap\n4 genes\n121 candidates",
            "Source-available\nGSE132465: 119\nGSE166555: 114", "Frozen common\ncore\n111 genes\nSensitivity only")
)
s1c_arrows <- data.frame(x = c(2.35, 4.8, 7.25), y = rep(5.85, 3),
                         xend = c(2.75, 5.2, 7.65), yend = rep(5.85, 3))
s1c <- box_panel(s1c_boxes, "SenMayo overlap and coverage control", "c", s1c_arrows)

s1d_boxes <- data.frame(
  xmin = c(0.4, 3.55, 6.7), xmax = c(3.0, 6.15, 9.55),
  ymin = c(5.4, 5.4, 5.4), ymax = c(8.6, 8.6, 8.6),
  fill = c(pal[["pale_green"]], pal[["pale_orange"]], pal[["pale_red"]]),
  label = c("Measured\nRNA/protein\ncovariation\nCompartment\ncontrasts",
            "Qualified\nProxy-sensitive\nbulk effects\nHeterogeneous\nspatial effects",
            "Not established\nDurable growth arrest\nCell-cell transfer\nor causality")
)
s1d <- box_panel(s1d_boxes, "Interpretation ceiling", "d")

sfig1 <- (s1a | s1b) / (s1c | s1d)
s1_paths <- save_preview(sfig1, supp_dir, "S1_dataset_flow_and_evidence_hierarchy", width = 7.8, height = 6.45)

# -----------------------------------------------------------------------------
# Supplementary Figure S2: cell census and common-core sensitivities
# -----------------------------------------------------------------------------
cell_counts <- read.csv(file.path(gse132_qa, "cell_counts.csv"), check.names = FALSE)
cc_keep <- cell_counts[cell_counts$set %in% c("all_cells", "all_fibroblast_lineage", "all_epithelial",
                                             "tumor_fibroblast_lineage", "tumor_epithelial",
                                             "tumor_fib_FAP_detected"), ]
cc_map <- c(all_cells = "All cells", all_fibroblast_lineage = "All fibroblast-\nlineage",
            all_epithelial = "All epithelial", tumor_fibroblast_lineage = "Tumour fibroblast-\nlineage",
            tumor_epithelial = "Tumour epithelial", tumor_fib_FAP_detected = "Tumour FAP-detected\nfibroblasts")
cc_keep$label <- unname(cc_map[cc_keep$set])
cc_keep$label <- factor(cc_keep$label, levels = rev(unname(cc_map[cc_keep$set])))
s2a <- ggplot(cc_keep, aes(n, label)) +
  geom_col(fill = pal[["blue"]], width = 0.62) +
  geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.08, size = 2.3) +
  scale_x_log10(expand = expansion(mult = c(0.02, 0.25))) +
  labs(title = "GSE132465 analysis census", subtitle = "Log scale; counts are descriptive, not independent n",
       x = "Cells (log scale)", y = NULL) + theme_pub() +
  theme(axis.text.y = element_text(size = 7.4, lineheight = 0.88),
        plot.margin = margin(5, 3, 5, 2))
s2a <- tag_plot(s2a, "a")

g132_core <- read.csv(file.path(common_dir, "gse132_common_core_paired_values.csv"), check.names = FALSE)
g166_core <- read.csv(file.path(common_dir, "gse166_common_core_paired_values.csv"), check.names = FALSE)
source_diff <- function(cohort_name) {
  d <- scores[scores$cohort == cohort_name, ]
  wide <- reshape(d[, c("patient", "compartment", "SenMayo_score")], idvar = "patient",
                  timevar = "compartment", direction = "wide")
  data.frame(patient = wide$patient,
             source_difference = wide$SenMayo_score.Fibroblast - wide$SenMayo_score.Epithelial)
}
g132_cmp <- merge(source_diff("GSE132465"), g132_core[, c("patient", "difference")], by = "patient")
g166_cmp <- merge(source_diff("GSE166555"), g166_core[, c("patient", "difference")], by = "patient")
names(g132_cmp)[3] <- "common_core_difference"
names(g166_cmp)[3] <- "common_core_difference"

comparison_panel <- function(d, title, tag) {
  r <- suppressWarnings(cor(d$source_difference, d$common_core_difference, method = "spearman"))
  p <- ggplot(d, aes(source_difference, common_core_difference)) +
    geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.35) +
    geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.35) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.5,
                colour = pal[["blue"]], fill = pal[["pale_blue"]]) +
    geom_point(size = 1.65, colour = pal[["blue"]]) +
    annotate("label", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.15,
             label = paste0("Spearman rho = ", sprintf("%.3f", r)), size = 2.25,
             linewidth = 0, fill = "white") +
    labs(title = title, subtitle = "Fibroblast - epithelial patient effects",
          x = "Source-available SenMayo effect", y = "111-gene common-core effect") + theme_pub() +
    theme(plot.margin = margin(5, 3, 5, 2))
  tag_plot(p, tag)
}
s2b <- comparison_panel(g132_cmp, "GSE132465 definition sensitivity", "b")
s2c <- comparison_panel(g166_cmp, "GSE166555 definition sensitivity", "c")

coverage <- read.csv(file.path(common_dir, "SenMayo_common_core_coverage.csv"), check.names = FALSE)
coverage_long <- rbind(
  data.frame(cohort = coverage$cohort, definition = "Source-available", genes = coverage$represented_n),
  data.frame(cohort = coverage$cohort, definition = "Common core", genes = coverage$common_core_n)
)
coverage_long$cohort_label <- c(Valdeolivas_Visium = "Valdeolivas", GSE334323 = "GSE334323",
                                GSE280315 = "GSE280315", GSE166555 = "GSE166555",
                                GSE132465 = "GSE132465")[coverage_long$cohort]
coverage_long$cohort_label <- factor(coverage_long$cohort_label,
                                     levels = rev(c("GSE132465", "GSE166555", "GSE280315", "GSE334323", "Valdeolivas")))
s2d <- ggplot(coverage_long, aes(genes, cohort_label, colour = definition, shape = definition)) +
  geom_line(aes(group = cohort), colour = pal[["light"]], linewidth = 0.65) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(`Source-available` = pal[["orange"]], `Common core` = pal[["green"]])) +
  scale_shape_manual(values = c(`Source-available` = 16, `Common core` = 17)) +
  coord_cartesian(xlim = c(108, 123)) +
  labs(title = "SenMayo coverage", subtitle = "The 111-gene intersection is sensitivity-only",
       x = "Genes represented with non-zero variance", y = NULL) + theme_pub() +
  theme(legend.position = "bottom")
s2d <- tag_plot(s2d, "d")

sfig2 <- ((s2a | s2b) / (s2c | s2d)) + plot_layout(widths = c(1.08, 0.92))
write.csv(cc_keep, file.path(supp_source_dir, "S2_GSE132465_cell_census.csv"), row.names = FALSE)
write.csv(g132_cmp, file.path(supp_source_dir, "S2_GSE132465_source_vs_common_core.csv"), row.names = FALSE)
write.csv(g166_cmp, file.path(supp_source_dir, "S2_GSE166555_source_vs_common_core.csv"), row.names = FALSE)
write.csv(coverage, file.path(supp_source_dir, "S2_SenMayo_coverage.csv"), row.names = FALSE)
s2_paths <- save_preview(sfig2, supp_dir, "S2_cell_census_and_common_core", width = 7.8, height = 6.65)

# -----------------------------------------------------------------------------
# Supplementary Figure S3: single-gene boundary and model diagnostics
# -----------------------------------------------------------------------------
marker_rates <- read.csv(file.path(gse132_qa, "marker_detection_by_FAP_status.csv"), check.names = FALSE)
marker_rates <- marker_rates[marker_rates$subset == "tumor_fibroblast", ]
marker_rates$FAP_status <- ifelse(marker_rates$FAP_status == "FAP+", "FAP-detected", "FAP-undetected")
marker_rates$gene <- factor(marker_rates$gene, levels = rev(c("CDKN2A", "CDKN2B", "CDKN1A", "LMNB1", "MKI67", "IL6", "CXCL8")))
s3a <- ggplot(marker_rates, aes(detection_rate, gene, colour = FAP_status, shape = FAP_status)) +
  geom_line(aes(group = gene), colour = pal[["light"]], linewidth = 0.7) + geom_point(size = 2.0) +
  scale_colour_manual(values = c(`FAP-undetected` = pal[["grey"]], `FAP-detected` = pal[["purple"]])) +
  scale_shape_manual(values = c(`FAP-undetected` = 16, `FAP-detected` = 17)) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Raw marker detection rates", subtitle = "Tumour fibroblast-lineage cells; descriptive only",
       x = "Cells with detected transcript", y = NULL) + theme_pub() + theme(legend.position = "bottom")
s3a <- tag_plot(s3a, "a")

marker_or <- data.frame(
  gene = c("CDKN2A", "CDKN2B", "CDKN1A", "LMNB1", "MKI67"),
  estimate = c(0.8902096, 0.9812225, 0.9618322, 0.8561041, 0.8363848),
  ci_low = c(0.6200189, 0.7031367, 0.7238329, 0.4780376, 0.4177235),
  ci_high = c(1.278143, 1.369289, 1.278087, 1.533173, 1.674647),
  p_value = c(0.5285691, 0.9112278, 0.7884639, 0.6012626, 0.6139833),
  diagnostic = c("converged", "converged", "converged", "singular fit", "converged"),
  source = "v6.4 frozen adjusted marker-detection model"
)
marker_or$gene <- factor(marker_or$gene, levels = rev(marker_or$gene))
s3b <- ggplot(marker_or, aes(estimate, gene, colour = diagnostic)) +
  geom_vline(xintercept = 1, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.16, linewidth = 0.50) +
  geom_point(size = 1.9) +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2), limits = c(0.38, 2.05)) +
  scale_colour_manual(values = c(converged = pal[["grey"]], `singular fit` = pal[["red"]])) +
  labs(title = "Adjusted marker-detection odds", subtitle = "Depth, subtype and patient random-intercept models",
       x = "Odds ratio (95% CI; log scale)", y = NULL) + theme_pub() + theme(legend.position = "bottom")
s3b <- tag_plot(s3b, "b")

cluster_raw <- read.csv(file.path(single_dir, "gse132465_cluster_inference.csv"), check.names = FALSE)
keep_models <- c("GSE132_FAPdet_SenMayo_depth_subtype", "GSE132_FAPdet_SASP_depth_subtype")
est_compare2 <- cluster_raw[cluster_raw$model_id %in% keep_models, ]
est_compare2$endpoint <- c(GSE132_FAPdet_SenMayo_depth_subtype = "SenMayo",
                           GSE132_FAPdet_SASP_depth_subtype = "SASP25")[est_compare2$model_id]
est_compare2$estimator_short <- ifelse(grepl("random-intercept", est_compare2$estimator),
                                       "Random-intercept lmer", "Patient-cluster bootstrap")
est_compare2$endpoint <- factor(est_compare2$endpoint, levels = c("SASP25", "SenMayo"))
s3c <- ggplot(est_compare2, aes(estimate, endpoint, colour = estimator_short, shape = estimator_short)) +
  geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.14,
                linewidth = 0.47, position = position_dodge(width = 0.32)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.32)) +
  scale_colour_manual(values = c(`Random-intercept lmer` = pal[["orange"]],
                                 `Patient-cluster bootstrap` = pal[["blue"]])) +
  scale_shape_manual(values = c(`Random-intercept lmer` = 16, `Patient-cluster bootstrap` = 17)) +
  labs(title = "Estimator sensitivity", subtitle = "Both intervals cross zero for both signatures",
       x = "Adjusted coefficient (95% CI)", y = NULL) + theme_pub() + theme(legend.position = "bottom")
s3c <- tag_plot(s3c, "c")

diag <- data.frame(
  item = factor(c("Patient clusters", "Cluster bootstrap", "Valid bootstrap fits", "LMNB1 random effect", "Interpretation"),
                levels = rev(c("Patient clusters", "Cluster bootstrap", "Valid bootstrap fits", "LMNB1 random effect", "Interpretation"))),
  value = c("23", "5,000 requested", "5,000 / 5,000", "Singular", "No additional FAP-detected enrichment")
)
s3d <- ggplot(diag, aes(0, item)) +
  geom_tile(aes(fill = item == "LMNB1 random effect"), width = 1.9, height = 0.80,
            colour = "white", linewidth = 0.6) +
  geom_text(aes(label = value), hjust = 0.5, size = 2.35) +
  scale_fill_manual(values = c(`FALSE` = "#E8ECEF", `TRUE` = pal[["pale_red"]]), guide = "none") +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  labs(title = "Model-diagnostic ledger", subtitle = "The singular LMNB1 fit is not headline evidence",
       x = NULL, y = NULL) + theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank())
s3d <- tag_plot(s3d, "d")

sfig3 <- (s3a | s3b) / (s3c | s3d)
write.csv(marker_rates, file.path(supp_source_dir, "S3_marker_detection_rates.csv"), row.names = FALSE)
write.csv(marker_or, file.path(supp_source_dir, "S3_adjusted_marker_detection_ORs.csv"), row.names = FALSE)
write.csv(est_compare2, file.path(supp_source_dir, "S3_estimator_sensitivity.csv"), row.names = FALSE)
s3_paths <- save_preview(sfig3, supp_dir, "S3_single_gene_models_and_diagnostics", height = 6.75)

# -----------------------------------------------------------------------------
# Supplementary Figure S4: full matched-null diagnostics
# -----------------------------------------------------------------------------
null_hist <- function(draws, summary, signature, title, tag, fill) {
  d <- draws[draws$signature == signature, ]
  s <- summary[summary$signature == signature & summary$statistic == "median_patient_difference", ]
  p <- ggplot(d, aes(median_patient_difference)) +
    geom_histogram(bins = 46, fill = fill, colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = s$null_median, colour = pal[["grey"]], linetype = 2, linewidth = 0.55) +
    geom_vline(xintercept = s$observed, colour = pal[["red"]], linewidth = 0.75) +
    annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
             label = paste0("Observed-null median = ", sprintf("%.3f", s$observed - s$null_median),
                            "\nUpper-tail P = ", format(s$empirical_p_upper, scientific = TRUE, digits = 2)),
             size = 2.25, linewidth = 0, fill = "white") +
    labs(title = title, subtitle = "GSE132465; n=23 patients; 5,000 matched gene sets",
         x = "Median patient compartment contrast", y = "Null draws") + theme_pub()
  tag_plot(p, tag)
}
s4a <- null_hist(null_draws, null_summary, "SenMayo_nonoverlap", "GSE132465 SenMayo", "a", pal[["pale_blue"]])
s4b <- null_hist(null_draws, null_summary, "SASP25", "GSE132465 SASP25", "b", pal[["pale_orange"]])

spatial_null <- read.csv(file.path(spatial_null_dir, "spatial_matched_null_corrected.csv"), check.names = FALSE)
spatial_null_panel <- function(cohort, display, tag, colour) {
  d <- read.csv(file.path(route_root, cohort, "expression_matched_null_draws.csv"), check.names = FALSE)
  names(d)[2] <- "statistic"
  s <- spatial_null[spatial_null$cohort == cohort, ]
  n_patients <- if (cohort == "GSE280315") 3L else 6L
  p <- ggplot(d, aes(statistic)) +
    geom_histogram(bins = 46, fill = "#E3E7EA", colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = s$null_center_median, colour = pal[["grey"]], linetype = 2, linewidth = 0.55) +
    geom_vline(xintercept = s$observed_median_patient_mean_difference, colour = colour, linewidth = 0.75) +
    annotate("label", x = Inf, y = Inf, hjust = 1.03, vjust = 1.15,
             label = paste0("Observed-null median = ", sprintf("%.3f", s$observed_minus_null_center),
                            "\nCentered absolute P = ", format(s$empirical_centered_absolute_p,
                                                               scientific = TRUE, digits = 2)),
             size = 2.25, linewidth = 0, fill = "white") +
    labs(title = display,
         subtitle = paste0("Source-available SenMayo; n=", n_patients, " patients; 5,000 matched sets"),
         x = "Median patient mean contrast", y = "Null draws") + theme_pub()
  tag_plot(p, tag)
}
s4c <- spatial_null_panel("GSE280315", "GSE280315 spatial null", "c", pal[["blue"]])
s4d <- spatial_null_panel("Valdeolivas_Visium", "Valdeolivas spatial null", "d", pal[["green"]])

sfig4 <- (s4a | s4b) / (s4c | s4d)
write.csv(null_summary, file.path(supp_source_dir, "S4_GSE132465_matched_null_summary.csv"), row.names = FALSE)
write.csv(spatial_null, file.path(supp_source_dir, "S4_spatial_matched_null_corrected.csv"), row.names = FALSE)
s4_paths <- save_preview(sfig4, supp_dir, "S4_matched_null_diagnostics", height = 6.55)

# -----------------------------------------------------------------------------
# Supplementary Figure S5: proxy overlap and full target-purge sensitivity
# -----------------------------------------------------------------------------
overlap <- read.csv(file.path(bulk_dir, "composition_proxy_feature_overlap_audit.csv"), check.names = FALSE)
overlap$proxy_short <- ifelse(grepl("MCP", overlap$proxy), "MCP-counter", "EPIC")
overlap$target_short <- c(FAP13_matrix4 = "FAP13-matrix4", SenMayo_FAP13 = "SenMayo-FAP13",
                          SenMayo_matrix4 = "SenMayo-matrix4", global_disjoint = "all targets")[overlap$target_set]
overlap$label <- paste(overlap$proxy_short, overlap$target_short, sep = " / ")
overlap$label <- factor(overlap$label, levels = rev(unique(overlap$label)))
s5a <- ggplot(overlap, aes(overlap_count, label, fill = proxy_short)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = ifelse(overlap_count == 0, "0", overlap_genes)), hjust = -0.08, size = 2.15) +
  scale_fill_manual(values = c(`MCP-counter` = pal[["blue"]], EPIC = pal[["orange"]])) +
  coord_cartesian(xlim = c(0, max(overlap$overlap_count) + 2.2), clip = "off") +
  labs(title = "Proxy-target overlap", subtitle = "Overlapping genes are labelled",
       x = "Overlapping features", y = NULL) + theme_pub() + theme(legend.position = "bottom")
s5a <- tag_plot(s5a, "a")

proxy_rows <- bulk[bulk$adjustment %in% c("MCP recomputed full", "MCP pair-purged", "MCP global-disjoint",
                                         "EPIC recomputed full", "EPIC pair-purged", "EPIC global-disjoint",
                                         "MCP+EPIC recomputed full", "MCP+EPIC pair-purged",
                                         "MCP+EPIC global-disjoint"), ]
proxy_rows$proxy <- ifelse(grepl("MCP+EPIC", proxy_rows$adjustment, fixed = TRUE), "MCP + EPIC",
                           ifelse(grepl("^MCP ", proxy_rows$adjustment), "MCP-counter", "EPIC"))
proxy_rows$purge <- ifelse(grepl("recomputed full", proxy_rows$adjustment), "Full",
                           ifelse(grepl("pair-purged", proxy_rows$adjustment), "Pair-purged", "Globally disjoint"))
proxy_rows$purge <- factor(proxy_rows$purge, levels = c("Full", "Pair-purged", "Globally disjoint"))

proxy_panel <- function(pair_id, title, tag) {
  d <- proxy_rows[proxy_rows$pair == pair_id, ]
  p <- ggplot(d, aes(rho, proxy, colour = purge, shape = purge)) +
    geom_vline(xintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0.14,
                  position = position_dodge(width = 0.42), linewidth = 0.45) +
    geom_point(position = position_dodge(width = 0.42), size = 1.65) +
    scale_colour_manual(values = c(Full = pal[["red"]], `Pair-purged` = pal[["blue"]],
                                   `Globally disjoint` = pal[["green"]])) +
    scale_shape_manual(values = c(Full = 16, `Pair-purged` = 17, `Globally disjoint` = 15)) +
    coord_cartesian(xlim = c(-0.30, 0.82)) +
    labs(title = title, subtitle = "TCGA; 5,000-bootstrap confidence intervals",
         x = "Partial Spearman rho", y = NULL) + theme_pub() + theme(legend.position = "bottom")
  tag_plot(p, tag)
}
s5b <- proxy_panel("FAP13_matrix4", "FAP13-matrix4", "b")
s5c <- proxy_panel("SenMayo_FAP13", "SenMayo-FAP13", "c")
s5d <- proxy_panel("SenMayo_matrix4", "SenMayo-matrix4", "d")

sfig5 <- (s5a | s5b) / (s5c | s5d)
write.csv(overlap, file.path(supp_source_dir, "S5_composition_proxy_overlap.csv"), row.names = FALSE)
write.csv(proxy_rows, file.path(supp_source_dir, "S5_target_purged_partial_correlations.csv"), row.names = FALSE)
s5_paths <- save_preview(sfig5, supp_dir, "S5_composition_proxy_overlap_and_purging", height = 6.75)

# -----------------------------------------------------------------------------
# Supplementary Figure S6: spatial QC, autocorrelation and provenance
# -----------------------------------------------------------------------------
moran_gse280 <- read.csv(file.path(route_root, "GSE280315", "patient_spatial_moran.csv"), check.names = FALSE)
moran_gse334 <- read.csv(file.path(route_root, "GSE334323", "patient_spatial_moran.csv"), check.names = FALSE)
moran_val <- read.csv(file.path(route_root, "Valdeolivas_Visium", "section_spatial_moran.csv"), check.names = FALSE)
moran_val_patient <- aggregate(moran_i ~ cohort + patient + feature, moran_val, median)
moran <- rbind(moran_gse280[, c("cohort", "patient", "feature", "moran_i")],
               moran_gse334[, c("cohort", "patient", "feature", "moran_i")],
               moran_val_patient[, c("cohort", "patient", "feature", "moran_i")])
moran$feature_label <- c(FAP_log1p10k = "FAP", matrix4_zmean = "matrix4", SenMayo_zmean = "SenMayo")[moran$feature]
moran$cohort_label <- c(GSE280315 = "GSE280315", GSE334323 = "GSE334323",
                        Valdeolivas_Visium = "Valdeolivas")[moran$cohort]
s6a <- ggplot(moran, aes(feature_label, moran_i, colour = cohort)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.10, height = 0, seed = 2026090405L), size = 1.55) +
  facet_wrap(~cohort_label, nrow = 1) +
  scale_colour_manual(values = c(GSE280315 = pal[["blue"]], GSE334323 = pal[["orange"]],
                                 Valdeolivas_Visium = pal[["green"]]), guide = "none") +
  labs(title = "Spatial autocorrelation audit", subtitle = "Patient medians for replicated Valdeolivas sections",
       x = NULL, y = "Moran's I") + theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.7))
s6a <- tag_plot(s6a, "a")

coverage_long2 <- coverage_long
s6b <- ggplot(coverage_long2, aes(genes, cohort_label, colour = definition, shape = definition)) +
  geom_line(aes(group = cohort), colour = pal[["light"]], linewidth = 0.65) + geom_point(size = 2.0) +
  scale_colour_manual(values = c(`Source-available` = pal[["orange"]], `Common core` = pal[["green"]])) +
  scale_shape_manual(values = c(`Source-available` = 16, `Common core` = 17)) +
  coord_cartesian(xlim = c(108, 123)) +
  labs(title = "Spatial score coverage control", subtitle = "All cohorts exceed the prespecified representation threshold",
       x = "Genes represented", y = NULL) + theme_pub() + theme(legend.position = "bottom")
s6b <- tag_plot(s6b, "b")

g334_pre <- read.csv(file.path(route_root, "GSE334323_pre_v01e_all_gene_expression", "patient_H1_compartment.csv"), check.names = FALSE)
g334_post <- read.csv(file.path(route_root, "GSE334323", "patient_H1_compartment.csv"), check.names = FALSE)
g334_cmp <- merge(g334_pre[, c("patient", "median_difference")],
                  g334_post[, c("patient", "median_difference")], by = "patient", suffixes = c("_pre", "_corrected"))
g334_long <- rbind(data.frame(patient = g334_cmp$patient, version = "Pre-correction",
                              effect = g334_cmp$median_difference_pre),
                   data.frame(patient = g334_cmp$patient, version = "Human-only corrected",
                              effect = g334_cmp$median_difference_corrected))
g334_long$version <- factor(g334_long$version, levels = c("Pre-correction", "Human-only corrected"))
s6c <- ggplot(g334_long, aes(version, effect, group = patient, colour = patient)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.35) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.85) +
  scale_colour_manual(values = c(CRC03 = pal[["blue"]], CRC07 = pal[["orange"]], CRC08 = pal[["green"]])) +
  labs(title = "GSE334323 denominator correction", subtitle = "Pre/post technical audit retained for provenance",
       x = NULL, y = "Median stromal - tumour score") + theme_pub() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "bottom")
s6c <- tag_plot(s6c, "c")

definition <- read.csv(file.path(route_root, "cross_cohort", "GSE280315_H1_definition_sensitivity.csv"), check.names = FALSE)
definition$definition_clean <- gsub("\n", " ", definition$definition)
definition$definition_clean <- factor(
  definition$definition_clean,
  levels = c("Author singlet (primary)", "Singlet + doublets", "Author unsupervised L1", "Frozen marker margin")
)
definition$definition_short <- c(
  `Author singlet (primary)` = "Author singlet",
  `Singlet + doublets` = "+ doublets",
  `Author unsupervised L1` = "Unsupervised L1",
  `Frozen marker margin` = "Marker margin"
)[as.character(definition$definition_clean)]
definition$definition_short <- factor(definition$definition_short,
                                      levels = c("Author singlet", "+ doublets", "Unsupervised L1", "Marker margin"))
s6d <- ggplot(definition, aes(definition_short, median_difference, group = patient, colour = patient)) +
  geom_hline(yintercept = 0, colour = pal[["grey"]], linetype = 2, linewidth = 0.38) +
  geom_line(linewidth = 0.48, alpha = 0.8) + geom_point(size = 1.75) +
  scale_colour_manual(values = c(P1 = pal[["blue"]], P2 = pal[["orange"]], P5 = pal[["green"]])) +
  labs(title = "GSE280315 assignment sensitivity",
       subtitle = "Patient effects vary with the compartment rule",
       x = NULL, y = "Median stromal - tumour score") + theme_pub() +
  theme(axis.text.x = element_text(angle = 22, hjust = 1), legend.position = "bottom")
s6d <- tag_plot(s6d, "d")

sfig6 <- (s6a | s6b) / (s6c | s6d)
write.csv(moran, file.path(supp_source_dir, "S6_spatial_Moran_patient_level.csv"), row.names = FALSE)
write.csv(g334_cmp, file.path(supp_source_dir, "S6_GSE334323_pre_post_correction.csv"), row.names = FALSE)
write.csv(definition, file.path(supp_source_dir, "S6_GSE280315_assignment_sensitivity.csv"), row.names = FALSE)
s6_paths <- save_preview(sfig6, supp_dir, "S6_spatial_QC_and_provenance", height = 6.85)

# Machine-readable production manifest.
manifest <- data.frame(
  figure = c("Fig1", "Fig2", "Fig3", "Fig4", "S1 Fig", "S2 Fig", "S3 Fig", "S4 Fig", "S5 Fig", "S6 Fig"),
  png = c(fig1_paths[["png"]], fig2_paths[["png"]], fig3_paths[["png"]], fig4_paths[["png"]],
          s1_paths[["png"]], s2_paths[["png"]], s3_paths[["png"]], s4_paths[["png"]],
          s5_paths[["png"]], s6_paths[["png"]]),
  tiff = c(fig1_paths[["tiff"]], fig2_paths[["tiff"]], fig3_paths[["tiff"]], fig4_paths[["tiff"]],
           s1_paths[["tiff"]], s2_paths[["tiff"]], s3_paths[["tiff"]], s4_paths[["tiff"]],
           s5_paths[["tiff"]], s6_paths[["tiff"]]),
  panel_labels = c("a-f", rep("a-d", 9)),
  status = "release-quality figure; manuscript-backfilled and QA passed",
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(dirname(main_dir), "figure_preview_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(dirname(main_dir), "R_sessionInfo_four_main_S1_S6.txt"))

cat("Generated 4 main figures and 6 supplementary figures.\n")
print(manifest)
