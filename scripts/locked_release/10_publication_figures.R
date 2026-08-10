options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ragg)
})

set.seed(2026080605)

runtime_dir <- file.path(root, "results", "analysis")
extended_dir <- runtime_dir
package_derived_dir <- runtime_dir
out_dir <- file.path(root, "results", "publication")
figure_dir <- file.path(out_dir, "figures")
source_dir <- file.path(out_dir, "source_data")
table_dir <- file.path(out_dir, "tables")
qa_dir <- file.path(out_dir, "qa")
for (path in c(out_dir, figure_dir, source_dir, table_dir, qa_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_csv <- function(path) {
  stopifnot(file.exists(path))
  read.csv(path, check.names = FALSE)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2),
         sprintf("%.3f", p))
}

bootstrap_spearman <- function(x, y, seed, reps = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  estimate <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  boot <- replicate(reps, {
    index <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(
    n = n,
    rho = estimate,
    ci_low = unname(quantile(boot, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(boot, 0.975, type = 6, na.rm = TRUE)),
    p_value = test$p.value
  )
}

palette <- c(
  navy = "#315B7D", blue = "#4F86A6", teal = "#2F8F83",
  gold = "#C08A35", red = "#B65B4B", grey = "#72777C",
  light_grey = "#D9DDDF", dark = "#25282A"
)

theme_pub <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4, colour = palette["dark"]),
      plot.title = element_text(size = base_size + 0.6, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2,
                                   colour = palette["grey"],
                                   margin = margin(b = 2)),
      plot.tag = element_text(size = base_size + 1.2, face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(5, 6, 5, 5)
    )
}

theme_set(theme_pub())

save_pub <- function(plot, stem, width_mm = 183, height_mm = 125,
                     tiff_dpi = 600L, png_dpi = 300L) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  png_file <- paste0(stem, ".png")
  tiff_file <- paste0(stem, ".tiff")
  svg_file <- paste0(stem, ".svg")
  pdf_file <- paste0(stem, ".pdf")

  agg_png(png_file, width = width_in, height = height_in, units = "in",
          res = png_dpi, background = "white")
  print(plot)
  dev.off()

  agg_tiff(tiff_file, width = width_in, height = height_in, units = "in",
           res = tiff_dpi, compression = "lzw", background = "white")
  print(plot)
  dev.off()

  grDevices::svg(svg_file, width = width_in, height = height_in,
                 family = "Arial", bg = "white", onefile = TRUE)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(pdf_file, width = width_in, height = height_in,
                       family = "Arial", bg = "white", onefile = TRUE)
  print(plot)
  dev.off()

  data.frame(
    file = normalizePath(c(png_file, tiff_file, svg_file, pdf_file),
                         winslash = "/"),
    format = c("PNG", "TIFF", "SVG", "PDF"),
    width_mm = width_mm,
    height_mm = height_mm,
    dpi = c(png_dpi, tiff_dpi, NA, NA),
    bytes = file.info(c(png_file, tiff_file, svg_file, pdf_file))$size
  )
}

scatter_panel <- function(data, x, y, title, subtitle, x_label, y_label,
                          colour) {
  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(shape = 21, size = 1.5, stroke = 0.18, alpha = 0.65,
               colour = "white", fill = colour) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = colour, fill = colour, alpha = 0.13,
                linewidth = 0.5) +
    labs(title = title, subtitle = subtitle, x = x_label, y = y_label) +
    theme_pub()
}

forest_panel <- function(data, title, subtitle = NULL, x_label,
                         colour_by = NULL, x_limits = NULL) {
  p <- ggplot(data, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["light_grey"])
  if (is.null(colour_by)) {
    p <- p +
      geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.15,
                    orientation = "y", linewidth = 0.5,
                    colour = palette["grey"]) +
      geom_point(size = 2, colour = palette["navy"])
  } else {
    p <- p +
      geom_errorbar(aes(xmin = ci_low, xmax = ci_high,
                        colour = .data[[colour_by]]),
                    width = 0.15, orientation = "y", linewidth = 0.55) +
      geom_point(aes(colour = .data[[colour_by]]), size = 2) +
      scale_colour_manual(values = c(
        "Matrix covariation" = unname(palette["teal"]),
        "Receptor comparator" = unname(palette["red"]),
        "Full classifier" = unname(palette["navy"]),
        "Overlap omitted" = unname(palette["gold"])
      ))
  }
  p <- p +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL,
         colour = NULL) +
    theme_pub() +
    theme(legend.position = "bottom")
  if (!is.null(x_limits)) p <- p + coord_cartesian(xlim = x_limits)
  p
}

# Figure 1: cross-platform matrix covariation and receptor boundaries.
bulk <- read_csv(file.path(runtime_dir, "Figure1_bulk_source_data.csv"))
cptac <- read_csv(file.path(runtime_dir, "Figure1_CPTAC_source_data.csv"))
tcga <- bulk[bulk$cohort == "TCGA-COADREAD", ]
gse <- bulk[bulk$cohort == "GSE39582", ]

fig1_stats <- rbind(
  cbind(cohort = "TCGA-COAD/READ", comparison = "FAP13 vs matrix4",
        family = "Matrix covariation",
        bootstrap_spearman(tcga$FAP13, tcga$matrix4, 2026081501L)),
  cbind(cohort = "TCGA-COAD/READ", comparison = "FAP13 vs receptor2",
        family = "Receptor comparator",
        bootstrap_spearman(tcga$FAP13, tcga$receptor2, 2026081502L)),
  cbind(cohort = "GSE39582", comparison = "FAP13 vs matrix4",
        family = "Matrix covariation",
        bootstrap_spearman(gse$FAP13, gse$matrix4, 2026081503L)),
  cbind(cohort = "GSE39582", comparison = "FAP13 vs receptor2",
        family = "Receptor comparator",
        bootstrap_spearman(gse$FAP13, gse$receptor2, 2026081504L)),
  cbind(cohort = "CPTAC-COAD", comparison = "FAP vs matrix3",
        family = "Matrix covariation",
        bootstrap_spearman(cptac$FAP, cptac$matrix3, 2026081505L)),
  cbind(cohort = "CPTAC-COAD", comparison = "matrix3 vs receptor2",
        family = "Receptor comparator",
        bootstrap_spearman(cptac$matrix3, cptac$receptor2, 2026081506L))
)
fig1_stats$p_value <- as.numeric(fig1_stats$p_value)
fig1_stats$fdr_bh_within_platform_family <- ave(
  fig1_stats$p_value, fig1_stats$cohort,
  FUN = function(x) p.adjust(x, method = "BH")
)
write.csv(fig1_stats, file.path(source_dir, "Figure1_statistics.csv"),
          row.names = FALSE)
write.csv(bulk, file.path(source_dir, "Figure1_bulk_scores.csv"),
          row.names = FALSE)
write.csv(cptac, file.path(source_dir, "Figure1_CPTAC_scores.csv"),
          row.names = FALSE)

stat_text <- function(cohort, comparison) {
  row <- fig1_stats[fig1_stats$cohort == cohort &
                      fig1_stats$comparison == comparison, ]
  sprintf("rho = %.2f (95%% CI %.2f-%.2f); n = %d",
          row$rho, row$ci_low, row$ci_high, row$n)
}

p1a <- scatter_panel(
  tcga, "FAP13", "matrix4", "TCGA-COAD/READ",
  stat_text("TCGA-COAD/READ", "FAP13 vs matrix4"),
  "FAP13 score (z mean)", "matrix4 score (z mean)", palette["navy"]
)
p1b <- scatter_panel(
  gse, "FAP13", "matrix4", "GSE39582",
  stat_text("GSE39582", "FAP13 vs matrix4"),
  "FAP13 score (z mean)", "matrix4 score (z mean)", palette["blue"]
)
p1c <- scatter_panel(
  cptac, "FAP", "matrix3", "CPTAC-COAD proteome",
  stat_text("CPTAC-COAD", "FAP vs matrix3"),
  "FAP protein (z)", "matrix3 protein score (z mean)", palette["teal"]
)
fig1_forest <- fig1_stats
fig1_forest$estimate <- fig1_forest$rho
fig1_forest$label <- paste(fig1_forest$cohort,
                           ifelse(fig1_forest$family == "Matrix covariation",
                                  "matrix", "SDC4/CD44"), sep = ": ")
fig1_forest$label <- factor(fig1_forest$label,
                            levels = rev(fig1_forest$label))
p1d <- forest_panel(
  fig1_forest, "Matrix signal and receptor comparator",
  "Spearman rho with bootstrap 95% CI",
  "Spearman rho (95% CI)", colour_by = "family", x_limits = c(-0.2, 1)
)
fig1 <- (p1a | p1b | p1c) / p1d +
  plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "a")

# Figure 2: anatomical/purity robustness and matched empirical nulls.
tcga_cor <- read_csv(file.path(runtime_dir, "TCGA_full_primary_correlations.csv"))
purity <- read_csv(file.path(package_derived_dir,
                             "TCGA_ABSOLUTE_purity_adjustment.csv"))
null_dist <- read_csv(file.path(runtime_dir,
                                "TCGA_matched_null_extended_distribution.csv"))
null_summary <- read_csv(file.path(runtime_dir,
                                   "TCGA_matched_null_extended_summary.csv"))

robustness <- tcga_cor[
  tcga_cor$comparison == "FAP13 vs matrix4" &
    tcga_cor$cohort %in% c("TCGA-COAD", "TCGA-READ"),
  c("cohort", "n", "rho", "ci_low", "ci_high")
]
purity_row <- purity[
  purity$population == "ABSOLUTE called" &
    purity$model == "Rank residuals adjusted for ABSOLUTE purity and project", ]
robustness <- rbind(
  robustness,
  data.frame(cohort = "Purity/project-adjusted", n = purity_row$n,
             rho = purity_row$rho, ci_low = purity_row$ci_low,
             ci_high = purity_row$ci_high)
)
robustness$estimate <- robustness$rho
robustness$label <- factor(
  paste0(robustness$cohort, " (n=", robustness$n, ")"),
  levels = rev(paste0(robustness$cohort, " (n=", robustness$n, ")"))
)
write.csv(robustness, file.path(source_dir, "Figure2_robustness.csv"),
          row.names = FALSE)
write.csv(null_dist, file.path(source_dir, "Figure2_empirical_null_draws.csv"),
          row.names = FALSE)
write.csv(null_summary, file.path(source_dir, "Figure2_empirical_null_summary.csv"),
          row.names = FALSE)

p2a <- forest_panel(
  robustness, "Anatomical and purity sensitivity",
  "FAP13-matrix4 rank association",
  "Spearman or adjusted rank-residual rho", x_limits = c(0.75, 1)
)
null_labels <- c(
  random_vs_random = "Random vs random",
  fixed_FAP13_vs_random_matrix4 = "FAP13 vs random matrix4",
  random_FAP13_vs_fixed_matrix4 = "Random FAP13-like vs matrix4"
)
null_dist$null_label <- factor(null_labels[null_dist$null_type],
                               levels = unname(null_labels))
null_summary$null_label <- factor(null_labels[null_summary$null_type],
                                  levels = unname(null_labels))
p2b <- ggplot(null_dist, aes(x = rho)) +
  geom_density(fill = palette["light_grey"], colour = palette["grey"],
               linewidth = 0.45, adjust = 1.05) +
  geom_vline(data = null_summary, aes(xintercept = observed_rho),
             colour = palette["red"], linewidth = 0.65) +
  geom_text(
    data = null_summary,
    aes(x = observed_rho, y = Inf,
        label = paste0("Observed rho = ", sprintf("%.2f", observed_rho),
                       "\nEmpirical P = ",
                       format(empirical_p_two_sided, scientific = TRUE,
                              digits = 1))),
    inherit.aes = FALSE, hjust = 1.03, vjust = 1.15,
    size = 2.25, family = "Arial", colour = palette["red"]
  ) +
  facet_wrap(~null_label, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(-0.4, 1), clip = "off") +
  labs(title = "Matched empirical gene-set nulls",
       subtitle = "10,000 draws per null; red line denotes observed rho",
       x = "Random-set Spearman rho", y = "Density") +
  theme_pub() +
  theme(strip.text = element_text(hjust = 0),
        plot.margin = margin(5, 8, 5, 5))
fig2 <- (p2a | p2b) + plot_layout(widths = c(0.9, 1.45)) +
  plot_annotation(tag_levels = "a")

# Figure 3: patient-level single-cell pseudobulk and paired support.
patient <- read_csv(file.path(runtime_dir, "Figure2_patient_source_data.csv"))
sc_cor <- read_csv(file.path(runtime_dir, "GSE132465",
                             "GSE132465_raw_UMI_correlations.csv"))
paired <- read_csv(file.path(runtime_dir, "GSE132465",
                             "GSE132465_paired_tumour_normal_scores.csv"))
paired_tests <- read_csv(file.path(runtime_dir, "GSE132465",
                                   "GSE132465_paired_tumour_normal_tests.csv"))
write.csv(patient, file.path(source_dir, "Figure3_patient_pseudobulk_scores.csv"),
          row.names = FALSE)
write.csv(sc_cor, file.path(source_dir, "Figure3_patient_correlations.csv"),
          row.names = FALSE)
write.csv(paired, file.path(source_dir, "Figure3_paired_scores.csv"),
          row.names = FALSE)
write.csv(paired_tests, file.path(source_dir, "Figure3_paired_tests.csv"),
          row.names = FALSE)

main_sc <- sc_cor[sc_cor$comparison ==
                    "Fibroblast-lineage FAP vs fibroblast-lineage matrix4", ]
p3a <- scatter_panel(
  patient, "fib_FAP_logCPM", "fib_matrix4_zmean",
  "Patient-level fibroblast pseudobulk",
  sprintf("rho %.2f (95%% CI %.2f-%.2f); FDR %.3f; n=%d",
          main_sc$rho, main_sc$ci_low, main_sc$ci_high,
          main_sc$fdr_bh_four, main_sc$n),
  "Fibroblast-lineage FAP (logCPM)", "Fibroblast matrix4 (z mean)",
  palette["teal"]
)

paired_long <- rbind(
  data.frame(patient = paired$patient, feature = "FAP",
             tissue = "Normal", value = paired$normal_FAP_logCPM),
  data.frame(patient = paired$patient, feature = "FAP",
             tissue = "Tumour", value = paired$tumour_FAP_logCPM),
  data.frame(patient = paired$patient, feature = "matrix4",
             tissue = "Normal", value = paired$normal_matrix4_mean_logCPM),
  data.frame(patient = paired$patient, feature = "matrix4",
             tissue = "Tumour", value = paired$tumour_matrix4_mean_logCPM)
)
paired_long$tissue <- factor(paired_long$tissue,
                             levels = c("Normal", "Tumour"))
write.csv(paired_long, file.path(source_dir, "Figure3_paired_long.csv"),
          row.names = FALSE)

paired_panel <- function(feature, colour, y_label) {
  test <- paired_tests[grepl(feature, paired_tests$feature,
                            ignore.case = TRUE), ][1, ]
  display_feature <- ifelse(feature == "matrix4", "Matrix4", feature)
  ggplot(paired_long[paired_long$feature == feature, ],
         aes(x = tissue, y = value, group = patient)) +
    geom_line(colour = palette["light_grey"], linewidth = 0.5) +
    geom_point(aes(fill = tissue), shape = 21, size = 2,
               stroke = 0.25, colour = "white") +
    scale_fill_manual(values = c(Normal = unname(palette["grey"]),
                                 Tumour = unname(colour))) +
    labs(title = paste(display_feature, "paired tumour-normal"),
         subtitle = sprintf("Median difference %.2f; FDR %.3f; %d pairs",
                            test$median_tumour_minus_normal,
                            test$fdr_bh_two, test$n_pairs),
         x = NULL, y = y_label) +
    theme_pub() + theme(legend.position = "none")
}
p3b <- paired_panel("FAP", palette["navy"], "FAP (logCPM)")
p3c <- paired_panel("matrix4", palette["teal"], "matrix4 (mean logCPM)")

boundary <- sc_cor[sc_cor$comparison !=
                     "Fibroblast-lineage FAP vs fibroblast-lineage matrix4", ]
boundary$estimate <- boundary$rho
boundary$label <- factor(
  c("FAP vs epithelial receptor2", "matrix4 vs epithelial receptor2",
    "FAP vs fibroblast receptor2"),
  levels = rev(c("FAP vs epithelial receptor2",
                 "matrix4 vs epithelial receptor2",
                 "FAP vs fibroblast receptor2"))
)
p3d <- forest_panel(
  boundary, "SDC4/CD44 comparator correlations",
  "Patient-level rho with bootstrap 95% CI; all BH FDR >= 0.248",
  "Spearman rho (95% CI)", x_limits = c(-1, 1)
)
fig3_design <- "
AAABBB
AAACCC
DDDDDD
"
fig3 <- p3a + p3b + p3c + p3d +
  plot_layout(design = fig3_design, heights = c(1, 1, 0.85)) +
  plot_annotation(tag_levels = "a")

# Figure 4: exploratory cross-compartment pathway prioritisation.
fgsea_unadj <- read_csv(file.path(
  extended_dir, "GSE132465_FAP_continuous_hallmark_fgsea.csv"
))
fgsea_adj <- read_csv(file.path(
  extended_dir, "GSE132465_FAP_cellcount_adjusted_hallmark_fgsea.csv"
))
partial <- read_csv(file.path(
  extended_dir, "GSE132465_fibFAP_epithelial_pathway_partial_correlations.csv"
))
gene_summary <- read_csv(file.path(
  extended_dir, "GSE132465_gene_model_audit_summary.csv"
))
candidate_summary <- read_csv(file.path(
  extended_dir, "GSE132465_candidate_screen_audit_summary.csv"
))

clean_pathway <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("Il6", "IL6", x, fixed = TRUE)
  x <- gsub("Jak", "JAK", x, fixed = TRUE)
  x <- gsub("Stat3", "STAT3", x, fixed = TRUE)
  x <- gsub("Kras", "KRAS", x, fixed = TRUE)
  x <- gsub("Pi3k", "PI3K", x, fixed = TRUE)
  x <- gsub("Akt", "AKT", x, fixed = TRUE)
  x <- gsub("Mtor", "mTOR", x, fixed = TRUE)
  x <- gsub("Tgf", "TGF", x, fixed = TRUE)
  x <- gsub("Tnfa", "TNFA", x, fixed = TRUE)
  x <- gsub("Nfkb", "NFKB", x, fixed = TRUE)
  x <- gsub("Wnt", "WNT", x, fixed = TRUE)
  x
}

fgsea_plot <- fgsea_adj[fgsea_adj$compartment == "epithelial", ]
fgsea_plot$omission <- ifelse(fgsea_plot$score_genes_omitted,
                              "Score genes omitted", "Full set")
fgsea_plot$pathway_label <- clean_pathway(fgsea_plot$pathway)
fgsea_plot$pathway_label <- factor(
  fgsea_plot$pathway_label,
  levels = rev(unique(clean_pathway(fgsea_plot$pathway)))
)
write.csv(fgsea_plot, file.path(source_dir, "Figure4_adjusted_fgsea.csv"),
          row.names = FALSE)
write.csv(partial, file.path(source_dir, "Figure4_partial_pathway_scores.csv"),
          row.names = FALSE)
write.csv(gene_summary, file.path(source_dir, "Figure4_gene_model_counts.csv"),
          row.names = FALSE)
write.csv(candidate_summary,
          file.path(source_dir, "Figure4_candidate_screen_counts.csv"),
          row.names = FALSE)

p4a <- ggplot(fgsea_plot, aes(x = NES, y = pathway_label,
                              shape = omission, fill = padj < 0.05)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["light_grey"]) +
  geom_point(size = 2.6, stroke = 0.5, colour = palette["dark"]) +
  scale_shape_manual(values = c("Full set" = 21,
                                "Score genes omitted" = 24)) +
  scale_fill_manual(values = c(`TRUE` = unname(palette["teal"]),
                               `FALSE` = "white"),
                    guide = "none") +
  labs(title = "Cell-count-adjusted epithelial enrichment",
       subtitle = "limma-ranked genes; filled symbols indicate BH FDR < 0.05",
       x = "Normalised enrichment score", y = NULL,
       shape = NULL, fill = "BH FDR < 0.05") +
  theme_pub() + theme(legend.position = "bottom")

emt_compare <- rbind(
  transform(fgsea_unadj[
    fgsea_unadj$compartment == "epithelial" &
      fgsea_unadj$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ],
    adjustment = "Unadjusted"),
  transform(fgsea_adj[
    fgsea_adj$compartment == "epithelial" &
      fgsea_adj$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ],
    adjustment = "Cell-count adjusted")
)
emt_compare$omission <- ifelse(emt_compare$score_genes_omitted,
                               "Genes omitted", "Full set")
write.csv(emt_compare, file.path(source_dir, "Figure4_EMT_sensitivity.csv"),
          row.names = FALSE)
p4b <- ggplot(emt_compare, aes(x = adjustment, y = NES, fill = omission)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62,
           colour = "white", linewidth = 0.2) +
  geom_text(aes(label = paste0("FDR ",
                               format(padj, scientific = TRUE, digits = 1))),
            position = position_dodge(width = 0.72), vjust = -0.5,
            size = 2.1, family = "Arial") +
  scale_fill_manual(values = c("Full set" = unname(palette["navy"]),
                               "Genes omitted" = unname(palette["teal"]))) +
  coord_cartesian(ylim = c(0, max(emt_compare$NES) * 1.18), clip = "off") +
  labs(title = "EMT ranking sensitivity",
       subtitle = "Positive ranking persists after both controls",
       x = NULL, y = "Normalised enrichment score", fill = NULL) +
  theme_pub() + theme(legend.position = "bottom")

partial_pick <- partial[
  (partial$method == "Epithelial Hallmark mean-z" &
     partial$pathway %in% c(
       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
       "HALLMARK_TGF_BETA_SIGNALING",
       "HALLMARK_WNT_BETA_CATENIN_SIGNALING")) |
    (partial$method == "PROGENy top500" &
       partial$pathway %in% c("Hypoxia", "TGFb", "WNT")), ]
partial_pick$estimate <- partial_pick$partial_rho
partial_pick$ci_low <- partial_pick$partial_ci_low
partial_pick$ci_high <- partial_pick$partial_ci_high
partial_pick$label <- ifelse(
  partial_pick$method == "PROGENy top500",
  paste0("PROGENy ", partial_pick$pathway),
  clean_pathway(partial_pick$pathway)
)
partial_pick$label <- factor(partial_pick$label,
                             levels = rev(partial_pick$label))
p4c <- forest_panel(
  partial_pick, "Patient-level scalar pathway scores",
  "Cell-count-adjusted partial rho; all BH FDR > 0.4",
  "Partial Spearman rho (95% CI)", x_limits = c(-1, 1)
)

screen_counts <- rbind(
  data.frame(stage = "Epi genes\nunadj.",
             tested = gene_summary$genes_tested[
               gene_summary$adjustment == "unadjusted" &
                 gene_summary$compartment == "epithelial"],
             retained = gene_summary$fdr_lt_0_05[
               gene_summary$adjustment == "unadjusted" &
                 gene_summary$compartment == "epithelial"]),
  data.frame(stage = "Epi genes\nadj.",
             tested = gene_summary$genes_tested[
               gene_summary$adjustment == "compartment_specific_cell_counts" &
                 gene_summary$compartment == "epithelial"],
             retained = gene_summary$fdr_lt_0_05[
               gene_summary$adjustment == "compartment_specific_cell_counts" &
                 gene_summary$compartment == "epithelial"]),
  data.frame(stage = "Candidates\nunadj.",
             tested = NA,
             retained = sum(candidate_summary$candidate_FDR)),
  data.frame(stage = "Candidates\nadj.",
             tested = NA,
             retained = sum(candidate_summary$partial_candidate_FDR))
)
screen_counts$stage <- factor(screen_counts$stage, levels = screen_counts$stage)
write.csv(screen_counts, file.path(source_dir, "Figure4_screen_attrition.csv"),
          row.names = FALSE)
p4d <- ggplot(screen_counts, aes(x = stage, y = retained)) +
  geom_col(width = 0.62, fill = palette["gold"], colour = "white") +
  geom_text(aes(label = ifelse(is.na(tested),
                               paste0(retained, " retained"),
                               paste0(retained, "/", tested))),
            vjust = -0.45, size = 2.25, family = "Arial") +
  coord_cartesian(ylim = c(0, max(screen_counts$retained) + 1.2), clip = "off") +
  labs(title = "Single-gene and candidate-screen attrition",
       subtitle = "Adjusted epithelial genes: 0; partial candidates: 0",
       x = NULL, y = "Results with FDR < 0.05") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.margin = margin(5, 8, 5, 5))
fig4 <- ((p4a | p4b) / (p4c | p4d)) +
  plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "a")

# Supplementary Figure S1: CMS sensitivity and interaction tests.
cms_cor <- read_csv(file.path(extended_dir,
                              "TCGA_CMS_stratified_correlations.csv"))
cms_interaction <- read_csv(file.path(extended_dir,
                                      "TCGA_CMS_global_interaction.csv"))
cms_cor$version <- ifelse(grepl("overlap-omitted", cms_cor$classifier),
                          "Overlap omitted", "Full classifier")
cms_cor$method <- ifelse(grepl("^RF", cms_cor$classifier), "RF", "SSP")
cms_cor$estimate <- cms_cor$rho
cms_cor$label <- factor(paste(cms_cor$method, cms_cor$CMS, sep = ": "),
                        levels = rev(unique(paste(cms_cor$method, cms_cor$CMS,
                                                  sep = ": "))))
write.csv(cms_cor, file.path(source_dir, "Supplementary_Figure_S1_CMS.csv"),
          row.names = FALSE)
write.csv(cms_interaction,
          file.path(source_dir, "Supplementary_Figure_S1_interactions.csv"),
          row.names = FALSE)
p_s1a <- forest_panel(
  cms_cor, "Within-CMS FAP13-matrix4 covariation",
  "RF and SSP high-confidence labels",
  "Spearman rho (95% CI)", colour_by = "version", x_limits = c(0.5, 1)
)
cms_interaction$version <- ifelse(grepl("overlap-omitted",
                                        cms_interaction$classifier),
                                  "Overlap omitted", "Full classifier")
cms_interaction$method <- ifelse(grepl("^RF", cms_interaction$classifier),
                                 "RF", "SSP")
p_s1b <- ggplot(cms_interaction,
                 aes(x = method, y = fdr_bh_four, fill = version)) +
  geom_hline(yintercept = 0.05, linetype = 2, linewidth = 0.45,
             colour = palette["red"]) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_text(aes(label = sprintf("%.3f", fdr_bh_four)),
            position = position_dodge(width = 0.72), vjust = -0.4,
            size = 2.3, family = "Arial") +
  scale_fill_manual(values = c("Full classifier" = unname(palette["navy"]),
                               "Overlap omitted" = unname(palette["gold"]))) +
  coord_cartesian(ylim = c(0, 0.62), clip = "off") +
  labs(title = "Global CMS interaction tests",
       subtitle = "20,000 blocked permutations; BH across 4 tests",
       x = "CMS classifier", y = "BH-adjusted permutation P", fill = NULL) +
  theme_pub() + theme(legend.position = "bottom")
fig_s1 <- (p_s1a | p_s1b) + plot_layout(widths = c(1.25, 0.75)) +
  plot_annotation(tag_levels = "a")

# Supplementary Figure S2: negative clinical boundaries.
t_trend <- read_csv(file.path(runtime_dir, "TCGA_full_primary_T_trends.csv"))
nodal <- read_csv(file.path(runtime_dir, "TCGA_full_primary_nodal_proxy.csv"))
survival <- read_csv(file.path(runtime_dir, "GSE39582_tumor_only_KM_final.csv"))
t_trend$estimate <- t_trend$rho
t_trend$label <- factor(t_trend$variable,
                        levels = rev(c("FAP", "FAP13", "matrix4", "receptor2")))
p_s2a <- forest_panel(
  t_trend, "Ordinal pathological T-category trends",
  "TCGA-COAD/READ; n=374",
  "Spearman rho (95% CI)", x_limits = c(-0.1, 0.35)
)
nodal_pick <- nodal[nodal$analysis != "unadjusted median split, Fisher exact", ]
nodal_pick$estimate <- log(nodal_pick$odds_ratio)
nodal_pick$ci_low <- log(nodal_pick$ci_low)
nodal_pick$ci_high <- log(nodal_pick$ci_high)
nodal_pick$label <- factor(
  c("Median-split proxy", "Continuous proxy"),
  levels = rev(c("Median-split proxy", "Continuous proxy"))
)
p_s2b <- forest_panel(
  nodal_pick, "T-adjusted nodal models",
  "TCGA-COAD/READ; both P > 0.05",
  "Log odds ratio (95% CI)", x_limits = c(-0.2, 1)
)
survival$estimate <- log(survival$HR_high_vs_low)
survival$ci_low <- log(survival$ci_low)
survival$ci_high <- log(survival$ci_high)
survival$label <- factor(survival$endpoint,
                         levels = rev(survival$endpoint))
p_s2c <- forest_panel(
  survival, "GSE39582 survival analyses",
  "Median split; OS n=556; RFS n=519",
  "Log hazard ratio (95% CI)", x_limits = c(-0.4, 0.8)
)
fig_s2 <- (p_s2a | p_s2b | p_s2c) + plot_annotation(tag_levels = "a")

export_manifest <- rbind(
  save_pub(fig1, file.path(figure_dir, "Figure1_cross_platform_covariation"),
           183, 145),
  save_pub(fig2, file.path(figure_dir, "Figure2_robustness_nulls"),
           183, 125),
  save_pub(fig3, file.path(figure_dir, "Figure3_patient_pseudobulk"),
           183, 135),
  save_pub(fig4, file.path(figure_dir, "Figure4_pathway_prioritisation"),
           183, 150),
  save_pub(fig_s1, file.path(figure_dir, "Supplementary_Figure_S1_CMS"),
           183, 115),
  save_pub(fig_s2,
           file.path(figure_dir, "Supplementary_Figure_S2_clinical_boundaries"),
           183, 85)
)
export_manifest_public <- export_manifest
export_manifest_public$file <- substring(export_manifest$file, nchar(root) + 2)
write.csv(export_manifest_public, file.path(qa_dir, "figure_export_manifest.csv"),
          row.names = FALSE)

# Publication tables. These are journal-neutral source tables for three-line
# Word formatting by the manuscript builder.
table1 <- data.frame(
  dataset = c("TCGA-COAD/READ", "GSE39582", "CPTAC-COAD", "GSE132465"),
  modality = c("Bulk RNA-seq", "Microarray", "Proteomics",
               "Single-cell RNA-seq pseudobulk"),
  analysis_unit = c("Primary tumour", "Tumour", "Primary tumour",
                    "Patient"),
  primary_n = c(380, 566, 97, 15),
  principal_role = c("Discovery and robustness", "Independent transcriptomic replication",
                     "Proteomic replication", "Cell-lineage patient-level support"),
  notes = c("COAD n=285; READ n=94; one project label missing",
            "Tumour-only analytical set", "FAP-free matrix3 protein score",
            "Seven paired tumour-normal patients")
)
write.csv(table1, file.path(table_dir, "Table1_dataset_overview.csv"),
          row.names = FALSE)

table2 <- fig1_stats
table2$effect_95CI <- sprintf("%.3f (%.3f to %.3f)", table2$rho,
                              table2$ci_low, table2$ci_high)
table2$BH_FDR <- fmt_p(table2$fdr_bh_within_platform_family)
table2 <- table2[c("cohort", "comparison", "n", "effect_95CI", "BH_FDR")]
write.csv(table2, file.path(table_dir, "Table2_primary_correlations.csv"),
          row.names = FALSE)

table3 <- rbind(
  data.frame(analysis = robustness$cohort, n = robustness$n,
             estimate = sprintf("rho %.3f (%.3f to %.3f)",
                                robustness$rho, robustness$ci_low,
                                robustness$ci_high),
             inference = "Sensitivity estimate"),
  data.frame(analysis = paste("Empirical null:", null_summary$null_type),
             n = null_summary$random_draws,
             estimate = sprintf("observed rho %.3f; null 97.5th percentile %.3f",
                                null_summary$observed_rho,
                                null_summary$null_q975),
             inference = paste("two-sided empirical P",
                               fmt_p(null_summary$empirical_p_two_sided)))
)
write.csv(table3, file.path(table_dir, "Table3_robustness_analyses.csv"),
          row.names = FALSE)

emt_adj <- fgsea_plot[fgsea_plot$pathway ==
                        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", ]
table4 <- data.frame(
  analysis = c("Adjusted epithelial EMT fgsea, full set",
               "Adjusted epithelial EMT fgsea, score genes omitted",
               "Patient-level epithelial EMT mean-z, count-adjusted",
               "Adjusted epithelial single-gene screen",
               "Partial ligand/receptor candidate screen"),
  n = c(15, 15, 15, 15, 15),
  result = c(
    sprintf("NES %.3f; BH FDR %s",
            emt_adj$NES[!emt_adj$score_genes_omitted],
            fmt_p(emt_adj$padj[!emt_adj$score_genes_omitted])),
    sprintf("NES %.3f; BH FDR %s",
            emt_adj$NES[emt_adj$score_genes_omitted],
            fmt_p(emt_adj$padj[emt_adj$score_genes_omitted])),
    sprintf("partial rho %.3f; BH FDR %.3f",
            partial$partial_rho[partial$method == "Epithelial Hallmark mean-z" &
                                  partial$pathway ==
                                    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"],
            partial$fdr_bh_within_method[
              partial$method == "Epithelial Hallmark mean-z" &
                partial$pathway ==
                  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]),
    paste0(sum(gene_summary$fdr_lt_0_05[
      gene_summary$adjustment == "log_fibroblast_cell_count" &
        gene_summary$compartment == "epithelial"]),
      " genome-wide FDR-significant genes"),
    paste0(sum(candidate_summary$partial_candidate_FDR),
           " candidate-family FDR-significant genes")
  ),
  interpretation = c("Exploratory ranked-gene enrichment",
                     "Not explained by score-gene overlap",
                     "No scalar-score replication after multiplicity correction",
                     "No stable downstream single-gene target",
                     "No supported ligand-receptor axis")
)
write.csv(table4, file.path(table_dir, "Table4_exploratory_pathways.csv"),
          row.names = FALSE)

write.csv(cms_cor, file.path(table_dir, "Supplementary_Table_S1_CMS.csv"),
          row.names = FALSE)
write.csv(cms_interaction,
          file.path(table_dir, "Supplementary_Table_S2_CMS_interactions.csv"),
          row.names = FALSE)
write.csv(fgsea_adj,
          file.path(table_dir, "Supplementary_Table_S3_adjusted_fgsea.csv"),
          row.names = FALSE)
write.csv(partial,
          file.path(table_dir, "Supplementary_Table_S4_pathway_scores.csv"),
          row.names = FALSE)
write.csv(read_csv(file.path(
  extended_dir, "GSE132465_candidate_ligand_receptor_screen.csv")),
  file.path(table_dir, "Supplementary_Table_S5_ligand_receptor_screen.csv"),
  row.names = FALSE)

# R-only raster QA: dimensions, dynamic range and non-white pixel fraction.
png_files <- export_manifest$file[export_manifest$format == "PNG"]
if (requireNamespace("png", quietly = TRUE)) {
  png_qa <- do.call(rbind, lapply(png_files, function(path) {
    image <- png::readPNG(path)
    rgb <- image[, , seq_len(min(3, dim(image)[3])), drop = FALSE]
    data.frame(
      file = path,
      width_px = dim(image)[2],
      height_px = dim(image)[1],
      channel_sd = stats::sd(as.numeric(rgb)),
      nonwhite_fraction = mean(rgb < 0.985),
      nonblank = stats::sd(as.numeric(rgb)) > 0.005 &&
        mean(rgb < 0.985) > 0.005
    )
  }))
} else {
  png_qa <- data.frame(file = png_files, width_px = NA, height_px = NA,
                       channel_sd = NA, nonwhite_fraction = NA,
                       nonblank = file.info(png_files)$size > 10000)
}
png_qa_public <- png_qa
png_qa_public$file <- substring(png_qa$file, nchar(root) + 2)
write.csv(png_qa_public, file.path(qa_dir, "figure_raster_QA.csv"),
          row.names = FALSE)
stopifnot(all(export_manifest$bytes > 1000), all(png_qa$nonblank))

all_outputs <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
all_outputs <- all_outputs[file.info(all_outputs)$isdir %in% FALSE]
sha <- tools::sha256sum(all_outputs)
checksums <- data.frame(
  file = substring(normalizePath(names(sha), winslash = "/"),
                   nchar(normalizePath(out_dir, winslash = "/")) + 2),
  sha256 = unname(sha)
)
write.csv(checksums, file.path(qa_dir, "SHA256SUMS.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(qa_dir, "publication_figure_sessionInfo.txt"))

cat("Publication figures and tables complete.\n")
cat("Figures:", length(png_files), "\n")
cat("All raster QA checks passed:", all(png_qa$nonblank), "\n")
