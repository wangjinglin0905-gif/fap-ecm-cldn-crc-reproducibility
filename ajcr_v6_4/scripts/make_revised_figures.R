options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript make_revised_figures.R <GSE132465_result_directory> <figure_output_directory>")
}
in_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols <- c(
  fibroblast = "#2F7F8F",
  epithelial = "#D9784A",
  fap = "#9B4D96",
  matrix = "#3A7D44",
  senescence = "#C4553D",
  neutral = "#6E7781",
  pale = "#D9E2E7"
)

base_theme <- theme_classic(base_family = "Arial", base_size = 8.2) +
  theme(
    plot.title = element_text(face = "bold", size = 8.7, margin = margin(b = 3)),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.2, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

fmt_p <- function(p) {
  if (p < 0.001) format(p, scientific = TRUE, digits = 2) else sprintf("%.3f", p)
}

save_pub <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot = plot,
         width = width_in, height = height_in, units = "in", dpi = 600,
         bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), plot = plot,
         width = width_in, height = height_in, units = "in", dpi = 600,
         compression = "lzw", bg = "white", limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot = plot,
         width = width_in, height = height_in, units = "in", device = cairo_pdf,
         bg = "white", limitsize = FALSE)
}

paired_panel <- function(file, labels, title, ylab, annotation) {
  d <- read.csv(file.path(in_dir, file), check.names = FALSE)
  long <- rbind(
    data.frame(Patient = d$id, group = labels[1], value = d$a),
    data.frame(Patient = d$id, group = labels[2], value = d$b)
  )
  long$group <- factor(long$group, levels = labels)
  ggplot(long, aes(group, value, group = Patient)) +
    geom_line(colour = "#AAB2B8", linewidth = 0.35, alpha = 0.72) +
    geom_point(aes(colour = group), size = 1.65, alpha = 0.95) +
    stat_summary(aes(group = group), fun = median, geom = "crossbar",
                 width = 0.42, linewidth = 0.55, colour = "black") +
    scale_colour_manual(values = setNames(c(cols["fibroblast"], cols["epithelial"]), labels)) +
    labs(title = title, x = NULL, y = ylab) +
    annotate("label", x = 1.5, y = Inf, vjust = 1.25, label = annotation,
             size = 2.25, label.size = 0, fill = "white", colour = "#30363D") +
    guides(colour = "none") +
    base_theme
}

# Figure 1: tumour-only, patient-level compartment contrasts.
p1a <- paired_panel(
  "tumor_compartment_SenMayo_patient_values.csv",
  c("Fibroblast", "Epithelial"),
  "SenMayo score",
  "Mean log-normalized score",
  "23/23 higher in fibroblasts\nexact paired P = 2.38e-7"
)
p1b <- paired_panel(
  "tumor_compartment_SASP_patient_values.csv",
  c("Fibroblast", "Epithelial"),
  "SASP score",
  "Mean log-normalized score",
  "22/23 higher in fibroblasts\nexact paired P = 4.77e-7"
)
p1c <- paired_panel(
  "tumor_compartment_MKI67_patient_values.csv",
  c("Fibroblast", "Epithelial"),
  "MKI67 expression",
  "Mean log-normalized expression",
  "21/23 lower in fibroblasts\nexact paired P = 1.81e-4"
)
figure1 <- (p1a | p1b | p1c) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", family = "Arial", size = 11)))
save_pub(figure1, "Figure1_tumor_only_compartment", 180, 72)

# Figure 2: FAP-specific inference and patient-level separation of axes.
score_fx <- data.frame(
  feature = c("SenMayo", "SASP"),
  estimate = c(0.02382830, 0.03225654),
  low = c(-0.0051, -0.0429),
  high = c(0.0527, 0.1074),
  p = c(0.1063647, 0.4003005)
)
score_fx$feature <- factor(score_fx$feature, levels = rev(score_fx$feature))
p2a <- ggplot(score_fx, aes(estimate, feature)) +
  geom_vline(xintercept = 0, colour = "#8C959F", linetype = 2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.16, linewidth = 0.55, colour = cols["fap"]) +
  geom_point(size = 2.0, colour = cols["fap"]) +
  geom_text(aes(label = paste0("P = ", sprintf("%.3f", p))), x = 0.135, hjust = 1, size = 2.35) +
  coord_cartesian(xlim = c(-0.055, 0.14), clip = "off") +
  labs(title = "Adjusted FAP effect on gene-set scores",
       subtitle = "log UMI + subtype + MKI67 + patient random intercept",
       x = "Adjusted beta (95% CI)", y = NULL) +
  base_theme + theme(plot.subtitle = element_text(size = 6.8, margin = margin(b = 3)))

marker_fx <- data.frame(
  feature = c("CDKN2A", "CDKN2B", "CDKN1A", "LMNB1", "MKI67"),
  OR = c(0.8902096, 0.9812225, 0.9618322, 0.8561041, 0.8363848),
  low = c(0.6200189, 0.7031367, 0.7238329, 0.4780376, 0.4177235),
  high = c(1.278143, 1.369289, 1.278087, 1.533173, 1.674647),
  p = c(0.5285691, 0.9112278, 0.7884639, 0.6012626, 0.6139833)
)
marker_fx$feature <- factor(marker_fx$feature, levels = rev(marker_fx$feature))
p2b <- ggplot(marker_fx, aes(OR, feature)) +
  geom_vline(xintercept = 1, colour = "#8C959F", linetype = 2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.16, linewidth = 0.5, colour = cols["neutral"]) +
  geom_point(size = 1.8, colour = cols["neutral"]) +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 2), limits = c(0.38, 2.05)) +
  labs(title = "Adjusted FAP effect on marker detection",
       subtitle = "log UMI + subtype + patient random intercept",
       x = "Odds ratio (95% CI; log scale)", y = NULL) +
  base_theme + theme(plot.subtitle = element_text(size = 6.8, margin = margin(b = 3)))

patient_scores <- read.csv(file.path(in_dir, "tumor_fibroblast_patient_scores.csv"), check.names = FALSE)
scatter_panel <- function(x, y, xlab, ylab, title, rho, p, colour) {
  ggplot(patient_scores, aes(.data[[x]], .data[[y]])) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = colour,
                fill = colour, alpha = 0.10, linewidth = 0.55) +
    geom_point(shape = 21, size = 2.1, stroke = 0.45, fill = "white", colour = colour) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.25,
             label = paste0("Spearman rho = ", sprintf("%.3f", rho), "\nP = ", fmt_p(p)),
             size = 2.25, label.size = 0, fill = "white") +
    labs(title = title, x = xlab, y = ylab) +
    base_theme
}
p2c <- scatter_panel("FAP_expr", "SenMayo_mean", "FAP", "SenMayo", "FAP versus SenMayo", 0.3310277, 0.1228564, cols["senescence"])
p2d <- scatter_panel("FAP_expr", "matrix4_mean", "FAP", "matrix4", "FAP versus matrix4", 0.6492095, 0.0008028, cols["matrix"])
p2e <- scatter_panel("SenMayo_mean", "matrix4_mean", "SenMayo", "matrix4", "SenMayo versus matrix4", -0.0484190, 0.8263478, cols["neutral"])
figure2 <- ((p2a | p2b) / (p2c | p2d | p2e)) +
  plot_layout(heights = c(0.92, 1.08)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", family = "Arial", size = 11)))
save_pub(figure2, "Figure2_FAP_and_matrix_separation", 180, 142)

# Figure 3: matched-signature nulls and bulk composition sensitivity.
null_draws <- read.csv(file.path(in_dir, "signature_matched_null_draws.csv"), check.names = FALSE)
null_summary <- read.csv(file.path(in_dir, "signature_matched_null_summary.csv"), check.names = FALSE)
null_hist <- function(signature, stat, title, fill) {
  d <- null_draws[null_draws$signature == signature, ]
  s <- null_summary[null_summary$signature == signature & null_summary$statistic == stat, ]
  xcol <- stat
  ggplot(d, aes(.data[[xcol]])) +
    geom_histogram(bins = 45, fill = fill, colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = s$observed, colour = "#B31B1B", linewidth = 0.8) +
    annotate("label", x = Inf, y = Inf, hjust = 1.04, vjust = 1.22,
             label = paste0("Observed = ", sprintf("%.3f", s$observed),
                            "\nNull 97.5% = ", sprintf("%.3f", s$null_q975),
                            "\nEmpirical P = ", format(as.numeric(s$empirical_p_upper), scientific = TRUE, digits = 2)),
             size = 2.2, label.size = 0, fill = "white") +
    labs(title = title, x = "Mean fibroblast - epithelial difference", y = "Matched-gene draws") +
    base_theme
}
p3a <- null_hist("SenMayo_nonoverlap", "mean_patient_difference", "SenMayo: expression-matched specificity", "#8DB9C4")
p3b <- null_hist("SASP25", "mean_patient_difference", "SASP: expression-matched specificity", "#E7A47C")

bulk <- data.frame(
  pair = rep(c("FAP13-matrix4", "SenMayo-FAP13", "SenMayo-matrix4"), each = 4),
  adjustment = rep(c("fib5", "MCP-counter", "EPIC", "MCP + EPIC"), 3),
  rho = c(0.650, 0.113, 0.288, 0.150, 0.388, 0.378, 0.575, 0.384, 0.133, -0.098, 0.082, -0.039),
  low = c(0.574, -0.015, 0.176, 0.040, 0.284, 0.270, 0.490, 0.273, 0.020, -0.221, -0.015, -0.137),
  high = c(0.715, 0.250, 0.398, 0.257, 0.485, 0.476, 0.648, 0.485, 0.248, 0.043, 0.188, 0.059)
)
bulk$adjustment <- factor(bulk$adjustment, levels = rev(c("fib5", "MCP-counter", "EPIC", "MCP + EPIC")))
bulk$pair <- factor(bulk$pair, levels = c("FAP13-matrix4", "SenMayo-FAP13", "SenMayo-matrix4"))
p3c <- ggplot(bulk, aes(rho, adjustment, colour = adjustment)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#8C959F", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.18, linewidth = 0.55) +
  geom_point(size = 2.0) +
  facet_wrap(~pair, nrow = 1) +
  scale_colour_manual(values = c("fib5" = "#7A5195", "MCP-counter" = "#2F7F8F", "EPIC" = "#D9784A", "MCP + EPIC" = "#3A7D44")) +
  coord_cartesian(xlim = c(-0.25, 0.75)) +
  labs(title = "TCGA composition-adjusted partial correlations", x = "Partial Spearman rho (95% bootstrap CI)", y = NULL) +
  guides(colour = "none") +
  base_theme + theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 7.4))

figure3 <- ((p3a | p3b) / p3c) +
  plot_layout(heights = c(1, 0.95)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", family = "Arial", size = 11)))
save_pub(figure3, "Figure3_specificity_and_composition", 180, 138)

# Figure 5: biologically cognate, tumour-only ligand-receptor screen.
cognate <- read.csv(file.path(in_dir, "tumor_only_cognate_ligand_receptor_screen.csv"), check.names = FALSE)
cognate$pair <- paste0(cognate$ligand, " → ", cognate$receptor)
cognate$immune_compartment <- factor(cognate$immune_compartment, levels = c("T_cells", "Myeloid"), labels = c("T cells", "Myeloid"))
pair_order <- unique(c("IL6 → IL6R", "IL6 → IL6ST", "CXCL8 → CXCR1", "CXCL8 → CXCR2", "CCL2 → CCR2", "TGFB1 → TGFBR1", "TGFB1 → TGFBR2"))
cognate$pair <- factor(cognate$pair, levels = rev(pair_order))
cognate$label <- sprintf("%.2f", cognate$rho)
figure5 <- ggplot(cognate, aes(immune_compartment, pair, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = label), size = 2.8, family = "Arial", colour = "black") +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C84C4C", midpoint = 0,
                       limits = c(-0.5, 0.5), oob = scales::squish, name = "Spearman\nrho") +
  labs(title = "Tumour-only patient-level screen of cognate SASP ligand-receptor pairs",
       subtitle = "FAP-detected fibroblast ligand expression versus immune-compartment receptor expression; n = 22 patients",
       x = NULL, y = NULL, caption = "No association survived BH correction (minimum q = 0.598). Correlations do not establish cell-cell communication.") +
  base_theme +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.subtitle = element_text(size = 7.2, margin = margin(b = 5)),
    plot.caption = element_text(size = 6.8, hjust = 0, colour = "#4A5560"),
    legend.position = "right"
  )
save_pub(figure5, "Figure5_cognate_SASP_immune_screen", 150, 92)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_figures.txt"))
cat("Revised figures written to", normalizePath(out_dir), "\n")
