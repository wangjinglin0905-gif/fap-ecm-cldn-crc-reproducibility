options(stringsAsFactors = FALSE, scipen = 999)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
user_library <- Sys.getenv("FAP_R_LIBRARY", unset = "")
if (nzchar(user_library)) .libPaths(c(user_library, .libPaths()))

required_packages <- c("ggplot2", "patchwork", "dplyr", "tidyr", "ggrepel", "ragg", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(ragg)
  library(scales)
})

result_root <- file.path(project_root, "work", "reproducibility", "results")
figure_root <- file.path(project_root, "outputs", "figures")
main_dir <- file.path(figure_root, "main")
supp_dir <- file.path(figure_root, "supplementary")
source_dir <- file.path(figure_root, "source_data")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

palette <- c(
  navy = "#264653", teal = "#2A9D8F", gold = "#E9C46A", orange = "#F4A261",
  red = "#D95D39", blue = "#457B9D", light_blue = "#A8DADC", purple = "#7B6D8D",
  grey = "#8A8A8A", light_grey = "#D9D9D9", dark = "#222222", white = "#FFFFFF"
)

theme_manuscript <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = palette["dark"]),
      axis.ticks = element_line(linewidth = 0.35, colour = palette["dark"]),
      axis.title = element_text(size = base_size, colour = palette["dark"]),
      axis.text = element_text(size = base_size - 0.5, colour = palette["dark"]),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "#F2F4F5", colour = NA),
      plot.title = element_text(size = base_size + 0.7, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#555555"),
      plot.tag = element_text(size = 9, face = "bold"),
      panel.grid = element_blank(),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(3.8, "mm"),
      plot.margin = margin(4, 5, 4, 5)
    )
}

theme_set(theme_manuscript())

save_figure <- function(plot, directory, stem, width_mm = 183, height_mm = 125) {
  png_path <- file.path(directory, paste0(stem, ".png"))
  tiff_path <- file.path(directory, paste0(stem, ".tiff"))
  ragg::agg_png(png_path, width = width_mm, height = height_mm, units = "mm", res = 300,
                background = "white", scaling = 1)
  print(plot)
  dev.off()
  ragg::agg_tiff(tiff_path, width = width_mm, height = height_mm, units = "mm", res = 600,
                 compression = "lzw", background = "white", scaling = 1)
  print(plot)
  dev.off()
  invisible(c(png_path, tiff_path))
}

read_result <- function(...) read.csv(file.path(result_root, ...), check.names = FALSE)

copy_source <- function(source_path, target_name = basename(source_path)) {
  file.copy(source_path, file.path(source_dir, target_name), overwrite = TRUE)
}

fisher_ci <- function(rho, n, level = 0.95) {
  z <- atanh(pmax(pmin(rho, 0.999999), -0.999999))
  se <- 1 / sqrt(n - 3)
  critical <- qnorm(1 - (1 - level) / 2)
  cbind(low = tanh(z - critical * se), high = tanh(z + critical * se))
}

format_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p)))
}

# Figure 1: bulk and protein evidence
tn <- read_result("L0_TCGA", "tumor_normal_gene_comparisons.csv")
scores <- read_result("L0_TCGA", "tcga_recomputed_scores.csv")
cms <- read_result("L0_TCGA", "cms_classifier_results.csv")
cptac <- read_result("CPTAC_protein", "cptac_prespecified_correlations.csv")
bulk_cldn <- read_result(
  "L0_TCGA_claudin_sensitivity",
  "reduced_fap_caf_individual_claudin_correlations.csv"
) %>%
  filter(scope %in% c("all_profiles", "primary_tumours")) %>%
  mutate(
    scope_label = recode(
      scope,
      all_profiles = "All profiles",
      primary_tumours = "Primary tumours"
    ),
    scope_label = factor(
      scope_label,
      levels = c("All profiles", "Primary tumours")
    ),
    outcome = factor(outcome, levels = rev(c("CLDN1", "CLDN2", "CLDN4", "CLDN_core")),
                     labels = rev(c("CLDN1", "CLDN2", "CLDN4", "CLDN core")))
  )

p1a <- ggplot(tn, aes(x = reorder(gene, median_difference), y = median_difference,
                      fill = median_difference > 0)) +
  geom_hline(yintercept = 0, linewidth = 0.35, colour = "#555555") +
  geom_col(width = 0.68) +
  geom_text(aes(y = ifelse(median_difference > 0, median_difference, 0.18),
                label = sprintf("%.2f", median_difference)),
            hjust = ifelse(tn$median_difference > 0, -0.1, 0),
            size = 2.2, family = "Arial") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = c(`TRUE` = unname(palette["teal"]), `FALSE` = unname(palette["orange"])), guide = "none") +
  expand_limits(y = c(min(tn$median_difference) - 0.6, max(tn$median_difference) + 1.0)) +
  labs(title = "Tumour–normal transcript differences", x = NULL, y = "Median difference (log2 expression)")

p1b <- ggplot(scores, aes(FAP_CAF, FAP_CAF_de_ligand)) +
  geom_point(size = 0.65, alpha = 0.45, colour = palette["blue"]) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.65, colour = palette["dark"]) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3,
           label = "Spearman rho = 0.987\nn = 434", size = 2.25, family = "Arial") +
  labs(title = "FAP–CAF score robustness", x = "Original FAP–CAF score", y = "Reduced FAP–CAF score")

p1c <- ggplot(bulk_cldn, aes(spearman_rho, outcome, colour = scope_label, shape = scope_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = palette["grey"]) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    width = 0.16,
    orientation = "y",
    position = position_dodge(width = 0.38),
    linewidth = 0.45
  ) +
  geom_point(position = position_dodge(width = 0.38), size = 1.7) +
  scale_colour_manual(
    values = c("All profiles" = unname(palette["teal"]),
               "Primary tumours" = unname(palette["purple"])),
    name = NULL
  ) +
  scale_shape_manual(values = c("All profiles" = 16, "Primary tumours" = 17), name = NULL) +
  coord_cartesian(xlim = c(-0.34, 0.34)) +
  labs(
    title = "Member-specific bulk associations",
    subtitle = "All profiles n = 434; primary tumours n = 380",
    x = "Spearman rho (95% CI)",
    y = NULL
  ) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.6),
    legend.box.spacing = grid::unit(0, "mm")
  )

cms_ssp <- cms %>% filter(SSP_predicted %in% paste0("CMS", 1:4)) %>%
  mutate(SSP_predicted = factor(SSP_predicted, levels = paste0("CMS", 1:4)))
p1d <- ggplot(cms_ssp, aes(SSP_predicted, FAP_CAF, fill = SSP_predicted)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.16, size = 0.45, alpha = 0.35, colour = palette["dark"]) +
  scale_fill_manual(values = c("#5B8FF9", "#61DDAA", "#F6BD16", "#E8684A"), guide = "none") +
  labs(title = "FAP–CAF enrichment in CMS4", subtitle = "SSP labels; Kruskal–Wallis P = 1.51 × 10⁻²⁸",
       x = NULL, y = "FAP–CAF score")

cptac <- cptac %>% mutate(ci = fisher_ci(rho, n)[, 1], ci_high = fisher_ci(rho, n)[, 2],
                          class = ifelse(grepl("CLDN", y), "CLDN endpoint", "ECM endpoint"),
                          x_label = recode(x, ECM_protein_score = "ECM protein score"),
                          label = paste(x_label, "vs", y))
p1e <- ggplot(filter(cptac, class == "ECM endpoint"), aes(rho, reorder(label, rho))) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = ci, xmax = ci_high), width = 0.16, orientation = "y",
                linewidth = 0.5, colour = palette["blue"]) +
  geom_point(size = 2.2, colour = palette["blue"]) +
  xlim(-0.1, 1) +
  labs(title = "CPTAC: FAP tracks ECM proteins", x = "Spearman rho (95% CI)", y = NULL)

p1f <- ggplot(filter(cptac, class == "CLDN endpoint"), aes(rho, reorder(label, rho))) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = ci, xmax = ci_high), width = 0.16, orientation = "y",
                linewidth = 0.5, colour = palette["purple"]) +
  geom_point(size = 2.2, colour = palette["purple"]) +
  xlim(-0.7, 0.5) +
  labs(title = "No positive CPTAC link", x = "Spearman rho (95% CI)", y = NULL)

fig1 <- (p1a | p1b | p1c) / (p1d | p1e | p1f) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(fig1, main_dir, "Figure_01", height_mm = 132)

# Figure 2: single-cell compartmentalisation and independent validation
umap_path <- file.path(result_root, "L1_Harmony", "harmony_umap_embeddings.csv.gz")
umap <- read.csv(gzfile(umap_path))
set.seed(20260714)
umap_plot <- umap %>% group_by(broad_group) %>% group_modify(~slice_sample(.x, n = min(nrow(.x), 7000))) %>% ungroup()
umap_colors <- c(B = "#4C78A8", Epithelial = "#E45756", Myeloid = "#F2CF5B", Stroma = "#72B7B2", T = "#B279A2")
p2a <- ggplot(umap_plot, aes(umap_1, umap_2, colour = broad_group)) +
  geom_point(size = 0.12, alpha = 0.55) +
  scale_colour_manual(values = umap_colors) +
  coord_equal() +
  labs(title = "Integrated single-cell atlas", subtitle = "117,792 cells; display downsampled within cell groups",
       x = "Harmony UMAP 1", y = "Harmony UMAP 2", colour = "Cell group") +
  theme(legend.position = "bottom")

marker_path <- file.path(result_root, "L1_single_cell", "gse132465_marker_expression_by_subtype.csv")
markers <- read.csv(marker_path)
selected_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2", "Stromal 3", "CMS1", "CMS2", "CMS3", "CMS4")
selected_genes <- c("FAP", "COL1A1", "COL1A2", "FN1", "CLDN1", "CLDN2", "CLDN4", "EPCAM")
dot <- markers %>% filter(class == "Tumor", cell_subtype %in% selected_subtypes, gene %in% selected_genes) %>%
  mutate(cell_subtype = factor(cell_subtype, levels = rev(selected_subtypes)), gene = factor(gene, levels = selected_genes))
p2b <- ggplot(dot, aes(gene, cell_subtype, size = positive_percent, colour = mean_log_normalized)) +
  geom_point() +
  scale_size(range = c(0.4, 4.2), breaks = c(25, 50, 75, 100)) +
  scale_colour_gradientn(colours = c("#F1F3F5", palette["gold"], palette["red"]), name = "Mean log\nexpression") +
  labs(title = "FAP/ECM and CLDN occupy different compartments", x = NULL, y = NULL, size = "% positive") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

paired <- read_result("L1_TISCH_GSE166555", "paired_tumor_normal_tests.csv") %>%
  mutate(analysis_group_label = recode(analysis_group, CAF_like = "CAF-like"),
         label = paste(analysis_group_label, gene, sep = ": "),
         label = factor(label, levels = rev(paste(analysis_group_label, gene, sep = ": "))))
p2c <- ggplot(paired, aes(median_tumor_minus_normal, label, fill = median_tumor_minus_normal > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["grey"]) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0("FDR=", format_p(fdr_bh))),
            hjust = ifelse(paired$median_tumor_minus_normal > 0, -0.06, 1.06), size = 2.05, family = "Arial") +
  scale_fill_manual(values = c(`TRUE` = unname(palette["teal"]), `FALSE` = unname(palette["orange"])), guide = "none") +
  expand_limits(x = c(-1.8, 4.4)) +
  labs(title = "GSE166555 paired tumour–normal validation", x = "Median tumour–normal difference (logCPM)", y = NULL)

patient_scores <- read_result("L1_TISCH_GSE166555", "patient_cross_compartment_scores.csv") %>% filter(caf_cells >= 20, epithelial_cells >= 100)
p2d <- ggplot(patient_scores, aes(FAP_CAF_de_ligand, epithelial_CLDN_core)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = palette["light_grey"]) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = palette["light_grey"]) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, linetype = 2, colour = palette["grey"]) +
  geom_point(size = 2.1, colour = palette["purple"]) +
  ggrepel::geom_text_repel(aes(label = patient), size = 2.0, family = "Arial", max.overlaps = Inf) +
  labs(title = "No positive patient-level association detected",
       subtitle = "Spearman rho = -0.382; P = 0.276; FDR = 0.750; n = 10 patients",
       x = "CAF reduced FAP–CAF score", y = "Epithelial CLDN core score")

fig2 <- (p2a | p2b) / (p2c | p2d) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(fig2, main_dir, "Figure_02", height_mm = 138)

# Figure 3: spatial organisation and independent raw-data validation
habitats <- read_result("L2_spatial", "spatial_habitat_summary.csv")
moran <- read_result("L2_spatial", "spatial_morans_i.csv")
cooccur <- read_result("L2_spatial", "spatial_fap_epithelial_cooccurrence.csv")
spatial_patient <- read_result("L2_Valdeolivas_spatial", "spatial_patient_averaged_metrics.csv")
spatial_samples <- read_result("L2_Valdeolivas_spatial", "spatial_sample_metrics.csv")

habitat_colors <- c("Epithelial/other" = "#D9D9D9", "FAP+CAF-dominant" = unname(palette["teal"]),
                    "FAP+CAF/SPP1+TAM" = unname(palette["purple"]), "SPP1+TAM-dominant" = unname(palette["gold"]),
                    "TLS-like immune" = unname(palette["blue"]))
p3a <- ggplot(habitats, aes(patient, percent, fill = habitat)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = habitat_colors) +
  labs(title = "Qi spatial habitats", x = NULL, y = "Spots (%)", fill = "Habitat") +
  theme(legend.position = "bottom", legend.key.width = grid::unit(4, "mm"))

p3b <- ggplot(moran, aes(patient, moran_i)) +
  geom_col(width = 0.62, fill = palette["teal"]) +
  geom_text(aes(label = sprintf("%.3f", moran_i)), vjust = -0.4, size = 2.2, family = "Arial") +
  ylim(0, 0.75) +
  labs(title = "FAP scores are spatially autocorrelated", subtitle = "All permutation P = 0.001 (999 permutations)",
       x = NULL, y = "Moran's I")

p3c <- ggplot(cooccur, aes(patient, percent_all_spots)) +
  geom_col(width = 0.62, fill = palette["orange"]) +
  geom_text(aes(label = sprintf("%.2f%%", percent_all_spots)), vjust = -0.4, size = 2.2, family = "Arial") +
  ylim(0, 4) +
  labs(title = "Rare FAP-high/epithelial-high spots", x = NULL, y = "All spots (%)")

p3d <- ggplot(spatial_patient, aes(adjusted_bivariate_moran, reorder(patient_id, adjusted_bivariate_moran))) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2, colour = palette["grey"]) +
  geom_segment(aes(x = 0, xend = adjusted_bivariate_moran, yend = patient_id), linewidth = 0.7, colour = palette["blue"]) +
  geom_point(size = 2.2, colour = palette["blue"]) +
  labs(title = "No positive adjusted spatial coupling", subtitle = "Patient median = -0.0256; P = 0.078; FDR = 0.117",
       x = "Adjusted bivariate Moran statistic", y = NULL)

p3e_data <- spatial_patient %>% filter(!is.na(high_high_enrichment)) %>% mutate(enrichment_minus_one = high_high_enrichment - 1)
p3e <- ggplot(p3e_data, aes(enrichment_minus_one, reorder(patient_id, enrichment_minus_one))) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = 2, colour = palette["grey"]) +
  geom_segment(aes(x = 0, xend = enrichment_minus_one, yend = patient_id), linewidth = 0.7, colour = palette["purple"]) +
  geom_point(size = 2.2, colour = palette["purple"]) +
  labs(title = "High–high adjacency below expectation", subtitle = "Patient median = -0.118; P = 0.063; FDR = 0.117",
       x = "Adjacency enrichment minus one", y = NULL)

replicate_wide <- spatial_samples %>% select(patient_id, replicate, raw_bivariate_moran) %>%
  pivot_wider(names_from = replicate, values_from = raw_bivariate_moran)
p3f <- ggplot(replicate_wide, aes(Rep1, Rep2)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = palette["light_grey"]) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = palette["light_grey"]) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, linetype = 2, colour = palette["grey"]) +
  geom_point(size = 2.2, colour = palette["red"]) +
  ggrepel::geom_text_repel(aes(label = patient_id), size = 1.9, family = "Arial", max.overlaps = Inf) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.2,
           label = "Spearman rho = 0.75\nP = 0.052", size = 2.15, family = "Arial") +
  labs(title = "Replicate concordance", x = "Replicate 1", y = "Replicate 2")

fig3 <- (p3a | p3b | p3c) / (p3d | p3e | p3f) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"), legend.position = "bottom")
save_figure(fig3, main_dir, "Figure_03", height_mm = 142)

# Figure 4: ECM-receptor communication reproducibility
scrna_pairs_path <- file.path(project_root, "work", "cellchat_reanalysis", "results",
                              "cellchat_FAPhigh_to_epithelial_significant_excluding_SMC20.csv")
scrna_pairs <- read.csv(scrna_pairs_path) %>%
  filter(interaction_name %in% c("COL1A1_SDC4", "COL1A2_SDC4", "FN1_SDC4", "FN1_CD44", "COL1A1_CD44", "COL1A2_CD44"))
prevalence <- read_result("L3_spatial_CellChat", "spatial_cellchat_pair_prevalence.csv")
patient_pairs <- read_result("L3_spatial_CellChat", "spatial_cellchat_patient_pair_summary.csv")

p4a <- ggplot(scrna_pairs, aes(prob, reorder(interaction_name_2, prob), fill = receptor)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(CD44 = unname(palette["purple"]), SDC4 = unname(palette["teal"]))) +
  labs(title = "scRNA CellChat: FAP-high CAF to epithelium", subtitle = "Six prespecified ECM–receptor interactions",
       x = "Communication probability", y = NULL, fill = "Receptor")

prev_long <- prevalence %>% select(interaction_name_2, eligible_patients, both_replicates_patients, any_replicate_patients) %>%
  pivot_longer(c(both_replicates_patients, any_replicate_patients), names_to = "criterion", values_to = "patients") %>%
  mutate(criterion = recode(criterion, both_replicates_patients = "Both replicates", any_replicate_patients = "At least one replicate"))
p4b <- ggplot(prev_long, aes(patients, reorder(interaction_name_2, patients), fill = criterion)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = paste0(patients, "/", eligible_patients)), position = position_dodge(width = 0.7),
            hjust = -0.08, size = 2.05, family = "Arial") +
  scale_fill_manual(values = c("Both replicates" = unname(palette["teal"]), "At least one replicate" = unname(palette["light_blue"]))) +
  xlim(0, 7) +
  labs(title = "Patient-level spatial reproducibility", x = "Patients", y = NULL, fill = NULL) +
  theme(legend.position = "bottom")

heat <- patient_pairs %>% mutate(patient_id = factor(patient_id),
                                 interaction_name_2 = factor(interaction_name_2, levels = prevalence$interaction_name_2))
p4c <- ggplot(heat, aes(interaction_name_2, patient_id, fill = significant_replicates)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = paste0(significant_replicates, "/", analyzed_replicates)), size = 2.0, family = "Arial") +
  scale_fill_gradient(low = "#F2F2F2", high = palette["teal"], limits = c(0, 2), breaks = 0:2) +
  labs(title = "Replicate-level support by patient", x = NULL, y = "Patient", fill = "Significant\nreplicates") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

summary_nodes <- data.frame(
  x = c(0.8, 2.0, 3.2, 4.4), y = 1,
  label = c("FAP-high\nmyofibroblasts", "COL1A1/2\nand FN1", "SDC4/CD44\nin silico support", "CLDN effect\nnot established"),
  status = c("Observed", "Supported", "Reproduced", "Unresolved")
)
p4d <- ggplot(summary_nodes, aes(x, y)) +
  geom_segment(data = data.frame(x = c(1.18, 2.38, 3.58), xend = c(1.62, 2.82, 4.02), y = 1, yend = 1),
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 0.7, colour = palette["dark"], inherit.aes = FALSE) +
  geom_label(aes(label = label, fill = status), size = 1.65, family = "Arial", linewidth = 0.35,
             label.padding = grid::unit(1.4, "mm"), colour = palette["dark"]) +
  scale_fill_manual(values = c(Observed = "#DCEAF3", Supported = "#CFE9E2", Reproduced = "#BFD8D2",
                               Unresolved = "#F6E1D8"), guide = "none") +
  xlim(0.2, 5.0) + ylim(0.65, 1.35) +
  labs(title = "Evidence boundary") + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold", margin = margin(b = 7)))

fig4 <- (p4a | p4b) / (p4c | p4d) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(fig4, main_dir, "Figure_04", height_mm = 132)

# Figure 5: immune context and limits of secondary mechanism inference
mcp <- read_result("L4_immune", "mcpcounter_high_vs_low.csv") %>% arrange(median_difference)
epic <- read_result("L4_immune", "epic_high_vs_low.csv") %>%
  mutate(cell_label = gsub("_", " ", cell_type),
         cell_label = recode(cell_label, otherCells = "Other cells", NKcells = "NK cells",
                             Bcells = "B cells", CAFs = "CAFs", Endothelial = "Endothelial cells",
                             Macrophages = "Macrophages")) %>% arrange(median_difference)
rppa <- read_result("RPPA", "rppa_tgf_beta_correlations_from_raw_api.csv") %>%
  mutate(ci_low = fisher_ci(rho, n)[, 1], ci_high = fisher_ci(rho, n)[, 2],
         protein_label = recode(protein, SMAD3_RPPA = "SMAD3", SMAD4_RPPA = "SMAD4", CDH1_RPPA = "E-cadherin"))
progeny <- read_result("L3_pathway_TF", "progeny_high_vs_low_and_burden.csv") %>%
  mutate(feature_label = recode(feature, TGFb = "TGF-β"))
dorothea <- read_result("L3_pathway_TF", "dorothea_selected_emt_tfs.csv")

p5a <- ggplot(mcp, aes(median_difference, reorder(cell_type, median_difference), fill = fdr_bh < 0.05)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["grey"]) +
  geom_col(width = 0.64) +
  scale_fill_manual(values = c(`TRUE` = unname(palette["blue"]), `FALSE` = unname(palette["light_grey"])), guide = "none") +
  labs(title = "MCP-counter: FAP-high minus FAP-low", x = "Median score difference", y = NULL)

p5b <- ggplot(epic, aes(median_difference, reorder(cell_label, median_difference), fill = fdr_bh < 0.05)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["grey"]) +
  geom_col(width = 0.64) +
  scale_fill_manual(values = c(`TRUE` = unname(palette["teal"]), `FALSE` = unname(palette["light_grey"])), guide = "none") +
  labs(title = "EPIC directionality check", subtitle = "Interpret with algorithm convergence warnings",
       x = "Median estimated fraction difference", y = NULL)

p5c <- ggplot(rppa, aes(rho, reorder(protein_label, rho))) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.16, orientation = "y",
                linewidth = 0.55, colour = palette["purple"]) +
  geom_point(size = 2.2, colour = palette["purple"]) +
  labs(title = "RPPA associations", x = "Spearman rho (95% CI)", y = NULL)

p5d <- ggplot(progeny, aes(mean_difference, reorder(feature_label, mean_difference), fill = wilcoxon_fdr_bh < 0.05)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["grey"]) +
  geom_col(width = 0.62) +
  scale_fill_manual(values = c(`TRUE` = unname(palette["red"]), `FALSE` = unname(palette["light_grey"])), guide = "none") +
  labs(title = "PROGENy patient-level contrasts", subtitle = "No pathway survives BH-FDR correction",
       x = "FAP-high minus FAP-low mean activity", y = NULL)

p5e <- ggplot(dorothea, aes(mean_difference, reorder(feature, mean_difference), fill = wilcoxon_fdr_bh < 0.05)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = palette["grey"]) +
  geom_col(width = 0.62) +
  scale_fill_manual(values = c(`TRUE` = unname(palette["red"]), `FALSE` = unname(palette["light_grey"])), guide = "none") +
  labs(title = "DoRothEA EMT-related TF contrasts", subtitle = "No selected TF survives BH-FDR correction",
       x = "FAP-high minus FAP-low mean activity", y = NULL)

nichenet_summary <- read_result("L3_NicheNet", "epithelial_DE_summary.csv")
p5f <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#F4F5F6", colour = "#B8B8B8", linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.72, label = "NicheNet receiver-gene audit", fontface = "bold", size = 3.0, family = "Arial") +
  annotate("text", x = 0.5, y = 0.48,
           label = paste0(format(nichenet_summary$tested_genes, big.mark = ",", scientific = FALSE), " genes tested\n0 genes at FDR < 0.10"),
           size = 3.2, family = "Arial") +
  annotate("text", x = 0.5, y = 0.20, label = "Ligand ranking is exploratory", colour = palette["red"],
           size = 2.4, family = "Arial") +
  xlim(0, 1) + ylim(0, 1) + theme_void()

fig5 <- (p5a | p5b | p5c) / (p5d | p5e | p5f) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(fig5, main_dir, "Figure_05", height_mm = 140)

# Figure 6: evidence synthesis with explicit claim boundaries
evidence <- data.frame(
  layer = factor(c("Bulk RNA", "Single-cell", "Protein", "Spatial", "Spatial communication", "Causal perturbation"),
                 levels = rev(c("Bulk RNA", "Single-cell", "Protein", "Spatial", "Spatial communication", "Causal perturbation"))),
  FAP_ECM = c("Supported", "Supported", "Supported", "Supported", "Supported", "Not tested"),
  FAP_CLDN = c("Not supported", "Not supported", "Not supported", "Not supported", "Not applicable", "Not tested"),
  ECM_receptor = c("Indirect", "Supported", "Not tested", "Indirect", "Supported", "Not tested"),
  CLDN_effect = c("Linkage untested", "Linkage untested", "Linkage untested", "Linkage untested", "Not established", "Not tested")
) %>% pivot_longer(-layer, names_to = "claim", values_to = "status") %>%
  mutate(claim = recode(claim, FAP_ECM = "FAP–ECM association", FAP_CLDN = "Consistent positive FAP–CLDN covariation",
                        ECM_receptor = "ECM–receptor interface", CLDN_effect = "Downstream CLDN effect"),
         claim = factor(claim, levels = c("FAP–ECM association", "Consistent positive FAP–CLDN covariation",
                                          "ECM–receptor interface", "Downstream CLDN effect")))

status_colors <- c("Supported" = "#5AB4AC", "Indirect" = "#A8DDB5", "Linkage untested" = "#FDD49E",
                   "Not supported" = "#D7301F", "Not established" = "#FC8D59", "Not applicable" = "#D9D9D9",
                   "Not tested" = "#F0F0F0")
p6a <- ggplot(evidence, aes(claim, layer, fill = status)) +
  geom_tile(colour = "white", linewidth = 1.0) +
  geom_text(aes(label = status), size = 2.45, family = "Arial", lineheight = 0.9) +
  scale_fill_manual(values = status_colors) +
  labs(title = "Cross-layer evidence matrix", x = NULL, y = NULL, fill = "Evidence status") +
  theme(axis.text.x = element_text(angle = 18, hjust = 1), legend.position = "bottom")

flow <- data.frame(
  x = c(1, 2.25, 3.5, 4.75), y = 1,
  label = c("FAP-associated\nCAF state", "Collagen/FN1-rich\nECM programme", "SDC4/CD44 interface\nreproduced in silico", "CLDN remodelling\nrequires perturbation"),
  fill = c("Observed", "Supported", "Reproduced", "Hypothesis")
)
p6b <- ggplot(flow, aes(x, y)) +
  geom_segment(data = data.frame(x = c(1.35, 2.60, 3.85), xend = c(1.90, 3.15, 4.40), y = 1, yend = 1),
               aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
               arrow = grid::arrow(length = grid::unit(2.5, "mm")), linewidth = 0.8, colour = palette["dark"]) +
  geom_label(aes(label = label, fill = fill), size = 2.7, family = "Arial", linewidth = 0.4,
             label.padding = grid::unit(2.8, "mm")) +
  scale_fill_manual(values = c(Observed = "#DCEAF3", Supported = "#CFE9E2", Reproduced = "#BFD8D2",
                               Hypothesis = "#F6E1D8"), guide = "none") +
  xlim(0.5, 5.25) + ylim(0.65, 1.35) +
  labs(title = "Defensible model and experimental boundary") + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8.5, face = "bold", margin = margin(b = 8)))

fig6 <- p6a / p6b + plot_layout(heights = c(2.2, 1)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold"))
save_figure(fig6, main_dir, "Figure_06", height_mm = 124)

# Supplementary Figure S1: UALCAN CPTAC display summaries
ualcan <- read_result("UALCAN_CPTAC", "ualcan_cptac_primary_vs_normal.csv")
ualcan_long <- bind_rows(
  ualcan %>% transmute(gene, group = "Normal", ymin = normal_low, lower = normal_q1, middle = normal_median, upper = normal_q3, ymax = normal_high, n = normal_n),
  ualcan %>% transmute(gene, group = "Primary tumour", ymin = tumour_low, lower = tumour_q1, middle = tumour_median, upper = tumour_q3, ymax = tumour_high, n = tumour_n)
) %>% mutate(group = factor(group, levels = c("Normal", "Primary tumour")))
ps1 <- ggplot(ualcan_long, aes(group, middle, fill = group)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = "identity", width = 0.58) +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Normal" = "#B9C6CF", "Primary tumour" = unname(palette["teal"])), guide = "none") +
  labs(title = "UALCAN display of CPTAC colon proteomics", subtitle = "Same CPTAC source; not counted as an independent proteomic cohort",
       x = NULL, y = "UALCAN protein Z-value") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_figure(ps1, supp_dir, "Figure_S01", height_mm = 122)

# Supplementary Figure S2: GSE166555 association sensitivity
assoc <- read_result("L1_TISCH_GSE166555", "patient_level_associations.csv") %>%
  filter(outcome == "epithelial_CLDN_core") %>%
  mutate(label = recode(predictor, FAP_log_cpm = "FAP logCPM", FAP_positive_fraction = "FAP-positive fraction",
                        FAP_CAF = "FAP–CAF score", FAP_CAF_de_ligand = "Reduced FAP–CAF score"),
         threshold = paste0("≥", minimum_caf_cells, " CAF-like cells"))
ps2 <- ggplot(assoc, aes(spearman_rho, reorder(label, spearman_rho), colour = threshold)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = 2, colour = palette["grey"]) +
  geom_errorbar(aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high), width = 0.16,
                orientation = "y", position = position_dodge(width = 0.45)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.45)) +
  facet_wrap(~threshold) +
  scale_x_continuous(limits = c(-1.05, 0.55), breaks = c(-1, -0.5, 0, 0.5)) +
  scale_colour_manual(values = c("≥20 CAF-like cells" = unname(palette["purple"]), "≥5 CAF-like cells" = unname(palette["blue"])), guide = "none") +
  labs(title = "GSE166555 patient-level sensitivity analyses", subtitle = "All BH-FDR values > 0.70",
       x = "Spearman rho with epithelial CLDN core (bootstrap 95% CI)", y = NULL) +
  theme(panel.spacing.x = grid::unit(8, "mm"))
save_figure(ps2, supp_dir, "Figure_S02", height_mm = 100)

# Supplementary Figure S3: Valdeolivas section-level metrics
spatial_long <- spatial_samples %>% select(sample_id, patient_id, replicate, raw_bivariate_moran, adjusted_bivariate_moran, high_high_enrichment) %>%
  pivot_longer(c(raw_bivariate_moran, adjusted_bivariate_moran, high_high_enrichment), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, raw_bivariate_moran = "Raw bivariate Moran", adjusted_bivariate_moran = "Adjusted bivariate Moran",
                         high_high_enrichment = "High–high adjacency enrichment"))
ps3 <- ggplot(spatial_long, aes(patient_id, value, colour = replicate, group = patient_id)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, colour = palette["grey"]) +
  geom_point(size = 1.8, position = position_dodge(width = 0.38), na.rm = TRUE) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c(Rep1 = unname(palette["blue"]), Rep2 = unname(palette["orange"]))) +
  labs(title = "Valdeolivas spatial validation by section", x = "Patient", y = "Metric", colour = "Section") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
save_figure(ps3, supp_dir, "Figure_S03", height_mm = 165)

# Supplementary Figure S4: spatial CellChat probabilities
prespecified <- read_result("L3_spatial_CellChat", "spatial_cellchat_prespecified_direction.csv")
ps4 <- ggplot(prespecified, aes(interaction_name_2, patient_id, fill = probability)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_wrap(~replicate, ncol = 1) +
  scale_fill_gradientn(colours = c("#F4F5F6", palette["light_blue"], palette["teal"]), na.value = "#EEEEEE") +
  labs(title = "Spatially constrained CellChat probabilities", subtitle = "Fibroblastic stroma → pathologist-annotated tumour; 200 µm range",
       x = NULL, y = "Patient", fill = "Probability") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_figure(ps4, supp_dir, "Figure_S04", height_mm = 128)

# Supplementary Figure S5: exploratory NicheNet ligand ranking
ligands <- read_result("L3_NicheNet", "nichenet_ligand_activity_dataset_specific.csv") %>%
  arrange(desc(pearson)) %>% slice_head(n = 20)
ps5 <- ggplot(ligands, aes(pearson, reorder(test_ligand, pearson), fill = test_ligand %in% c("COL1A1", "COL1A2", "FN1", "TGFB1", "CCL2"))) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(`TRUE` = unname(palette["orange"]), `FALSE` = unname(palette["light_grey"])), guide = "none") +
  labs(title = "Exploratory NicheNet ligand activity", subtitle = "Receiver-gene differential expression yielded 0 genes at FDR < 0.10",
       x = "Pearson ligand activity", y = NULL)
save_figure(ps5, supp_dir, "Figure_S05", height_mm = 125)

# Copy source tables used by the figures
source_files <- c(
  file.path(result_root, "L0_TCGA", "tumor_normal_gene_comparisons.csv"),
  file.path(result_root, "L0_TCGA", "tcga_recomputed_scores.csv"),
  file.path(result_root, "L0_TCGA", "cms_classifier_results.csv"),
  file.path(result_root, "L0_TCGA_claudin_sensitivity", "reduced_fap_caf_individual_claudin_correlations.csv"),
  file.path(result_root, "L0_TCGA_claudin_sensitivity", "matched_primary_normal_tests.csv"),
  file.path(result_root, "CPTAC_protein", "cptac_prespecified_correlations.csv"),
  file.path(result_root, "L1_single_cell", "gse132465_marker_expression_by_subtype.csv"),
  file.path(result_root, "L1_TISCH_GSE166555", "paired_tumor_normal_tests.csv"),
  file.path(result_root, "L1_TISCH_GSE166555", "patient_cross_compartment_scores.csv"),
  file.path(result_root, "L2_spatial", "spatial_habitat_summary.csv"),
  file.path(result_root, "L2_spatial", "spatial_morans_i.csv"),
  file.path(result_root, "L2_spatial", "spatial_fap_epithelial_cooccurrence.csv"),
  file.path(result_root, "L2_Valdeolivas_spatial", "spatial_patient_averaged_metrics.csv"),
  file.path(result_root, "L2_Valdeolivas_spatial", "spatial_sample_metrics.csv"),
  scrna_pairs_path,
  file.path(result_root, "L3_spatial_CellChat", "spatial_cellchat_pair_prevalence.csv"),
  file.path(result_root, "L3_spatial_CellChat", "spatial_cellchat_patient_pair_summary.csv"),
  file.path(result_root, "L4_immune", "mcpcounter_high_vs_low.csv"),
  file.path(result_root, "L4_immune", "epic_high_vs_low.csv"),
  file.path(result_root, "RPPA", "rppa_tgf_beta_correlations_from_raw_api.csv"),
  file.path(result_root, "L3_pathway_TF", "progeny_high_vs_low_and_burden.csv"),
  file.path(result_root, "L3_pathway_TF", "dorothea_selected_emt_tfs.csv"),
  file.path(result_root, "L3_NicheNet", "epithelial_DE_summary.csv"),
  file.path(result_root, "UALCAN_CPTAC", "ualcan_cptac_primary_vs_normal.csv")
)
invisible(lapply(source_files[file.exists(source_files)], copy_source))
write.csv(evidence, file.path(source_dir, "Figure_06_evidence_matrix.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(figure_root, "figure_sessionInfo.txt"))

main_png_paths <- file.path(main_dir, paste0(sprintf("Figure_%02d", 1:6), ".png"))
supp_png_paths <- file.path(supp_dir, paste0(sprintf("Figure_S%02d", 1:5), ".png"))
main_tiff_paths <- file.path(main_dir, paste0(sprintf("Figure_%02d", 1:6), ".tiff"))
supp_tiff_paths <- file.path(supp_dir, paste0(sprintf("Figure_S%02d", 1:5), ".tiff"))
all_png_paths <- c(main_png_paths, supp_png_paths)
all_tiff_paths <- c(main_tiff_paths, supp_tiff_paths)
relative_to_project <- function(path) {
  sub(paste0(project_root, "/"), "", gsub("\\\\", "/", path), fixed = TRUE)
}

figure_manifest <- data.frame(
  figure = c(sprintf("Figure_%02d", 1:6), sprintf("Figure_S%02d", 1:5)),
  png = vapply(all_png_paths, relative_to_project, character(1)),
  tiff = vapply(all_tiff_paths, relative_to_project, character(1)),
  png_dpi = 300,
  tiff_dpi = 600
)
figure_manifest$png_exists <- file.exists(all_png_paths)
figure_manifest$tiff_exists <- file.exists(all_tiff_paths)
write.csv(figure_manifest, file.path(figure_root, "figure_manifest.csv"), row.names = FALSE)

cat("Generated", nrow(figure_manifest), "figures in", figure_root, "\n")
