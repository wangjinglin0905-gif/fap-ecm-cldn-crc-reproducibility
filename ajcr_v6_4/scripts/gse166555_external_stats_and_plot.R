options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript gse166555_external_stats_and_plot.R",
    "<patient_compartment_scores.csv> <analysis_ledger.json>",
    "<analysis_output_directory> <figure_output_directory>"
  ))
}
in_file <- normalizePath(args[[1]], mustWork = TRUE)
ledger_file <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
figure_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(in_file, check.names = FALSE)
ledger <- fromJSON(ledger_file)

paired_values <- function(feature) {
  a <- d[d$compartment == "Fibroblast", c("patient", feature)]
  b <- d[d$compartment == "Epithelial", c("patient", feature)]
  names(a)[2] <- "fibroblast"
  names(b)[2] <- "epithelial"
  merge(a, b, by = "patient", all = FALSE)
}

summarize_paired <- function(feature, seed_offset = 0L) {
  x <- paired_values(feature)
  delta <- x$fibroblast - x$epithelial
  test <- suppressWarnings(wilcox.test(
    x$fibroblast, x$epithelial, paired = TRUE, exact = TRUE,
    alternative = "two.sided", correct = FALSE
  ))
  set.seed(20260814L + seed_offset)
  boot <- replicate(10000L, median(sample(delta, length(delta), replace = TRUE)))
  list(
    n = nrow(x),
    concordant = sum(delta > 0),
    median_difference = unname(median(delta)),
    mean_difference = unname(mean(delta)),
    ci_low = unname(quantile(boot, 0.025, names = FALSE)),
    ci_high = unname(quantile(boot, 0.975, names = FALSE)),
    p_exact = unname(test$p.value)
  )
}

summary <- list(
  dataset = "GSE166555",
  source = "GEO official supplementary count matrices and metadata",
  eligible_patients = ledger$eligible_patients,
  SenMayo_target_genes = 119L,
  SenMayo_available_genes = length(ledger$SenMayo119_available),
  SASP_available_genes = length(ledger$SASP25_available),
  SenMayo119 = summarize_paired("SenMayo119", 1L),
  SASP25 = summarize_paired("SASP25", 2L),
  MKI67 = summarize_paired("MKI67", 3L)
)

fib <- d[d$compartment == "Fibroblast", ]
cor_summary <- list(
  FAP13_matrix4 = list(
    rho = unname(cor(fib$FAP13, fib$matrix4, method = "spearman")),
    p = unname(cor.test(fib$FAP13, fib$matrix4, method = "spearman", exact = FALSE)$p.value)
  ),
  FAP13_SenMayo = list(
    rho = unname(cor(fib$FAP13, fib$SenMayo119, method = "spearman")),
    p = unname(cor.test(fib$FAP13, fib$SenMayo119, method = "spearman", exact = FALSE)$p.value)
  ),
  SenMayo_matrix4 = list(
    rho = unname(cor(fib$SenMayo119, fib$matrix4, method = "spearman")),
    p = unname(cor.test(fib$SenMayo119, fib$matrix4, method = "spearman", exact = FALSE)$p.value)
  )
)
summary$fibroblast_patient_correlations <- cor_summary

write_json(summary, file.path(out_dir, "validation_summary.json"), pretty = TRUE, auto_unbox = TRUE, digits = 12)

cols <- c(Fibroblast = "#2F7F8F", Epithelial = "#D9784A")
base_theme <- theme_classic(base_family = "Arial", base_size = 8.2) +
  theme(
    plot.title = element_text(face = "bold", size = 8.8, margin = margin(b = 3)),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.2, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    plot.tag = element_text(face = "bold", size = 11),
    plot.tag.position = c(0, 1)
  )

paired_panel <- function(feature, title, ylab, result, direction = "higher", invert_count = FALSE) {
  x <- paired_values(feature)
  long <- rbind(
    data.frame(patient = x$patient, compartment = "Fibroblast", value = x$fibroblast),
    data.frame(patient = x$patient, compartment = "Epithelial", value = x$epithelial)
  )
  long$compartment <- factor(long$compartment, levels = c("Fibroblast", "Epithelial"))
  shown_count <- if (invert_count) result$n - result$concordant else result$concordant
  annotation <- paste0(
    shown_count, "/", result$n, " ", direction, " in fibroblasts\n",
    "exact paired P = ", format(result$p_exact, digits = 3, scientific = TRUE)
  )
  ggplot(long, aes(compartment, value, group = patient)) +
    geom_line(colour = "#AAB2B8", linewidth = 0.38, alpha = 0.75) +
    geom_point(aes(colour = compartment), size = 1.9, alpha = 0.96) +
    stat_summary(aes(group = compartment), fun = median, geom = "crossbar",
                 width = 0.42, linewidth = 0.58, colour = "black") +
    scale_colour_manual(values = cols) +
    annotate("label", x = 1.5, y = Inf, vjust = 1.25, label = annotation,
             size = 2.25, fill = "white", colour = "#30363D") +
    guides(colour = "none") +
    labs(title = title, x = NULL, y = ylab) +
    base_theme
}

p1 <- paired_panel(
  "SenMayo119", "SenMayo score (114/119 genes)",
  "Patient-compartment z-mean", summary$SenMayo119, "higher"
)
p2 <- paired_panel(
  "SASP25", "SASP score", "Patient-compartment z-mean",
  summary$SASP25, "higher"
)
p3 <- paired_panel(
  "MKI67", "MKI67 expression", "Pseudobulk log1p(CPM)",
  summary$MKI67, "lower", TRUE
)

figure <- (p1 | p2 | p3) +
  plot_annotation(
    title = "External tumour-compartment evaluation in GSE166555",
    subtitle = "Tumour-only patient pseudobulks; n = 10 patients meeting frozen cell-count thresholds",
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 11),
      plot.subtitle = element_text(family = "Arial", size = 8),
      plot.tag = element_text(family = "Arial", face = "bold", size = 11)
    )
  )

width_in <- 180 / 25.4
height_in <- 82 / 25.4
ggsave(file.path(figure_dir, "Figure4_GSE166555_external_validation.png"), figure,
       width = width_in, height = height_in, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Figure4_GSE166555_external_validation.tiff"), figure,
       width = width_in, height = height_in, dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(figure_dir, "Figure4_GSE166555_external_validation.pdf"), figure,
       width = width_in, height = height_in, device = cairo_pdf, bg = "white")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE, digits = 8), "\n")
