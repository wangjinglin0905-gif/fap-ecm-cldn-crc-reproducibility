#!/usr/bin/env Rscript
# Layout-only edit of the original quantitative grid, using the frozen source tables.
options(stringsAsFactors = FALSE, scipen = 999)
if (nzchar(Sys.getenv("FAP_R_LIBRARY"))) .libPaths(c(Sys.getenv("FAP_R_LIBRARY"), .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(grid) })
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
root <- normalizePath(args[1], mustWork = TRUE)
public_package <- dir.exists(file.path(root, "figures/main/source_data"))
out <- file.path(root, if (public_package) "figures/main" else "outputs/v7.3_2026-09-09/figures/main")
qa <- file.path(root, if (public_package) "qa/fig3_layout_v7.3" else "qa/v7.3_2026-09-09/fig3")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(qa, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(root, if (public_package) "figures/main/source_data" else "figures/preview_gate_v02_2026-09-04/main/source_data")
source_script <- file.path(root, "scripts/21_build_four_main_and_supplement_previews_v7_1.R")
expressions <- parse(source_script)
assignment_name <- function(e) if (is.call(e) && identical(e[[1]], as.name("<-")) && is.symbol(e[[2]])) as.character(e[[2]]) else ""
for (e in expressions) if (assignment_name(e) %in% c("pal", "theme_pub", "tag_plot")) eval(e)
marginal <- read.csv(file.path(source_dir, "Fig3_bulk_marginal_correlations.csv"), check.names = FALSE)
marginal$label <- factor(marginal$cohort, levels = rev(c("TCGA-COAD/READ", "GSE39582")))
bulk_main <- read.csv(file.path(source_dir, "Fig3_TCGA_composition_sensitivity.csv"), check.names = FALSE)
bulk_main$label <- factor(bulk_main$label, levels = rev(c("Unadjusted", "fib5", "MCP + EPIC full", "MCP + EPIC pair-purged", "MCP + EPIC globally disjoint")))
global_disjoint <- read.csv(file.path(source_dir, "Fig3_TCGA_global_disjoint_correlations.csv"), check.names = FALSE)
global_disjoint$label <- factor(global_disjoint$label, levels = rev(c("FAP13-matrix4", "SenMayo-FAP13", "SenMayo-matrix4")))
cptac <- read.csv(file.path(source_dir, "Fig3_CPTAC_protein_correlations.csv"), check.names = FALSE)
cptac$short <- factor(cptac$short, levels = rev(cptac$short))
for (e in expressions) if (assignment_name(e) %in% c("f3a", "f3b", "f3c", "f3d", "fig3")) eval(e)
ragg::agg_capture(width = 2250, height = 1995, res = 300)
baseline <- patchworkGrob(fig3)
dev.off()
revised <- baseline
write.csv(baseline$layout, file.path(qa, "gtable_layout_before.csv"), row.names = FALSE)
for (row in c(1, 2)) {
  table_index <- which(revised$layout$name == paste0("patchwork-table-", row))
  tab <- revised$grobs[[table_index]]
  panel_index <- which(tab$layout$name == "panel-2")
  stopifnot(length(panel_index) == 1)
  left <- tab$layout$l[panel_index]
  for (kind in c("title", "subtitle", "tag")) {
    index <- which(tab$layout$name == paste0(kind, "-2"))
    stopifnot(length(index) == 1)
    tab$layout$l[index] <- left
  }
  stopifnot(identical(tab$grobs, baseline$grobs[[table_index]]$grobs))
  revised$grobs[[table_index]] <- tab
  write.csv(tab$layout, file.path(qa, paste0("row_", row, "_layout_after.csv")), row.names = FALSE)
}
# Layout columns/rows and quantitative grobs are preserved; only the title group anchors change.
stopifnot(identical(baseline$widths, revised$widths), identical(baseline$heights, revised$heights))
write.csv(revised$layout, file.path(qa, "gtable_layout_after.csv"), row.names = FALSE)
ggsave(file.path(qa, "Fig3_reconstructed_baseline.png"), baseline, width = 7.5, height = 6.65, dpi = 300, bg = "white")
ggsave(file.path(out, if (public_package) "Fig3_bulk_composition_and_proteomics_preview.png" else "Fig3.png"), revised, width = 7.5, height = 6.65, dpi = 300, bg = "white")
ggsave(file.path(out, if (public_package) "Fig3_bulk_composition_and_proteomics_review_600dpi.tiff" else "Fig3.tiff"), revised, width = 7.5, height = 6.65, dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(out, "Fig3.pdf"), revised, width = 7.5, height = 6.65, device = cairo_pdf, bg = "white")
saveRDS(revised, file.path(out, "Fig3_editable_gtable.rds"))
capture.output(sessionInfo(), file = file.path(qa, "R_sessionInfo.txt"))
cat("FIG3_EXPORT_OK\n")
