# =============================================================================
# A5: Monocle3 Pseudotime Trajectory of FAP+ CAF Activation
# Prerequisite: A2 CAF subset H5AD loaded in Seurat
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(SeuratWrappers)
  library(dplyr); library(ggplot2); library(ggpubr)
})

PROJ_ROOT <- "."
OUT_DIR   <- file.path(PROJ_ROOT, "output"); FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR   <- file.path(OUT_DIR, "tables"); DATA_DIR <- file.path(PROJ_ROOT, "data")

ECM_SDC4_CD44_GENES <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")

message("=== A5: FAP+ CAF Pseudotime Trajectory ===")

# ---- Step 1: Convert Seurat CAF object to Monocle3 CDS ----
message("[A5-1] Converting Seurat CAF subset to CDS...")

seurat_path <- file.path(DATA_DIR, "A2_CAF_seurat.rds")
if (!file.exists(seurat_path)) {
  # Convert from H5AD
  message("  Converting H5AD -> Seurat -> CDS...")
  if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    install.packages("remotes")
    remotes::install_github("mojaveazure/seurat-disk")
  }
  library(SeuratDisk)
  Convert(file.path(DATA_DIR, "A2_GSE132465_CAF_subset.h5ad"),
          file.path(DATA_DIR, "A2_CAF_subset.h5Seurat"), overwrite = TRUE)
  caf <- LoadH5Seurat(file.path(DATA_DIR, "A2_CAF_subset.h5Seurat"))
  saveRDS(caf, seurat_path)
} else {
  caf <- readRDS(seurat_path)
}

# Convert to CDS
cds <- as.cell_data_set(caf)
cds <- cluster_cells(cds, resolution = 1e-5)
cds <- learn_graph(cds)

# ---- Step 2: Order cells by FAP expression gradient ----
message("[A5-2] Ordering cells by FAP+ activation gradient...")

# Set root cells: FAP expression = 0 (least activated)
fap_expr <- FetchData(caf, "FAP")[,1]
root_cells <- colnames(caf)[fap_expr == 0]
if (length(root_cells) == 0) {
  root_cells <- colnames(caf)[which.min(fap_expr)]
}
cds <- order_cells(cds, root_cells = root_cells)

# ---- Step 3: Pseudotime vs ECM-SDC4/CD44 score ----
message("[A5-3] Mapping ECM-SDC4/CD44 along pseudotime...")

# Score ECM genes
caf <- AddModuleScore(caf, features = list(ECM_SDC4_CD44_GENES), name = "ECM_SDC4_CD44")
cds$ECM_score <- caf$ECM_SDC4_CD441

# Plot pseudotime vs ECM score
p_pseudo <- plot_cells(cds, color_cells_by = "pseudotime", 
                       label_cell_groups = FALSE, label_leaves = FALSE,
                       label_branch_points = FALSE) +
  labs(title = "FAP+ CAF Pseudotime Trajectory")
ggsave(file.path(FIG_DIR, "A5_pseudotime_trajectory.png"), p_pseudo, width = 8, height = 6, dpi = 300)

# Correlation: pseudotime ~ ECM score
cor_test <- cor.test(cds$pseudotime, cds$ECM_score, method = "spearman")
message("  Pseudotime vs ECM-SDC4/CD44 score: rho = ", signif(cor_test$estimate, 3),
        ", P = ", signif(cor_test$p.value, 3))

# ---- Step 4: Gene dynamics along pseudotime ----
message("[A5-4] Identifying pseudotime-dependent genes...")

# Test for genes that vary with pseudotime
pr_test <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
pr_test <- pr_test[order(pr_test$morans_I, decreasing = TRUE), ]
write.csv(pr_test, file.path(TAB_DIR, "A5_pseudotime_genes.csv"))

# Plot key ECM genes along pseudotime
ecm_in_data <- ECM_SDC4_CD44_GENES[ECM_SDC4_CD44_GENES %in% rownames(cds)]

for (gene in ecm_in_data) {
  p <- plot_cells(cds, genes = gene, show_trajectory_graph = FALSE,
                  label_cell_groups = FALSE, label_leaves = FALSE) +
    labs(title = paste(gene, "along Pseudotime"))
  ggsave(file.path(FIG_DIR, paste0("A5_pseudotime_", gene, ".png")),
         p, width = 8, height = 6, dpi = 300)
}

# ---- Step 5: Summary ----
message("\n=== A5 Complete ===")
message("Key pseudotime genes (top 10 Morans I):")
print(head(pr_test, 10))

message("Pseudotime ~ ECM-SDC4/CD44 correlation: rho = ", 
        signif(cor_test$estimate, 3), ", P = ", signif(cor_test$p.value, 3))
message("Output: ", FIG_DIR, "/A5_pseudotime_*.png")
