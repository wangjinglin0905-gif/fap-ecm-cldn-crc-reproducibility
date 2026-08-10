#!/usr/bin/env Rscript
# =====================================================================
# 08_A5_monocle3_pseudotime.R  (v2 - inject Seurat UMAP, fix root)
# A5 升级：Monocle3 拟时序，注入 Seurat UMAP 保证轨迹方向正确
# 根 = FAP 最低的细胞（静止成纤维端）-> FAP+ 活化 CAF 端
# =====================================================================
.libPaths(c('library', 'library', .libPaths()))
suppressPackageStartupMessages({
  library(Seurat); library(monocle3); library(ggplot2); library(dplyr)
})
set.seed(20260802)

SC_RDS <- 'data/GSE132465/integrated_seurat_Harmony.rds'
OUTD <- './output'
FIGD <- file.path(OUTD, 'figures'); TABD <- file.path(OUTD, 'tables')
dir.create(FIGD, recursive = TRUE, showWarnings = FALSE)
dir.create(TABD, recursive = TRUE, showWarnings = FALSE)

cat('=== 1. Load & subset ===\n')
sc <- readRDS(SC_RDS)
stromal <- subset(sc, subset = Cell_type == 'Stromal cells')
fap_expr <- GetAssayData(stromal, assay = 'RNA', layer = 'data')['FAP', ]
stromal$FAP_status <- ifelse(fap_expr > 0, 'FAP+', 'FAP-')
cat('Stromal:', ncol(stromal), '| FAP+:', sum(stromal$FAP_status=='FAP+'), '\n')

cat('=== 2. Build cds & inject Seurat UMAP ===\n')
expr_mat <- GetAssayData(stromal, assay = 'RNA', layer = 'counts')
cell_meta <- stromal@meta.data
gene_meta <- data.frame(gene_short_name = rownames(expr_mat), row.names = rownames(expr_mat))
cds <- new_cell_data_set(expr_mat, cell_metadata = cell_meta, gene_metadata = gene_meta)

# Inject Seurat UMAP (cells in same order as expr_mat columns)
seurat_umap <- Embeddings(stromal, 'umap')
reducedDims(cds)[['UMAP']] <- as.matrix(seurat_umap[colnames(cds), , drop = FALSE])
cat('UMAP injected:', nrow(reducedDims(cds)[['UMAP']]), 'cells x', ncol(reducedDims(cds)[['UMAP']]), 'dims\n')

cat('=== 3. Preprocess (PCA for gene space) + cluster on injected UMAP ===\n')
cds <- preprocess_cds(cds, num_dim = 25)
cds <- cluster_cells(cds, reduction_method = 'UMAP')
cat('clusters:', length(unique(clusters(cds))), '\n')

cat('=== 4. Learn graph on injected UMAP ===\n')
cds <- learn_graph(cds, use_partition = FALSE)

cat('=== 5. Root = FAP-lowest cells (quiescent end) ===\n')
# FAP expression vector aligned to cds cell order
fap_cds <- fap_expr[colnames(cds)]
n_root <- min(100, sum(fap_cds == 0))
root_cells <- names(sort(fap_cds))[1:n_root]
cds <- order_cells(cds, root_cells = root_cells)
cat('root cells:', n_root, '| pseudotime range:', range(pseudotime(cds)), '\n')

cat('=== 6. Correlations ===\n')
# ECM-SDC4/CD44 z-mean score (independent of AddModuleScore column naming)
score_zmean <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  z <- t(scale(t(as.matrix(mat[genes, , drop = FALSE]))))
  colMeans(z, na.rm = TRUE)
}
ecm_genes <- c('COL1A1', 'COL1A2', 'COL3A1', 'FN1', 'SDC4', 'CD44')
ecm_score <- score_zmean(GetAssayData(stromal, assay='RNA', layer='data'), ecm_genes)
pseudo <- pseudotime(cds)
common <- intersect(names(pseudo), names(ecm_score))
pseudo_vec <- pseudo[common]; ecm_vec <- ecm_score[common]; fap_vec <- fap_expr[common]

cor_ecm <- cor.test(pseudo_vec, ecm_vec, method='spearman', exact=FALSE)
cor_fap <- cor.test(pseudo_vec, fap_vec, method='spearman', exact=FALSE)
cat(sprintf('pseudotime vs ECM-SDC4/CD44: rho=%.3f, P=%.2e (n=%d)\n', cor_ecm$estimate, cor_ecm$p.value, length(common)))
cat(sprintf('pseudotime vs FAP:           rho=%.3f, P=%.2e\n', cor_fap$estimate, cor_fap$p.value))

res_df <- data.frame(comparison=c('pseudotime_vs_ECM_SDC4_CD44','pseudotime_vs_FAP'),
                     rho=c(cor_ecm$estimate, cor_fap$estimate),
                     p=c(cor_ecm$p.value, cor_fap$p.value), n=length(common))
write.csv(res_df, file.path(TABD, 'A5_pseudotime_correlations.csv'), row.names=FALSE)

pt_df <- data.frame(cell=common, pseudotime=pseudo_vec, ECM_SDC4_CD44=ecm_vec,
                    FAP=fap_vec, FAP_status=stromal$FAP_status[common])
write.csv(pt_df, file.path(TABD, 'A5_pseudotime_cell_table.csv'), row.names=FALSE)

cat('=== 7. graph_test ===\n')
pr_test <- graph_test(cds, neighbor_graph='principal_graph', cores=4)
pr_test <- pr_test[order(pr_test$morans_I, decreasing=TRUE), ]
write.csv(pr_test, file.path(TABD, 'A5_pseudotime_graph_test.csv'), row.names=FALSE)
cat('Top 10:\n'); print(head(pr_test[, c('gene_short_name','morans_I','q_value')], 10))
ecm_rank <- pr_test[pr_test$gene_short_name %in% c(ecm_genes, 'FAP'), c('gene_short_name','morans_I','q_value')]
cat('\nECM/FAP in graph_test:\n'); print(ecm_rank)

cat('=== 8. Figures ===\n')
p1 <- plot_cells(cds, color_cells_by='pseudotime', label_groups_by_cluster=FALSE,
                 label_leaves=FALSE, label_branch_points=FALSE, cell_size=0.7) +
  ggtitle('Stromal cell pseudotime trajectory (GSE132465)') + theme(legend.position='right')
ggsave(file.path(FIGD,'A5_trajectory_pseudotime.png'), p1, width=7, height=5.5, dpi=300)
ggsave(file.path(FIGD,'A5_trajectory_pseudotime.tiff'), p1, width=7, height=5.5, dpi=600, compression='lzw')

p2 <- ggplot(pt_df, aes(pseudotime, ECM_SDC4_CD44)) +
  geom_point(aes(color=FAP_status), alpha=0.5, size=1) +
  geom_smooth(method='loess', color='#E64B35', se=TRUE) +
  scale_color_manual(values=c('FAP-'='#3C5488','FAP+'='#E64B35')) +
  labs(title=sprintf('ECM-SDC4/CD44 along pseudotime (rho=%.3f, P<0.001)', cor_ecm$estimate),
       x='Pseudotime', y='ECM-SDC4/CD44 score') + theme_classic()
ggsave(file.path(FIGD,'A5_ECM_vs_pseudotime.png'), p2, width=6.5, height=5, dpi=300)
ggsave(file.path(FIGD,'A5_ECM_vs_pseudotime.tiff'), p2, width=6.5, height=5, dpi=600, compression='lzw')

pt_sorted <- pt_df[order(pt_df$pseudotime), ]
gene_dyn <- lapply(c('FAP','COL1A1','FN1','SDC4','CD44'), function(g) {
  gv <- GetAssayData(stromal, assay='RNA', layer='data')[g, pt_sorted$cell]
  data.frame(pseudotime=pt_sorted$pseudotime, gene=g, expr=as.numeric(gv))
})
gd <- do.call(rbind, gene_dyn)
p3 <- ggplot(gd, aes(pseudotime, expr)) +
  geom_smooth(method='loess', se=TRUE, aes(color=gene), linewidth=0.9) +
  scale_color_manual(values=c('FAP'='#E64B35','COL1A1'='#3C5488','FN1'='#00A087','SDC4'='#F39B7F','CD44'='#8491B4')) +
  labs(title='Key gene expression along pseudotime', x='Pseudotime', y='Normalized expression') +
  theme_classic() + theme(legend.position='right')
ggsave(file.path(FIGD,'A5_gene_dynamics.png'), p3, width=7, height=5, dpi=300)
ggsave(file.path(FIGD,'A5_gene_dynamics.tiff'), p3, width=7, height=5, dpi=600, compression='lzw')

saveRDS(cds, file.path(OUTD, 'A5_monocle3_cds.rds'))
cat('\n[DONE] A5 Monocle3 pseudotime complete (v2, injected UMAP)\n')
