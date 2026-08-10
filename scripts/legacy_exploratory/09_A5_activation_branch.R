.libPaths(c('library', 'library', .libPaths()))
suppressPackageStartupMessages({library(monocle3); library(Seurat); library(ggplot2); library(dplyr)})
OUTD <- './output'
FIGD <- file.path(OUTD, 'figures')

cds <- readRDS(file.path(OUTD, 'A5_monocle3_cds.rds'))
sc <- readRDS('data/GSE132465/integrated_seurat_Harmony.rds')
stromal <- subset(sc, subset = Cell_type == 'Stromal cells')
fap <- GetAssayData(stromal, assay='RNA', layer='data')['FAP', ]
score_zmean <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  z <- t(scale(t(as.matrix(mat[genes, , drop=FALSE])))); colMeans(z, na.rm=TRUE)
}
ecm <- score_zmean(GetAssayData(stromal, assay='RNA', layer='data'),
                   c('COL1A1','COL1A2','COL3A1','FN1','SDC4','CD44'))

m <- data.frame(cell = colnames(cds), pseudotime = pseudotime(cds),
                cluster = as.character(clusters(cds)))
m$FAP <- fap[m$cell]; m$ECM <- ecm[m$cell]
m$FAP_status <- ifelse(m$FAP > 0, 'FAP+', 'FAP-')

# Activation branch = cluster 1 (FAP+ activated CAF) vs cluster 2 (quiescent root)
m$branch <- ifelse(m$cluster == '1', 'Activated CAF (FAP+)',
            ifelse(m$cluster == '2', 'Quiescent fibroblast (FAP-)', 'Other stromal'))

# ---- Fig: ECM score by branch (boxplot + stats) ----
p <- ggplot(m[m$branch != 'Other stromal', ], aes(branch, ECM, fill = branch)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.85) +
  scale_fill_manual(values = c('Activated CAF (FAP+)' = '#E64B35', 'Quiescent fibroblast (FAP-)' = '#3C5488')) +
  labs(title = 'ECM-SDC4/CD44 module score: activated vs quiescent fibroblasts',
       x = '', y = 'ECM-SDC4/CD44 score (z-mean of COL1A1/COL1A2/COL3A1/FN1/SDC4/CD44)') +
  theme_classic() + theme(legend.position = 'none')
ggsave(file.path(FIGD, 'A5_activation_branch_ECM.png'), p, width = 5.5, height = 5, dpi = 300)
ggsave(file.path(FIGD, 'A5_activation_branch_ECM.tiff'), p, width = 5.5, height = 5, dpi = 600, compression = 'lzw')

# stats
a <- m$ECM[m$cluster == '1']; q <- m$ECM[m$cluster == '2']
wt <- wilcox.test(a, q)
cat(sprintf('Activated vs Quiescent ECM: mean %.2f vs %.2f, Wilcoxon P=%.2e\n',
            mean(a), mean(q), wt$p.value))

# ---- Fig: trajectory with activation branch highlighted ----
umap <- reducedDims(cds)[['UMAP']]
colnames(umap) <- c('UMAP1', 'UMAP2')
um_df <- data.frame(cell = rownames(umap), umap)
um_df$branch <- m$branch[match(um_df$cell, m$cell)]
um_df$pseudotime <- m$pseudotime[match(um_df$cell, m$cell)]
um_df$FAP <- m$FAP[match(um_df$cell, m$cell)]

p2 <- ggplot(um_df, aes(UMAP1, UMAP2, color = branch)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c('Activated CAF (FAP+)' = '#E64B35',
                                'Quiescent fibroblast (FAP-)' = '#3C5488',
                                'Other stromal' = 'grey80')) +
  labs(title = 'Stromal trajectory: FAP+ activated CAF branch',
       color = '') + theme_classic() + theme(legend.position = 'right')
ggsave(file.path(FIGD, 'A5_trajectory_branch.png'), p2, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(FIGD, 'A5_trajectory_branch.tiff'), p2, width = 7, height = 5.5, dpi = 600, compression = 'lzw')

# Save branch table
write.csv(m, file.path(OUTD, 'tables', 'A5_pseudotime_branch_table.csv'), row.names = FALSE)
cat('[DONE] activation branch figures\n')
