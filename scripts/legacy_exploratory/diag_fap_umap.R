.libPaths(c('library', 'library', .libPaths()))
suppressPackageStartupMessages({library(Seurat); library(monocle3)})
set.seed(20260802)

sc <- readRDS('data/GSE132465/integrated_seurat_Harmony.rds')
stromal <- subset(sc, subset = Cell_type == 'Stromal cells')
fap <- GetAssayData(stromal, assay='RNA', layer='data')['FAP', ]
umap <- as.data.frame(Embeddings(stromal, 'umap'))
colnames(umap) <- c('UMAP1', 'UMAP2')
umap$FAP <- fap
umap$FAP_status <- ifelse(fap > 0, 'FAP+', 'FAP-')

cat('=== FAP distribution on UMAP ===\n')
cat('FAP+ cells:', sum(fap > 0), '(', round(mean(fap > 0)*100,1), '%)\n')
cat('FAP max:', max(fap), '| median (all):', median(fap), '| median (FAP+ only):', median(fap[fap>0]), '\n')

# Where are FAP+ cells on UMAP? Median UMAP coords by status
agg <- aggregate(cbind(UMAP1, UMAP2) ~ FAP_status, data = umap, FUN = median)
print(agg)

# Correlation of FAP with UMAP axes
cat('\nFAP vs UMAP1 rho:', cor(fap, umap$UMAP1, method='spearman'), '\n')
cat('FAP vs UMAP2 rho:', cor(fap, umap$UMAP2, method='spearman'), '\n')

# Quick visual check: FAP-high cells (top decile) centroid
top10 <- quantile(fap[fap>0], 0.5)  # median of expressing
high <- umap[fap > top10, ]
cat('\nTop-expressing FAP cells (FAP > median of FAP+): n =', nrow(high), '\n')
cat('  median UMAP1:', median(high$UMAP1), ' UMAP2:', median(high$UMAP2), '\n')
cat('  FAP-low cells (FAP = 0) median UMAP1:', median(umap$UMAP1[umap$FAP == 0]), ' UMAP2:', median(umap$UMAP2[umap$FAP == 0]), '\n')

# Save UMAP with FAP for plotting
write.csv(umap, './output/tables/A5_umap_fap_diagnostic.csv', row.names = TRUE)
cat('\nDiagnostic saved\n')
