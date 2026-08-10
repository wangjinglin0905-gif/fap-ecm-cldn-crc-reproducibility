.libPaths(c('library', 'library', .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SingleCellExperiment)
})

cat('=== Monocle3 conversion path ===\n')
# 1. does monocle3 provide as.cell_data_set method for SCE?
m <- methods('as.cell_data_set')
cat('monocle3 as.cell_data_set methods:', paste(m, collapse = '; '), '\n')
# 2. Seurat has as.SingleCellExperiment
cat('Seurat as.SingleCellExperiment:', exists('as.SingleCellExperiment', where = asNamespace('Seurat')), '\n')

# 3. Full round-trip smoke test on tiny data
set.seed(1)
pbmc <- CreateSeuratObject(counts = matrix(rpois(200 * 50, 1), nrow = 200, ncol = 50), project = 'test')
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 10)
pbmc <- FindNeighbors(pbmc, dims = 1:5)
pbmc <- FindClusters(pbmc, resolution = 0.1)

sce <- as.SingleCellExperiment(pbmc)
cds <- as.cell_data_set(sce)
cat('\n=== Smoke test OK ===\n')
cat('Seurat clusters:', nlevels(Idents(pbmc)), '\n')
cat('cds cells:', ncol(cds), '| genes:', nrow(cds), '\n')
cat('cds class OK:', class(cds)[1], '\n')

# 4. Verify core analysis functions exist
cat('\n=== Key functions ===\n')
for (fn in c('learn_graph', 'order_cells', 'cluster_cells', 'plot_cells', 'graph_test')) {
  cat(sprintf('  monocle3::%s : %s\n', fn, exists(fn, where = asNamespace('monocle3'))))
}
cat('\nALL_CORE_READY\n')
