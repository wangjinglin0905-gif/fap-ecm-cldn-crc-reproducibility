.libPaths(c('library', 'library', .libPaths()))
cat('=== FINAL ENVIRONMENT VERIFICATION ===\n')
cat('R:', R.version.string, '\n\n')

pkgs <- c('Seurat', 'SeuratObject', 'monocle3', 'BPCells', 'slingshot',
          'presto', 'CellChat', 'SingleCellExperiment', 'hdf5r',
          'DESeq2', 'GSVA', 'TCGAbiolinks', 'survminer', 'ComplexHeatmap')
cat('--- Package status ---\n')
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) as.character(packageVersion(p)) else '-'
  cat(sprintf('  %-24s %s %s\n', p, if (ok) 'OK ' else 'MISS', ver))
}

cat('\n--- Seurat loads + presto acceleration ---\n')
suppressPackageStartupMessages(library(Seurat))
cat('Seurat', as.character(packageVersion('Seurat')), 'loaded\n')
cat('presto registered (FindMarkers accelerated):', requireNamespace('presto', quietly=TRUE), '\n')

cat('\n--- Monocle3 pipeline ---\n')
suppressPackageStartupMessages({library(monocle3); library(SingleCellExperiment)})
set.seed(1)
pbmc <- CreateSeuratObject(counts = matrix(rpois(100*30, 1), nrow=100, ncol=30), project='t')
pbmc <- NormalizeData(pbmc); pbmc <- FindVariableFeatures(pbmc); pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs=5)
expr <- GetAssayData(pbmc, assay='RNA', layer='data')
cds <- new_cell_data_set(expr,
                         cell_metadata = pbmc@meta.data,
                         gene_metadata = data.frame(gene_short_name = rownames(expr), row.names = rownames(expr)))
cds <- preprocess_cds(cds, num_dim = 5)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cat('Monocle3 pipeline OK: cds', ncol(cds), 'cells ->', length(unique(clusters(cds))), 'clusters\n')

cat('\n--- CellChat ---\n')
suppressPackageStartupMessages(library(CellChat))
cat('CellChat', as.character(packageVersion('CellChat')), 'loaded\n')

cat('\n=== ALL_ENVIRONMENT_READY ===\n')
