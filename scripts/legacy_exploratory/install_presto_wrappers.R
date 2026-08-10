.libPaths(c('library', 'library', .libPaths()))
options(timeout = 900)
options(repos = c(CRAN = 'https://mirrors.sjtug.sjtu.edu.cn/cran'))

cat('=== 1. Install presto (GitHub, C++ accelerated Wilcoxon) ===\n')
remotes::install_github('immunogenomics/presto', upgrade = 'never', quiet = FALSE)
cat('presto installed:', requireNamespace('presto', quietly = TRUE), '\n')

cat('\n=== 2. Install SeuratWrappers (GitHub, Seurat <-> Monocle3 bridge) ===\n')
remotes::install_github('satijalab/seurat-wrappers', upgrade = 'never', quiet = FALSE)
cat('SeuratWrappers installed:', requireNamespace('SeuratWrappers', quietly = TRUE), '\n')

cat('\n=== Final verification ===\n')
for (p in c('Seurat', 'monocle3', 'BPCells', 'slingshot', 'presto', 'SeuratWrappers', 'CellChat')) {
  ok <- requireNamespace(p, quietly = TRUE)
  ver <- if (ok) as.character(packageVersion(p)) else '-'
  cat(sprintf('  %-18s %s %s\n', p, if (ok) 'OK ' else 'MISS', ver))
}
cat('\nALL_DONE\n')
