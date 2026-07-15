local_lib <- normalizePath("work/Rlib", mustWork = FALSE)
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

options(
  repos = c(
    jinworks = "https://jinworks.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  ),
  timeout = 1200
)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_lib, type = "binary")
}

bioc_dependencies <- c("BiocGenerics", "ComplexHeatmap", "BiocNeighbors")
missing_bioc <- bioc_dependencies[!vapply(bioc_dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, lib = local_lib, ask = FALSE, update = FALSE, type = "binary")
}

cran_dependencies <- c(
  "data.table", "future", "future.apply", "pbapply", "dplyr", "ggplot2",
  "ggalluvial", "svglite", "circlize", "igraph", "NMF", "RcppEigen",
  "Rtsne", "cowplot", "reshape2", "reticulate", "patchwork", "colorspace",
  "stringr", "FNN", "irlba", "RSpectra", "uwot", "clue", "scam", "sna",
  "foreach", "doParallel"
)
missing_cran <- cran_dependencies[!vapply(cran_dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, lib = local_lib, type = "binary", Ncpus = 2)
}

if (!requireNamespace("CellChat", quietly = TRUE)) {
  install.packages("CellChat", lib = local_lib, type = "binary", dependencies = FALSE)
}

packages <- c("CellChat", "data.table", "Matrix", "future")
for (package in packages) {
  cat(package, "=", as.character(packageVersion(package)), "\n")
}
