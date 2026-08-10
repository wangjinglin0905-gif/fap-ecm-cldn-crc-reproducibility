.libPaths(c('library', 'library', .libPaths()))
options(timeout = 900)
options(repos = c(CRAN = 'https://mirrors.sjtug.sjtu.edu.cn/cran'))

# SeuratWrappers: small package, mostly R code with as.cell_data_set conversion
# Try pak (better GitHub installer) or fallback to direct tarball + R CMD INSTALL
cat('=== Trying SeuratWrappers via remotes with forced build tools ===\n')
Sys.setenv(RTOOLS45_HOME = 'C:/rtools45')
# remotes checks pkgbuild::has_build_tools() which fails due to sandbox safe-delete on tmp.def
# Workaround: skip build-tools check by installing source tarball manually
cat('Downloading SeuratWrappers tarball...\n')
dl <- tryCatch(download.file(
  'https://github.com/satijalab/seurat-wrappers/archive/refs/heads/master.tar.gz',
  destfile = './seurat-wrappers.tar.gz',
  method = 'libcurl', quiet = FALSE), error = function(e) {cat('DL ERR:', conditionMessage(e), '\n'); NA})
cat('download status:', dl, '\n')
