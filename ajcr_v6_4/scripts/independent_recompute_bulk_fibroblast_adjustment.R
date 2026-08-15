options(width = 180, digits = 8)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(paste(
    "Usage: Rscript independent_recompute_bulk_fibroblast_adjustment.R",
    "<bulk_input_directory> <deconvolution_directory> <senmayo_genes.txt>"
  ))
}
base_dir <- normalizePath(args[[1]], mustWork = TRUE)
deconv_dir <- normalizePath(args[[2]], mustWork = TRUE)
senmayo_file <- normalizePath(args[[3]], mustWork = TRUE)

coad <- read.delim(gzfile(file.path(base_dir, "COAD_HiSeqV2.gz")), row.names = 1, check.names = FALSE)
read <- read.delim(gzfile(file.path(base_dir, "READ_HiSeqV2.gz")), row.names = 1, check.names = FALSE)
coad <- coad[, grepl("-01$", colnames(coad)), drop = FALSE]
read <- read[, grepl("-01$", colnames(read)), drop = FALSE]
common_genes <- intersect(rownames(coad), rownames(read))
expr <- cbind(coad[common_genes, , drop = FALSE], read[common_genes, , drop = FALSE])

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
MATRIX4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
FIB5 <- c("PDGFRA", "PDGFRB", "LUM", "DCN", "COL14A1")

zmean <- function(mat, genes) {
  sub <- mat[genes, , drop = FALSE]
  z <- t(scale(t(sub), center = TRUE, scale = TRUE))
  colMeans(z)
}

fap <- zmean(expr, FAP13)
matrix4 <- zmean(expr, MATRIX4)
fib5 <- zmean(expr, FIB5)

sen <- unique(trimws(readLines(senmayo_file)))
sen <- sen[nzchar(sen)]
sen <- setdiff(sen, intersect(sen, FAP13))
sen[sen == "CXCL8"] <- "IL8"
sen <- unique(sen)
sen_use <- intersect(sen, rownames(expr))
senmayo <- zmean(expr, sen_use)

rank_residual <- function(x, covars) {
  xrank <- rank(x, ties.method = "average")
  cov_rank <- apply(as.matrix(covars), 2, rank, ties.method = "average")
  if (is.null(dim(cov_rank))) cov_rank <- matrix(cov_rank, ncol = 1)
  resid(lm(xrank ~ cov_rank))
}

partial_standard <- function(x, y, covars) {
  cor(rank_residual(x, covars), rank_residual(y, covars), method = "pearson")
}

partial_double_rank <- function(x, y, covars) {
  cor(rank_residual(x, covars), rank_residual(y, covars), method = "spearman")
}

boot_ci <- function(x, y, covars, n_boot = 5000, seed = 20260810) {
  set.seed(seed)
  n <- length(x)
  vals <- replicate(n_boot, {
    ii <- sample.int(n, n, replace = TRUE)
    partial_standard(x[ii], y[ii], as.matrix(covars)[ii, , drop = FALSE])
  })
  quantile(vals, c(0.025, 0.975), na.rm = TRUE)
}

partial_p <- function(r, n, k) {
  tval <- r * sqrt((n - k - 2) / (1 - r^2))
  2 * pt(-abs(tval), df = n - k - 2)
}

mcp <- read.csv(file.path(deconv_dir, "mcpcounter_tumor_scores.csv"), row.names = 1, check.names = FALSE)
epic <- read.csv(file.path(deconv_dir, "epic_tumor_cell_fractions.csv"), row.names = 1, check.names = FALSE)
samples <- Reduce(intersect, list(names(fap), colnames(mcp), colnames(epic)))

cat("DATA\n")
cat("COAD primary:", ncol(coad), "READ primary:", ncol(read), "combined:", ncol(expr), "deconv-aligned:", length(samples), "\n")
cat("SenMayo overlap-removed represented:", length(sen_use), "missing:", paste(setdiff(sen, rownames(expr)), collapse = ","), "\n")
cat("marginal FAP13-matrix4 rho:", cor(fap, matrix4, method = "spearman"), "\n\n")

cat("FIB5 IMPLEMENTATION CHECK\n")
r_double <- partial_double_rank(fap, matrix4, fib5)
r_standard <- partial_standard(fap, matrix4, fib5)
cat("double-ranked residual Spearman (v5.5 script):", r_double, "\n")
cat("standard partial Spearman = Pearson(rank residuals):", r_standard, "CI", paste(boot_ci(fap, matrix4, fib5), collapse = " to "), "P", partial_p(r_standard, length(fap), 1), "\n\n")

specs <- list(
  fib5 = matrix(fib5[samples], ncol = 1, dimnames = list(samples, "fib5")),
  MCPcounter_Fibroblasts = matrix(as.numeric(mcp["Fibroblasts", samples]), ncol = 1, dimnames = list(samples, "MCP")),
  EPIC_CAFs = matrix(as.numeric(epic["CAFs", samples]), ncol = 1, dimnames = list(samples, "EPIC")),
  MCP_EPIC = cbind(MCP = as.numeric(mcp["Fibroblasts", samples]), EPIC = as.numeric(epic["CAFs", samples]))
)

cat("DECONVOLUTION-ADJUSTED RESULTS\n")
targets <- list(
  FAP13_matrix4 = list(x = fap[samples], y = matrix4[samples]),
  SenMayo_FAP13 = list(x = senmayo[samples], y = fap[samples]),
  SenMayo_matrix4 = list(x = senmayo[samples], y = matrix4[samples])
)

for (target_name in names(targets)) {
  cat(target_name, "\n")
  target <- targets[[target_name]]
  for (spec_name in names(specs)) {
    covars <- specs[[spec_name]]
    r <- partial_standard(target$x, target$y, covars)
    ci <- boot_ci(target$x, target$y, covars)
    p <- partial_p(r, length(target$x), ncol(covars))
    cat("  ", spec_name, ": rho=", r, " CI=", ci[1], " to ", ci[2], " P=", p, "\n", sep = "")
  }
}

cat("DONE\n")
