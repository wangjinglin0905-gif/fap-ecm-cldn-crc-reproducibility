# =============================================================================
# REVIEW_R2: Patient-level single-cell re-analysis (GSE132465)
# Verify v3.4/codex numbers:
#   - Stromal FAP vs stromal matrix4:        rho = 0.857 (FDR < 0.001)
#   - Stromal FAP vs epithelial receptor2:   rho = 0.593 (FDR = 0.026)
#   - Stromal matrix4 vs epithelial receptor2: rho = 0.604 (FDR = 0.026)
#   - Stromal FAP vs stromal receptor2:      rho = -0.075 (FDR = 0.791)
# Patient-level: 15 patients with >=20 cells per compartment.
# =============================================================================
suppressPackageStartupMessages({
  library(monocle3); library(SummarizedExperiment); library(dplyr)
})

PROJ <- "."
OUT  <- file.path(PROJ, "output", "review_r2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(file.path(PROJ, "output", "A5_monocle3_cds.rds"))
md  <- as.data.frame(colData(cds))
expr <- assay(cds)

MATRIX4   <- c("COL1A1","COL1A2","COL3A1","FN1")
RECEPTOR2 <- c("SDC4","CD44")
FAP13     <- c("FAP","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")

# z-score mean across cells for a given cell subset
zmean_vec <- function(expr_mat, genes, cells) {
  g <- genes[genes %in% rownames(expr_mat)]
  if (length(g) == 0) return(rep(NA, length(cells)))
  m <- as.matrix(expr_mat[g, cells, drop = FALSE])
  z <- t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
  colMeans(z, na.rm = TRUE)
}

# ---- Compartment definitions ----
fib_lineage <- md$Class == "Tumor" & md$Cell_subtype %in% c("Myofibroblasts","Stromal 1","Stromal 2","Stromal 3")
epi_cells <- md$Class == "Tumor" & grepl("Epithelial", md$Cell_subtype)
cat("tumor fib-lineage cells:", sum(fib_lineage), "| tumor epithelial cells:", sum(epi_cells), "\n")

# Per-patient aggregate scores
patients <- unique(md$Patient[fib_lineage])
res <- data.frame(patient = patients, n_fib = NA, n_epi = NA,
                  fib_FAP = NA, fib_matrix4 = NA, fib_receptor2 = NA,
                  epi_receptor2 = NA, stringsAsFactors = FALSE)
for (i in seq_len(nrow(res))) {
  p <- res$patient[i]
  fc <- md$Patient == p & fib_lineage
  ec <- md$Patient == p & epi_cells
  res$n_fib[i] <- sum(fc); res$n_epi[i] <- sum(ec)
  if (sum(fc) >= 20) {
    res$fib_FAP[i]       <- mean(expr["FAP", md$Patient == p & fib_lineage])
    res$fib_matrix4[i]   <- zmean_vec(expr, MATRIX4, md$Patient == p & fib_lineage)
    res$fib_receptor2[i] <- zmean_vec(expr, RECEPTOR2, md$Patient == p & fib_lineage)
  }
  if (sum(ec) >= 20) {
    res$epi_receptor2[i] <- zmean_vec(expr, RECEPTOR2, md$Patient == p & epi_cells)
  }
}
cat("patients with >=20 fib cells:", sum(res$n_fib >= 20), "\n")
cat("patients with >=20 epi cells:", sum(res$n_epi >= 20), "\n")

# ---- Patient-level correlations (complete pairs) ----
valid <- res[!is.na(res$fib_FAP) & !is.na(res$fib_matrix4) & !is.na(res$fib_receptor2) & !is.na(res$epi_receptor2), ]
cat("\ncomplete-pair patients:", nrow(valid), "\n")

c1 <- cor.test(valid$fib_FAP, valid$fib_matrix4, method = "spearman")
c2 <- cor.test(valid$fib_FAP, valid$epi_receptor2, method = "spearman")
c3 <- cor.test(valid$fib_matrix4, valid$epi_receptor2, method = "spearman")
c4 <- cor.test(valid$fib_FAP, valid$fib_receptor2, method = "spearman")
cat(sprintf("fib FAP vs fib matrix4:     rho=%.3f P=%.4f\n", c1$estimate, c1$p.value))
cat(sprintf("fib FAP vs epi receptor2:   rho=%.3f P=%.4f\n", c2$estimate, c2$p.value))
cat(sprintf("fib matrix4 vs epi receptor2: rho=%.3f P=%.4f\n", c3$estimate, c3$p.value))
cat(sprintf("fib FAP vs fib receptor2:   rho=%.3f P=%.4f\n", c4$estimate, c4$p.value))

# BH across 4 tests
pvals <- c(c1$p.value, c2$p.value, c3$p.value, c4$p.value)
fdr <- p.adjust(pvals, method = "BH")
cat("FDR:", round(fdr, 3), "\n")

sum_tab <- data.frame(
  comparison = c("Stromal FAP vs stromal matrix4",
                 "Stromal FAP vs epithelial receptor2",
                 "Stromal matrix4 vs epithelial receptor2",
                 "Stromal FAP vs stromal receptor2"),
  rho = round(c(c1$estimate, c2$estimate, c3$estimate, c4$estimate), 3),
  P = signif(pvals, 3),
  FDR = signif(fdr, 3),
  stringsAsFactors = FALSE)
write.csv(sum_tab, file.path(OUT, "review_patient_level.csv"), row.names = FALSE)
print(sum_tab)
cat("\n=== [DONE] REVIEW_R2 patient-level ===\n")
