# =============================================================================
# REVIEW_R3: Fibroblast-lineage state ordering (light verification)
# Verify direction of codex claims: activation rho=0.580, matrix4 rho=0.422,
# SDC4/CD44 rho=-0.308 (patient-level, signed-rank, BH).
# Light version: use the A5 monocle3 object restricted to tumor fibroblast-lineage
# cells; compute a quiescence score (PI16/COL14A1/CFD/DCN/DPT/C7) as ordering axis;
# correlate patient-level mean scores with ordering axis.
# =============================================================================
suppressPackageStartupMessages({
  library(monocle3); library(SummarizedExperiment); library(dplyr)
})
PROJ <- "."
OUT  <- file.path(PROJ, "output", "review_r3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cds <- readRDS(file.path(PROJ, "output", "A5_monocle3_cds.rds"))
md  <- as.data.frame(colData(cds))
expr <- assay(cds)

FIB_LINEAGE <- c("Myofibroblasts","Stromal 1","Stromal 2","Stromal 3")
QUIESCENCE  <- c("PI16","COL14A1","CFD","DCN","DPT","C7")
ACTIVATION  <- c("FAP","POSTN","THY1","TAGLN","ACTA2","COL1A1","COL1A2","COL3A1","FN1")
MATRIX4     <- c("COL1A1","COL1A2","COL3A1","FN1")
RECEPTOR2   <- c("SDC4","CD44")

# restrict to tumor fibroblast-lineage cells
keep <- md$Class == "Tumor" & md$Cell_subtype %in% FIB_LINEAGE
cat("tumor fib-lineage cells:", sum(keep), "\n")
mdk <- md[keep, ]
all_genes_needed <- unique(c(QUIESCENCE, ACTIVATION, RECEPTOR2))
present <- all_genes_needed[all_genes_needed %in% rownames(expr)]
cat("genes present:", length(present), "/", length(all_genes_needed), "\n")
cat("missing:", setdiff(all_genes_needed, rownames(expr)), "\n")

# normalize full subset: log1p(CPM), guard zero-total columns
sub <- as.matrix(expr[present, keep, drop = FALSE])
tot <- colSums(sub)
tot[tot == 0] <- NA
norm <- log1p(t(sub) / tot * 10000)
norm <- t(norm)
keep2 <- !is.na(colSums(norm))
norm <- norm[, keep2]
mdk <- mdk[keep2, ]
cat("cells after zero-total guard:", ncol(norm), "\n")
cat("assay max after norm:", max(norm), "\n")

zmean_cells <- function(mat, genes, cells) {
  g <- genes[genes %in% rownames(mat)]
  if (length(g) == 0) return(rep(NA, length(cells)))
  m <- as.matrix(mat[g, cells, drop = FALSE])
  z <- t(apply(m, 1, function(r) (r - mean(r)) / sd(r)))
  colMeans(z, na.rm = TRUE)
}

cells_idx <- which(keep)
qui <- zmean_cells(norm, QUIESCENCE, seq_len(ncol(norm)))
act <- zmean_cells(norm, ACTIVATION, seq_len(ncol(norm)))
m4  <- zmean_cells(norm, MATRIX4, seq_len(ncol(norm)))
r2  <- zmean_cells(norm, RECEPTOR2, seq_len(ncol(norm)))

# ordering axis: quiescence score (lower = more activated)
df <- data.frame(cell = colnames(norm), Patient = mdk$Patient,
                 qui = qui, act = act, m4 = m4, r2 = r2)
df <- df[!is.na(df$qui) & !is.na(df$act) & !is.na(df$m4) & !is.na(df$r2), ]
cat("cells with all scores:", nrow(df), "| patients:", length(unique(df$Patient)), "\n")

# patient-level: correlate mean score with mean quiescence across patients
pat_res <- df %>% group_by(Patient) %>%
  summarise(qui = mean(qui), act = mean(act), m4 = mean(m4), r2 = mean(r2), n = n()) %>%
  filter(n >= 20)
cat("patients >=20 cells:", nrow(pat_res), "\n")

c_act <- cor.test(pat_res$act, pat_res$qui, method = "spearman")
c_m4  <- cor.test(pat_res$m4,  pat_res$qui, method = "spearman")
c_r2  <- cor.test(pat_res$r2,  pat_res$qui, method = "spearman")
cat(sprintf("activation vs quiescence: rho=%.3f P=%.4f\n", c_act$estimate, c_act$p.value))
cat(sprintf("matrix4 vs quiescence:    rho=%.3f P=%.4f\n", c_m4$estimate, c_m4$p.value))
cat(sprintf("receptor2 vs quiescence:  rho=%.3f P=%.4f\n", c_r2$estimate, c_r2$p.value))
# Note: quiescence HIGH = resting; so activation should be NEGATIVELY correlated with quiescence.
# codex reports "activation rho=0.580" likely vs pseudotime (activated end). We report raw sign here.

pvals <- c(c_act$p.value, c_m4$p.value, c_r2$p.value)
fdr <- p.adjust(pvals, method = "BH")
cat("FDR:", round(fdr, 3), "\n")

sum_tab <- data.frame(
  comparison = c("activation vs quiescence", "matrix4 vs quiescence", "receptor2 vs quiescence"),
  rho = round(c(c_act$estimate, c_m4$estimate, c_r2$estimate), 3),
  P = signif(pvals, 3), FDR = signif(fdr, 3), stringsAsFactors = FALSE)
write.csv(sum_tab, file.path(OUT, "review_state_ordering.csv"), row.names = FALSE)
print(sum_tab)
cat("\n=== [DONE] REVIEW_R3 ===\n")
