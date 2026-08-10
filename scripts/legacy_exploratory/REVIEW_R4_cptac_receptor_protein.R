# =============================================================================
# REVIEW_R4: CPTAC protein-level validation of the receptor-uncoupling hypothesis
# Question: does SDC4/CD44 protein track the FAP/ECM program in CPTAC-COAD?
# Hypothesis (from v3.5): matrix couples with FAP (positive), receptors do NOT
# co-induce (null or negative). AJCR version only tested FAP-ECM-CLDN proteins;
# this adds SDC4/CD44 receptor proteins.
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})
PROJ <- "."
OUT  <- file.path(PROJ, "output", "cptac")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cpt <- read.csv(file.path(OUT, "cptac_protein_sdc4_cd44.csv"), check.names = FALSE)
cat("CPTAC samples:", nrow(cpt), "\n")
cat("columns:", colnames(cpt), "\n")

# ---- ECM protein score (FAP + COL1A1 + COL1A2 + FN1, z-mean) ----
ecm_genes <- intersect(c("FAP", "COL1A1", "COL1A2", "FN1"), colnames(cpt))
zmean <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
for (g in ecm_genes) cpt[[paste0(g, "_z")]] <- zmean(cpt[[g]])
cpt$ECM_protein <- rowMeans(cpt[, paste0(ecm_genes, "_z")], na.rm = TRUE)
# Receptor protein score
rec_genes <- intersect(c("SDC4", "CD44"), colnames(cpt))
for (g in rec_genes) cpt[[paste0(g, "_z")]] <- zmean(cpt[[g]])
cpt$Receptor_protein <- rowMeans(cpt[, paste0(rec_genes, "_z")], na.rm = TRUE)

# ---- Correlations ----
pairs <- list(
  c("FAP", "SDC4"), c("FAP", "CD44"),
  c("FAP", "ECM_protein"),
  c("ECM_protein", "Receptor_protein"),
  c("ECM_protein", "SDC4"), c("ECM_protein", "CD44"),
  c("SDC4", "CD44"),
  c("ECM_protein", "CLDN4"), c("FAP", "CLDN4")
)
res <- data.frame(x = character(), y = character(), n = numeric(),
                  rho = numeric(), P = numeric(), stringsAsFactors = FALSE)
for (pr in pairs) {
  x <- cpt[[pr[1]]]; y <- cpt[[pr[2]]]
  ok <- complete.cases(x, y)
  if (sum(ok) >= 5) {
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    res <- rbind(res, data.frame(x = pr[1], y = pr[2], n = sum(ok),
                                 rho = unname(ct$estimate), P = ct$p.value))
  } else {
    res <- rbind(res, data.frame(x = pr[1], y = pr[2], n = sum(ok), rho = NA, P = NA))
  }
}
res$FDR <- p.adjust(res$P, method = "BH")
print(res, row.names = FALSE)
write.csv(res, file.path(OUT, "cptac_receptor_correlations.csv"), row.names = FALSE)

# ---- FAP-high vs FAP-low receptor protein (median split) ----
cpt$FAP_grp <- ifelse(cpt$FAP > median(cpt$FAP, na.rm = TRUE), "High", "Low")
for (g in c("SDC4", "CD44", "Receptor_protein", "ECM_protein", "CLDN4")) {
  if (!g %in% colnames(cpt)) next
  hi <- cpt[[g]][cpt$FAP_grp == "High"]; lo <- cpt[[g]][cpt$FAP_grp == "Low"]
  w <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
  cat(sprintf("%s: FAP-high median=%.3f vs FAP-low median=%.3f | P=%.4f\n",
      g, median(hi, na.rm = TRUE), median(lo, na.rm = TRUE), w$p.value))
}
cat("\n=== [DONE] REVIEW_R4 CPTAC receptor protein ===\n")
