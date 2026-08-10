# REVIEW_R6b: left/right colon stratification - FAP13-matrix4 correlations
suppressPackageStartupMessages({library(dplyr)})
PROJ <- "."
TCGA_EXPR <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"
FAP13 <- c("FAP","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
MATRIX4 <- c("COL1A1","COL1A2","COL3A1","FN1")
RECEPTOR2 <- c("SDC4","CD44")

zmean <- function(df, genes) {
  avail <- genes[genes %in% colnames(df)]
  sub <- as.data.frame(lapply(avail, function(g) (df[[g]] - mean(df[[g]], na.rm=TRUE))/sd(df[[g]], na.rm=TRUE)))
  rowMeans(sub, na.rm = TRUE)
}
expr_raw <- read.table(gzfile(TCGA_EXPR), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
tcga <- read.csv(file.path(PROJ, "data", "A1_tcga_coad_merged.csv"))
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) vals[i] <- expr_mat[gene, sid]
    else { pid <- substr(sid, 1, 12); m <- grep(pid, colnames(expr_mat), value=TRUE, fixed=TRUE); if (length(m)>0) vals[i] <- expr_mat[gene, m[1]] }
  }
  vals
}
for (g in unique(c(FAP13, MATRIX4, RECEPTOR2))) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga <- tcga[tcga$sample_scope == "Tumor", ]
tcga$FAP13 <- zmean(tcga, FAP13)
tcga$matrix4 <- zmean(tcga, MATRIX4)
tcga$receptor2 <- zmean(tcga, RECEPTOR2)
site <- tcga$site_of_resection_or_biopsy
tcga$side <- ifelse(site %in% c("Cecum","Ascending colon","Hepatic flexure of colon","Transverse colon"), "Right",
             ifelse(site %in% c("Splenic flexure of colon","Descending colon","Sigmoid colon"), "Left", NA))
cat("n by side:", table(tcga$side, useNA="ifany"), "\n\n")
for (sd in c("Left", "Right")) {
  sub <- tcga[tcga$side == sd, ]
  cat(sprintf("=== %s colon (n=%d) ===\n", sd, nrow(sub)))
  cc <- cor.test(sub$FAP13, sub$matrix4, method = "spearman")
  cat(sprintf("FAP13-matrix4: rho=%.3f P=%.4f\n", cc$estimate, cc$p.value))
  cr <- cor.test(sub$FAP13, sub$receptor2, method = "spearman")
  cat(sprintf("FAP13-receptor2: rho=%.3f P=%.4f\n", cr$estimate, cr$p.value))
  cm <- cor.test(sub$matrix4, sub$receptor2, method = "spearman")
  cat(sprintf("matrix4-receptor2: rho=%.3f P=%.4f\n", cm$estimate, cm$p.value))
}
# overall for reference
cat("\n=== Overall (n=89) ===\n")
cat(sprintf("FAP13-matrix4 rho=%.3f\n", cor.test(tcga$FAP13, tcga$matrix4, method="spearman")$estimate))
write.csv(tcga[, c("sample_id","side","FAP13","matrix4","receptor2")], file.path(PROJ,"output","review_r6","tcga_side_scores.csv"), row.names=FALSE)
cat("\n[DONE] R6b left/right\n")
