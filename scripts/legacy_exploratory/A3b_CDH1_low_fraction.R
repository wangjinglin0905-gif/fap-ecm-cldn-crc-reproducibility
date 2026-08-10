# A3b: CDH1-low threshold analysis (mirror IHC negativity logic)
suppressPackageStartupMessages({library(dplyr)})

tcga <- read.csv("./data/A1_tcga_coad_merged.csv")
expr_file <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"
expr_raw <- read.table(gzfile(expr_file), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)

extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) vals[i] <- expr_mat[gene, sid]
    else {
      pid <- substr(sid, 1, 12)
      m <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(m) > 0) vals[i] <- expr_mat[gene, m[1]]
    }
  }
  vals
}

tcga$CDH1 <- extract_gene("CDH1", expr_raw, tcga)
tcga$T_grad <- ifelse(grepl("^T1$", tcga$ajcc_pathologic_t), "T1",
                ifelse(grepl("^T2$", tcga$ajcc_pathologic_t), "T2",
                ifelse(grepl("^T3$", tcga$ajcc_pathologic_t), "T3",
                ifelse(grepl("^T4", tcga$ajcc_pathologic_t), "T4", NA))))
tcga <- tcga[!is.na(tcga$T_grad) & tcga$sample_scope == "Tumor" & !is.na(tcga$CDH1), ]

# Dichotomize at cohort median
tcga$CDH1_low <- tcga$CDH1 < median(tcga$CDH1)
cat("=== TCGA CDH1-low fraction by T stage ===\n")
print(tcga %>% group_by(T_grad) %>% summarise(n = n(), CDH1_low_pct = round(100 * mean(CDH1_low), 1)))

t_ord <- as.numeric(factor(tcga$T_grad, levels = c("T1","T2","T3","T4")))
glm1 <- glm(CDH1_low ~ t_ord, family = binomial, data = tcga)
cat(sprintf("Logistic trend (T ordinal -> CDH1-low): OR=%.2f, P=%.4f\n",
    exp(coef(glm1)[2]), summary(glm1)$coefficients[2, 4]))

cat("T1-2 low%:", round(100 * mean(tcga$CDH1_low[tcga$T_grad %in% c("T1","T2")]), 1),
    "| T3-4 low%:", round(100 * mean(tcga$CDH1_low[tcga$T_grad %in% c("T3","T4")]), 1), "\n")
ft <- fisher.test(table(tcga$CDH1_low, tcga$T_grad %in% c("T1","T2")))
cat("Fisher T3-4 vs T1-2 P:", ft$p.value, "\n")

# Also compare with per-T distribution of the actual CDH1 values
cat("\nTCGA CDH1 summary by T:\n")
print(tcga %>% group_by(T_grad) %>% summarise(n = n(), median = round(median(CDH1), 2), mean = round(mean(CDH1), 2), sd = round(sd(CDH1), 2)))
