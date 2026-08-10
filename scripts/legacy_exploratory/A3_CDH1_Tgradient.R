# =============================================================================
# A3: CDH1 (E-cadherin) expression gradient across T stages
# Hypothesis from thesis [9b]: E-cadherin IHC negativity increases with invasion
# depth (HGIN+pT1a 6% -> pT1b 8% -> pT2 36% -> pT3 42% -> pT4 56%).
# Transcriptomic validation: CDH1 mRNA should decline from T1 to T4 in both
# TCGA-COAD and GSE39582.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(ggpubr)
})

PROJ_ROOT <- "."
FIG_DIR   <- file.path(PROJ_ROOT, "output", "figures")
TAB_DIR   <- file.path(PROJ_ROOT, "output", "tables")
DATA_DIR  <- file.path(PROJ_ROOT, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

CDH1_PROBES <- c("201130_s_at", "201131_s_at")

# =============================================================================
# PART 1: TCGA-COAD CDH1 by T stage
# =============================================================================
cat("=== PART 1: TCGA-COAD CDH1 gradient ===\n")
expr_file <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"
expr_raw <- read.table(gzfile(expr_file), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
cat("Expression:", nrow(expr_raw), "genes x", ncol(expr_raw), "samples\n")

tcga <- read.csv(file.path(DATA_DIR, "A1_tcga_coad_merged.csv"))
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) vals[i] <- expr_mat[gene, sid]
    else {
      pid <- substr(sid, 1, 12)
      matches <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(matches) > 0) vals[i] <- expr_mat[gene, matches[1]]
    }
  }
  vals
}
tcga$CDH1 <- extract_gene("CDH1", expr_raw, tcga)

# T grouping
tcga$T_grad <- ifelse(grepl("^T1$", tcga$ajcc_pathologic_t), "T1",
                ifelse(grepl("^T2$", tcga$ajcc_pathologic_t), "T2",
                ifelse(grepl("^T3$", tcga$ajcc_pathologic_t), "T3",
                ifelse(grepl("^T4", tcga$ajcc_pathologic_t), "T4", NA))))
tcga <- tcga[!is.na(tcga$T_grad) & tcga$sample_scope == "Tumor" & !is.na(tcga$CDH1), ]
cat("TCGA tumors with CDH1:", nrow(tcga), "\n")
print(table(tcga$T_grad))

# Kruskal-Wallis across T1-T4
kw <- kruskal.test(CDH1 ~ factor(T_grad, levels = c("T1","T2","T3","T4")), data = tcga)
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw$p.value))
# Trend: Spearman CDH1 vs T ordinal
t_ord <- as.numeric(factor(tcga$T_grad, levels = c("T1","T2","T3","T4")))
tr <- cor.test(tcga$CDH1, t_ord, method = "spearman")
cat(sprintf("CDH1 vs T-ordinal Spearman rho = %.3f, P = %.4f\n", tr$estimate, tr$p.value))

# T1-2 vs T3-4
t12 <- tcga$CDH1[tcga$T_grad %in% c("T1","T2")]
t34 <- tcga$CDH1[tcga$T_grad %in% c("T3","T4")]
wt <- wilcox.test(t34, t12)
cat(sprintf("CDH1 T3-4 vs T1-2: P = %.4f (median %.3f vs %.3f)\n", wt$p.value,
    median(t34), median(t12)))

# Box plot TCGA
tcga$T_grad_f <- factor(tcga$T_grad, levels = c("T1","T2","T3","T4"))
p1 <- ggplot(tcga, aes(x = T_grad_f, y = CDH1, fill = T_grad_f)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.35) +
  scale_fill_manual(values = c("T1" = "#4DBBD5", "T2" = "#00A087", "T3" = "#E64B35", "T4" = "#B2182B")) +
  labs(x = "T stage (TCGA-COAD)", y = "CDH1 expression (log-normalized)",
       title = sprintf("CDH1 (E-cadherin) Across T Stages - TCGA-COAD (KW P=%.4f, rho=%.3f)",
                       kw$p.value, tr$estimate)) +
  theme_classic(base_size = 13) + theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("T1","T4")), method = "wilcox.test", size = 3)
ggsave(file.path(FIG_DIR, "A3_CDH1_TCGA_Tgradient.png"), p1, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A3_CDH1_TCGA_Tgradient.tiff"), p1, width = 7.5, height = 6, dpi = 600, compression = "lzw")

# =============================================================================
# PART 2: GSE39582 CDH1 by T stage
# =============================================================================
cat("\n=== PART 2: GSE39582 CDH1 gradient ===\n")
# GSE39582_expr.csv only contains probes for the 22 previously targeted genes;
# CDH1 was not among them. Extract CDH1 probes directly from the series matrix.
sm_file <- file.path(DATA_DIR, "GSE39582", "GSE39582_series_matrix.txt.gz")
con <- gzfile(sm_file, "rt")
expr_rows <- list(); sample_names <- NULL; in_table <- FALSE
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_geo_accession")) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar("!Sample_geo_accession") + 1, nchar(line)), "\t")[[1]])
    sample_names <- sample_names[sample_names != ""]
  } else if (grepl("series_matrix_table_begin", line)) {
    in_table <- TRUE
  } else if (in_table && startsWith(line, "!series_matrix_table_end")) {
    in_table <- FALSE
  } else if (in_table && !startsWith(line, "!Series_")) {
    # expression row: first col = probe ID (quoted), rest = sample values
    cells <- strsplit(line, "\t")[[1]]
    probe_id <- gsub('"', '', trimws(cells[1]))
    if (probe_id %in% CDH1_PROBES) {
      vals <- suppressWarnings(as.numeric(cells[-1]))
      expr_rows[[probe_id]] <- vals
    }
  }
}
close(con)
cat("CDH1 probes extracted:", names(expr_rows), "\n")
if (length(expr_rows) > 0) {
  n <- length(sample_names)
  mat <- do.call(rbind, lapply(expr_rows, function(v) {
    if (length(v) == n) v else rep(NA, n)
  }))
  gse_cdh1 <- colMeans(mat, na.rm = TRUE)
} else {
  gse_cdh1 <- rep(NA, length(sample_names))
}

# Re-parse tnm.t from series matrix (same logic as A1T2)
sm_file <- file.path(DATA_DIR, "GSE39582", "GSE39582_series_matrix.txt.gz")
con <- gzfile(sm_file, "rt")
field_rows <- list(); sample_names <- NULL
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_characteristics_ch1")) {
    cells <- strsplit(substr(line, nchar("!Sample_characteristics_ch1") + 1, nchar(line)), "\t")[[1]]
    cells <- gsub('"', '', trimws(cells))
    field_rows[[length(field_rows) + 1]] <- cells
  } else if (startsWith(line, "!Sample_geo_accession")) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar("!Sample_geo_accession") + 1, nchar(line)), "\t")[[1]])
    sample_names <- sample_names[sample_names != ""]
  }
}
close(con)
field_names <- sapply(field_rows, function(r) {
  hit <- which(grepl(":", r))[1]
  if (!is.na(hit)) sub(":.*", "", r[hit]) else ""
})
get_field <- function(fname) {
  idx <- which(field_names == fname)
  if (length(idx) == 0) return(NULL)
  vals <- field_rows[[idx[1]]]
  vals <- vals[vals != ""]
  vals <- sub(paste0(fname, ":"), "", vals, fixed = TRUE)
  trimws(vals)
}
tnm_t <- get_field("tnm.t")
n <- length(sample_names)
gse <- data.frame(sample = sample_names,
                  tnm_t = if (length(tnm_t) == n) tnm_t else NA,
                  CDH1 = gse_cdh1, stringsAsFactors = FALSE)
gse$T_grad <- ifelse(grepl("^T1$", gse$tnm_t), "T1",
               ifelse(grepl("^T2$", gse$tnm_t), "T2",
               ifelse(grepl("^T3$", gse$tnm_t), "T3",
               ifelse(grepl("^T4", gse$tnm_t), "T4", NA))))
gse <- gse[!is.na(gse$T_grad) & !is.na(gse$CDH1), ]
cat("GSE39582 samples with CDH1:", nrow(gse), "\n")
print(table(gse$T_grad))

kw2 <- kruskal.test(CDH1 ~ factor(T_grad, levels = c("T1","T2","T3","T4")), data = gse)
cat(sprintf("Kruskal-Wallis P = %.4f\n", kw2$p.value))
t_ord2 <- as.numeric(factor(gse$T_grad, levels = c("T1","T2","T3","T4")))
tr2 <- cor.test(gse$CDH1, t_ord2, method = "spearman")
cat(sprintf("CDH1 vs T-ordinal Spearman rho = %.3f, P = %.4f\n", tr2$estimate, tr2$p.value))
g12 <- gse$CDH1[gse$T_grad %in% c("T1","T2")]
g34 <- gse$CDH1[gse$T_grad %in% c("T3","T4")]
wt2 <- wilcox.test(g34, g12)
cat(sprintf("CDH1 T3-4 vs T1-2: P = %.4f (median %.3f vs %.3f)\n", wt2$p.value,
    median(g34), median(g12)))

gse$T_grad_f <- factor(gse$T_grad, levels = c("T1","T2","T3","T4"))
p2 <- ggplot(gse, aes(x = T_grad_f, y = CDH1, fill = T_grad_f)) +
  geom_boxplot(outlier.size = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
  scale_fill_manual(values = c("T1" = "#4DBBD5", "T2" = "#00A087", "T3" = "#E64B35", "T4" = "#B2182B")) +
  labs(x = "T stage (GSE39582)", y = "CDH1 expression (microarray)",
       title = sprintf("CDH1 (E-cadherin) Across T Stages - GSE39582 (KW P=%.4f, rho=%.3f)",
                       kw2$p.value, tr2$estimate)) +
  theme_classic(base_size = 13) + theme(legend.position = "none") +
  stat_compare_means(comparisons = list(c("T1","T4")), method = "wilcox.test", size = 3)
ggsave(file.path(FIG_DIR, "A3_CDH1_GSE39582_Tgradient.png"), p2, width = 7.5, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A3_CDH1_GSE39582_Tgradient.tiff"), p2, width = 7.5, height = 6, dpi = 600, compression = "lzw")

# =============================================================================
# PART 3: Summary + thesis cross-check
# =============================================================================
cat("\n=== PART 3: Summary ===\n")
# Thesis E-cadherin negativity gradient (IHC, from thesis Table 6 [9b]):
# HGIN+pT1a 2/33 (6%) -> pT1b 1/13 (8%) -> pT2 4/11 (36%) -> pT3 28/66 (42%) -> pT4 14/25 (56%)
thesis_ecad_neg <- c(T1 = 6, T2 = 36, T3 = 42, T4 = 56)  # % negative by IHC
# Transcriptomic: mean CDH1 by T (z-score within cohort for comparability)
tcga_sum <- tcga %>% group_by(T_grad) %>% summarise(mean_CDH1 = mean(CDH1, na.rm = TRUE), .groups = "drop")
gse_sum <- gse %>% group_by(T_grad) %>% summarise(mean_CDH1 = mean(CDH1, na.rm = TRUE), .groups = "drop")
cat("TCGA mean CDH1 by T:\n"); print(tcga_sum)
cat("\nGSE39582 mean CDH1 by T:\n"); print(gse_sum)
cat("\nThesis E-cadherin IHC negativity % (from [9b]):\n"); print(thesis_ecad_neg)

sum_tab <- data.frame(
  cohort = c("TCGA-COAD", "GSE39582"),
  KW_P = c(signif(kw$p.value, 3), signif(kw2$p.value, 3)),
  rho_Tordinal = c(round(tr$estimate, 3), round(tr2$estimate, 3)),
  rho_P = c(signif(tr$p.value, 3), signif(tr2$p.value, 3)),
  T34_vs_T12_P = c(signif(wt$p.value, 3), signif(wt2$p.value, 3)),
  median_T12 = c(round(median(t12), 3), round(median(g12), 3)),
  median_T34 = c(round(median(t34), 3), round(median(g34), 3)),
  stringsAsFactors = FALSE)
write.csv(sum_tab, file.path(TAB_DIR, "A3_CDH1_summary.csv"), row.names = FALSE)
print(sum_tab)
cat("\n=== [DONE] A3 CDH1 gradient complete ===\n")
