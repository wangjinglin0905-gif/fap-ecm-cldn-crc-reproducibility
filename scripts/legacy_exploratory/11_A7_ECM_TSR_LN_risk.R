# =============================================================================
# A7: ECM score as transcriptomic proxy for pathological TSR (<1 vs >=1)
# Thesis finding: TSR<1 (stroma-rich) + FAP IHC are LN-metastasis risk factors
# (FAP OR=7.836). Here we use ECM_SDC4_CD44 z-mean (stromal program) dichotomized
# at median as a TSR proxy, and test LN-risk stratification in TCGA & GSE39582.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(pROC)
})

PROJ_ROOT <- "."
FIG_DIR   <- file.path(PROJ_ROOT, "output", "figures")
TAB_DIR   <- file.path(PROJ_ROOT, "output", "tables")
DATA_DIR  <- file.path(PROJ_ROOT, "data")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

ECM_SDC4_CD44 <- c("COL1A1","COL1A2","COL3A1","FN1","SDC4","CD44")
CLDN_CORE     <- c("CLDN1","CLDN2","CLDN4")

compute_zmean <- function(df, genes) {
  avail <- genes[genes %in% colnames(df)]
  sub <- as.data.frame(lapply(avail, function(g) {
    v <- df[[g]]
    if (all(is.na(v))) return(v)
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }))
  names(sub) <- avail
  rowMeans(sub, na.rm = TRUE)
}

cat("=== A7: ECM score as TSR proxy for LN-risk stratification ===\n")

# =============================================================================
# PART A: TCGA-COAD
# =============================================================================
cat("\n--- TCGA-COAD ---\n")
tcga <- read.csv(file.path(DATA_DIR, "A1_tcga_coad_merged.csv"))
tcga$N_pos <- ifelse(grepl("N[12]", tcga$ajcc_pathologic_n), 1,
              ifelse(grepl("N0", tcga$ajcc_pathologic_n), 0, NA))
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
      matches <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
      if (length(matches) > 0) vals[i] <- expr_mat[gene, matches[1]]
    }
  }
  vals
}
for (g in c(ECM_SDC4_CD44, CLDN_CORE)) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga$ECM_score <- compute_zmean(tcga, ECM_SDC4_CD44)
tcga$CLDN_score <- compute_zmean(tcga, CLDN_CORE)
tcga <- tcga[!is.na(tcga$N_pos) & tcga$sample_scope == "Tumor", ]
cat("n =", nrow(tcga), "(N0:", sum(tcga$N_pos == 0), "| N+:", sum(tcga$N_pos == 1), ")\n")

# Dichotomize ECM score at median (TSR proxy: high = stroma-rich = "TSR<1-like")
med_t <- median(tcga$ECM_score, na.rm = TRUE)
tcga$ECM_high <- ifelse(tcga$ECM_score >= med_t, 1, 0)
cat("ECM median =", round(med_t, 3), "\n")
cat("ECM-high:", sum(tcga$ECM_high == 1), "| ECM-low:", sum(tcga$ECM_high == 0), "\n")

# 2x2 table: ECM-high x N+
tb_t <- table(ECM_high = tcga$ECM_high, N_pos = tcga$N_pos)
print(tb_t)
f_t <- fisher.test(tb_t)
or_t <- (tb_t[2,2] * tb_t[1,1]) / (tb_t[1,2] * tb_t[2,1])  # OR for N+ in ECM-high
cat(sprintf("ECM-high vs low: N+ rate %.1f%% vs %.1f%% | OR=%.2f (%.2f-%.2f) P=%.4f\n",
    100*tb_t[2,2]/sum(tb_t[2,]), 100*tb_t[1,2]/sum(tb_t[1,]),
    or_t, f_t$conf.int[1], f_t$conf.int[2], f_t$p.value))

# Same for FAP_CAF (thesis: FAP OR=7.84)
med_f <- median(tcga$FAP_CAF, na.rm = TRUE)
tcga$FAP_high <- ifelse(tcga$FAP_CAF >= med_f, 1, 0)
tb_f <- table(FAP_high = tcga$FAP_high, N_pos = tcga$N_pos)
print(tb_f)
f_f <- fisher.test(tb_f)
or_f <- (tb_f[2,2] * tb_f[1,1]) / (tb_f[1,2] * tb_f[2,1])
cat(sprintf("FAP-high vs low: N+ rate %.1f%% vs %.1f%% | OR=%.2f (%.2f-%.2f) P=%.4f\n",
    100*tb_f[2,2]/sum(tb_f[2,]), 100*tb_f[1,2]/sum(tb_f[1,]),
    or_f, f_f$conf.int[1], f_f$conf.int[2], f_f$p.value))

# Multivariable: adjust for T stage
tcga$T_group <- ifelse(grepl("^T[12]$", tcga$ajcc_pathologic_t), "T1-2",
                ifelse(grepl("^T[34]", tcga$ajcc_pathologic_t), "T3-4", NA))
tcga$T34 <- ifelse(tcga$T_group == "T3-4", 1, ifelse(tcga$T_group == "T1-2", 0, NA))
log_mv_t <- glm(N_pos ~ ECM_high + T34, data = tcga, family = binomial)
cat("Multivariable (ECM-high + T3-4) TCGA:\n")
print(summary(log_mv_t)$coefficients)

# AUC of continuous ECM score
roc_t <- roc(tcga$N_pos, tcga$ECM_score, quiet = TRUE)
cat(sprintf("AUC(ECM continuous) = %.3f (%.3f-%.3f)\n", auc(roc_t), ci.auc(roc_t)[1], ci.auc(roc_t)[3]))

# =============================================================================
# PART B: GSE39582
# =============================================================================
cat("\n--- GSE39582 ---\n")
expr <- read.csv(file.path(DATA_DIR, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(DATA_DIR, "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}
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
tnm_t <- get_field("tnm.t"); tnm_n <- get_field("tnm.n")
n <- length(sample_names)
gse <- data.frame(sample = sample_names,
                  tnm_t = if (length(tnm_t) == n) tnm_t else NA,
                  tnm_n = if (length(tnm_n) == n) tnm_n else NA,
                  stringsAsFactors = FALSE)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
gse$FAP_CAF <- compute_zmean(gse, c("FAP","COL1A1","COL1A2","COL3A1","FN1","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A"))
gse$ECM_score <- compute_zmean(gse, ECM_SDC4_CD44)
gse$N_pos <- ifelse(grepl("N[12]", gse$tnm_n), 1, ifelse(grepl("N0", gse$tnm_n), 0, NA))
gse <- gse[!is.na(gse$N_pos), ]
cat("n =", nrow(gse), "(N0:", sum(gse$N_pos == 0), "| N+:", sum(gse$N_pos == 1), ")\n")

med_g <- median(gse$ECM_score, na.rm = TRUE)
gse$ECM_high <- ifelse(gse$ECM_score >= med_g, 1, 0)
tb_g <- table(ECM_high = gse$ECM_high, N_pos = gse$N_pos)
print(tb_g)
f_g <- fisher.test(tb_g)
or_g <- (tb_g[2,2] * tb_g[1,1]) / (tb_g[1,2] * tb_g[2,1])
cat(sprintf("ECM-high vs low: N+ rate %.1f%% vs %.1f%% | OR=%.2f (%.2f-%.2f) P=%.4f\n",
    100*tb_g[2,2]/sum(tb_g[2,]), 100*tb_g[1,2]/sum(tb_g[1,]),
    or_g, f_g$conf.int[1], f_g$conf.int[2], f_g$p.value))

med_fg <- median(gse$FAP_CAF, na.rm = TRUE)
gse$FAP_high <- ifelse(gse$FAP_CAF >= med_fg, 1, 0)
tb_fg <- table(FAP_high = gse$FAP_high, N_pos = gse$N_pos)
print(tb_fg)
f_fg <- fisher.test(tb_fg)
or_fg <- (tb_fg[2,2] * tb_fg[1,1]) / (tb_fg[1,2] * tb_fg[2,1])
cat(sprintf("FAP-high vs low: N+ rate %.1f%% vs %.1f%% | OR=%.2f (%.2f-%.2f) P=%.4f\n",
    100*tb_fg[2,2]/sum(tb_fg[2,]), 100*tb_fg[1,2]/sum(tb_fg[1,]),
    or_fg, f_fg$conf.int[1], f_fg$conf.int[2], f_fg$p.value))

gse$T_group <- ifelse(grepl("^T[12]$", gse$tnm_t), "T1-2",
               ifelse(grepl("^T[34]", gse$tnm_t), "T3-4", NA))
gse$T34 <- ifelse(gse$T_group == "T3-4", 1, ifelse(gse$T_group == "T1-2", 0, NA))
log_mv_g <- glm(N_pos ~ ECM_high + T34, data = gse, family = binomial)
cat("Multivariable (ECM-high + T3-4) GSE39582:\n")
print(summary(log_mv_g)$coefficients)

roc_g <- roc(gse$N_pos, gse$ECM_score, quiet = TRUE)
cat(sprintf("AUC(ECM continuous) = %.3f (%.3f-%.3f)\n", auc(roc_g), ci.auc(roc_g)[1], ci.auc(roc_g)[3]))

# =============================================================================
# PART C: Combined plot (ECM-high vs low N+ rate, both cohorts)
# =============================================================================
comb <- data.frame(
  cohort = c("TCGA-COAD", "TCGA-COAD", "GSE39582", "GSE39582"),
  ECM_group = rep(c("ECM-low", "ECM-high"), 2),
  N_plus_rate = c(100*tb_t[1,2]/sum(tb_t[1,]), 100*tb_t[2,2]/sum(tb_t[2,]),
                  100*tb_g[1,2]/sum(tb_g[1,]), 100*tb_g[2,2]/sum(tb_g[2,])),
  n_Nplus = c(tb_t[1,2], tb_t[2,2], tb_g[1,2], tb_g[2,2]),
  n_total = c(sum(tb_t[1,]), sum(tb_t[2,]), sum(tb_g[1,]), sum(tb_g[2,])),
  stringsAsFactors = FALSE
)
comb$ECM_group <- factor(comb$ECM_group, levels = c("ECM-low", "ECM-high"))
p1 <- ggplot(comb, aes(x = cohort, y = N_plus_rate, fill = ECM_group)) +
  geom_bar(stat = "identity", position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", N_plus_rate, n_Nplus, n_total)),
            position = position_dodge(0.75), vjust = -0.25, size = 3.2) +
  scale_fill_manual(values = c("ECM-low" = "#4DBBD5", "ECM-high" = "#E64B35")) +
  labs(x = "", y = "Lymph node metastasis rate (%)",
       title = "LN+ Rate by ECM Score (transcriptomic TSR proxy)",
       subtitle = "ECM-high mimics pathological TSR<1 (stroma-rich); thesis FAP OR=7.84") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  ylim(0, max(comb$N_plus_rate) * 1.25)
ggsave(file.path(FIG_DIR, "A7_ECM_TSR_LN_rate.png"), p1, width = 7, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "A7_ECM_TSR_LN_rate.tiff"), p1, width = 7, height = 6, dpi = 600, compression = "lzw")

# Summary table
sum_tab <- data.frame(
  analysis = c("TCGA ECM-high vs low OR",
               "TCGA FAP-high vs low OR",
               "TCGA AUC(ECM continuous)",
               "TCGA MV ECM-high adj T",
               "GSE39582 ECM-high vs low OR",
               "GSE39582 FAP-high vs low OR",
               "GSE39582 AUC(ECM continuous)",
               "GSE39582 MV ECM-high adj T"),
  value = c(sprintf("OR=%.2f (%.2f-%.2f) P=%.4f", or_t, f_t$conf.int[1], f_t$conf.int[2], f_t$p.value),
            sprintf("OR=%.2f (%.2f-%.2f) P=%.4f", or_f, f_f$conf.int[1], f_f$conf.int[2], f_f$p.value),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_t), ci.auc(roc_t)[1], ci.auc(roc_t)[3]),
            sprintf("OR=%.2f P=%.4f", exp(coef(log_mv_t)["ECM_high"]), summary(log_mv_t)$coefficients["ECM_high", 4]),
            sprintf("OR=%.2f (%.2f-%.2f) P=%.4f", or_g, f_g$conf.int[1], f_g$conf.int[2], f_g$p.value),
            sprintf("OR=%.2f (%.2f-%.2f) P=%.4f", or_fg, f_fg$conf.int[1], f_fg$conf.int[2], f_fg$p.value),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_g), ci.auc(roc_g)[1], ci.auc(roc_g)[3]),
            sprintf("OR=%.2f P=%.4f", exp(coef(log_mv_g)["ECM_high"]), summary(log_mv_g)$coefficients["ECM_high", 4])),
  stringsAsFactors = FALSE
)
write.csv(sum_tab, file.path(TAB_DIR, "A7_ECM_TSR_LN_summary.csv"), row.names = FALSE)
cat("\n=== [DONE] A7 ECM-as-TSR LN-risk complete ===\n")
print(sum_tab)
