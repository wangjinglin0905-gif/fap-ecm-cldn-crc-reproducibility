# =============================================================================
# REVIEW_R5: Stromal-content deconvolution proxy for TSR - lymph node risk
# Strategy: (1) use MCP-counter Fibroblast score from AJCR package (validated
# stromal marker score); (2) compute an ESTIMATE-style stromal signature score
# locally from published ESTIMATE stromal genes available in TCGA matrix.
# Test whether a comprehensive stromal score (vs the 6-gene ECM proxy) can
# stratify lymph node metastasis (N+ vs N0).
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(pROC)})

PROJ <- "."
OUT  <- file.path(PROJ, "output", "review_r5")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TCGA_EXPR <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"

# ---- ESTIMATE stromal signature (from Yoshihara 2013, curated core genes
#      available in TCGA) ----
ESTIMATE_STROMAL_CORE <- c(
  # collagens & ECM (stromal core, non-overlapping with FAP13/matrix4)
  "COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COL11A1","COL12A1",
  "COL15A1","COL16A1","COL18A1","COL1A1","COL1A2","COL3A1","COL4A1","COL4A2",
  "FN1","SPARC","DCN","LUM","BGN","POSTN","THY1","PDPN","TAGLN","ACTA2",
  "FAP","MMP2","MMP11","MMP14","TIMP1","TIMP2","TIMP3","LOXL1","LOXL2",
  "LOX","PLOD1","PLOD2","PLOD3","VCAN","FBN1","FBLN1","FBLN2","ELN","MFAP4",
  "IGFBP2","IGFBP3","IGFBP4","IGFBP5","IGFBP7","TGFBI","CTGF","CYR61",
  "PDGFRB","PDGFRA","CXCL12","SDC1","SDC2","SDC4","CD44","VIM","S100A4",
  "ITGA1","ITGA2","ITGA5","ITGB1","ITGB3","ITGB5","EMP3","MGP","CNN1",
  "CALD1","LMOD1","MYLK","TPM2","FLNA","FLNB","ACTN1","ACTG2","DES","VWF"
)
ESTIMATE_STROMAL_CORE <- unique(ESTIMATE_STROMAL_CORE)

zmean <- function(df, genes) {
  avail <- genes[genes %in% colnames(df)]
  if (length(avail) < 5) return(rep(NA, nrow(df)))
  sub <- as.data.frame(lapply(avail, function(g) {
    v <- df[[g]]
    if (all(is.na(v))) return(v)
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }))
  names(sub) <- avail
  rowMeans(sub, na.rm = TRUE)
}

# ---- TCGA ----
cat("=== TCGA-COAD ===\n")
expr_raw <- read.table(gzfile(TCGA_EXPR), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
tcga <- read.csv(file.path(PROJ, "data", "A1_tcga_coad_merged.csv"))
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
for (g in ESTIMATE_STROMAL_CORE) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga <- tcga[tcga$sample_scope == "Tumor", ]
tcga$StromalScore <- zmean(tcga, ESTIMATE_STROMAL_CORE)
tcga$N_pos <- ifelse(grepl("N[12]", tcga$ajcc_pathologic_n), 1,
              ifelse(grepl("N0", tcga$ajcc_pathologic_n), 0, NA))
tcga <- tcga[!is.na(tcga$N_pos), ]
cat("n =", nrow(tcga), "(N0:", sum(tcga$N_pos == 0), "N+:", sum(tcga$N_pos == 1), ")\n")
cat("stromal genes available:", sum(ESTIMATE_STROMAL_CORE %in% colnames(tcga)), "/", length(ESTIMATE_STROMAL_CORE), "\n")

# median split
med <- median(tcga$StromalScore, na.rm = TRUE)
tcga$Stromal_high <- ifelse(tcga$StromalScore >= med, 1, 0)
tb <- table(Stromal_high = tcga$Stromal_high, N_pos = tcga$N_pos)
print(tb)
f <- fisher.test(tb)
or_ <- (tb[2,2]*tb[1,1])/(tb[1,2]*tb[2,1])
cat(sprintf("StromalScore-high vs low: N+ rate %.1f%% vs %.1f%% | OR=%.2f (%.2f-%.2f) P=%.4f\n",
    100*tb[2,2]/sum(tb[2,]), 100*tb[1,2]/sum(tb[1,]),
    or_, f$conf.int[1], f$conf.int[2], f$p.value))
roc_t <- roc(tcga$N_pos, tcga$StromalScore, quiet = TRUE)
cat(sprintf("AUC(StromalScore continuous) = %.3f (%.3f-%.3f)\n", auc(roc_t), ci.auc(roc_t)[1], ci.auc(roc_t)[3]))

# multivariable with T
tcga$T34 <- ifelse(grepl("^T[34]", tcga$ajcc_pathologic_t), 1,
            ifelse(grepl("^T[12]$", tcga$ajcc_pathologic_t), 0, NA))
m <- glm(N_pos ~ Stromal_high + T34, data = tcga, family = binomial)
cat("MV (Stromal_high + T34):\n")
print(summary(m)$coefficients)

# ---- GSE39582 ----
cat("\n=== GSE39582 ===\n")
expr <- read.csv(file.path(PROJ, "data", "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(PROJ, "data", "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}
sm_file <- file.path(PROJ, "data", "GSE39582", "GSE39582_series_matrix.txt.gz")
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
n <- length(sample_names)
gse <- data.frame(sample = sample_names,
                  dataset = if (length(get_field("dataset")) == n) get_field("dataset") else NA,
                  tnm_t = if (length(get_field("tnm.t")) == n) get_field("tnm.t") else NA,
                  tnm_n = if (length(get_field("tnm.n")) == n) get_field("tnm.n") else NA,
                  stringsAsFactors = FALSE)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
gse <- gse[!is.na(gse$dataset) & gse$dataset != "Non Tumoral", ]
gse$StromalScore <- zmean(gse, ESTIMATE_STROMAL_CORE)
gse$N_pos <- ifelse(grepl("N[12]", gse$tnm_n), 1, ifelse(grepl("N0", gse$tnm_n), 0, NA))
gse <- gse[!is.na(gse$N_pos), ]
cat("n =", nrow(gse), "(N0:", sum(gse$N_pos == 0), "N+:", sum(gse$N_pos == 1), ")\n")
cat("stromal genes available:", sum(ESTIMATE_STROMAL_CORE %in% colnames(gse)), "/", length(ESTIMATE_STROMAL_CORE), "\n")

med_g <- median(gse$StromalScore, na.rm = TRUE)
gse$Stromal_high <- ifelse(gse$StromalScore >= med_g, 1, 0)
tb_g <- table(Stromal_high = gse$Stromal_high, N_pos = gse$N_pos)
print(tb_g)
f_g <- fisher.test(tb_g)
or_g <- (tb_g[2,2]*tb_g[1,1])/(tb_g[1,2]*tb_g[2,1])
cat(sprintf("GSE StromalScore-high vs low: N+ rate %.1f%% vs %.1f%% | OR=%.2f (%.2f-%.2f) P=%.4f\n",
    100*tb_g[2,2]/sum(tb_g[2,]), 100*tb_g[1,2]/sum(tb_g[1,]),
    or_g, f_g$conf.int[1], f_g$conf.int[2], f_g$p.value))
roc_g <- roc(gse$N_pos, gse$StromalScore, quiet = TRUE)
cat(sprintf("GSE AUC(StromalScore) = %.3f (%.3f-%.3f)\n", auc(roc_g), ci.auc(roc_g)[1], ci.auc(roc_g)[3]))

gse$T34 <- ifelse(grepl("^T[34]", gse$tnm_t), 1, ifelse(grepl("^T[12]$", gse$tnm_t), 0, NA))
m_g <- glm(N_pos ~ Stromal_high + T34, data = gse, family = binomial)
cat("GSE MV (Stromal_high + T34):\n")
print(summary(m_g)$coefficients)

# Summary
sum_tab <- data.frame(
  analysis = c("TCGA StromalScore-high OR", "TCGA AUC", "TCGA MV Stromal adj T",
               "GSE StromalScore-high OR", "GSE AUC", "GSE MV Stromal adj T"),
  value = c(sprintf("OR=%.2f (%.2f-%.2f) P=%.4f", or_, f$conf.int[1], f$conf.int[2], f$p.value),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_t), ci.auc(roc_t)[1], ci.auc(roc_t)[3]),
            sprintf("OR=%.2f P=%.4f", exp(coef(m)["Stromal_high"]), summary(m)$coefficients["Stromal_high", 4]),
            sprintf("OR=%.2f (%.2f-%.2f) P=%.4f", or_g, f_g$conf.int[1], f_g$conf.int[2], f_g$p.value),
            sprintf("%.3f (%.3f-%.3f)", auc(roc_g), ci.auc(roc_g)[1], ci.auc(roc_g)[3]),
            sprintf("OR=%.2f P=%.4f", exp(coef(m_g)["Stromal_high"]), summary(m_g)$coefficients["Stromal_high", 4])),
  stringsAsFactors = FALSE)
write.csv(sum_tab, file.path(OUT, "review_r5_summary.csv"), row.names = FALSE)
print(sum_tab)
cat("\n=== [DONE] REVIEW_R5 StromalScore TSR proxy ===\n")
