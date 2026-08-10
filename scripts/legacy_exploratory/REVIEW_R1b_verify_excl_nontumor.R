# =============================================================================
# REVIEW_R1b: Re-analysis v2 - exclude GSE39582 non-tumoral samples (19)
# and reproduce v3.4/codex numbers precisely.
# Fixes from R1: (1) exclude Non Tumoral; (2) exclude stage-0 & N/A stage;
# (3) multivariable OS with factor stage (n=557, 190 deaths per codex);
# (4) CLDN: extract all 10 genes from expression matrix.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(survival)
})

PROJ   <- "."
DATA   <- file.path(PROJ, "data")
OUT    <- file.path(PROJ, "output", "review_r1b")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TCGA_EXPR <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"

FAP13    <- c("FAP","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
MATRIX4  <- c("COL1A1","COL1A2","COL3A1","FN1")
RECEPTOR2<- c("SDC4","CD44")
CLDN_CORE<- c("CLDN1","CLDN2","CLDN4")
CLDN10   <- c("CLDN1","CLDN2","CLDN3","CLDN4","CLDN7","CLDN8","CLDN12","CLDN15","CLDN18","CLDN23")

zmean <- function(df, genes) {
  avail <- genes[genes %in% colnames(df)]
  sub <- as.data.frame(lapply(avail, function(g) {
    v <- df[[g]]
    if (all(is.na(v))) return(v)
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }))
  names(sub) <- avail
  rowMeans(sub, na.rm = TRUE)
}

# =============================================================================
# PART 1: TCGA (reuse R1 approach, all 10 CLDN)
# =============================================================================
cat("=== PART 1: TCGA ===\n")
expr_raw <- read.table(gzfile(TCGA_EXPR), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
tcga <- read.csv(file.path(DATA, "A1_tcga_coad_merged.csv"))
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
for (g in unique(c(FAP13, MATRIX4, RECEPTOR2, CLDN10, "FAP"))) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga <- tcga[tcga$sample_scope == "Tumor", ]
tcga$FAP13     <- zmean(tcga, FAP13)
tcga$matrix4   <- zmean(tcga, MATRIX4)
tcga$receptor2 <- zmean(tcga, RECEPTOR2)
tcga$FAP_CAF   <- if ("FAP_CAF" %in% colnames(tcga)) tcga$FAP_CAF else zmean(tcga, c(FAP13, MATRIX4))
cat("TCGA tumors:", nrow(tcga), "\n")
cat(sprintf("TCGA FAP13-matrix4 rho=%.3f | FAP13-receptor2 rho=%.3f P=%.3f\n",
    cor.test(tcga$FAP13, tcga$matrix4, method="spearman")$estimate,
    cor.test(tcga$FAP13, tcga$receptor2, method="spearman")$estimate,
    cor.test(tcga$FAP13, tcga$receptor2, method="spearman")$p.value))

# CLDN all 10, stage-stratified, BH
tcga$stage_grp <- ifelse(grepl("Stage I", tcga$ajcc_pathologic_stage), "Early", "Late")
cld <- data.frame(gene = CLDN10, rho_e = NA, P_e = NA, rho_l = NA, P_l = NA)
for (i in seq_along(CLDN10)) {
  g <- CLDN10[i]
  if (!g %in% colnames(tcga)) next
  e <- tcga[tcga$stage_grp == "Early", ]; l <- tcga[tcga$stage_grp == "Late", ]
  ce <- tryCatch(cor.test(e[[g]], e$FAP_CAF, method = "spearman"), error = function(x) NULL)
  cl <- tryCatch(cor.test(l[[g]], l$FAP_CAF, method = "spearman"), error = function(x) NULL)
  cld$rho_e[i] <- if (!is.null(ce)) ce$estimate else NA
  cld$P_e[i]   <- if (!is.null(ce)) ce$p.value else NA
  cld$rho_l[i] <- if (!is.null(cl)) cl$estimate else NA
  cld$P_l[i]   <- if (!is.null(cl)) cl$p.value else NA
}
cld$FDR_e <- p.adjust(cld$P_e, method = "BH")
cld$FDR_l <- p.adjust(cld$P_l, method = "BH")
cat("\nCLDN BH (TCGA):\n"); print(cld, row.names = FALSE)
cat("n significant after BH:", sum(cld$FDR_e < 0.05, na.rm = TRUE) + sum(cld$FDR_l < 0.05, na.rm = TRUE), "\n")
write.csv(cld, file.path(OUT, "review_CLDN_BH_full.csv"), row.names = FALSE)

# =============================================================================
# PART 2: GSE39582 - EXCLUDE non-tumoral (19) & stage 0 (4) & N/A (2)
# =============================================================================
cat("\n=== PART 2: GSE39582 (tumor-only) ===\n")
expr <- read.csv(file.path(DATA, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(DATA, "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}
sm_file <- file.path(DATA, "GSE39582", "GSE39582_series_matrix.txt.gz")
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
gse <- data.frame(
  sample = sample_names,
  dataset = if (length(get_field("dataset")) == n) get_field("dataset") else NA,
  tnm_t = if (length(get_field("tnm.t")) == n) get_field("tnm.t") else NA,
  tnm_n = if (length(get_field("tnm.n")) == n) get_field("tnm.n") else NA,
  tnm_m = if (length(get_field("tnm.m")) == n) get_field("tnm.m") else NA,
  stage = if (length(get_field("tnm.stage")) == n) get_field("tnm.stage") else NA,
  os_event = as.numeric(if (length(get_field("os.event")) == n) get_field("os.event") else NA),
  os_delay = as.numeric(if (length(get_field("os.delay (months)")) == n) get_field("os.delay (months)") else NA),
  age = as.numeric(if (length(get_field("age.at.diagnosis (year)")) == n) get_field("age.at.diagnosis (year)") else NA),
  Sex = if (length(get_field("Sex")) == n) get_field("Sex") else NA,
  mmr = if (length(get_field("mmr.status")) == n) get_field("mmr.status") else NA,
  stringsAsFactors = FALSE
)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
# EXCLUDE non-tumoral
gse_t <- gse[!is.na(gse$dataset) & gse$dataset != "Non Tumoral", ]
cat("Total:", nrow(gse), "| after excluding Non Tumoral:", nrow(gse_t), "\n")

gse_t$FAP13     <- zmean(gse_t, FAP13)
gse_t$matrix4   <- zmean(gse_t, MATRIX4)
gse_t$receptor2 <- zmean(gse_t, RECEPTOR2)
gse_t$FAP_CAF   <- zmean(gse_t, c(FAP13, MATRIX4))
gse_t$ECM6      <- zmean(gse_t, c(MATRIX4, RECEPTOR2))
gse_t$CLDN_core <- zmean(gse_t, CLDN_CORE)

c1 <- cor.test(gse_t$FAP13, gse_t$matrix4, method = "spearman")
c2 <- cor.test(gse_t$FAP13, gse_t$receptor2, method = "spearman")
cat(sprintf("GSE tumor-only FAP13-matrix4: rho=%.3f P=%.4f\n", c1$estimate, c1$p.value))
cat(sprintf("GSE tumor-only FAP13-receptor2: rho=%.3f P=%.4f\n", c2$estimate, c2$p.value))

# T-stratified
gse_t$T_grp <- ifelse(grepl("^T[12]$", gse_t$tnm_t), "T1-2", ifelse(grepl("^T[34]", gse_t$tnm_t), "T3-4", NA))
gs <- gse_t[!is.na(gse_t$T_grp), ]
ct12 <- cor.test(gs$FAP13[gs$T_grp == "T1-2"], gs$matrix4[gs$T_grp == "T1-2"], method = "spearman")
ct34 <- cor.test(gs$FAP13[gs$T_grp == "T3-4"], gs$matrix4[gs$T_grp == "T3-4"], method = "spearman")
cat(sprintf("GSE FAP13-matrix4 T1-2 rho=%.3f (n=%d) | T3-4 rho=%.3f (n=%d)\n",
    ct12$estimate, sum(gs$T_grp == "T1-2"), ct34$estimate, sum(gs$T_grp == "T3-4")))

# Multivariable OS: n=557, factor stage (codex claims n=557, 190 deaths)
gse_os <- gse_t[!is.na(gse_t$os_delay) & gse_t$os_delay > 0 & !is.na(gse_t$os_event) &
                !is.na(gse_t$FAP13) & !is.na(gse_t$age) & !is.na(gse_t$Sex) &
                !is.na(gse_t$stage) & gse_t$stage %in% c("1","2","3","4"), ]
gse_os$stageF <- factor(gse_os$stage, levels = c("1","2","3","4"))
cox_g <- coxph(Surv(os_delay, os_event) ~ FAP13 + age + Sex + stageF, data = gse_os)
cat(sprintf("\nGSE OS model (tumor-only, stage factor): n=%d, events=%d\n", nrow(gse_os), sum(gse_os$os_event)))
cat(sprintf("FAP13 HR=%.2f (%.2f-%.2f) P=%.3f | PH global P=%.3f\n",
    exp(coef(cox_g)["FAP13"]), exp(confint(cox_g)["FAP13",1]), exp(confint(cox_g)["FAP13",2]),
    summary(cox_g)$coefficients["FAP13",5], cox.zph(cox_g)$table[1,3]))

# MMR
gsem <- gse_t[!is.na(gse_t$mmr) & gse_t$mmr %in% c("dMMR","pMMR"), ]
cd <- cor.test(gsem$FAP_CAF[gsem$mmr=="dMMR"], gsem$ECM6[gsem$mmr=="dMMR"], method="spearman")
cp <- cor.test(gsem$FAP_CAF[gsem$mmr=="pMMR"], gsem$ECM6[gsem$mmr=="pMMR"], method="spearman")
cat(sprintf("\nMMR: dMMR rho=%.3f (n=%d) | pMMR rho=%.3f (n=%d)\n",
    cd$estimate, sum(gsem$mmr=="dMMR"), cp$estimate, sum(gsem$mmr=="pMMR")))

# Summary CSV
sum_tab <- data.frame(
  analysis = c("TCGA FAP13-matrix4", "TCGA FAP13-receptor2",
               "GSE tumor-only FAP13-matrix4", "GSE tumor-only FAP13-receptor2",
               "GSE FAP13-matrix4 T1-2", "GSE FAP13-matrix4 T3-4",
               "GSE OS n/events", "GSE OS FAP13 HR",
               "MMR dMMR/pMMR rho",
               "CLDN n_sig after BH"),
  value = c(sprintf("rho=%.3f", cor.test(tcga$FAP13, tcga$matrix4, method="spearman")$estimate),
            sprintf("rho=%.3f P=%.3f", cor.test(tcga$FAP13, tcga$receptor2, method="spearman")$estimate, cor.test(tcga$FAP13, tcga$receptor2, method="spearman")$p.value),
            sprintf("rho=%.3f P=%.4f", c1$estimate, c1$p.value),
            sprintf("rho=%.3f P=%.4f", c2$estimate, c2$p.value),
            sprintf("rho=%.3f n=%d", ct12$estimate, sum(gs$T_grp=="T1-2")),
            sprintf("rho=%.3f n=%d", ct34$estimate, sum(gs$T_grp=="T3-4")),
            sprintf("n=%d events=%d", nrow(gse_os), sum(gse_os$os_event)),
            sprintf("HR=%.2f P=%.3f", exp(coef(cox_g)["FAP13"]), summary(cox_g)$coefficients["FAP13",5]),
            sprintf("%.3f / %.3f", cd$estimate, cp$estimate),
            sprintf("%d", sum(cld$FDR_e < 0.05, na.rm=TRUE) + sum(cld$FDR_l < 0.05, na.rm=TRUE))),
  stringsAsFactors = FALSE
)
write.csv(sum_tab, file.path(OUT, "review_r1b_summary.csv"), row.names = FALSE)
print(sum_tab)
cat("\n=== [DONE] REVIEW_R1b ===\n")
