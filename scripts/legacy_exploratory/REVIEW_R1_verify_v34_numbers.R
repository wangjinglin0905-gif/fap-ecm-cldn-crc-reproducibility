# =============================================================================
# REVIEW_R1: Fourth-round independent re-analysis of core numbers cited in v3.4
# Purpose: verify (1) non-overlapping score correlations FAP13/matrix4/receptor2
# in TCGA-COAD & GSE39582; (2) ordinal T-category trends with patient-level
# bootstrap + BH correction; (3) GSE39582 multivariable OS (FAP13 HR);
# (4) CLDN member BH-corrected results; (5) MMR-stratified correlations.
# R 4.6.1 | outputs go to output/review_r1/
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(survival)
})

PROJ   <- "."
DATA   <- file.path(PROJ, "data")
OUT    <- file.path(PROJ, "output", "review_r1")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

TCGA_EXPR <- "data/TCGA/TCGA_COADREAD_expression.txt.gz"

# ---- Score definitions (v3.4 primary, non-overlapping) ----
FAP13    <- c("FAP","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
MATRIX4  <- c("COL1A1","COL1A2","COL3A1","FN1")
RECEPTOR2<- c("SDC4","CD44")
CLDN_CORE<- c("CLDN1","CLDN2","CLDN4")

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

res <- list()

# =============================================================================
# PART 1: TCGA-COAD non-overlapping score correlations
# =============================================================================
cat("=== PART 1: TCGA-COAD non-overlapping scores ===\n")
expr_raw <- read.table(gzfile(TCGA_EXPR), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
tcga <- read.csv(file.path(DATA, "A1_tcga_coad_merged.csv"))
# map samples to expression columns
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
for (g in unique(c(FAP13, MATRIX4, RECEPTOR2, CLDN_CORE))) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
tcga <- tcga[tcga$sample_scope == "Tumor", ]
tcga$FAP13     <- zmean(tcga, FAP13)
tcga$matrix4   <- zmean(tcga, MATRIX4)
tcga$receptor2 <- zmean(tcga, RECEPTOR2)
tcga$FAP_CAF   <- if ("FAP_CAF" %in% colnames(tcga)) tcga$FAP_CAF else zmean(tcga, c(FAP13, MATRIX4))
tcga$ECM6      <- zmean(tcga, c(MATRIX4, RECEPTOR2))
cat("TCGA tumors:", nrow(tcga), "\n")

c_f13m4 <- cor.test(tcga$FAP13, tcga$matrix4, method = "spearman")
c_f13r2 <- cor.test(tcga$FAP13, tcga$receptor2, method = "spearman")
c_ecm6  <- cor.test(tcga$FAP_CAF, tcga$ECM6, method = "spearman")
cat(sprintf("TCGA FAP13-matrix4: rho=%.3f P=%.4f\n", c_f13m4$estimate, c_f13m4$p.value))
cat(sprintf("TCGA FAP13-receptor2: rho=%.3f P=%.4f\n", c_f13r2$estimate, c_f13r2$p.value))
cat(sprintf("TCGA FAP_CAF(17g)-ECM6(overlap): rho=%.3f P=%.4f\n", c_ecm6$estimate, c_ecm6$p.value))
res$tcga_f13m4 <- c(rho = unname(c_f13m4$estimate), P = c_f13m4$p.value)
res$tcga_f13r2 <- c(rho = unname(c_f13r2$estimate), P = c_f13r2$p.value)
res$tcga_overlap <- c(rho = unname(c_ecm6$estimate), P = c_ecm6$p.value)

# Stage-stratified (Stage I vs II-IV)
tcga$stage_grp <- ifelse(tcga$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB"), "I", "II-IV")
for (sg in c("I", "II-IV")) {
  sub <- tcga[tcga$stage_grp == sg, ]
  if (nrow(sub) >= 5) {
    cc <- cor.test(sub$FAP13, sub$matrix4, method = "spearman")
    cat(sprintf("TCGA FAP13-matrix4 Stage %s: rho=%.3f P=%.4f (n=%d)\n", sg, cc$estimate, cc$p.value, nrow(sub)))
  }
}

# =============================================================================
# PART 2: Ordinal T-category trends (TCGA, patient-level bootstrap + BH)
# =============================================================================
cat("\n=== PART 2: TCGA ordinal T-category trends ===\n")
tnum <- ifelse(grepl("^T1$", tcga$ajcc_pathologic_t), 1,
        ifelse(grepl("^T2$", tcga$ajcc_pathologic_t), 2,
        ifelse(grepl("^T3$", tcga$ajcc_pathologic_t), 3,
        ifelse(grepl("^T4", tcga$ajcc_pathologic_t), 4, NA))))
tcga$T_ord <- tnum
tt <- tcga[!is.na(tcga$T_ord), ]
cat("T counts:", table(tt$T_ord), "\n")

trend_test <- function(score, Tord) {
  cor.test(score, Tord, method = "spearman")$estimate
}
scores_t <- list(FAP_expr = tt$FAP, FAP13 = tt$FAP13, matrix4 = tt$matrix4, receptor2 = tt$receptor2)
trends <- sapply(scores_t, function(s) cor.test(s, tt$T_ord, method = "spearman")$estimate)
cat("Point trends:\n"); print(round(trends, 3))

# Patient-level bootstrap CI (5000 resamples)
set.seed(42)
n <- nrow(tt)
boot <- function(score) {
  replicate(5000, {
    idx <- sample(seq_len(n), n, replace = TRUE)
    cor.test(score[idx], tt$T_ord[idx], method = "spearman")$estimate
  })
}
boot_ci <- function(b) quantile(b, c(0.025, 0.975))
boot_fap <- boot(tt$FAP); boot_f13 <- boot(tt$FAP13)
boot_m4  <- boot(tt$matrix4); boot_r2 <- boot(tt$receptor2)
ci_fap <- boot_ci(boot_fap); ci_f13 <- boot_ci(boot_f13)
ci_m4  <- boot_ci(boot_m4);  ci_r2  <- boot_ci(boot_r2)
# BH across 4 prespecified tests: use bootstrap p as approximate (two-sided via z)
p_trend <- sapply(scores_t, function(s) cor.test(s, tt$T_ord, method = "spearman")$p.value)
fdr <- p.adjust(p_trend, method = "BH")
cat("Trend table (v3.4 claims FAP 0.336/F13 0.290/m4 0.280/r2 0.001; FDR 0.008/0.016/0.016/0.995):\n")
tab2 <- data.frame(score = names(trends), rho = round(trends, 3),
                   boot_low = round(c(ci_fap[1], ci_f13[1], ci_m4[1], ci_r2[1]), 3),
                   boot_high = round(c(ci_fap[2], ci_f13[2], ci_m4[2], ci_r2[2]), 3),
                   P = signif(p_trend, 3), FDR = signif(fdr, 3))
print(tab2)
write.csv(tab2, file.path(OUT, "review_T_trends.csv"), row.names = FALSE)
res$t_trends <- tab2

# =============================================================================
# PART 3: GSE39582 non-overlapping scores + multivariable OS (FAP13)
# =============================================================================
cat("\n=== PART 3: GSE39582 ===\n")
expr <- read.csv(file.path(DATA, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(DATA, "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}
# Parse clinical fields
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
  tnm_t = if (length(get_field("tnm.t")) == n) get_field("tnm.t") else NA,
  tnm_n = if (length(get_field("tnm.n")) == n) get_field("tnm.n") else NA,
  os_event = as.numeric(if (length(get_field("os.event")) == n) get_field("os.event") else NA),
  os_delay = as.numeric(if (length(get_field("os.delay (months)")) == n) get_field("os.delay (months)") else NA),
  age = as.numeric(if (length(get_field("age.at.diagnosis (year)")) == n) get_field("age.at.diagnosis (year)") else NA),
  Sex = if (length(get_field("Sex")) == n) get_field("Sex") else NA,
  mmr = if (length(get_field("mmr.status")) == n) get_field("mmr.status") else NA,
  stringsAsFactors = FALSE
)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
gse$FAP13     <- zmean(gse, FAP13)
gse$matrix4   <- zmean(gse, MATRIX4)
gse$receptor2 <- zmean(gse, RECEPTOR2)
gse$FAP_CAF   <- zmean(gse, c(FAP13, MATRIX4))
gse$ECM6      <- zmean(gse, c(MATRIX4, RECEPTOR2))
gse$CLDN_core <- zmean(gse, CLDN_CORE)
cat("GSE39582 samples:", nrow(gse), "\n")

cg_f13m4 <- cor.test(gse$FAP13, gse$matrix4, method = "spearman")
cg_f13r2 <- cor.test(gse$FAP13, gse$receptor2, method = "spearman")
cg_ov    <- cor.test(gse$FAP_CAF, gse$ECM6, method = "spearman")
cat(sprintf("GSE FAP13-matrix4: rho=%.3f P=%.4f\n", cg_f13m4$estimate, cg_f13m4$p.value))
cat(sprintf("GSE FAP13-receptor2: rho=%.3f P=%.4f\n", cg_f13r2$estimate, cg_f13r2$p.value))
cat(sprintf("GSE FAP_CAF-ECM6(overlap): rho=%.3f P=%.4f\n", cg_ov$estimate, cg_ov$p.value))
res$gse_f13m4 <- c(rho = unname(cg_f13m4$estimate), P = cg_f13m4$p.value)
res$gse_f13r2 <- c(rho = unname(cg_f13r2$estimate), P = cg_f13r2$p.value)
res$gse_overlap <- c(rho = unname(cg_ov$estimate), P = cg_ov$p.value)

# T-stratified
gse$T_grp <- ifelse(grepl("^T[12]$", gse$tnm_t), "T1-2", ifelse(grepl("^T[34]", gse$tnm_t), "T3-4", NA))
gse_s <- gse[!is.na(gse$T_grp), ]
c12 <- cor.test(gse_s$FAP13[gse_s$T_grp == "T1-2"], gse_s$matrix4[gse_s$T_grp == "T1-2"], method = "spearman")
c34 <- cor.test(gse_s$FAP13[gse_s$T_grp == "T3-4"], gse_s$matrix4[gse_s$T_grp == "T3-4"], method = "spearman")
cat(sprintf("GSE FAP13-matrix4 T1-2: rho=%.3f (n=%d) | T3-4: rho=%.3f (n=%d)\n",
    c12$estimate, sum(gse_s$T_grp == "T1-2"), c34$estimate, sum(gse_s$T_grp == "T3-4")))

# Multivariable OS (FAP13) - replicate v3.4 claim: n=557, 190 deaths, HR=1.08
gse_os <- gse[!is.na(gse$os_delay) & gse$os_delay > 0 & !is.na(gse$os_event) & !is.na(gse$FAP13) &
              !is.na(gse$age) & !is.na(gse$Sex) & !is.na(gse$tnm_t), ]
# stage factor from tnm.t + tnm.n + tnm.m (simplify: use T group + N group as factor)
gse_os$Tg <- gse_os$tnm_t
gse_os$Ng <- ifelse(grepl("N0", gse_os$tnm_n), "N0", ifelse(grepl("N[12]", gse_os$tnm_n), "N+", "NA"))
gse_os <- gse_os[gse_os$Tg %in% c("T1","T2","T3","T4") & gse_os$Ng %in% c("N0","N+"), ]
cox_g <- coxph(Surv(os_delay, os_event) ~ FAP13 + age + Sex + factor(Tg) + factor(Ng), data = gse_os)
cat(sprintf("GSE OS model: n=%d, events=%d | FAP13 HR=%.2f (%.2f-%.2f) P=%.3f\n",
    nrow(gse_os), sum(gse_os$os_event),
    exp(coef(cox_g)["FAP13"]), exp(confint(cox_g)["FAP13", 1]), exp(confint(cox_g)["FAP13", 2]),
    summary(cox_g)$coefficients["FAP13", 5]))
ph_test <- cox.zph(cox_g)
cat(sprintf("PH global P=%.3f\n", ph_test$table[1, 3]))
res$gse_os <- c(n = nrow(gse_os), events = sum(gse_os$os_event),
                HR = exp(coef(cox_g)["FAP13"]), P = summary(cox_g)$coefficients["FAP13", 5])

# =============================================================================
# PART 4: CLDN members - BH correction check (v3.4 claims: none survive BH)
# =============================================================================
cat("\n=== PART 4: CLDN BH check (TCGA) ===\n")
CLDN10 <- c("CLDN1","CLDN2","CLDN3","CLDN4","CLDN7","CLDN8","CLDN12","CLDN15","CLDN18","CLDN23")
tcga$stage_grp2 <- ifelse(grepl("Stage I", tcga$ajcc_pathologic_stage), "Early", "Late")
cld_res <- data.frame(gene = character(), rho_early = numeric(), P_early = numeric(),
                      rho_late = numeric(), P_late = numeric(), stringsAsFactors = FALSE)
for (g in CLDN10) {
  if (!g %in% colnames(tcga)) { cld_res <- rbind(cld_res, data.frame(gene = g, rho_early = NA, P_early = NA, rho_late = NA, P_late = NA)); next }
  e <- tcga[tcga$stage_grp2 == "Early", ]; l <- tcga[tcga$stage_grp2 == "Late", ]
  ce <- tryCatch(cor.test(e[[g]], e$FAP_CAF, method = "spearman"), error = function(x) NULL)
  cl <- tryCatch(cor.test(l[[g]], l$FAP_CAF, method = "spearman"), error = function(x) NULL)
  cld_res <- rbind(cld_res, data.frame(
    gene = g,
    rho_early = ifelse(is.null(ce), NA, ce$estimate),
    P_early = ifelse(is.null(ce), NA, ce$p.value),
    rho_late = ifelse(is.null(cl), NA, cl$estimate),
    P_late = ifelse(is.null(cl), NA, cl$p.value)))
}
cld_res$FDR_early <- p.adjust(cld_res$P_early, method = "BH")
cld_res$FDR_late <- p.adjust(cld_res$P_late, method = "BH")
cat("CLDN BH check (any FDR < 0.05?):\n")
print(cld_res)
cat("n_sig_after_BH:", sum(cld_res$FDR_early < 0.05, na.rm = TRUE) + sum(cld_res$FDR_late < 0.05, na.rm = TRUE), "\n")
write.csv(cld_res, file.path(OUT, "review_CLDN_BH.csv"), row.names = FALSE)

# =============================================================================
# PART 5: MMR-stratified correlations (v3.4 claims dMMR 0.887 / pMMR 0.917)
# =============================================================================
cat("\n=== PART 5: MMR stratification (GSE39582) ===\n")
gsem <- gse[!is.na(gse$mmr) & gse$mmr %in% c("dMMR", "pMMR"), ]
cd <- cor.test(gsem$FAP_CAF[gsem$mmr == "dMMR"], gsem$ECM6[gsem$mmr == "dMMR"], method = "spearman")
cp <- cor.test(gsem$FAP_CAF[gsem$mmr == "pMMR"], gsem$ECM6[gsem$mmr == "pMMR"], method = "spearman")
fz <- function(r1, n1, r2, n2) { z <- (atanh(r1) - atanh(r2)) / sqrt(1/(n1-3) + 1/(n2-3)); 2 * pnorm(-abs(z)) }
pz <- fz(cd$estimate, sum(gsem$mmr == "dMMR"), cp$estimate, sum(gsem$mmr == "pMMR"))
cat(sprintf("MMR FAP-ECM: dMMR rho=%.3f (n=%d) | pMMR rho=%.3f (n=%d) | Fisher z P=%.3f\n",
    cd$estimate, sum(gsem$mmr == "dMMR"), cp$estimate, sum(gsem$mmr == "pMMR"), pz))
res$mmr <- c(dMMR = cd$estimate, pMMR = cp$estimate, fisher_z_P = pz)

# Save summary
sum_df <- data.frame(analysis = names(res), stringsAsFactors = FALSE)
for (nm in names(res)) {
  cat(sprintf("\n%s: %s\n", nm, paste(names(res[[nm]]), res[[nm]], sep = "=", collapse = " | ")))
}
saveRDS(res, file.path(OUT, "review_results.rds"))
cat("\n=== [DONE] REVIEW_R1 complete ===\n")
