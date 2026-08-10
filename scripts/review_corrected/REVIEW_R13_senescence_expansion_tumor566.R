# R13 v2: P0 Dual entry + P1 confounder/multi-signature/prognosis
# Unified format: samples x genes (columns = gene symbols)
# Run from the repository root or set FAP_REPO_ROOT and FAP_DATA_ROOT.
repo_root <- normalizePath(Sys.getenv('FAP_REPO_ROOT', unset = getwd()), mustWork = TRUE)
data_root <- normalizePath(Sys.getenv('FAP_DATA_ROOT', unset = file.path(repo_root, 'data')),
                           mustWork = FALSE)
extra_lib <- Sys.getenv('FAP_R_LIB', unset = '')
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
suppressPackageStartupMessages({library(dplyr)})
set.seed(20260807)
OUT <- file.path(repo_root, 'derived_results', 'senescence_corrected_rerun')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------- Gene sets ----------
FAP13 <- c('FAP','POSTN','THY1','PDPN','TAGLN','ACTA2','MMP2','MMP9','CXCL12','TGFB1','INHBA','WNT2','WNT5A')
matrix4 <- c('COL1A1','COL1A2','COL3A1','FN1')
receptor2 <- c('SDC4','CD44')
SASP25 <- c('IL6','CXCL8','IL1A','IL1B','CCL2','CCL5','CXCL1','CXCL2','CXCL3','CXCL10',
            'MMP1','MMP3','MMP9','MMP10','MMP13','SERPINE1','PLAU','TIMP2','VEGFA','GDF15',
            'IGFBP3','TNF','CSF2','HGF','FAS')
colonSASP7 <- c('GDF15','MMP3','CXCL8','CXCL1','CXCL5','STC1','CCL5')
read_genes <- function(p) {
  g <- character(0)
  for (line in readLines(p, warn = FALSE)) {
    line <- trimws(line)
    if (nchar(line) == 0 || startsWith(line, '#')) next
    g <- c(g, strsplit(line, '\\s+')[[1]])
  }
  unique(g)
}
senmayo <- read_genes(file.path(repo_root, 'config', 'gene_sets', 'SenMayo_125_genes.txt'))
cellage <- read_genes(file.path(repo_root, 'config', 'gene_sets', 'CellAge_inducing_genes.txt'))
score_sets <- list(SenMayo = senmayo, SASP25 = SASP25,
                   colonSASP = colonSASP7, CellAge = cellage)
overlap_background <- unique(c(FAP13, matrix4))
clean_sets <- lapply(score_sets, setdiff, y = overlap_background)
cat('SenMayo genes:', length(senmayo), '| CellAge inducing genes:', length(cellage), '\n')
for (nm in names(score_sets)) {
  ov <- intersect(score_sets[[nm]], overlap_background)
  cat(nm, 'overlap with FAP13/matrix4:',
      if (length(ov)) paste(ov, collapse = ',') else 'none', '\n')
}

# ---------- Helpers ----------
zmean <- function(m, genes) {          # m: samples  genes
  avail <- intersect(genes, colnames(m))
  if (length(avail) < 3) return(rep(NA_real_, nrow(m)))
  s <- apply(m[, avail, drop = FALSE], 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  rowMeans(s, na.rm = TRUE)
}
single_col <- function(m, gene) {
  g <- intersect(gene, colnames(m))
  if (length(g) == 0) return(rep(NA_real_, nrow(m)))
  as.numeric(m[, g[1]])
}
spear <- function(x, y) {
  k <- complete.cases(x, y)
  if (sum(k) < 10) return(c(n = sum(k), rho = NA, p = NA))
  ct <- suppressWarnings(cor.test(x[k], y[k], method = 'spearman', exact = FALSE))
  c(n = sum(k), rho = unname(ct$estimate), p = unname(ct$p.value))
}
partial_spear <- function(x, y, Z) {   # Z: matrix/df of confounders
  rx <- rank(x, na.last = 'keep'); ry <- rank(y, na.last = 'keep')
  Z <- as.data.frame(Z)
  rz <- do.call(data.frame, lapply(Z, function(z) rank(z, na.last = 'keep')))
  ok <- complete.cases(rx, ry, rz)
  if (sum(ok) < 10) return(c(n = sum(ok), rho = NA, p = NA))
  ex <- residuals(lm(rx[ok] ~ ., data = rz[ok, , drop = FALSE]))
  ey <- residuals(lm(ry[ok] ~ ., data = rz[ok, , drop = FALSE]))
  ct <- suppressWarnings(cor.test(ex, ey, method = 'spearman', exact = FALSE))
  c(n = sum(ok), rho = unname(ct$estimate), p = unname(ct$p.value))
}

# ---------- Data loading ----------
# TCGA380 (Xena, genes x samples) -> samples x genes
expr380 <- as.matrix(read.table(gzfile(file.path(data_root, 'TCGA', 'TCGA_COADREAD_expression.txt.gz')),
                                header = TRUE, row.names = 1, check.names = FALSE))
is_tumor <- grepl('^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-0[1-9]', colnames(expr380))
tc <- colnames(expr380)[is_tumor]; tc <- tc[!duplicated(substr(tc, 1, 12))]
m380 <- t(expr380[, tc, drop = FALSE])
cat('TCGA380:', nrow(m380), 'samples ', ncol(m380), 'genes\n')

# TCGA592 (cBioPortal CSV, samples x genes)
m592 <- NULL
f592 <- file.path(data_root, 'TCGA592', 'TCGA592_target_genes_expression.csv')
if (file.exists(f592)) {
  d <- read.csv(f592, check.names = FALSE, stringsAsFactors = FALSE)
  m592 <- as.data.frame(lapply(d[, -1], as.numeric)); rownames(m592) <- d$sample
  cat('TCGA592:', nrow(m592), 'samples ', ncol(m592), 'genes\n')
} else cat('WARN: n=592 \n')

# GSE39582 (series matrix full probe parsing -> samples x gene)
suppressPackageStartupMessages({library(hgu133plus2.db); library(AnnotationDbi)})
probe_anno <- suppressMessages(select(hgu133plus2.db, keys = keys(hgu133plus2.db),
                                      columns = 'SYMBOL', keytype = 'PROBEID'))
probe_anno <- probe_anno[!is.na(probe_anno$SYMBOL), ]
probe2gene <- setNames(probe_anno$SYMBOL, probe_anno$PROBEID)
gene2probes <- split(names(probe2gene), probe2gene)
gse_pheno <- read.csv(file.path(data_root, 'GSE39582', 'GSE39582_pheno.csv'),
                      check.names = FALSE, stringsAsFactors = FALSE)
target <- unique(c(FAP13, matrix4, receptor2, unlist(score_sets),
                   'MKI67','CDKN2A','CDKN2B'))
need_probes <- unique(unlist(gene2probes[intersect(names(gene2probes), target)]))
cat('GSE39582 target genes:', length(target), '| probes:', length(need_probes), '\n')
sm_file <- file.path(data_root, 'GSE39582', 'GSE39582_series_matrix.txt.gz')
con <- gzfile(sm_file, 'rt')
sample_names <- NULL; in_table <- FALSE
gene_vals <- list(); hdr <- NULL
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE); if (length(line) == 0) break
  if (startsWith(line, '!Sample_geo_accession')) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar('!Sample_geo_accession') + 1, nchar(line)), '\t')[[1]])
    sample_names <- sample_names[sample_names != '']
  } else if (startsWith(line, '!series_matrix_table_begin')) {
    in_table <- TRUE; next
  } else if (startsWith(line, '!series_matrix_table_end')) {
    break
  } else if (in_table) {
    cells <- strsplit(line, '\t')[[1]]
    probe <- gsub('"', '', cells[1])
    if (probe == 'ID_REF') { hdr <- cells; next }
    if (probe %in% need_probes) {
      vals <- suppressWarnings(as.numeric(gsub('"', '', cells[-1])))
      gene <- probe2gene[probe]
      if (!is.na(gene)) {
        if (is.null(gene_vals[[gene]])) gene_vals[[gene]] <- list()
        gene_vals[[gene]][[length(gene_vals[[gene]]) + 1]] <- vals
      }
    }
  }
}
close(con)
m_gse <- as.data.frame(lapply(gene_vals, function(v) if (length(v) == 1) v[[1]] else colMeans(do.call(rbind, v), na.rm = TRUE)))
rownames(m_gse) <- sample_names
# 仅保留肿瘤样本（pheno 匹配）
sample_map <- read.csv(file.path(repo_root, 'config', 'gse39582_sample_dataset.csv'),
                       check.names = FALSE, stringsAsFactors = FALSE)
tumour_ids <- sample_map$sample[sample_map$dataset != 'Non Tumoral']
m_gse <- m_gse[rownames(m_gse) %in% tumour_ids, , drop = FALSE]
gse_pheno <- gse_pheno[gse_pheno$sample %in% tumour_ids, , drop = FALSE]
stopifnot(nrow(m_gse) == 566L, nrow(gse_pheno) == 566L)
cat('GSE39582:', nrow(m_gse), 'samples; genes:', ncol(m_gse), '\n')
cat('GSE target hit:', sum(target %in% colnames(m_gse)), '/', length(target), '\n')

# Clinical (TCGA)
clin <- read.delim(file.path(data_root, 'TCGA', 'COADREAD_clinicalMatrix.txt'),
                   check.names = FALSE, stringsAsFactors = FALSE)
colnames(clin)[1] <- 'sample'
age_col <- intersect(c('age_at_initial_pathologic_diagnosis','age_at_diagnosis','age'), colnames(clin))
clin$age <- if (length(age_col)) suppressWarnings(as.numeric(clin[[age_col[1]]])) else NA
clin$sample_norm <- sub('\\-01.*$', '', clin$sample)
# OS from vital_status + days_to_death / days_to_last_followup
clin$os_event <- if ('vital_status' %in% colnames(clin)) {
  ifelse(clin$vital_status %in% c('DECEASED','Dead','1'), 1,
         ifelse(clin$vital_status %in% c('ALIVE','Live','0'), 0, NA))
} else NA
clin$os_days <- NA_real_
if ('days_to_death' %in% colnames(clin) && 'days_to_last_followup' %in% colnames(clin)) {
  dd <- suppressWarnings(as.numeric(clin$days_to_death))
  dl <- suppressWarnings(as.numeric(clin$days_to_last_followup))
  clin$os_days <- ifelse(clin$os_event == 1, dd, dl)
}
cat('clinical: age col=', age_col, '; OS events=', sum(clin$os_event == 1, na.rm = TRUE),
    '; OS time non-NA=', sum(!is.na(clin$os_days)), '\n')

# CAF abundance + EPIC purity
mcp <- read.csv(file.path(repo_root, 'derived_inputs', 'immune_deconvolution', 'mcpcounter_tumor_scores.csv'),
                check.names = FALSE, stringsAsFactors = FALSE)
mcp_fib <- data.frame(sample = names(mcp)[-1], caf = as.numeric(unlist(mcp[1, -1])))
epic <- read.csv(file.path(repo_root, 'derived_inputs', 'immune_deconvolution', 'epic_tumor_cell_fractions.csv'),
                 check.names = FALSE, stringsAsFactors = FALSE)
epic_tum <- data.frame(sample = names(epic)[-1], epic_tumor = as.numeric(unlist(epic[1, -1])))

# ---------- Scoring ----------
mk_scores <- function(m) data.frame(
  sample = rownames(m),
  FAP13 = zmean(m, FAP13), matrix4 = zmean(m, matrix4), receptor2 = zmean(m, receptor2),
  SenMayo = zmean(m, senmayo), SASP25 = zmean(m, SASP25), colonSASP = zmean(m, colonSASP7),
  CellAge = zmean(m, cellage),
  SenMayo_clean = zmean(m, clean_sets$SenMayo),
  SASP25_clean = zmean(m, clean_sets$SASP25),
  colonSASP_clean = zmean(m, clean_sets$colonSASP),
  CellAge_clean = zmean(m, clean_sets$CellAge),
  MKI67 = single_col(m, 'MKI67'), CDKN2A = single_col(m, 'CDKN2A'), CDKN2B = single_col(m, 'CDKN2B'),
  stringsAsFactors = FALSE)
sc380 <- mk_scores(m380)
sc592 <- if (!is.null(m592)) mk_scores(m592) else NULL
sc_gse <- mk_scores(m_gse)

coverage_one <- function(m, label) do.call(rbind, lapply(names(score_sets), function(nm) {
  full <- score_sets[[nm]]; clean <- clean_sets[[nm]]
  data.frame(dataset = label, signature = nm, source_genes = length(full),
             genes_after_overlap_removal = length(clean),
             observed_full = length(intersect(full, colnames(m))),
             observed_clean = length(intersect(clean, colnames(m))),
             removed_overlap = paste(intersect(full, overlap_background), collapse = ';'),
             stringsAsFactors = FALSE)
}))
coverage <- rbind(coverage_one(m380, 'TCGA380'),
                  if (!is.null(m592)) coverage_one(m592, 'TCGA592'),
                  coverage_one(m_gse, 'GSE39582_tumor566'))
write.csv(coverage, file.path(OUT, 'signature_coverage_and_overlap.csv'), row.names = FALSE)

# ============================================================
# P0: Dual entry SenMayoFAP13/matrix4 + MKI67  + CDKN2A/B
# ============================================================
run_p0 <- function(sc, label) {
  out <- list()
  for (g in c('FAP13','matrix4','receptor2')) {
    out[[paste0('SenMayo_', g)]] <- spear(sc$SenMayo, sc[[g]])
    out[[paste0('SASP25_', g)]] <- spear(sc$SASP25, sc[[g]])
    out[[paste0('colonSASP_', g)]] <- spear(sc$colonSASP, sc[[g]])
  }
  out[['SenMayo_FAP13_MKI67adj']] <- partial_spear(sc$SenMayo, sc$FAP13, data.frame(mk = sc$MKI67))
  out[['SenMayo_matrix4_MKI67adj']] <- partial_spear(sc$SenMayo, sc$matrix4, data.frame(mk = sc$MKI67))
  med <- median(sc$FAP13, na.rm = TRUE); hi <- !is.na(sc$FAP13) & sc$FAP13 >= med; lo <- !is.na(sc$FAP13) & sc$FAP13 < med
  out[['CDKN2A_p']] <- if (sum(hi) >= 5 && sum(lo) >= 5 && sum(!is.na(sc$CDKN2A[hi])) >= 3 &&
                            sum(!is.na(sc$CDKN2A[lo])) >= 3)
    suppressWarnings(wilcox.test(sc$CDKN2A[hi], sc$CDKN2A[lo])$p.value) else NA
  out[['CDKN2B_p']] <- if (sum(hi) >= 5 && sum(lo) >= 5 && sum(!is.na(sc$CDKN2B[hi])) >= 3 &&
                            sum(!is.na(sc$CDKN2B[lo])) >= 3)
    suppressWarnings(wilcox.test(sc$CDKN2B[hi], sc$CDKN2B[lo])$p.value) else NA
  out[['n']] <- nrow(sc)
  out
}
r380 <- run_p0(sc380, 'TCGA380'); r592 <- if (!is.null(sc592)) run_p0(sc592, 'TCGA592') else NULL
r_gse <- run_p0(sc_gse, 'GSE39582')

fmtp0 <- function(r, k) sprintf('%s: rho=%.3f (n=%d, P=%s)', k, r[[k]]['rho'], r[[k]]['n'],
                                format(r[[k]]['p'], digits = 3, scientific = TRUE))
cat('\n===== P0 Dual entry =====\n')
cat('[TCGA n=380] ', fmtp0(r380, 'SenMayo_FAP13'), ' | ', fmtp0(r380, 'SenMayo_matrix4'), '\n', sep = '')
cat('[TCGA n=380] ', fmtp0(r380, 'SenMayo_FAP13_MKI67adj'), ' (MKI67 adj)\n', sep = '')
if (!is.null(r592)) {
  cat('[TCGA n=592] ', fmtp0(r592, 'SenMayo_FAP13'), ' | ', fmtp0(r592, 'SenMayo_matrix4'), '\n', sep = '')
  cat('[TCGA n=592] ', fmtp0(r592, 'SenMayo_FAP13_MKI67adj'), ' (MKI67 adj)\n', sep = '')
}
cat('[GSE39582]  ', fmtp0(r_gse, 'SenMayo_FAP13'), ' | ', fmtp0(r_gse, 'SenMayo_matrix4'), '\n', sep = '')
cat('[GSE39582]  CDKN2A p=', format(r_gse$CDKN2A_p, digits=3, scientific=TRUE),
    '; CDKN2B p=', format(r_gse$CDKN2B_p, digits=3, scientific=TRUE), '\n', sep='')

summ <- function(r, label) data.frame(entry = label,
  SenMayo_FAP13 = r$SenMayo_FAP13['rho'], SenMayo_FAP13_p = r$SenMayo_FAP13['p'],
  SenMayo_matrix4 = r$SenMayo_matrix4['rho'], SenMayo_matrix4_p = r$SenMayo_matrix4['p'],
  SenMayo_FAP13_MKI67adj = r$SenMayo_FAP13_MKI67adj['rho'],
  SASP25_FAP13 = r$SASP25_FAP13['rho'], colonSASP_FAP13 = r$colonSASP_FAP13['rho'],
  CDKN2A_p = r$CDKN2A_p, CDKN2B_p = r$CDKN2B_p, n = r$n)
p0_tab <- rbind(summ(r380, 'TCGA380'),
                if (!is.null(r592)) summ(r592, 'TCGA592'),
                summ(r_gse, 'GSE39582'))
write.csv(p0_tab, file.path(OUT, 'P0_dual_entry_senescence.csv'), row.names = FALSE)
cat('P0 \n')

# ============================================================
# P1a: Confounder adjustmentTCGA 380
# ============================================================
sc380m <- merge(sc380, mcp_fib, by = 'sample', all.x = TRUE)
sc380m <- merge(sc380m, epic_tum, by = 'sample', all.x = TRUE)
sc380m$sample_norm <- sub('\\-01.*$', '', sc380m$sample)
sc380m <- merge(sc380m, clin[, c('sample_norm','age')], by = 'sample_norm', all.x = TRUE)
conf <- list(
  raw = spear(sc380m$SenMayo, sc380m$FAP13),
  adj_CAF = partial_spear(sc380m$SenMayo, sc380m$FAP13, data.frame(caf = sc380m$caf)),
  adj_purity = partial_spear(sc380m$SenMayo, sc380m$FAP13, data.frame(pur = sc380m$epic_tumor)),
  adj_age = partial_spear(sc380m$SenMayo, sc380m$FAP13, data.frame(age = sc380m$age)),
  adj_MKI67 = partial_spear(sc380m$SenMayo, sc380m$FAP13, data.frame(mk = sc380m$MKI67)),
  adj_CAF_age = partial_spear(sc380m$SenMayo, sc380m$FAP13, data.frame(caf = sc380m$caf, age = sc380m$age)),
  adj_CAF_purity_MKI67 = partial_spear(sc380m$SenMayo, sc380m$FAP13,
    data.frame(caf = sc380m$caf, pur = sc380m$epic_tumor, mk = sc380m$MKI67)))
cat('\n===== P1a Confounder adjustmentTCGA 380, SenMayo~FAP13=====\n')
for (nm in names(conf)) cat(sprintf('%-24s rho=%.3f (n=%d, P=%s)\n', nm, conf[[nm]]['rho'],
    conf[[nm]]['n'], format(conf[[nm]]['p'], digits = 3, scientific = TRUE)))
conf_tab <- do.call(rbind, lapply(names(conf), function(nm)
  data.frame(comparison = nm, rho = conf[[nm]]['rho'], n = conf[[nm]]['n'], p = conf[[nm]]['p'])))
write.csv(conf_tab, file.path(OUT, 'P1a_confounder_adjustment.csv'), row.names = FALSE)

# ============================================================
# P1b: TCGA 380 + GSE39582
# ============================================================
cat('\n===== P1b Multi-signature consistency =====\n')
sig_names <- c('SenMayo', 'SenMayo_clean', 'SASP25', 'SASP25_clean',
               'colonSASP', 'colonSASP_clean', 'CellAge', 'CellAge_clean')
multi <- list()
for (sn in sig_names) {
  for (g in c('FAP13','matrix4')) {
    r1 <- spear(sc380m[[sn]], sc380m[[g]])
    r2 <- spear(sc_gse[[sn]], sc_gse[[g]])
    multi[[paste0(sn,'_',g,'_TCGA380')]] <- r1
    multi[[paste0(sn,'_',g,'_GSE39582')]] <- r2
    cat(sprintf('%-26s TCGA rho=%.3f (P=%s) | GSE rho=%.3f (P=%s)\n', paste0(sn,'',g),
        r1['rho'], format(r1['p'], digits=3, scientific=TRUE),
        r2['rho'], format(r2['p'], digits=3, scientific=TRUE)))
  }
}
multi_tab <- do.call(rbind, lapply(names(multi), function(nm)
  data.frame(comparison = nm, rho = multi[[nm]]['rho'], n = multi[[nm]]['n'], p = multi[[nm]]['p'])))
write.csv(multi_tab, file.path(OUT, 'P1b_multisignature.csv'), row.names = FALSE)

# ============================================================
# P1c: Scoring Cox
# ============================================================
suppressPackageStartupMessages(library(survival))
cat('\n===== P1c  =====\n')
cox_out <- list()
# --- TCGA ---
if ('os_event' %in% colnames(clin)) {
  tcga_surv <- data.frame(sample_norm = clin$sample_norm,
                          os_event = clin$os_event,
                          os_time = clin$os_days / 365.25)
  surv_df <- merge(sc380m, tcga_surv, by = 'sample_norm')
  surv_df <- surv_df[complete.cases(surv_df[, c('os_event','os_time','FAP13','matrix4','SenMayo','SASP25')]), ]
  cat('TCGA Cox samples:', nrow(surv_df), '\n')
  for (g in c('FAP13','matrix4','SenMayo','SASP25','colonSASP')) {
    if (sum(!is.na(surv_df[[g]])) < 30) next
    m <- coxph(as.formula(paste0('Surv(os_time, os_event) ~ scale(', g, ')')), data = surv_df)
    s <- summary(m)$coefficients[1, ]
    cox_out[[paste0('TCGA_', g)]] <- c(HR = unname(exp(s[1])),
                                       ci_low = unname(exp(s[1]-1.96*s[3])),
                                       ci_high = unname(exp(s[1]+1.96*s[3])),
                                       P = unname(s[5]), n = nrow(surv_df))
    cat(sprintf('TCGA %-10s HR=%.2f (95%%CI %.2f-%.2f, P=%s)\n', g,
        exp(s[1]), exp(s[1]-1.96*s[3]), exp(s[1]+1.96*s[3]),
        format(s[5], digits=3, scientific=TRUE)))
  }
}
# --- GSE39582os_delay_months / 12  ---
if ('os_event' %in% colnames(gse_pheno) && 'os_delay_months' %in% colnames(gse_pheno)) {
  gse_surv <- data.frame(sample = gse_pheno$sample,
                         os_event = as.numeric(gse_pheno$os_event),
                         os_time = as.numeric(gse_pheno$os_delay_months) / 12)
  surv_gse <- merge(sc_gse, gse_surv, by = 'sample')
  surv_gse <- surv_gse[complete.cases(surv_gse[, c('os_event','os_time','FAP13','matrix4','SenMayo')]), ]
  cat('GSE Cox samples:', nrow(surv_gse), '\n')
  for (g in c('FAP13','matrix4','SenMayo','SASP25','colonSASP')) {
    if (sum(!is.na(surv_gse[[g]])) < 30) next
    m <- coxph(as.formula(paste0('Surv(os_time, os_event) ~ scale(', g, ')')), data = surv_gse)
    s <- summary(m)$coefficients[1, ]
    cox_out[[paste0('GSE_', g)]] <- c(HR = unname(exp(s[1])),
                                      ci_low = unname(exp(s[1]-1.96*s[3])),
                                      ci_high = unname(exp(s[1]+1.96*s[3])),
                                      P = unname(s[5]), n = nrow(surv_gse))
    cat(sprintf('GSE %-10s HR=%.2f (95%%CI %.2f-%.2f, P=%s)\n', g,
        exp(s[1]), exp(s[1]-1.96*s[3]), exp(s[1]+1.96*s[3]),
        format(s[5], digits=3, scientific=TRUE)))
  }
}
if (length(cox_out)) {
  cox_tab <- do.call(rbind, lapply(names(cox_out), function(nm) {
    v <- cox_out[[nm]]
    if (is.na(v['HR'])) return(NULL)
    data.frame(comparison = nm, HR = v['HR'], ci_low = v['ci_low'],
               ci_high = v['ci_high'], P = v['P'], n = v['n'])
  }))
  write.csv(cox_tab, file.path(OUT, 'P1c_prognosis_cox.csv'), row.names = FALSE)
  cat('P1c CSV saved\n')
}
cat('\nR13 done:', OUT, '\n')
