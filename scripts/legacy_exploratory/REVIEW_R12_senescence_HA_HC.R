# -*- coding: utf-8 -*-
# H-A / H-C 衰老分析 v2：SenMayo / SASP × FAP13/matrix4/receptor2 + MKI67 校正
# 修复：cor.test n 用 length；GSE39582 用 hgu133plus2.db 全基因映射
suppressPackageStartupMessages({library(dplyr); library(hgu133plus2.db); library(AnnotationDbi)})

set.seed(20260806)
PROJ <- '.'
OUT  <- file.path(PROJ, 'output', 'review_r12_senescence')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------- 基因集 ----------
senmayo_raw <- readLines(file.path(PROJ, 'data', 'senmayo', 'SenMayo_125_genes.txt'))
senmayo <- unique(unlist(strsplit(senmayo_raw, '\\s+', perl = TRUE)))
senmayo <- senmayo[!grepl('^#', senmayo) & nchar(senmayo) > 0]
FAP13 <- c('FAP','POSTN','THY1','PDPN','TAGLN','ACTA2','MMP2','MMP9','CXCL12','TGFB1','INHBA','WNT2','WNT5A')
matrix4 <- c('COL1A1','COL1A2','COL3A1','FN1')
receptor2 <- c('SDC4','CD44')
SASP <- c('IL6','CXCL8','IL1A','IL1B','CCL2','CCL5','CXCL1','CXCL2','CXCL3','CXCL10',
          'MMP1','MMP3','MMP9','MMP10','MMP13','SERPINE1','PLAU','TIMP2','VEGFA','GDF15',
          'IGFBP3','TNF','CSF2','HGF','FAS')

zmean <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) < 2) return(rep(NA_real_, ncol(mat)))
  s <- apply(mat[genes, , drop = FALSE], 1, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  if (is.null(dim(s))) return(s)
  rowMeans(s, na.rm = TRUE)
}
spear2 <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 10) return(c(rho = NA, p = NA, n = sum(ok)))
  ct <- cor.test(x[ok], y[ok], method = 'spearman')
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
partial_spear <- function(x, y, z) {
  rx <- rank(x, na.last = 'keep'); ry <- rank(y, na.last = 'keep'); rz <- rank(z, na.last = 'keep')
  ok <- complete.cases(rx, ry, rz)
  if (sum(ok) < 10) return(c(rho = NA, p = NA, n = sum(ok)))
  ex <- residuals(lm(rx[ok] ~ rz[ok])); ey <- residuals(lm(ry[ok] ~ rz[ok]))
  ct <- cor.test(ex, ey, method = 'spearman')
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
report <- function(tag, cohort, rho, p, n, note = '') {
  r <- ifelse(is.na(rho), NA, round(rho, 3))
  cat(sprintf('%-34s %-9s rho=%7s  P=%10s  n=%d  %s\n', tag, cohort,
              ifelse(is.na(r), 'NA', format(r)), format(p, digits = 3, scientific = TRUE), n, note))
}

results <- list()

# ============================================================
# 1. TCGA-COADREAD（Xena HiSeqV2）
# ============================================================
cat('===== TCGA-COADREAD =====\n')
expr_file <- 'data/TCGA/TCGA_COADREAD_expression.txt.gz'
expr_raw <- read.table(gzfile(expr_file), header = TRUE, row.names = 1, check.names = FALSE)
expr_raw <- as.matrix(expr_raw)
is_tumor <- grepl('^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-0[1-9]', colnames(expr_raw))
tumor_cols <- colnames(expr_raw)[is_tumor]
tumor_cols <- tumor_cols[!duplicated(substr(tumor_cols, 1, 12))]
tcga <- expr_raw[, tumor_cols, drop = FALSE]
cat('TCGA unique tumors:', ncol(tcga), '\n')

scores_t <- list(
  FAP13 = zmean(tcga, FAP13), matrix4 = zmean(tcga, matrix4),
  receptor2 = zmean(tcga, receptor2), SenMayo = zmean(tcga, senmayo), SASP = zmean(tcga, SASP))
scores_t$MKI67 <- if ('MKI67' %in% rownames(tcga)) (tcga['MKI67', ] - mean(tcga['MKI67', ])) / sd(tcga['MKI67', ]) else rep(NA, ncol(tcga))
cat('SenMayo detected in TCGA:', sum(senmayo %in% rownames(tcga)), '/', length(senmayo), '\n')
cat('SASP detected in TCGA:', sum(SASP %in% rownames(tcga)), '/', length(SASP), '\n')

pairs_t <- list(c('SenMayo','FAP13'), c('SenMayo','matrix4'), c('SenMayo','receptor2'),
                c('SASP','FAP13'), c('SASP','matrix4'), c('SASP','receptor2'))
for (pair in pairs_t) {
  r <- spear2(scores_t[[pair[1]]], scores_t[[pair[2]]])
  report(paste(pair[1], '×', pair[2]), 'TCGA', r['rho'], r['p'], r['n'])
  results[[paste0('TCGA_', pair[1], '_x_', pair[2])]] <- r
}
# MKI67 部分相关
r1 <- partial_spear(scores_t$FAP13, scores_t$SenMayo, scores_t$MKI67)
report('FAP13 ~ SenMayo | MKI67', 'TCGA', r1['rho'], r1['p'], r1['n'], '(partial rank)')
results[['TCGA_FAP13_SenMayo_partial_MKI67']] <- r1
r2 <- partial_spear(scores_t$matrix4, scores_t$SenMayo, scores_t$MKI67)
report('matrix4 ~ SenMayo | MKI67', 'TCGA', r2['rho'], r2['p'], r2['n'], '(partial rank)')
results[['TCGA_matrix4_SenMayo_partial_MKI67']] <- r2
# H-C：SASP×matrix4 共变 vs SASP×receptor2 解耦
r3 <- spear2(scores_t$SASP, scores_t$matrix4)
report('SASP × matrix4 (ECM 共变)', 'TCGA', r3['rho'], r3['p'], r3['n'])
results[['TCGA_SASP_x_matrix4']] <- r3
r4 <- spear2(scores_t$SASP, scores_t$receptor2)
report('SASP × receptor2 (解耦)', 'TCGA', r4['rho'], r4['p'], r4['n'])
results[['TCGA_SASP_x_receptor2']] <- r4
# CDKN2A/2B
med <- median(scores_t$FAP13, na.rm = TRUE)
hi <- scores_t$FAP13 >= med; lo <- scores_t$FAP13 < med
for (g in c('CDKN2A', 'CDKN2B')) {
  if (g %in% rownames(tcga)) {
    v <- tcga[g, ]
    w <- wilcox.test(v[hi], v[lo])
    report(paste(g, 'FAP13-hi vs lo'), 'TCGA', NA, w$p.value, sum(hi) + sum(lo),
           sprintf('med hi %.2f / lo %.2f', median(v[hi]), median(v[lo])))
    results[[paste0('TCGA_', g, '_FAP13_hi_lo')]] <- c(rho = NA, p = w$p.value, n = sum(hi) + sum(lo))
  }
}

# ============================================================
# 2. GSE39582（全基因组探针 → 基因）
# ============================================================
cat('\n===== GSE39582 =====\n')
# 探针→基因映射（hgu133plus2.db）
probe_anno <- select(hgu133plus2.db, keys = keys(hgu133plus2.db), columns = 'SYMBOL', keytype = 'PROBEID')
probe_anno <- probe_anno[!is.na(probe_anno$SYMBOL), ]
probe2gene <- setNames(probe_anno$SYMBOL, probe_anno$PROBEID)
gene2probes <- split(names(probe2gene), probe2gene)

# 从 series matrix 全量解析目标基因表达
sm_file <- file.path(PROJ, 'data', 'GSE39582', 'GSE39582_series_matrix.txt.gz')
target_genes <- unique(c(senmayo, SASP, FAP13, matrix4, receptor2, 'MKI67', 'CDKN2A', 'CDKN2B'))
need_probes <- unique(unlist(gene2probes[intersect(names(gene2probes), target_genes)]))
cat('GSE39582 目标基因数:', length(target_genes), '| 对应探针数:', length(need_probes), '\n')

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
cat('解析到基因数:', length(gene_vals), '\n')

# 汇总基因表达（多探针取均值）
gse_mat <- do.call(rbind, lapply(names(gene_vals), function(g) {
  m <- do.call(rbind, gene_vals[[g]])
  colMeans(m, na.rm = TRUE)
}))
rownames(gse_mat) <- names(gene_vals)
colnames(gse_mat) <- sample_names
cat('GSE39582 基因矩阵:', nrow(gse_mat), 'x', ncol(gse_mat), '\n')

# 排除非肿瘤
con2 <- gzfile(sm_file, 'rt'); field_rows <- list()
while (TRUE) {
  line <- readLines(con2, n = 1, warn = FALSE); if (length(line) == 0) break
  if (startsWith(line, '!Sample_characteristics_ch1')) {
    cells <- strsplit(substr(line, nchar('!Sample_characteristics_ch1') + 1, nchar(line)), '\t')[[1]]
    field_rows[[length(field_rows) + 1]] <- gsub('"', '', trimws(cells))
  }
}
close(con2)
field_names <- sapply(field_rows, function(r) {
  hit <- which(grepl(':', r))[1]; if (!is.na(hit)) sub(':.*', '', r[hit]) else ''
})
dataset <- field_rows[[which(field_names == 'dataset')[1]]]
dataset <- dataset[dataset != '']; dataset <- trimws(sub('dataset:', '', dataset, fixed = TRUE))
tumor_mask <- dataset != 'Non Tumoral'
tumor_mask[is.na(tumor_mask)] <- TRUE
gse_t <- gse_mat[, tumor_mask, drop = FALSE]
cat('GSE39582 tumor samples:', ncol(gse_t), '\n')

scores_g <- list(
  FAP13 = zmean(gse_t, FAP13), matrix4 = zmean(gse_t, matrix4),
  receptor2 = zmean(gse_t, receptor2), SenMayo = zmean(gse_t, senmayo), SASP = zmean(gse_t, SASP))
scores_g$MKI67 <- if ('MKI67' %in% rownames(gse_t)) (gse_t['MKI67', ] - mean(gse_t['MKI67', ])) / sd(gse_t['MKI67', ]) else rep(NA, ncol(gse_t))
cat('SenMayo detected in GSE39582:', sum(senmayo %in% rownames(gse_t)), '/', length(senmayo), '\n')
cat('SASP detected in GSE39582:', sum(SASP %in% rownames(gse_t)), '/', length(SASP), '\n')

for (pair in pairs_t) {
  r <- spear2(scores_g[[pair[1]]], scores_g[[pair[2]]])
  report(paste(pair[1], '×', pair[2]), 'GSE39582', r['rho'], r['p'], r['n'])
  results[[paste0('GSE_', pair[1], '_x_', pair[2])]] <- r
}
r1 <- partial_spear(scores_g$FAP13, scores_g$SenMayo, scores_g$MKI67)
report('FAP13 ~ SenMayo | MKI67', 'GSE39582', r1['rho'], r1['p'], r1['n'], '(partial rank)')
results[['GSE_FAP13_SenMayo_partial_MKI67']] <- r1
r3 <- spear2(scores_g$SASP, scores_g$matrix4)
report('SASP × matrix4 (ECM 共变)', 'GSE39582', r3['rho'], r3['p'], r3['n'])
results[['GSE_SASP_x_matrix4']] <- r3
r4 <- spear2(scores_g$SASP, scores_g$receptor2)
report('SASP × receptor2 (解耦)', 'GSE39582', r4['rho'], r4['p'], r4['n'])
results[['GSE_SASP_x_receptor2']] <- r4

# ============================================================
# 3. 保存结果
# ============================================================
res_df <- do.call(rbind, lapply(names(results), function(k) {
  v <- results[[k]]
  data.frame(test = k, rho = v['rho'], p = v['p'], n = v['n'], stringsAsFactors = FALSE)
}))
write.csv(res_df, file.path(OUT, 'R12_senescence_results.csv'), row.names = FALSE)
cat('\nDONE. Saved:', file.path(OUT, 'R12_senescence_results.csv'), '\n')
