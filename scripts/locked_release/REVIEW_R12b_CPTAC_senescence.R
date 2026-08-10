# -*- coding: utf-8 -*-
# H-C 蛋白层验证：CPTAC SASP 相关蛋白 × FAP-ECM / 受体蛋白
# 目的：检验"SASP 与 FAP-ECM 共变、与受体解耦"在蛋白层的表现
suppressPackageStartupMessages(library(dplyr))

PROJ <- 'C:/Users/Spica/WorkBuddy/Claw/deliverables/fap_early_crc'
OUT  <- file.path(PROJ, 'output', 'review_r12_senescence')

# CPTAC 蛋白数据（已有：FAP, COL1A1, COL1A2, FN1, SDC4, CD44, CLDN1, CLDN4）
cptac <- read.csv(file.path(PROJ, 'output', 'cptac', 'cptac_protein_sdc4_cd44.csv'), check.names = FALSE)
cat('CPTAC n =', nrow(cptac), '| cols:', paste(colnames(cptac), collapse = ', '), '\n')

# CPTAC 蛋白组可测的 SASP 相关蛋白（在现有 8 蛋白内）：
# COL1A1/COL1A2/FN1 是 ECM（也是 SASP 组分）；SDC4/CD44 是受体
# 用 ECM 蛋白评分（COL1A1/COL1A2/FN1，不含 FAP 避免循环）作为"间质衰老输出"代理
ECM_PROT <- c('COL1A1', 'COL1A2', 'FN1')
REC_PROT <- c('SDC4', 'CD44')

zmean_col <- function(df, cols) {
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) < 2) return(rep(NA_real_, nrow(df)))
  s <- as.data.frame(lapply(cols, function(g) {
    v <- as.numeric(df[[g]])
    (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
  }))
  rowMeans(s, na.rm = TRUE)
}

cptac$ECM_prot  <- zmean_col(cptac, ECM_PROT)
cptac$REC_prot  <- zmean_col(cptac, REC_PROT)
cptac$FAP_prot  <- as.numeric(cptac$FAP)

spear <- function(x, y) {
  ok <- complete.cases(x, y)
  ct <- cor.test(x[ok], y[ok], method = 'spearman')
  c(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}
report <- function(tag, rho, p, n) {
  cat(sprintf('%-38s rho=%7.3f  P=%10s  n=%d\n', tag, rho, format(p, digits = 3, scientific = TRUE), n))
}

cat('\n===== CPTAC 蛋白层 =====\n')
# 间质程序（FAP/ECM）vs 受体：解耦检验
r <- spear(cptac$FAP_prot, cptac$ECM_prot)
report('FAP蛋白 × ECM蛋白 (共变)', r['rho'], r['p'], r['n'])
r <- spear(cptac$FAP_prot, cptac$REC_prot)
report('FAP蛋白 × 受体蛋白 (解耦?)', r['rho'], r['p'], r['n'])
r <- spear(cptac$ECM_prot, cptac$REC_prot)
report('ECM蛋白 × 受体蛋白 (解耦?)', r['rho'], r['p'], r['n'])
# FAP-high vs low 的受体蛋白差异
med <- median(cptac$FAP_prot, na.rm = TRUE)
hi <- cptac$FAP_prot >= med; lo <- cptac$FAP_prot < med
for (g in c('SDC4', 'CD44')) {
  v <- as.numeric(cptac[[g]])
  w <- wilcox.test(v[hi], v[lo])
  report(paste0(g, ' 蛋白 FAP-hi vs lo'), NA, w$p.value, sum(hi) + sum(lo))
}

# 保存
res <- data.frame(
  test = c('FAP_prot_x_ECM_prot', 'FAP_prot_x_REC_prot', 'ECM_prot_x_REC_prot',
           'SDC4_FAP_hi_lo', 'CD44_FAP_hi_lo'),
  rho = c(spear(cptac$FAP_prot, cptac$ECM_prot)['rho'],
          spear(cptac$FAP_prot, cptac$REC_prot)['rho'],
          spear(cptac$ECM_prot, cptac$REC_prot)['rho'], NA, NA),
  p = c(spear(cptac$FAP_prot, cptac$ECM_prot)['p'],
        spear(cptac$FAP_prot, cptac$REC_prot)['p'],
        spear(cptac$ECM_prot, cptac$REC_prot)['p'],
        wilcox.test(as.numeric(cptac$SDC4[hi]), as.numeric(cptac$SDC4[lo]))$p.value,
        wilcox.test(as.numeric(cptac$CD44[hi]), as.numeric(cptac$CD44[lo]))$p.value),
  n = c(spear(cptac$FAP_prot, cptac$ECM_prot)['n'],
        spear(cptac$FAP_prot, cptac$REC_prot)['n'],
        spear(cptac$ECM_prot, cptac$REC_prot)['n'], sum(hi) + sum(lo), sum(hi) + sum(lo))
)
write.csv(res, file.path(OUT, 'R12_CPTAC_protein_layer.csv'), row.names = FALSE)
cat('\nDONE. Saved:', file.path(OUT, 'R12_CPTAC_protein_layer.csv'), '\n')
