.libPaths(c('library', 'library', .libPaths()))
suppressPackageStartupMessages({library(monocle3); library(Seurat); library(dplyr)})
OUTD <- './output'

cds <- readRDS(file.path(OUTD, 'A5_monocle3_cds.rds'))
pt <- data.frame(cell = colnames(cds), pseudotime = pseudotime(cds),
                 cluster = as.character(clusters(cds)))
sc <- readRDS('data/GSE132465/integrated_seurat_Harmony.rds')
stromal <- subset(sc, subset = Cell_type == 'Stromal cells')
fap <- GetAssayData(stromal, assay='RNA', layer='data')['FAP', ]
score_zmean <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  z <- t(scale(t(as.matrix(mat[genes, , drop=FALSE])))); colMeans(z, na.rm=TRUE)
}
ecm <- score_zmean(GetAssayData(stromal, assay='RNA', layer='data'),
                   c('COL1A1','COL1A2','COL3A1','FN1','SDC4','CD44'))

m <- merge(pt, data.frame(cell=names(fap), FAP=fap), by='cell', all.x=TRUE)
m <- merge(m, data.frame(cell=names(ecm), ECM=ecm), by='cell', all.x=TRUE)
m$FAP_status <- ifelse(m$FAP > 0, 'FAP+', 'FAP-')

cat('=== Per-cluster profile ===\n')
prof <- m %>% group_by(cluster) %>%
  summarise(n=n(), FAP_med=median(FAP), ECM_med=median(ECM),
            pt_med=median(pseudotime), FAP_pos_pct=mean(FAP_status=='FAP+')*100, .groups='drop')
print(as.data.frame(prof), row.names=FALSE)

cat('\n=== FAP+ vs FAP- pseudotime ===\n')
print(tapply(m$pseudotime, m$FAP_status, summary))
wt <- wilcox.test(m$pseudotime[m$FAP_status=='FAP+'], m$pseudotime[m$FAP_status=='FAP-'])
cat('Wilcoxon P:', format(wt$p.value, digits=3), '\n')

cat('\n=== Pseudotime deciles vs ECM (nonlinear check) ===\n')
m$pt_decile <- cut(m$pseudotime, breaks=10, labels=FALSE)
dec <- m %>% group_by(pt_decile) %>% summarise(ECM_mean=mean(ECM), FAP_mean=mean(FAP), n=n(), .groups='drop')
print(as.data.frame(dec), row.names=FALSE)

cat('\n=== Correlation in FAP+ subset only ===\n')
fap_pos <- m[m$FAP_status=='FAP+', ]
cat('n(FAP+):', nrow(fap_pos), '\n')
ct <- cor.test(fap_pos$pseudotime, fap_pos$ECM, method='spearman', exact=FALSE)
cat(sprintf('FAP+ only: pseudotime vs ECM rho=%.3f P=%.2e\n', ct$estimate, ct$p.value))

write.csv(m, file.path(OUTD, 'tables', 'A5_pseudotime_diag_full.csv'), row.names=FALSE)
cat('\n[DONE] diagnostic\n')
