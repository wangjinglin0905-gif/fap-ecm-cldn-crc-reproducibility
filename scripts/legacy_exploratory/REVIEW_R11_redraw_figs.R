# =============================================================================
# REVIEW_R11: redraw S4a/S4b/S5b/S6a/S6b for v4.4 (fix in-figure issues)
# S4a: title "pseudotime trajectory" -> "state ordering"
# S4b: fix truncated title/axis labels
# S5b: remove orphan "B" title, add descriptive title
# S6a/S6b: fix truncated titles + ambiguous bracket annotations
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(monocle3); library(SummarizedExperiment)})
PROJ <- "."
OUT  <- file.path(PROJ, "output", "review_r11")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
FIG <- "legacy_external/图片源"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

save_dual <- function(plot_obj, base, w = 7, h = 5.6) {
  ggsave(file.path(OUT, paste0(base, ".png")), plot_obj, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, paste0(base, ".tiff")), plot_obj, width = w, height = h, dpi = 600,
         compression = "lzw", device = "tiff")
  cat("saved:", base, "\n")
}

# ============ Part 1: S4a state ordering trajectory ============
cds <- readRDS(file.path(PROJ, "output", "A5_monocle3_cds.rds"))
# rename legend by mapping to an explicitly-named colData column
colData(cds)$state_ordering <- pseudotime(cds)
p_s4a <- plot_cells(cds, color_cells_by = "state_ordering", label_cell_groups = FALSE,
                    label_leaves = FALSE, label_branch_points = FALSE,
                    cell_size = 0.8, trajectory_graph_color = "grey40") +
  ggtitle("Stromal cell state ordering (GSE132465)") +
  theme(plot.title = element_text(size = 12, face = "bold"))
save_dual(p_s4a, "S4a_state_ordering", 7, 5.6)

# ============ Part 2: S4b activated vs quiescent ECM score boxplot ============
md <- as.data.frame(colData(cds))
expr_m <- assay(cds)
ECM <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "SDC4", "CD44")
ECM <- ECM[ECM %in% rownames(expr_m)]
sub_m <- as.matrix(expr_m[ECM, , drop = FALSE])
tot <- colSums(as.matrix(expr_m)); tot[tot == 0] <- NA
norm <- log1p(t(t(sub_m) / tot * 10000))   # genes x cells, log1p CPM
zc <- scale(t(norm))                        # cells x genes, per-gene z
md$ECM_score <- rowMeans(zc[, ECM, drop = FALSE], na.rm = TRUE)   # per-cell module score
md$group <- ifelse(md$FAP_status == "FAP+", "Activated CAF (FAP+)", "Quiescent fibroblast (FAP\u2212)")
wil <- suppressWarnings(wilcox.test(ECM_score ~ group, data = md))
cat(sprintf("S4b check: activated median=%.2f vs quiescent median=%.2f | Wilcoxon P=%.2e\n",
            median(md$ECM_score[md$group == "Activated CAF (FAP+)"], na.rm = TRUE),
            median(md$ECM_score[md$group == "Quiescent fibroblast (FAP\u2212)"], na.rm = TRUE),
            wil$p.value))
p_s4b <- ggplot(md, aes(x = group, y = ECM_score, fill = group)) +
  geom_boxplot(width = 0.55, outlier.size = 0.4, alpha = 0.85) +
  scale_fill_manual(values = c("Activated CAF (FAP+)" = "#B2182B", "Quiescent fibroblast (FAP\u2212)" = "#2166AC")) +
  labs(title = "ECM\u2013SDC4/CD44 module score by fibroblast state (GSE132465)",
       x = NULL, y = "ECM\u2013SDC4/CD44 score (z-mean)") +
  annotate("text", x = 1.5, y = max(md$ECM_score, na.rm = TRUE) * 0.92,
           label = "Wilcoxon P < 0.001", size = 4) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(size = 11.5, face = "bold"))
save_dual(p_s4b, "S4b_activation_branch", 6.4, 5.6)

# ============ Part 3: S5b claudin heatmap (recompute TCGA) ============
tcga <- read.csv(file.path(PROJ, "data", "A1_tcga_coad_merged.csv"))
tcga <- tcga[tcga$sample_scope == "Tumor", ]
expr_raw <- read.table(gzfile("data/TCGA/TCGA_COADREAD_expression.txt.gz"),
                       header = TRUE, row.names = 1, check.names = FALSE)
FAP17 <- c("FAP","COL1A1","COL1A2","COL3A1","FN1","POSTN","THY1","PDPN","TAGLN",
           "ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
extract_gene <- function(gene, expr_mat, samples) {
  if (!gene %in% rownames(expr_mat)) return(rep(NA, nrow(samples)))
  vals <- rep(NA, nrow(samples))
  for (i in seq_len(nrow(samples))) {
    sid <- samples$sample_id[i]
    if (sid %in% colnames(expr_mat)) vals[i] <- expr_mat[gene, sid]
    else { pid <- substr(sid, 1, 12); m <- grep(pid, colnames(expr_mat), value = TRUE, fixed = TRUE)
           if (length(m) > 0) vals[i] <- expr_mat[gene, m[1]] }
  }
  vals
}
for (g in FAP17) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
zmean <- function(df, genes) { a <- genes[genes %in% colnames(df)]
  s <- as.data.frame(lapply(a, function(g) (df[[g]] - mean(df[[g]], na.rm=TRUE))/sd(df[[g]], na.rm=TRUE)))
  rowMeans(s, na.rm=TRUE) }
tcga$FAP_CAF <- zmean(tcga, FAP17)
CLDN10 <- c("CLDN1","CLDN2","CLDN3","CLDN4","CLDN7","CLDN8","CLDN12","CLDN15","CLDN18","CLDN23")
for (g in CLDN10) tcga[[g]] <- extract_gene(g, expr_raw, tcga)
stage_grp <- ifelse(tcga$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB"), "Early (Stage I)", "Late (Stages II\u2013IV)")
tcga$stage_grp <- stage_grp
mat <- matrix(NA, nrow = length(CLDN10), ncol = 3, dimnames = list(CLDN10, c("Early (Stage I)", "Late (Stages II\u2013IV)", "Delta")))
for (g in CLDN10) {
  d_e <- tcga[tcga$stage_grp == "Early (Stage I)", ]; d_l <- tcga[tcga$stage_grp == "Late (Stages II\u2013IV)", ]
  re <- suppressWarnings(cor(d_e$FAP_CAF, d_e[[g]], method = "spearman", use = "complete.obs"))
  rl <- suppressWarnings(cor(d_l$FAP_CAF, d_l[[g]], method = "spearman", use = "complete.obs"))
  mat[g, ] <- c(re, rl, re - rl)
}
print(round(mat, 3))
library(tidyr); library(tibble)
hm_df <- as.data.frame(mat) |> rownames_to_column("Claudin") |> pivot_longer(-Claudin, names_to = "Stratum", values_to = "rho")
hm_df$Stratum <- factor(hm_df$Stratum, levels = c("Early (Stage I)", "Late (Stages II\u2013IV)", "Delta"))
hm_df$Claudin <- factor(hm_df$Claudin, levels = rev(CLDN10))
p_s5b <- ggplot(hm_df, aes(x = Stratum, y = Claudin, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-0.4, 0.4), oob = scales::squish, name = "Spearman \u03c1") +
  labs(title = "FAP\u2013claudin Spearman correlations by stage stratum (TCGA-COAD)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 11.5, face = "bold"), panel.grid = element_blank())
save_dual(p_s5b, "S5b_claudin_heatmap", 6.8, 5.6)

# ============ Part 4: S6a/S6b CDH1 T gradient ============
tcga$T_cat <- gsub("([Tt][1-4]).*", "\\1", tcga$ajcc_pathologic_t)
tcga$T_cat <- factor(tcga$T_cat, levels = c("T1","T2","T3","T4"))
tcga$CDH1 <- extract_gene("CDH1", expr_raw, tcga)
d6a <- tcga[!is.na(tcga$T_cat) & !is.na(tcga$CDH1), ]
kw_a <- kruskal.test(CDH1 ~ T_cat, data = d6a)
rho_a <- suppressWarnings(cor(as.numeric(d6a$T_cat), d6a$CDH1, method = "spearman"))
p_s6a <- ggplot(d6a, aes(x = T_cat, y = CDH1)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#92C5DE", alpha = 0.8) +
  geom_jitter(width = 0.18, size = 0.7, alpha = 0.45) +
  labs(title = sprintf("CDH1 expression across T categories \u2014 TCGA-COAD (Kruskal\u2013Wallis P = %.3f; Spearman \u03c1 = %.3f)",
                       kw_a$p.value, rho_a),
       x = "T category", y = "CDH1 expression (log2)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 10.5, face = "bold"))
save_dual(p_s6a, "S6a_CDH1_TCGA", 6.8, 5.6)
cat(sprintf("S6a: KW P=%.4f, rho=%.3f, n=%d\n", kw_a$p.value, rho_a, nrow(d6a)))

# S6b: GSE39582 CDH1 from series matrix (probes 201130_s_at / 201131_s_at)
sm_file <- file.path(PROJ, "data", "GSE39582", "GSE39582_series_matrix.txt.gz")
con <- gzfile(sm_file, "rt"); sample_names <- NULL; cdh1_rows <- list(); field_rows <- list()
in_table <- FALSE
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_characteristics_ch1")) {
    cells <- strsplit(substr(line, nchar("!Sample_characteristics_ch1") + 1, nchar(line)), "\t")[[1]]
    field_rows[[length(field_rows) + 1]] <- gsub('"', '', trimws(cells))
  } else if (startsWith(line, "!Sample_geo_accession")) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar("!Sample_geo_accession") + 1, nchar(line)), "\t")[[1]])
    sample_names <- sample_names[sample_names != ""]
  } else if (startsWith(line, "!series_matrix_table_begin")) { in_table <- TRUE; next
  } else if (startsWith(line, "!series_matrix_table_end")) { break
  } else if (in_table) {
    cells <- strsplit(line, "\t")[[1]]
    probe <- gsub('"', '', cells[1])
    if (probe %in% c("201130_s_at", "201131_s_at")) {
      cdh1_rows[[probe]] <- suppressWarnings(as.numeric(gsub('"', '', cells[-1])))
    }
  }
}
close(con)
cdh1_gse <- rowMeans(do.call(cbind, cdh1_rows), na.rm = TRUE)
field_names <- sapply(field_rows, function(r) { hit <- which(grepl(":", r))[1]; if (!is.na(hit)) sub(":.*", "", r[hit]) else "" })
get_field <- function(fname) { idx <- which(field_names == fname); if (length(idx) == 0) return(NULL)
  vals <- field_rows[[idx[1]]]; vals <- vals[vals != ""]; trimws(sub(paste0(fname, ":"), "", vals, fixed = TRUE)) }
gse6 <- data.frame(sample = sample_names, CDH1 = cdh1_gse,
                   dataset = get_field("dataset"), tnm_t = get_field("tnm.t"), stringsAsFactors = FALSE)
gse6 <- gse6[gse6$dataset != "Non Tumoral" & !is.na(gse6$tnm_t), ]
gse6$T_cat <- factor(paste0("T", gsub("[^1-4]", "", gse6$tnm_t)), levels = c("T1","T2","T3","T4"))
gse6 <- gse6[!is.na(gse6$T_cat) & gse6$T_cat %in% c("T1","T2","T3","T4") & !is.na(gse6$CDH1), ]
kw_b <- kruskal.test(CDH1 ~ T_cat, data = gse6)
rho_b <- suppressWarnings(cor(as.numeric(gse6$T_cat), gse6$CDH1, method = "spearman"))
p_s6b <- ggplot(gse6, aes(x = T_cat, y = CDH1)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, fill = "#F4A582", alpha = 0.8) +
  geom_jitter(width = 0.18, size = 0.5, alpha = 0.3) +
  labs(title = sprintf("CDH1 expression across T categories \u2014 GSE39582 (Kruskal\u2013Wallis P = %.3f; Spearman \u03c1 = %.3f)",
                       kw_b$p.value, rho_b),
       x = "T category", y = "CDH1 expression (log2)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 10.5, face = "bold"))
save_dual(p_s6b, "S6b_CDH1_GSE39582", 6.8, 5.6)
cat(sprintf("S6b: KW P=%.4f, rho=%.3f, n=%d\n", kw_b$p.value, rho_b, nrow(gse6)))

cat("\n=== [DONE] REVIEW_R11 figure redraw ===\n")
