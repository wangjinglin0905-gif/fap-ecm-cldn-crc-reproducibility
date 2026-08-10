options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/")
r_user <- file.path(root, "runtime", "r_user")
dir.create(r_user, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_USER = r_user, XDG_CACHE_HOME = file.path(r_user, "cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(progeny)
  library(msigdbr)
})

set.seed(2026080604)

default_source <- file.path(root, "data", "raw", "GSE132465")
source_root <- Sys.getenv("GSE132465_SOURCE_ROOT", unset = default_source)
seurat_file <- file.path(source_root, "results", "GSE132465_seurat.rds")
annotation_file <- file.path(
  source_root, "data", "GSE132465",
  "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
)
stopifnot(all(file.exists(c(seurat_file, annotation_file))))

out_dir <- file.path(root, "results", "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
score_genes <- unique(c(FAP13, matrix4))
fibroblast_subtypes <- c("Myofibroblasts", "Stromal 1", "Stromal 2",
                         "Stromal 3")

candidate_ligands <- c("TGFB1", "TGFB2", "TGFB3", "INHBA", "CXCL12", "IL6",
                       "IL11", "TNF", "WNT2", "WNT5A", "HGF", "VEGFA")
candidate_receptors <- c("TGFBR1", "TGFBR2", "ACVR1B", "ACVR2A", "CXCR4",
                         "IL6R", "IL6ST", "IL11RA", "TNFRSF1A", "TNFRSF1B",
                         "FZD2", "FZD6", "FZD7", "MET", "KDR")

bootstrap_spearman <- function(x, y, seed, reps = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(seed)
  estimates <- replicate(reps, {
    index <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(n = n, rho = observed,
    ci_low = unname(quantile(estimates, 0.025, type = 6, na.rm = TRUE)),
    ci_high = unname(quantile(estimates, 0.975, type = 6, na.rm = TRUE)),
    p_value = test$p.value)
}

bootstrap_partial_spearman <- function(x, y, covariate, seed, reps = 5000L) {
  covariate <- as.matrix(covariate)
  keep <- is.finite(x) & is.finite(y) &
    apply(covariate, 1, function(row) all(is.finite(row)))
  x <- x[keep]
  y <- y[keep]
  covariate <- covariate[keep, , drop = FALSE]
  n <- length(x)
  estimate <- function(index) {
    rx <- rank(x[index])
    ry <- rank(y[index])
    ranked_covariates <- apply(covariate[index, , drop = FALSE], 2, rank)
    ranked_covariates <- as.data.frame(ranked_covariates)
    suppressWarnings(cor(residuals(lm(rx ~ ., data = ranked_covariates)),
                         residuals(lm(ry ~ ., data = ranked_covariates))))
  }
  observed <- estimate(seq_len(n))
  set.seed(seed)
  estimates <- replicate(reps, estimate(sample.int(n, n, replace = TRUE)))
  degrees_freedom <- n - ncol(covariate) - 2L
  statistic <- observed * sqrt(degrees_freedom / (1 - observed^2))
  p_value <- 2 * pt(abs(statistic), df = degrees_freedom,
                    lower.tail = FALSE)
  c(partial_n = n, partial_rho = observed,
    partial_ci_low = unname(quantile(estimates, 0.025, type = 6,
                                     na.rm = TRUE)),
    partial_ci_high = unname(quantile(estimates, 0.975, type = 6,
                                      na.rm = TRUE)),
    partial_p_value = p_value)
}

normalise_cpm <- function(counts) {
  dge <- normLibSizes(DGEList(counts = counts), method = "TMM")
  cpm(dge, log = TRUE, prior.count = 0.25)
}

score_zmean <- function(expression, genes) {
  genes <- intersect(genes, rownames(expression))
  if (length(genes) < 3L) return(rep(NA_real_, ncol(expression)))
  z <- t(scale(t(expression[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

make_indicator <- function(group, levels) {
  keep <- !is.na(group)
  sparseMatrix(i = which(keep), j = match(group[keep], levels), x = 1,
               dims = c(length(group), length(levels)),
               dimnames = list(NULL, levels))
}

make_design <- function(predictor, covariate = NULL) {
  if (is.null(covariate)) {
    model.matrix(~ scale(predictor))
  } else {
    covariate <- as.data.frame(covariate)
    model.matrix(~ . + scale(predictor), data = covariate)
  }
}

filter_gene_universe <- function(counts, predictor, covariate) {
  dge <- DGEList(counts = counts)
  design <- make_design(predictor, covariate)
  keep <- filterByExpr(dge, design = design, min.count = 5,
                       min.total.count = 15)
  rownames(counts)[keep]
}

fit_gene_models <- function(counts, predictor, label, covariate = NULL,
                            keep_genes = NULL) {
  dge <- DGEList(counts = counts)
  design <- make_design(predictor, covariate)
  if (is.null(keep_genes)) {
    keep <- filterByExpr(dge, design = design, min.count = 5,
                         min.total.count = 15)
  } else {
    keep <- rownames(dge) %in% keep_genes
  }
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- normLibSizes(dge, method = "TMM")
  voom_fit <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(voom_fit, design), robust = TRUE)
  coefficient <- ncol(design)
  table <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  table$gene <- rownames(table)
  table$compartment <- label
  table <- table[, c("compartment", "gene", "logFC", "AveExpr", "t", "P.Value",
                     "adj.P.Val", "B")]
  list(table = table, ranks = setNames(table$t, table$gene),
       logcpm = voom_fit$E, genes_tested = nrow(table))
}

object <- readRDS(seurat_file)
counts <- LayerData(object[["RNA"]], layer = "counts")
metadata <- object[[]]
annotation <- read.delim(gzfile(annotation_file), check.names = FALSE)
stopifnot(inherits(counts, "sparseMatrix"))
stopifnot(identical(colnames(counts), rownames(metadata)))
stopifnot(identical(colnames(counts), annotation$Index))

is_fibroblast <- metadata$Cell_type == "Stromal cells" &
  metadata$Cell_subtype %in% fibroblast_subtypes
is_epithelial <- metadata$Cell_type == "Epithelial cells"
compartment <- ifelse(is_fibroblast, "fibroblast_lineage",
                      ifelse(is_epithelial, "epithelial", NA_character_))
group <- ifelse(is.na(compartment), NA_character_,
                paste(metadata$Class, compartment, metadata$Patient, sep = "::"))
group_levels <- unique(group[!is.na(group)])
indicator <- make_indicator(group, group_levels)
pseudobulk_counts <- counts %*% indicator

parts <- do.call(rbind, strsplit(group_levels, "::", fixed = TRUE))
group_manifest <- data.frame(group = group_levels, class = parts[, 1],
                             compartment = parts[, 2], patient = parts[, 3],
                             cells = as.integer(Matrix::colSums(indicator)),
                             library_size = as.numeric(Matrix::colSums(
                               pseudobulk_counts)))
tumour <- group_manifest[group_manifest$class == "Tumor", , drop = FALSE]
wide_counts <- reshape(tumour[, c("patient", "compartment", "cells")],
                       idvar = "patient", timevar = "compartment",
                       direction = "wide")
names(wide_counts) <- sub("^cells\\.", "n_", names(wide_counts))
eligible <- wide_counts$patient[
  wide_counts$n_fibroblast_lineage >= 20L & wide_counts$n_epithelial >= 20L
]
stopifnot(length(eligible) == 15L)

fib_groups <- paste("Tumor", "fibroblast_lineage", eligible, sep = "::")
epi_groups <- paste("Tumor", "epithelial", eligible, sep = "::")
fib_counts <- pseudobulk_counts[, fib_groups, drop = FALSE]
epi_counts <- pseudobulk_counts[, epi_groups, drop = FALSE]
colnames(fib_counts) <- eligible
colnames(epi_counts) <- eligible
fib_logcpm <- normalise_cpm(fib_counts)
epi_logcpm <- normalise_cpm(epi_counts)

fib_fap <- as.numeric(fib_logcpm["FAP", eligible])
patient_table <- data.frame(
  patient = eligible,
  n_fib = wide_counts$n_fibroblast_lineage[match(eligible, wide_counts$patient)],
  n_epi = wide_counts$n_epithelial[match(eligible, wide_counts$patient)],
  fib_FAP_logCPM = fib_fap,
  stringsAsFactors = FALSE
)
fib_cell_covariate <- data.frame(log_n_fib = log1p(patient_table$n_fib))
epi_cell_covariates <- data.frame(
  log_n_fib = log1p(patient_table$n_fib),
  log_n_epi = log1p(patient_table$n_epi)
)

fib_gene_universe <- filter_gene_universe(
  fib_counts, fib_fap, fib_cell_covariate
)
epi_gene_universe <- filter_gene_universe(
  epi_counts, fib_fap, epi_cell_covariates
)

fib_model <- fit_gene_models(
  fib_counts, fib_fap, "fibroblast_lineage", keep_genes = fib_gene_universe
)
epi_model <- fit_gene_models(
  epi_counts, fib_fap, "epithelial", keep_genes = epi_gene_universe
)
fib_model_adjusted <- fit_gene_models(
  fib_counts, fib_fap, "fibroblast_lineage", fib_cell_covariate,
  keep_genes = fib_gene_universe
)
epi_model_adjusted <- fit_gene_models(
  epi_counts, fib_fap, "epithelial", epi_cell_covariates,
  keep_genes = epi_gene_universe
)
gene_models <- rbind(fib_model$table, epi_model$table)
write.csv(gene_models,
          file.path(out_dir, "GSE132465_FAP_continuous_gene_models.csv"),
          row.names = FALSE)
gene_models_adjusted <- rbind(fib_model_adjusted$table,
                              epi_model_adjusted$table)
write.csv(
  gene_models_adjusted,
  file.path(out_dir, "GSE132465_FAP_cellcount_adjusted_gene_models.csv"),
  row.names = FALSE
)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
prespecified_names <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
)
hallmark <- hallmark[hallmark$gs_name %in% prespecified_names, , drop = FALSE]
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
pathways <- lapply(pathways, unique)
stopifnot(all(prespecified_names %in% names(pathways)))

run_fgsea <- function(ranks, compartment, omit_score_genes = FALSE) {
  ranks <- ranks[is.finite(ranks)]
  if (omit_score_genes) ranks <- ranks[!names(ranks) %in% score_genes]
  ranks <- sort(ranks, decreasing = TRUE)
  result <- fgseaMultilevel(pathways = pathways, stats = ranks,
                            minSize = 10, maxSize = 500,
                            eps = 0, nproc = 1)
  result <- as.data.frame(result)
  result$leadingEdge <- vapply(result$leadingEdge, paste, collapse = ";",
                               FUN.VALUE = character(1))
  result$compartment <- compartment
  result$score_genes_omitted <- omit_score_genes
  result
}

enrichment <- rbind(
  run_fgsea(fib_model$ranks, "fibroblast_lineage", FALSE),
  run_fgsea(fib_model$ranks, "fibroblast_lineage", TRUE),
  run_fgsea(epi_model$ranks, "epithelial", FALSE),
  run_fgsea(epi_model$ranks, "epithelial", TRUE)
)
write.csv(enrichment,
          file.path(out_dir, "GSE132465_FAP_continuous_hallmark_fgsea.csv"),
          row.names = FALSE)

enrichment_adjusted <- rbind(
  run_fgsea(fib_model_adjusted$ranks, "fibroblast_lineage", FALSE),
  run_fgsea(fib_model_adjusted$ranks, "fibroblast_lineage", TRUE),
  run_fgsea(epi_model_adjusted$ranks, "epithelial", FALSE),
  run_fgsea(epi_model_adjusted$ranks, "epithelial", TRUE)
)
write.csv(
  enrichment_adjusted,
  file.path(out_dir, "GSE132465_FAP_cellcount_adjusted_hallmark_fgsea.csv"),
  row.names = FALSE
)

hallmark_scores <- do.call(cbind, lapply(pathways, function(genes) {
  score_zmean(epi_logcpm, genes)
}))
colnames(hallmark_scores) <- names(pathways)
rownames(hallmark_scores) <- eligible

hallmark_correlations <- do.call(rbind, lapply(
  seq_len(ncol(hallmark_scores)), function(i) {
    data.frame(
      method = "Epithelial Hallmark mean-z",
      pathway = colnames(hallmark_scores)[i],
      t(bootstrap_spearman(fib_fap, hallmark_scores[, i],
                           seed = 2026080700L + i))
    )
  }
))
hallmark_correlations$fdr_bh_family <- p.adjust(hallmark_correlations$p_value,
                                                method = "BH")

progeny_primary <- progeny(epi_logcpm, scale = TRUE, organism = "Human",
                           top = 500, perm = 1, verbose = FALSE,
                           z_scores = FALSE)
progeny_100 <- progeny(epi_logcpm, scale = TRUE, organism = "Human",
                       top = 100, perm = 1, verbose = FALSE,
                       z_scores = FALSE)
progeny_1000 <- progeny(epi_logcpm, scale = TRUE, organism = "Human",
                        top = 1000, perm = 1, verbose = FALSE,
                        z_scores = FALSE)

score_progeny <- function(scores, top, seed_base) {
  scores <- as.matrix(scores)
  if (nrow(scores) != length(eligible) && ncol(scores) == length(eligible)) {
    scores <- t(scores)
  }
  stopifnot(nrow(scores) == length(eligible))
  stopifnot(setequal(rownames(scores), eligible))
  scores <- scores[eligible, , drop = FALSE]
  do.call(rbind, lapply(seq_len(ncol(scores)), function(i) {
    data.frame(method = paste0("PROGENy top", top),
               pathway = colnames(scores)[i],
               t(bootstrap_spearman(fib_fap, scores[, i],
                                    seed = seed_base + i)))
  }))
}

progeny_correlations <- rbind(
  score_progeny(progeny_primary, 500, 2026080800L),
  score_progeny(progeny_100, 100, 2026080900L),
  score_progeny(progeny_1000, 1000, 2026081000L)
)
progeny_correlations$fdr_bh_within_top <- ave(
  progeny_correlations$p_value, progeny_correlations$method,
  FUN = function(x) p.adjust(x, method = "BH")
)
progeny_correlations$fdr_bh_family <- progeny_correlations$fdr_bh_within_top
hallmark_correlations$fdr_bh_within_top <- NA_real_

pathway_correlations <- rbind(
  transform(hallmark_correlations, sensitivity = "primary"),
  transform(progeny_correlations,
            sensitivity = ifelse(method == "PROGENy top500", "primary",
                                 "footprint-size sensitivity"))
)
write.csv(pathway_correlations,
          file.path(out_dir, "GSE132465_fibFAP_epithelial_pathway_correlations.csv"),
          row.names = FALSE)

primary_progeny_for_partial <- as.matrix(progeny_primary)
if (nrow(primary_progeny_for_partial) != length(eligible)) {
  primary_progeny_for_partial <- t(primary_progeny_for_partial)
}
stopifnot(setequal(rownames(primary_progeny_for_partial), eligible))
primary_progeny_for_partial <- primary_progeny_for_partial[
  eligible, , drop = FALSE
]

partial_pathway_correlations <- rbind(
  do.call(rbind, lapply(seq_len(ncol(hallmark_scores)), function(i) {
    data.frame(
      method = "Epithelial Hallmark mean-z",
      pathway = colnames(hallmark_scores)[i],
      t(bootstrap_partial_spearman(
        fib_fap, hallmark_scores[, i], epi_cell_covariates,
        seed = 2026081300L + i
      ))
    )
  })),
  do.call(rbind, lapply(seq_len(ncol(primary_progeny_for_partial)),
                        function(i) {
    data.frame(
      method = "PROGENy top500",
      pathway = colnames(primary_progeny_for_partial)[i],
      t(bootstrap_partial_spearman(
        fib_fap, primary_progeny_for_partial[, i], epi_cell_covariates,
        seed = 2026081400L + i
      ))
    )
  }))
)
partial_pathway_correlations$fdr_bh_within_method <- ave(
  partial_pathway_correlations$partial_p_value,
  partial_pathway_correlations$method,
  FUN = function(x) p.adjust(x, method = "BH")
)
write.csv(
  partial_pathway_correlations,
  file.path(out_dir,
            "GSE132465_fibFAP_epithelial_pathway_partial_correlations.csv"),
  row.names = FALSE
)

patient_scores <- cbind(patient_table, hallmark_scores)
primary_progeny <- as.matrix(progeny_primary)
if (nrow(primary_progeny) != length(eligible)) primary_progeny <- t(primary_progeny)
stopifnot(setequal(rownames(primary_progeny), eligible))
primary_progeny <- primary_progeny[eligible, , drop = FALSE]
colnames(primary_progeny) <- paste0("PROGENy_", colnames(primary_progeny))
patient_scores <- cbind(patient_scores, primary_progeny)
write.csv(patient_scores,
          file.path(out_dir, "GSE132465_cross_compartment_patient_scores.csv"),
          row.names = FALSE)

extract_candidates <- function(model_table, adjusted_model_table, genes, family,
                               seed_base, expression, predictor, covariate) {
  present <- intersect(genes, rownames(expression))
  rows <- do.call(rbind, lapply(seq_along(present), function(i) {
    gene <- present[i]
    model_row <- model_table[match(gene, model_table$gene), , drop = FALSE]
    adjusted_row <- adjusted_model_table[
      match(gene, adjusted_model_table$gene), , drop = FALSE
    ]
    data.frame(family = family, gene = gene,
               member_of_FAP13_or_matrix4 = gene %in% score_genes,
               model_logFC_per_1SD_FAP = model_row$logFC,
               model_t = model_row$t,
               model_p = model_row$P.Value,
               model_fdr_genomewide = model_row$adj.P.Val,
               adjusted_model_logFC_per_1SD_FAP = adjusted_row$logFC,
               adjusted_model_t = adjusted_row$t,
               adjusted_model_p = adjusted_row$P.Value,
               adjusted_model_fdr_genomewide = adjusted_row$adj.P.Val,
               t(bootstrap_spearman(predictor,
                                    as.numeric(expression[gene, ]),
                                    seed = seed_base + i)),
               t(bootstrap_partial_spearman(
                 predictor, as.numeric(expression[gene, ]), covariate,
                 seed = seed_base + 100L + i
               )))
  }))
  rows$fdr_bh_candidate_family <- p.adjust(rows$p_value, method = "BH")
  rows$partial_fdr_bh_candidate_family <- p.adjust(
    rows$partial_p_value, method = "BH"
  )
  rows
}

candidates <- rbind(
  extract_candidates(fib_model$table, fib_model_adjusted$table,
                     candidate_ligands,
                     "fibroblast candidate ligand", 2026081100L,
                     fib_logcpm, fib_fap, fib_cell_covariate),
  extract_candidates(epi_model$table, epi_model_adjusted$table,
                     candidate_receptors,
                     "epithelial candidate receptor", 2026081200L,
                     epi_logcpm, fib_fap, epi_cell_covariates)
)
write.csv(candidates,
          file.path(out_dir, "GSE132465_candidate_ligand_receptor_screen.csv"),
          row.names = FALSE)

subtypes <- as.data.frame(table(Class = metadata$Class,
                                Cell_type = metadata$Cell_type,
                                Cell_subtype = metadata$Cell_subtype),
                          stringsAsFactors = FALSE)
subtypes <- subtypes[subtypes$Freq > 0, , drop = FALSE]
write.csv(subtypes,
          file.path(out_dir, "GSE132465_cell_subtype_counts.csv"),
          row.names = FALSE)

analysis_manifest <- data.frame(
  item = c("patients", "minimum_fibroblast_cells", "minimum_epithelial_cells",
           "fibroblast_gene_models", "epithelial_gene_models",
           "cellcount_adjusted_fibroblast_gene_models",
           "cellcount_adjusted_epithelial_gene_models",
           "gene_filtering", "epithelial_cellcount_covariates",
           "hallmark_sets", "PROGENy_primary_top", "score_gene_omission",
           "replication_unit", "causal_boundary", "MSigDB_release",
           "msigdbr_version", "MSigDB_archive_MD5", "R_version",
           "Seurat_RDS_MD5", "annotation_MD5"),
  value = c(length(eligible), 20, 20, fib_model$genes_tested,
            epi_model$genes_tested, fib_model_adjusted$genes_tested,
            epi_model_adjusted$genes_tested,
            "Common adjusted-design gene universe within each compartment",
            "log1p fibroblast and epithelial cell counts", length(pathways), 500,
            "Full and FAP13/matrix4-omitted fgsea results reported",
            "patient", "association and pathway prioritization only",
            "2026.1", as.character(packageVersion("msigdbr")),
            unname(tools::md5sum(file.path(
              tools::R_user_dir("msigdbr", "cache"),
              "msigdb.2026.1.zip"
            ))),
            R.version.string, unname(tools::md5sum(seurat_file)),
            unname(tools::md5sum(annotation_file)))
)
write.csv(analysis_manifest,
          file.path(out_dir, "GSE132465_cross_compartment_manifest.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "GSE132465_cross_compartment_sessionInfo.txt"))

cat("GSE132465 cross-compartment pathway screen complete.\n")
cat("Eligible patients:", length(eligible), "\n")
cat("\nPrimary pathway correlations:\n")
print(subset(pathway_correlations, sensitivity == "primary"), row.names = FALSE)
cat("\nPrimary epithelial fgsea results:\n")
print(subset(enrichment, compartment == "epithelial" &
               !score_genes_omitted)[, c("pathway", "NES", "pval", "padj")],
      row.names = FALSE)
cat("\nCandidate ligands and receptors:\n")
print(candidates[, c("family", "gene", "rho", "ci_low", "ci_high",
                     "fdr_bh_candidate_family", "partial_rho",
                     "partial_fdr_bh_candidate_family")], row.names = FALSE)
