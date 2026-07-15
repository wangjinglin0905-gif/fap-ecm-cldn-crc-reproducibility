suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(future)
})

set.seed(42)
plan("sequential")

input_dir <- "work/cellchat_reanalysis/input"
output_dir <- "work/cellchat_reanalysis/results"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readMM(file.path(input_dir, "counts.mtx"))
counts <- as(counts, "dgCMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), encoding = "UTF-8")
meta <- read.delim(
  file.path(input_dir, "metadata.tsv"),
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(counts) <- genes
colnames(counts) <- rownames(meta)
meta <- meta[colnames(counts), , drop = FALSE]
stopifnot(identical(colnames(counts), rownames(meta)))
meta$samples <- meta$patient

library_sizes <- Matrix::colSums(counts)
if (any(library_sizes == 0)) stop("Zero-library cells were detected")
normalized <- counts %*% Diagonal(x = 10000 / library_sizes)
normalized@x <- log1p(normalized@x)

data(CellChatDB.human)
base_object <- createCellChat(
  object = normalized,
  meta = meta,
  group.by = "group"
)
base_object@DB <- CellChatDB.human
base_object <- subsetData(base_object)
base_object <- identifyOverExpressedGenes(base_object, do.fast = FALSE)
base_object <- identifyOverExpressedInteractions(base_object)

run_inference <- function(object, population_size, label) {
  object <- computeCommunProb(
    object,
    type = "triMean",
    population.size = population_size,
    raw.use = TRUE,
    nboot = 100,
    seed.use = 42
  )
  object <- filterCommunication(object, min.cells = 10)
  object <- computeCommunProbPathway(object)
  object <- aggregateNet(object)

  all_pairs <- subsetCommunication(object, thresh = 1)
  direction_pairs <- subsetCommunication(
    object,
    sources.use = "FAP_high_myofibroblast",
    targets.use = "Tumor_epithelial",
    thresh = 1
  )
  direction_significant <- subsetCommunication(
    object,
    sources.use = "FAP_high_myofibroblast",
    targets.use = "Tumor_epithelial",
    thresh = 0.05
  )
  all_pairs <- all_pairs[order(all_pairs$prob, decreasing = TRUE), ]
  direction_pairs <- direction_pairs[order(direction_pairs$prob, decreasing = TRUE), ]
  direction_significant <- direction_significant[order(direction_significant$prob, decreasing = TRUE), ]
  write.csv(all_pairs, file.path(output_dir, paste0("cellchat_all_pairs_", label, ".csv")), row.names = FALSE)
  write.csv(direction_pairs, file.path(output_dir, paste0("cellchat_FAPhigh_to_epithelial_", label, ".csv")), row.names = FALSE)
  write.csv(direction_significant, file.path(output_dir, paste0("cellchat_FAPhigh_to_epithelial_significant_", label, ".csv")), row.names = FALSE)
  saveRDS(object, file.path(output_dir, paste0("cellchat_object_", label, ".rds")), compress = "xz")
  object
}

unweighted <- run_inference(base_object, FALSE, "unweighted")
weighted <- run_inference(base_object, TRUE, "population_weighted")

summarize_pairs <- function(path, label) {
  pairs <- read.csv(path, check.names = FALSE)
  pairs$rank <- seq_len(nrow(pairs))
  collagen <- grepl("^COL", pairs$ligand) | pairs$pathway_name == "COLLAGEN"
  fibronectin <- pairs$ligand == "FN1" | pairs$pathway_name == "FN1"
  tgfb <- grepl("^TGFB|^INHBA$", pairs$ligand) | pairs$pathway_name %in% c("TGFb", "ACTIVIN")
  ccl2 <- pairs$ligand == "CCL2"
  data.frame(
    analysis = label,
    pair_count = nrow(pairs),
    collagen_pair_count = sum(collagen),
    fibronectin_pair_count = sum(fibronectin),
    tgfb_or_inhba_pair_count = sum(tgfb),
    ccl2_pair_count = sum(ccl2),
    top_pair = if (nrow(pairs)) pairs$interaction_name[1] else NA_character_,
    top_probability = if (nrow(pairs)) pairs$prob[1] else NA_real_,
    stringsAsFactors = FALSE
  )
}

summary_table <- rbind(
  summarize_pairs(file.path(output_dir, "cellchat_FAPhigh_to_epithelial_significant_unweighted.csv"), "unweighted"),
  summarize_pairs(file.path(output_dir, "cellchat_FAPhigh_to_epithelial_significant_population_weighted.csv"), "population_weighted")
)
write.csv(summary_table, file.path(output_dir, "cellchat_reanalysis_summary.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
print(summary_table)
