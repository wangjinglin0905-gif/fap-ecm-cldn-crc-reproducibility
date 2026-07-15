suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(future)
})

set.seed(42)
plan("sequential")

input_dir <- "work/cellchat_reanalysis/input"
output_dir <- "work/cellchat_reanalysis/results"

counts <- as(readMM(file.path(input_dir, "counts.mtx")), "dgCMatrix")
genes <- readLines(file.path(input_dir, "genes.tsv"), encoding = "UTF-8")
meta <- read.delim(file.path(input_dir, "metadata.tsv"), row.names = 1, check.names = FALSE)
rownames(counts) <- genes
colnames(counts) <- rownames(meta)

keep <- meta$patient != "SMC20"
counts <- counts[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]
meta <- meta[colnames(counts), , drop = FALSE]
meta$samples <- factor(meta$patient)

library_sizes <- Matrix::colSums(counts)
normalized <- counts %*% Diagonal(x = 10000 / library_sizes)
normalized@x <- log1p(normalized@x)

data(CellChatDB.human)
object <- createCellChat(normalized, meta = meta, group.by = "group")
object@DB <- CellChatDB.human
object <- subsetData(object)
object <- identifyOverExpressedGenes(object, do.fast = FALSE)
object <- identifyOverExpressedInteractions(object)
object <- computeCommunProb(
  object,
  type = "triMean",
  population.size = FALSE,
  raw.use = TRUE,
  nboot = 100,
  seed.use = 42
)
object <- filterCommunication(object, min.cells = 10)
object <- computeCommunProbPathway(object)
object <- aggregateNet(object)

pairs <- subsetCommunication(
  object,
  sources.use = "FAP_high_myofibroblast",
  targets.use = "Tumor_epithelial",
  thresh = 0.05
)
pairs <- pairs[order(pairs$prob, decreasing = TRUE), ]
write.csv(pairs, file.path(output_dir, "cellchat_FAPhigh_to_epithelial_significant_excluding_SMC20.csv"), row.names = FALSE)
saveRDS(object, file.path(output_dir, "cellchat_object_excluding_SMC20.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_excluding_SMC20.txt"))

cat("Cells retained:", ncol(counts), "\n")
cat("Sender cells retained:", sum(meta$group == "FAP_high_myofibroblast"), "\n")
cat("Significant pairs:", nrow(pairs), "\n")
print(head(pairs[, c("ligand", "receptor", "prob", "pval", "pathway_name")], 15))
