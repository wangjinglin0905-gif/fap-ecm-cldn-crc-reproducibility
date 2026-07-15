#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260714)

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(readxl)
  library(Seurat)
  library(harmony)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
input_dir <- file.path(task_root, "inputs", "harmony")
output_dir <- file.path(task_root, "results", "L1_Harmony")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gse_matrix_file <- file.path(input_dir, "gse_harmony_counts.mtx")
gse_gene_file <- file.path(input_dir, "gse_harmony_genes.tsv")
gse_cell_file <- file.path(input_dir, "gse_harmony_cells.tsv")
qi_matrix_file <- Sys.getenv("FAP_QI_MATRIX", file.path(task_root, "inputs", "41467_2022_29366_MOESM6_ESM.gz"))
qi_metadata_file <- Sys.getenv("FAP_QI_METADATA", file.path(task_root, "inputs", "41467_2022_29366_MOESM7_ESM.xlsx"))
gse_metadata_file <- Sys.getenv("FAP_GSE132465_METADATA", file.path(task_root, "inputs", "GSE132465_cell_annotation.txt.gz"))
stopifnot(file.exists(gse_matrix_file), file.exists(qi_matrix_file), file.exists(qi_metadata_file), file.exists(gse_metadata_file))

message("Loading selected GSE132465 matrix")
gse_counts <- readMM(gse_matrix_file)
rownames(gse_counts) <- readLines(gse_gene_file)
colnames(gse_counts) <- readLines(gse_cell_file)
gse_counts <- as(gse_counts, "dgCMatrix")

gse_meta_raw <- fread(gse_metadata_file)
gse_meta_raw <- gse_meta_raw[match(colnames(gse_counts), Index)]
stopifnot(identical(gse_meta_raw$Index, colnames(gse_counts)))
gse_broad_map <- c(
  "B cells" = "B",
  "Epithelial cells" = "Epithelial",
  "Mast cells" = "Myeloid",
  "Myeloids" = "Myeloid",
  "Stromal cells" = "Stroma",
  "T cells" = "T"
)
gse_meta <- data.frame(
  row.names = gse_meta_raw$Index,
  dataset = "GSE132465",
  patient = gse_meta_raw$Patient,
  tissue = gse_meta_raw$Class,
  broad_group = unname(gse_broad_map[gse_meta_raw$Cell_type]),
  source_annotation = gse_meta_raw$Cell_subtype
)

message("Loading Qi et al. matrix and metadata")
qi_counts_full <- readRDS(gzfile(qi_matrix_file))
common_genes <- intersect(rownames(gse_counts), rownames(qi_counts_full))
gse_counts <- gse_counts[common_genes, , drop = FALSE]
qi_counts <- qi_counts_full[common_genes, , drop = FALSE]
rm(qi_counts_full)
gc()

qi_meta_raw <- as.data.table(read_excel(qi_metadata_file))
qi_meta_raw <- qi_meta_raw[match(colnames(qi_counts), Barcode)]
stopifnot(identical(qi_meta_raw$Barcode, colnames(qi_counts)))
qi_meta <- data.frame(
  row.names = qi_meta_raw$Barcode,
  dataset = "Qi2022",
  patient = qi_meta_raw$PatientID,
  tissue = qi_meta_raw$Tissues,
  broad_group = qi_meta_raw$MainTypes,
  source_annotation = qi_meta_raw$`Cell Types`
)

message("Creating merged Seurat object: ", ncol(gse_counts) + ncol(qi_counts), " cells")
gse_object <- CreateSeuratObject(counts = gse_counts, meta.data = gse_meta, project = "GSE132465")
qi_object <- CreateSeuratObject(counts = qi_counts, meta.data = qi_meta, project = "Qi2022")
combined <- merge(gse_object, y = qi_object, merge.data = FALSE)
rm(gse_object, qi_object, gse_counts, qi_counts)
gc()

combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 30, verbose = FALSE)

message("Running Harmony with theta=2")
combined <- RunHarmony(
  combined,
  group.by.vars = "dataset",
  reduction.use = "pca",
  dims.use = 1:30,
  theta = 2,
  max_iter = 20,
  early_stop = TRUE,
  ncores = 4,
  verbose = TRUE
)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, seed.use = 20260714, verbose = FALSE)

embedding <- as.data.table(Embeddings(combined, reduction = "umap"), keep.rownames = "cell_id")
metadata <- as.data.table(combined[[]], keep.rownames = "cell_id")
embedding <- merge(embedding, metadata[, .(cell_id, dataset, patient, tissue, broad_group, source_annotation)],
                   by = "cell_id", all.x = TRUE)
fwrite(embedding, file.path(output_dir, "harmony_umap_embeddings.csv.gz"))

dataset_counts <- embedding[, .N, by = .(dataset, broad_group)]
fwrite(dataset_counts, file.path(output_dir, "harmony_cell_counts.csv"))
qc <- data.table(
  total_cells = nrow(embedding),
  gse132465_cells = sum(embedding$dataset == "GSE132465"),
  qi2022_cells = sum(embedding$dataset == "Qi2022"),
  common_genes = length(common_genes),
  variable_features = length(VariableFeatures(combined)),
  harmony_dimensions = ncol(Embeddings(combined, reduction = "harmony"))
)
fwrite(qc, file.path(output_dir, "harmony_qc.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Harmony integration complete: ", output_dir)
