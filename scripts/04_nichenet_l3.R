#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260714)

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(nichenetr)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
input_dir <- file.path(task_root, "results", "L1_single_cell")
output_dir <- file.path(task_root, "results", "L3_NicheNet")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

pseudobulk_file <- file.path(input_dir, "epithelial_pseudobulk_counts.tsv.gz")
sender_file <- file.path(input_dir, "strict_sender_gene_expression.tsv.gz")
burden_file <- file.path(input_dir, "patient_sender_burden.csv")
stopifnot(file.exists(pseudobulk_file), file.exists(sender_file), file.exists(burden_file))

resource_root <- Sys.getenv("FAP_NICHENET_RESOURCES", file.path(task_root, "inputs", "nichenet_resources"))
ligand_target_file <- file.path(resource_root, "ligand_target_matrix_nsga2r_final.rds")
lr_network_file <- file.path(resource_root, "lr_network_human_21122021.rds")
stopifnot(file.exists(ligand_target_file), file.exists(lr_network_file))

message("Loading epithelial pseudobulk counts")
count_table <- fread(pseudobulk_file)
genes <- count_table[[1]]
counts <- as.matrix(count_table[, -1])
storage.mode(counts) <- "integer"
rownames(counts) <- genes

burden <- fread(burden_file)
eligible <- burden[tumor_myofibroblasts >= 20]
eligible[, group := ifelse(strict_sender_fraction > median(strict_sender_fraction), "High", "Low")]
eligible[, group := factor(group, levels = c("Low", "High"))]
if (nrow(eligible) < 10 || length(unique(eligible$group)) != 2) stop("Insufficient eligible patients")
fwrite(eligible, file.path(output_dir, "nichenet_patient_groups.csv"))

counts <- counts[, eligible$patient, drop = FALSE]
dge <- DGEList(counts = counts, group = eligible$group)
keep <- filterByExpr(dge, group = eligible$group, min.count = 10)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
design <- model.matrix(~ group, data = eligible)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
test <- glmQLFTest(fit, coef = "groupHigh")
de_table <- as.data.table(topTags(test, n = Inf, sort.by = "PValue")$table, keep.rownames = "gene")
fwrite(de_table, file.path(output_dir, "epithelial_pseudobulk_high_vs_low_DE.csv"))

positive_de <- de_table[logFC > 0]
geneset_oi <- head(positive_de[order(PValue)]$gene, 200)
background_expressed_genes <- rownames(dge)
de_summary <- data.table(
  eligible_patients = nrow(eligible),
  high_patients = sum(eligible$group == "High"),
  low_patients = sum(eligible$group == "Low"),
  tested_genes = nrow(de_table),
  fdr_lt_0.05_positive = sum(de_table$FDR < 0.05 & de_table$logFC > 0),
  fdr_lt_0.10_positive = sum(de_table$FDR < 0.10 & de_table$logFC > 0),
  nichenet_top_positive_genes = length(geneset_oi)
)
fwrite(de_summary, file.path(output_dir, "epithelial_DE_summary.csv"))
fwrite(data.table(gene = geneset_oi), file.path(output_dir, "nichenet_receiver_gene_set_top200.csv"))

message("Loading NicheNet prior networks")
ligand_target_matrix <- readRDS(ligand_target_file)
lr_network <- as.data.table(readRDS(lr_network_file))
sender_expression <- fread(sender_file)

expressed_sender_genes <- sender_expression[positive_percent >= 10, gene]
expressed_sender_ligands <- intersect(unique(lr_network$from), expressed_sender_genes)
expressed_receiver_receptors <- intersect(unique(lr_network$to), background_expressed_genes)
potential_ligands <- unique(
  lr_network[from %in% expressed_sender_ligands & to %in% expressed_receiver_receptors, from]
)
potential_ligands <- intersect(potential_ligands, colnames(ligand_target_matrix))
if (length(potential_ligands) < 10) stop("Too few candidate ligands for NicheNet")

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
)
ligand_activities <- as.data.table(ligand_activities)
setorder(ligand_activities, -pearson)
ligand_activities[, rank := seq_len(.N)]
ligand_activities <- merge(
  ligand_activities,
  sender_expression[, .(test_ligand = gene, sender_mean_count = mean_count_per_cell,
                         sender_positive_percent = positive_percent)],
  by = "test_ligand",
  all.x = TRUE
)
setorder(ligand_activities, rank)

cellchat_file <- file.path(getwd(), "work", "cellchat_reanalysis", "results",
                           "cellchat_FAPhigh_to_epithelial_unweighted.csv")
if (file.exists(cellchat_file)) {
  cellchat <- fread(cellchat_file)
  cellchat_ligand_column <- if ("ligand" %in% names(cellchat)) "ligand" else "interaction_name"
  if (cellchat_ligand_column == "ligand") {
    direct_ligands <- unique(cellchat$ligand)
  } else {
    direct_ligands <- unique(sub("_.*$", "", cellchat$interaction_name))
  }
  ligand_activities[, cellchat_direct_support := test_ligand %in% direct_ligands]
}

fwrite(ligand_activities, file.path(output_dir, "nichenet_ligand_activity_dataset_specific.csv"))
fwrite(data.table(ligand = potential_ligands), file.path(output_dir, "nichenet_potential_ligands.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Dataset-specific NicheNet analysis complete: ", output_dir)
