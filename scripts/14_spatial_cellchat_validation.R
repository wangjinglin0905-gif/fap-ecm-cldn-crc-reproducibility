options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(CellChat)
  library(data.table)
  library(jsonlite)
  library(Matrix)
})

input_dir <- file.path(
  "work", "reproducibility", "inputs", "valdeolivas_spatial", "extracted"
)
annotation_dir <- file.path(input_dir, "Pathology_SpotAnnotations")
output_dir <- file.path(
  "work", "reproducibility", "results", "L3_spatial_CellChat"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- c(
  "SN048_A121573_Rep1",
  "SN048_A121573_Rep2",
  "SN048_A416371_Rep1",
  "SN048_A416371_Rep2",
  "SN123_A551763_Rep1",
  "SN124_A551763_Rep2",
  "SN123_A595688_Rep1",
  "SN124_A595688_Rep2",
  "SN123_A798015_Rep1",
  "SN124_A798015_Rep2",
  "SN123_A938797_Rep1_X",
  "SN124_A938797_Rep2",
  "SN84_A120838_Rep1",
  "SN84_A120838_Rep2"
)

focused_interactions <- c(
  "FN1_CD44",
  "COL1A1_CD44",
  "COL1A2_CD44",
  "COL1A1_SDC4",
  "COL1A2_SDC4",
  "FN1_SDC4"
)

classify_pathology <- function(annotation) {
  normalized <- tolower(trimws(annotation))
  category <- rep("Other_tissue", length(normalized))
  category[normalized == "" | is.na(normalized) | grepl("exclude", normalized)] <- "Exclude"
  category[grepl("^tumor$", normalized)] <- "Tumor"
  category[grepl("tumor&stroma", normalized)] <- "Tumor_stroma"
  category[grepl(
    "stroma_fibroblastic|connective tissue_[23]_fibroblastic",
    normalized
  )] <- "Fibroblastic_stroma"
  category
}

read_spatial_sample <- function(sample_id) {
  sample_dir <- file.path(input_dir, sample_id)
  matrix_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  features <- fread(
    file.path(matrix_dir, "features.tsv.gz"),
    header = FALSE,
    col.names = c("ensembl_id", "gene", "feature_type")
  )
  barcodes <- fread(
    file.path(matrix_dir, "barcodes.tsv.gz"),
    header = FALSE,
    col.names = "barcode"
  )$barcode
  counts <- readMM(gzfile(file.path(matrix_dir, "matrix.mtx.gz")))
  rownames(counts) <- make.unique(features$gene)
  colnames(counts) <- barcodes

  positions <- fread(
    file.path(sample_dir, "spatial", "tissue_positions_list.csv"),
    header = FALSE,
    col.names = c(
      "barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"
    )
  )
  positions <- positions[in_tissue == 1L & barcode %in% barcodes]
  setkey(positions, barcode)
  positions <- positions[barcodes]

  annotation_path <- file.path(
    annotation_dir,
    paste0("Pathologist_Annotations_", sample_id, ".csv")
  )
  annotations <- fread(annotation_path)
  setnames(annotations, names(annotations)[1:2], c("barcode", "annotation"))
  setkey(annotations, barcode)
  positions[, annotation := annotations[barcode, annotation]]
  positions[, pathology_category := classify_pathology(annotation)]

  retained_groups <- c("Fibroblastic_stroma", "Tumor", "Tumor_stroma")
  retained <- positions$pathology_category %in% retained_groups
  positions <- positions[retained]
  counts <- counts[, retained, drop = FALSE]
  total_counts <- Matrix::colSums(counts)
  normalized <- log1p(Matrix::t(Matrix::t(counts) / total_counts) * 1e4)

  scalefactors <- fromJSON(
    file.path(sample_dir, "spatial", "scalefactors_json.json")
  )
  ratio <- 55 / scalefactors$spot_diameter_fullres
  spatial_factors <- data.frame(ratio = ratio, tol = 27.5)
  rownames(spatial_factors) <- sample_id

  metadata <- data.frame(
    labels = factor(positions$pathology_category, levels = retained_groups),
    row.names = positions$barcode
  )
  coordinates <- as.matrix(positions[, .(pixel_col, pixel_row)])
  rownames(coordinates) <- positions$barcode

  list(
    normalized = normalized,
    metadata = metadata,
    coordinates = coordinates,
    spatial_factors = spatial_factors,
    group_counts = table(metadata$labels)
  )
}

database <- CellChatDB.human
interaction_table <- database[["interaction"]]
focused_table <- interaction_table[
  interaction_table[["interaction_name"]] %in% focused_interactions,
  ,
  drop = FALSE
]

communication_results <- list()
status_results <- list()
result_index <- 1L

for (sample_id in sample_ids) {
  sample_result <- tryCatch({
    sample_data <- read_spatial_sample(sample_id)
    group_counts <- sample_data$group_counts
    sender_spots <- unname(group_counts["Fibroblastic_stroma"])
    receiver_spots <- unname(group_counts["Tumor"])
    if (is.na(sender_spots)) sender_spots <- 0L
    if (is.na(receiver_spots)) receiver_spots <- 0L
    if (sender_spots < 20L || receiver_spots < 20L) {
      stop(
        "Insufficient prespecified sender or receiver spots: ",
        sender_spots,
        " / ",
        receiver_spots
      )
    }

    cellchat <- createCellChat(
      object = sample_data$normalized,
      meta = sample_data$metadata,
      group.by = "labels",
      datatype = "spatial",
      coordinates = sample_data$coordinates,
      spatial.factors = sample_data$spatial_factors,
      do.sparse = TRUE
    )
    cellchat@DB <- database
    cellchat <- subsetData(cellchat)
    cellchat@LR$LRsig <- focused_table
    cellchat <- computeCommunProb(
      cellchat,
      type = "truncatedMean",
      trim = 0.1,
      LR.use = focused_table,
      raw.use = TRUE,
      population.size = FALSE,
      distance.use = TRUE,
      interaction.range = 200,
      scale.distance = 0.012,
      k.min = 5,
      contact.dependent = TRUE,
      contact.range = 110,
      nboot = 100,
      seed.use = 20260714L + result_index
    )
    cellchat <- filterCommunication(cellchat, min.cells = 20)
    communication <- as.data.table(subsetCommunication(cellchat, thresh = 1))
    communication[, sample_id := sample_id]
    communication[, patient_id := sub(".*_(A[0-9]+)_.*", "\\1", sample_id)]
    communication[, replicate := sub(".*_(Rep[12]).*", "\\1", sample_id)]
    communication[, sender_spots := sender_spots]
    communication[, receiver_spots := receiver_spots]
    communication_results[[result_index]] <- communication
    status_results[[result_index]] <- data.table(
      sample_id = sample_id,
      patient_id = sub(".*_(A[0-9]+)_.*", "\\1", sample_id),
      replicate = sub(".*_(Rep[12]).*", "\\1", sample_id),
      sender_spots = sender_spots,
      receiver_spots = receiver_spots,
      status = "completed",
      error = NA_character_
    )
    message("Completed spatial CellChat: ", sample_id)
    TRUE
  }, error = function(condition) {
    status_results[[result_index]] <<- data.table(
      sample_id = sample_id,
      patient_id = sub(".*_(A[0-9]+)_.*", "\\1", sample_id),
      replicate = sub(".*_(Rep[12]).*", "\\1", sample_id),
      sender_spots = NA_integer_,
      receiver_spots = NA_integer_,
      status = "not_analyzed",
      error = conditionMessage(condition)
    )
    message("Skipped spatial CellChat: ", sample_id, " - ", conditionMessage(condition))
    FALSE
  })
  result_index <- result_index + 1L
}

status <- rbindlist(status_results, fill = TRUE)
fwrite(status, file.path(output_dir, "spatial_cellchat_status.csv"))

if (length(communication_results) > 0L) {
  communications <- rbindlist(communication_results, fill = TRUE)
  communications[, prespecified_direction :=
    source == "Fibroblastic_stroma" & target == "Tumor"]
  fwrite(
    communications,
    file.path(output_dir, "spatial_cellchat_focused_pairs_all.csv")
  )

  prespecified <- communications[prespecified_direction == TRUE]
  prespecified[, fdr_bh := p.adjust(pval, method = "BH"), by = sample_id]
  prespecified[, significant := fdr_bh < 0.05 & prob > 0]
  sample_pair_summary <- prespecified[, .(
    probability = max(prob, na.rm = TRUE),
    p_value = min(pval, na.rm = TRUE),
    fdr_bh = min(fdr_bh, na.rm = TRUE),
    significant = any(significant)
  ), by = .(
    sample_id,
    patient_id,
    replicate,
    interaction_name,
    interaction_name_2,
    pathway_name
  )]
  fwrite(
    sample_pair_summary,
    file.path(output_dir, "spatial_cellchat_prespecified_direction.csv")
  )

  patient_pair_summary <- sample_pair_summary[, .(
    analyzed_replicates = .N,
    significant_replicates = sum(significant),
    both_replicates_significant = .N == 2L && all(significant),
    any_replicate_significant = any(significant),
    mean_probability = mean(probability)
  ), by = .(patient_id, interaction_name, interaction_name_2, pathway_name)]
  fwrite(
    patient_pair_summary,
    file.path(output_dir, "spatial_cellchat_patient_pair_summary.csv")
  )

  prevalence <- patient_pair_summary[, .(
    eligible_patients = uniqueN(patient_id),
    both_replicates_patients = sum(both_replicates_significant),
    any_replicate_patients = sum(any_replicate_significant),
    median_probability = median(mean_probability)
  ), by = .(interaction_name, interaction_name_2, pathway_name)]
  fwrite(
    prevalence,
    file.path(output_dir, "spatial_cellchat_pair_prevalence.csv")
  )
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
