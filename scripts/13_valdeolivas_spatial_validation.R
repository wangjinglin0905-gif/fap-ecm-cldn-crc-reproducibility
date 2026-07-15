options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

input_dir <- file.path(
  "work", "reproducibility", "inputs", "valdeolivas_spatial", "extracted"
)
annotation_dir <- file.path(input_dir, "Pathology_SpotAnnotations")
output_dir <- file.path(
  "work", "reproducibility", "results", "L2_Valdeolivas_spatial"
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

fap_caf_genes <- c(
  "FAP", "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "THY1", "PDPN",
  "TAGLN", "ACTA2", "MMP2", "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A"
)
de_ligand_genes <- setdiff(fap_caf_genes, c("TGFB1", "INHBA", "WNT2", "WNT5A"))
cldn_core_genes <- c("CLDN1", "CLDN2", "CLDN4")
epithelial_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19")
target_genes <- unique(c(fap_caf_genes, cldn_core_genes, epithelial_genes))

score_zmean <- function(expression_matrix, genes) {
  available <- intersect(genes, rownames(expression_matrix))
  if (length(available) == 0L) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }
  gene_matrix <- expression_matrix[available, , drop = FALSE]
  gene_sd <- apply(gene_matrix, 1L, sd, na.rm = TRUE)
  variable <- is.finite(gene_sd) & gene_sd > 0
  if (!any(variable)) {
    return(rep(NA_real_, ncol(expression_matrix)))
  }
  scaled <- t(scale(t(gene_matrix[variable, , drop = FALSE])))
  colMeans(scaled, na.rm = TRUE)
}

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
  category[grepl(
    "epithelium|glandular|non neo",
    normalized
  )] <- "Non_neoplastic_epithelium"
  category
}

make_hex_edges <- function(array_row, array_col) {
  coordinate_key <- paste(array_row, array_col, sep = ":")
  index_lookup <- setNames(seq_along(coordinate_key), coordinate_key)
  offsets <- matrix(
    c(-1L, -1L, -1L, 1L, 0L, -2L, 0L, 2L, 1L, -1L, 1L, 1L),
    ncol = 2L,
    byrow = TRUE
  )
  from <- integer()
  to <- integer()
  for (offset_index in seq_len(nrow(offsets))) {
    neighbor_key <- paste(
      array_row + offsets[offset_index, 1L],
      array_col + offsets[offset_index, 2L],
      sep = ":"
    )
    neighbor <- unname(index_lookup[neighbor_key])
    present <- !is.na(neighbor)
    from <- c(from, which(present))
    to <- c(to, neighbor[present])
  }
  data.table(from = from, to = as.integer(to))
}

standardize <- function(values) {
  deviation <- sd(values, na.rm = TRUE)
  if (!is.finite(deviation) || deviation == 0) {
    return(rep(NA_real_, length(values)))
  }
  (values - mean(values, na.rm = TRUE)) / deviation
}

bivariate_moran <- function(predictor, outcome, edges) {
  standardized_predictor <- standardize(predictor)
  standardized_outcome <- standardize(outcome)
  numerator <- sum(
    standardized_predictor[edges$from] * standardized_outcome[edges$to]
  )
  denominator <- sum(standardized_predictor^2)
  length(predictor) / nrow(edges) * numerator / denominator
}

permute_within_strata <- function(values, strata) {
  permuted <- values
  for (stratum in unique(strata)) {
    indices <- which(strata == stratum)
    permuted[indices] <- sample(values[indices], length(indices), replace = FALSE)
  }
  permuted
}

permutation_test <- function(
  predictor,
  outcome,
  edges,
  strata,
  permutations = 999L,
  seed = 20260714L
) {
  set.seed(seed)
  observed <- bivariate_moran(predictor, outcome, edges)
  null <- replicate(permutations, {
    permuted_outcome <- permute_within_strata(outcome, strata)
    bivariate_moran(predictor, permuted_outcome, edges)
  })
  p_value <- (sum(abs(null) >= abs(observed)) + 1) / (permutations + 1)
  c(
    statistic = observed,
    permutation_p = p_value,
    null_mean = mean(null),
    null_sd = sd(null)
  )
}

high_high_test <- function(
  fap_score,
  cldn_score,
  pathology_category,
  edges,
  permutations = 999L,
  seed = 20260714L
) {
  stromal_reference <- pathology_category %in% c("Fibroblastic_stroma", "Tumor_stroma")
  tumor_reference <- pathology_category %in% c("Tumor", "Tumor_stroma")
  if (sum(stromal_reference) < 20L || sum(tumor_reference) < 20L) {
    return(c(observed_edges = NA, expected_edges = NA, enrichment = NA, permutation_p = NA))
  }
  fap_threshold <- quantile(fap_score[stromal_reference], 0.75, na.rm = TRUE)
  cldn_threshold <- quantile(cldn_score[tumor_reference], 0.75, na.rm = TRUE)
  fap_high <- stromal_reference & fap_score >= fap_threshold
  cldn_high <- tumor_reference & cldn_score >= cldn_threshold
  observed <- sum(fap_high[edges$from] & cldn_high[edges$to])
  set.seed(seed)
  null <- replicate(permutations, {
    permuted_cldn <- permute_within_strata(cldn_high, pathology_category)
    sum(fap_high[edges$from] & permuted_cldn[edges$to])
  })
  expected <- mean(null)
  c(
    observed_edges = observed,
    expected_edges = expected,
    enrichment = if (expected > 0) observed / expected else NA_real_,
    permutation_p = (sum(null >= observed) + 1) / (permutations + 1)
  )
}

read_sample <- function(sample_id) {
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
  count_matrix <- readMM(gzfile(file.path(matrix_dir, "matrix.mtx.gz")))
  colnames(count_matrix) <- barcodes

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
  if (anyNA(positions$barcode)) {
    stop("Missing tissue positions for sample ", sample_id)
  }

  annotation_path <- file.path(
    annotation_dir,
    paste0("Pathologist_Annotations_", sample_id, ".csv")
  )
  annotations <- fread(annotation_path)
  setnames(annotations, names(annotations)[1:2], c("barcode", "annotation"))
  setkey(annotations, barcode)
  positions[, annotation := annotations[barcode, annotation]]
  positions[, pathology_category := classify_pathology(annotation)]

  gene_indices <- split(seq_len(nrow(features)), features$gene)
  selected_counts <- matrix(
    0,
    nrow = length(target_genes),
    ncol = ncol(count_matrix),
    dimnames = list(target_genes, barcodes)
  )
  for (gene in target_genes) {
    indices <- gene_indices[[gene]]
    if (!is.null(indices)) {
      selected_counts[gene, ] <- Matrix::colSums(count_matrix[indices, , drop = FALSE])
    }
  }
  total_umi <- Matrix::colSums(count_matrix)
  log_cpm <- log1p(sweep(selected_counts, 2L, total_umi, "/") * 1e6)
  list(positions = positions, log_cpm = log_cpm, total_umi = total_umi)
}

sample_metrics <- list()
spot_tables <- list()
metric_index <- 1L

for (sample_id in sample_ids) {
  sample_data <- read_sample(sample_id)
  positions <- sample_data$positions
  log_cpm <- sample_data$log_cpm
  positions[, log_total_umi := log1p(sample_data$total_umi)]
  positions[, FAP_log_cpm := as.numeric(log_cpm["FAP", ])]
  positions[, FAP_CAF := score_zmean(log_cpm, fap_caf_genes)]
  positions[, FAP_CAF_de_ligand := score_zmean(log_cpm, de_ligand_genes)]
  positions[, CLDN_core := score_zmean(log_cpm, cldn_core_genes)]
  positions[, epithelial_score := score_zmean(log_cpm, epithelial_genes)]
  positions[, sample_id := sample_id]

  analysis_spots <- positions[pathology_category != "Exclude"]
  edges <- make_hex_edges(analysis_spots$array_row, analysis_spots$array_col)
  if (nrow(edges) == 0L) {
    stop("No spatial edges for sample ", sample_id)
  }

  raw_test <- permutation_test(
    analysis_spots$FAP_CAF_de_ligand,
    analysis_spots$CLDN_core,
    edges,
    analysis_spots$pathology_category,
    seed = 20260714L + metric_index
  )

  adjusted_data <- data.frame(
    FAP_CAF_de_ligand = analysis_spots$FAP_CAF_de_ligand,
    CLDN_core = analysis_spots$CLDN_core,
    log_total_umi = analysis_spots$log_total_umi,
    epithelial_score = analysis_spots$epithelial_score,
    pathology_category = factor(analysis_spots$pathology_category)
  )
  predictor_model <- lm(
    FAP_CAF_de_ligand ~ log_total_umi + epithelial_score + pathology_category,
    data = adjusted_data
  )
  outcome_model <- lm(
    CLDN_core ~ log_total_umi + epithelial_score + pathology_category,
    data = adjusted_data
  )
  adjusted_test <- permutation_test(
    residuals(predictor_model),
    residuals(outcome_model),
    edges,
    analysis_spots$pathology_category,
    seed = 20261714L + metric_index
  )

  high_high <- high_high_test(
    analysis_spots$FAP_CAF_de_ligand,
    analysis_spots$CLDN_core,
    analysis_spots$pathology_category,
    edges,
    seed = 20262714L + metric_index
  )

  sample_metrics[[metric_index]] <- data.table(
    sample_id = sample_id,
    patient_id = sub(".*_(A[0-9]+)_.*", "\\1", sample_id),
    replicate = sub(".*_(Rep[12]).*", "\\1", sample_id),
    analyzed_spots = nrow(analysis_spots),
    directed_edges = nrow(edges),
    raw_bivariate_moran = raw_test["statistic"],
    raw_permutation_p = raw_test["permutation_p"],
    adjusted_bivariate_moran = adjusted_test["statistic"],
    adjusted_permutation_p = adjusted_test["permutation_p"],
    high_high_observed_edges = high_high["observed_edges"],
    high_high_expected_edges = high_high["expected_edges"],
    high_high_enrichment = high_high["enrichment"],
    high_high_permutation_p = high_high["permutation_p"]
  )
  spot_tables[[metric_index]] <- positions[, .(
    sample_id,
    barcode,
    array_row,
    array_col,
    pixel_row,
    pixel_col,
    annotation,
    pathology_category,
    log_total_umi,
    FAP_log_cpm,
    FAP_CAF,
    FAP_CAF_de_ligand,
    CLDN_core,
    epithelial_score
  )]
  message("Processed ", sample_id, ": ", nrow(analysis_spots), " annotated spots")
  metric_index <- metric_index + 1L
}

metrics <- rbindlist(sample_metrics)
metrics[, raw_fdr_bh := p.adjust(raw_permutation_p, method = "BH")]
metrics[, adjusted_fdr_bh := p.adjust(adjusted_permutation_p, method = "BH")]
metrics[, high_high_fdr_bh := p.adjust(high_high_permutation_p, method = "BH")]
fwrite(metrics, file.path(output_dir, "spatial_sample_metrics.csv"))
fwrite(
  rbindlist(spot_tables),
  file.path(output_dir, "spatial_spot_scores.csv.gz"),
  compress = "gzip"
)

metric_columns <- c(
  "raw_bivariate_moran",
  "adjusted_bivariate_moran",
  "high_high_enrichment"
)
patient_metrics <- metrics[, lapply(.SD, function(values) {
  finite <- values[is.finite(values)]
  if (length(finite) == 0L) NA_real_ else mean(finite)
}), by = patient_id, .SDcols = metric_columns]
fwrite(
  patient_metrics,
  file.path(output_dir, "spatial_patient_averaged_metrics.csv")
)

replicate_concordance <- rbindlist(lapply(metric_columns, function(metric) {
  wide <- dcast(
    metrics[, .(patient_id, replicate, value = get(metric))],
    patient_id ~ replicate,
    value.var = "value"
  )
  complete <- wide[complete.cases(Rep1, Rep2)]
  test <- if (nrow(complete) >= 3L) {
    suppressWarnings(cor.test(complete$Rep1, complete$Rep2, method = "spearman", exact = FALSE))
  } else {
    NULL
  }
  data.table(
    metric = metric,
    patients = nrow(complete),
    spearman_rho = if (is.null(test)) NA_real_ else unname(test$estimate),
    p_value = if (is.null(test)) NA_real_ else test$p.value
  )
}))
replicate_concordance[, fdr_bh := p.adjust(p_value, method = "BH")]
fwrite(
  replicate_concordance,
  file.path(output_dir, "spatial_replicate_concordance.csv")
)

summarize_patient_metric <- function(values, metric) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(data.table(
      metric = metric,
      patients = 0L,
      median = NA_real_,
      minimum = NA_real_,
      maximum = NA_real_,
      positive_patients = 0L,
      negative_patients = 0L,
      wilcoxon_p = NA_real_
    ))
  }
  test <- suppressWarnings(wilcox.test(values, mu = 0, exact = TRUE))
  data.table(
    metric = metric,
    patients = length(values),
    median = median(values),
    minimum = min(values),
    maximum = max(values),
    positive_patients = sum(values > 0),
    negative_patients = sum(values < 0),
    wilcoxon_p = test$p.value
  )
}

patient_summary <- rbind(
  summarize_patient_metric(patient_metrics$raw_bivariate_moran, "raw_bivariate_moran"),
  summarize_patient_metric(patient_metrics$adjusted_bivariate_moran, "adjusted_bivariate_moran"),
  summarize_patient_metric(
    patient_metrics$high_high_enrichment - 1,
    "high_high_enrichment_minus_one"
  )
)
patient_summary[, fdr_bh := p.adjust(wilcoxon_p, method = "BH")]
fwrite(patient_summary, file.path(output_dir, "spatial_patient_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
