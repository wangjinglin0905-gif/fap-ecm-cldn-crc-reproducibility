#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260714)

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(spdep)
})

task_root <- normalizePath(file.path(getwd(), "work", "reproducibility"), mustWork = FALSE)
input_file <- file.path(task_root, "inputs", "41467_2022_29366_MOESM9_ESM.xlsx")
output_dir <- file.path(task_root, "results", "L2_spatial")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))

sheet_map <- list(
  P6 = c(scores = "Supplementary Figure 9k", immune = "Supplementary Figure 9s", types = "Figure 6b"),
  P8 = c(scores = "Supplementary Figure 9l", immune = "Supplementary Figure 9t", types = "Supplementary Figure 9b"),
  P9 = c(scores = "Supplementary Figure 9m", immune = "Supplementary Figure 9u", types = "Supplementary Figure 9e")
)

make_neighbors <- function(coords, k) {
  # Resolve equal-distance neighbours deterministically by original row index.
  # knearneigh() can choose different members of a distance tie across package
  # builds, which previously altered a small number of boundary-spot values.
  n <- nrow(coords)
  if (n <= k) stop("The number of spots must exceed k")
  index <- seq_len(n)
  neighbours <- lapply(index, function(i) {
    squared_distance <-
      (coords[, 1L] - coords[i, 1L])^2 +
      (coords[, 2L] - coords[i, 2L])^2
    ordered <- order(squared_distance, index, method = "radix")
    head(ordered[ordered != i], k)
  })
  attr(neighbours, "region.id") <- as.character(index)
  attr(neighbours, "call") <- match.call()
  attr(neighbours, "type") <- "deterministic k-nearest neighbours"
  attr(neighbours, "sym") <- FALSE
  class(neighbours) <- "nb"
  neighbours
}

read_patient <- function(patient, sheets) {
  scores <- as.data.table(read_excel(input_file, sheet = sheets[["scores"]]))
  immune <- as.data.table(read_excel(input_file, sheet = sheets[["immune"]]))
  types <- as.data.table(read_excel(input_file, sheet = sheets[["types"]]))
  stopifnot(nrow(scores) == nrow(immune), nrow(scores) == nrow(types))
  stopifnot(isTRUE(all.equal(scores$imagerow, immune$imagerow)), isTRUE(all.equal(scores$imagecol, immune$imagecol)))
  stopifnot(isTRUE(all.equal(scores$imagerow, types$imagerow)), isTRUE(all.equal(scores$imagecol, types$imagecol)))

  immune_genes <- c("CD3D", "CD8A", "CD4", "CD19", "MS4A1", "GZMA")
  data <- data.table(
    patient = patient,
    row = scores$imagerow,
    col = scores$imagecol,
    DefineTypes = types$DefineTypes,
    FAP_fibroblasts = scores$`FAP+.fibroblasts`,
    SPP1_macrophages = scores$`SPP1+.macrophages`,
    immune_score = rowMeans(as.matrix(immune[, ..immune_genes]), na.rm = TRUE),
    epithelial_flag = as.numeric(grepl("epithelial", types$DefineTypes, ignore.case = TRUE))
  )

  coords <- as.matrix(data[, .(row, col)])
  neighbors6 <- make_neighbors(coords, k = 6)
  weights6 <- nb2listw(neighbors6, style = "W", zero.policy = TRUE)
  for (score_name in c("FAP_fibroblasts", "SPP1_macrophages", "immune_score", "epithelial_flag")) {
    neighbor_mean <- lag.listw(weights6, data[[score_name]], zero.policy = TRUE)
    data[[paste0("nhood_", score_name)]] <- (data[[score_name]] + 6 * neighbor_mean) / 7
  }

  fap_threshold <- quantile(data$nhood_FAP_fibroblasts, 0.75, na.rm = TRUE)
  spp1_threshold <- quantile(data$nhood_SPP1_macrophages, 0.75, na.rm = TRUE)
  immune_threshold <- quantile(data$nhood_immune_score, 0.85, na.rm = TRUE)
  fap_high <- data$nhood_FAP_fibroblasts > fap_threshold
  spp1_high <- data$nhood_SPP1_macrophages > spp1_threshold
  immune_high <- data$nhood_immune_score > immune_threshold
  habitat <- fifelse(
    fap_high & spp1_high,
    "FAP+CAF/SPP1+TAM",
    fifelse(
      fap_high,
      "FAP+CAF-dominant",
      fifelse(spp1_high, "SPP1+TAM-dominant", fifelse(immune_high, "TLS-like immune", "Epithelial/other"))
    )
  )
  set(data, j = "habitat", value = habitat)

  thresholds <- data.table(
    patient = patient,
    fap_q75 = fap_threshold,
    spp1_q75 = spp1_threshold,
    immune_q85 = immune_threshold
  )

  neighbors8 <- make_neighbors(coords, k = 8)
  weights8 <- nb2listw(neighbors8, style = "W", zero.policy = TRUE)
  moran_analytic <- moran.test(data$FAP_fibroblasts, weights8, alternative = "greater", zero.policy = TRUE)
  moran_permutation <- moran.mc(data$FAP_fibroblasts, weights8, nsim = 999, alternative = "greater", zero.policy = TRUE)
  spp1_lag <- lag.listw(weights8, data$SPP1_macrophages, zero.policy = TRUE)
  cross_test <- suppressWarnings(cor.test(data$FAP_fibroblasts, spp1_lag, method = "spearman", exact = FALSE))

  moran <- data.table(
    patient = patient,
    n_spots = nrow(data),
    moran_i = unname(moran_analytic$estimate[[1]]),
    expected_i = unname(moran_analytic$estimate[[2]]),
    analytic_p = moran_analytic$p.value,
    permutation_p = moran_permutation$p.value,
    permutations = 999,
    fap_neighbor_spp1_rho = unname(cross_test$estimate),
    fap_neighbor_spp1_p = cross_test$p.value
  )

  fap_spatial_high <- data$nhood_FAP_fibroblasts > quantile(data$nhood_FAP_fibroblasts, 0.75, na.rm = TRUE)
  epithelial_spatial_high <- data$nhood_epithelial_flag > quantile(data$nhood_epithelial_flag, 0.75, na.rm = TRUE)
  cooccurrence <- data.table(
    patient = patient,
    n_spots = nrow(data),
    n_fap_high_epithelial_high = sum(fap_spatial_high & epithelial_spatial_high),
    percent_all_spots = 100 * mean(fap_spatial_high & epithelial_spatial_high),
    fap_q75 = quantile(data$nhood_FAP_fibroblasts, 0.75, na.rm = TRUE),
    epithelial_q75 = quantile(data$nhood_epithelial_flag, 0.75, na.rm = TRUE)
  )

  list(data = data, thresholds = thresholds, moran = moran, cooccurrence = cooccurrence)
}

patient_results <- lapply(names(sheet_map), function(patient) read_patient(patient, sheet_map[[patient]]))
names(patient_results) <- names(sheet_map)

spot_table <- rbindlist(lapply(patient_results, `[[`, "data"))
threshold_table <- rbindlist(lapply(patient_results, `[[`, "thresholds"))
moran_table <- rbindlist(lapply(patient_results, `[[`, "moran"))
cooccurrence_table <- rbindlist(lapply(patient_results, `[[`, "cooccurrence"))
habitat_summary <- spot_table[, .(n_spots = .N), by = .(patient, habitat)]
habitat_summary[, percent := 100 * n_spots / sum(n_spots), by = patient]

fwrite(spot_table, file.path(output_dir, "spatial_spot_habitats_recomputed.csv"))
fwrite(threshold_table, file.path(output_dir, "spatial_habitat_thresholds.csv"))
fwrite(habitat_summary, file.path(output_dir, "spatial_habitat_summary.csv"))
fwrite(moran_table, file.path(output_dir, "spatial_morans_i.csv"))
fwrite(cooccurrence_table, file.path(output_dir, "spatial_fap_epithelial_cooccurrence.csv"))

legacy_file <- Sys.getenv("FAP_LEGACY_HABITAT", file.path(task_root, "inputs", "CANVAS_like_CRC_FAP_SPP1_spot_habitats.csv"))
if (file.exists(legacy_file)) {
  legacy <- fread(legacy_file)
  comparison <- merge(
    spot_table[, .(patient, row, col, recomputed_habitat = habitat)],
    legacy[, .(patient, row, col, legacy_habitat = habitat)],
    by = c("patient", "row", "col"),
    all = TRUE
  )
  audit <- data.table(
    n_recomputed = nrow(spot_table),
    n_legacy = nrow(legacy),
    n_matched = sum(complete.cases(comparison$recomputed_habitat, comparison$legacy_habitat)),
    habitat_agreement = mean(comparison$recomputed_habitat == comparison$legacy_habitat, na.rm = TRUE)
  )
  fwrite(audit, file.path(output_dir, "legacy_habitat_agreement.csv"))
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("L2 spatial analysis complete: ", output_dir)
