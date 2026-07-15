options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(data.table))

output_dir <- file.path(
  "work", "reproducibility", "results", "L3_spatial_CellChat"
)
communications <- fread(
  file.path(output_dir, "spatial_cellchat_focused_pairs_all.csv")
)
prespecified <- communications[
  source == "Fibroblastic_stroma" & target == "Tumor"
]
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
