options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260805)

suppressPackageStartupMessages(library(jsonlite))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

root <- normalizePath(getwd(), winslash = "/")
input_file <- file.path(root, "data", "public",
                        "TCGA_COADREAD_expression.txt.gz")
out_dir <- file.path(root, "results", "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(input_file))

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
receptor2 <- c("SDC4", "CD44")
target_genes <- unique(c(FAP13, matrix4, receptor2))

score_zmean <- function(x, genes) {
  available <- intersect(genes, colnames(x))
  stopifnot(length(available) == length(genes))
  rowMeans(scale(x[, available, drop = FALSE]), na.rm = TRUE)
}

bootstrap_cor <- function(x, y, B = 5000L) {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  observed <- suppressWarnings(cor(x, y, method = "spearman"))
  boot <- replicate(B, {
    idx <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(n = length(x), rho = observed,
    ci_low = unname(quantile(boot, 0.025, na.rm = TRUE)),
    ci_high = unname(quantile(boot, 0.975, na.rm = TRUE)),
    p_value = test$p.value)
}

extract_diagnosis <- function(hit) {
  diagnoses <- hit$diagnoses
  if (is.null(diagnoses) || length(diagnoses) == 0) return(NULL)
  primary <- which(vapply(diagnoses, function(x) {
    isTRUE(x$diagnosis_is_primary_disease)
  }, logical(1)))
  diagnoses[[if (length(primary)) primary[1] else 1]]
}

filters <- toJSON(list(op = "in",
                       content = list(field = "project.project_id",
                                      value = c("TCGA-COAD", "TCGA-READ"))),
                  auto_unbox = TRUE)
fields <- paste(c("submitter_id", "project.project_id",
                  "diagnoses.ajcc_pathologic_stage",
                  "diagnoses.ajcc_pathologic_t",
                  "diagnoses.ajcc_pathologic_n",
                  "diagnoses.diagnosis_is_primary_disease"), collapse = ",")
api_url <- paste0(
  "https://api.gdc.cancer.gov/cases?filters=",
  URLencode(filters, reserved = TRUE),
  "&fields=", URLencode(fields, reserved = TRUE),
  "&expand=diagnoses&format=JSON&size=2000"
)

raw_json <- file.path(root, "data", "public",
                       "TCGA_GDC_clinical_query_raw.json")
if (!file.exists(raw_json) || file.info(raw_json)$size == 0) {
  download.file(api_url, raw_json, mode = "wb", method = "libcurl",
                quiet = TRUE)
}
gdc <- fromJSON(raw_json, simplifyVector = FALSE)
hits <- gdc$data$hits
stopifnot(length(hits) > 0)

clinical_rows <- lapply(hits, function(hit) {
  diagnosis <- extract_diagnosis(hit)
  data.frame(
    patient_id = hit$submitter_id %||% NA_character_,
    project_id = hit$project$project_id %||% NA_character_,
    pathologic_stage = if (is.null(diagnosis)) NA_character_ else
      diagnosis$ajcc_pathologic_stage %||% NA_character_,
    pathologic_t = if (is.null(diagnosis)) NA_character_ else
      diagnosis$ajcc_pathologic_t %||% NA_character_,
    pathologic_n = if (is.null(diagnosis)) NA_character_ else
      diagnosis$ajcc_pathologic_n %||% NA_character_
  )
})
clinical <- do.call(rbind, clinical_rows)

expr_all <- read.delim(gzfile(input_file), row.names = 1, check.names = FALSE)
missing_genes <- setdiff(target_genes, rownames(expr_all))
stopifnot(length(missing_genes) == 0)
sample_ids <- colnames(expr_all)
sample_type <- substr(sample_ids, 14, 15)
primary_ids <- sample_ids[sample_type == "01"]
stopifnot(length(primary_ids) > 0)

patient_ids <- substr(primary_ids, 1, 12)
keep_unique <- !duplicated(patient_ids)
primary_ids <- primary_ids[keep_unique]
patient_ids <- patient_ids[keep_unique]
matched <- match(patient_ids, clinical$patient_id)

expr <- as.data.frame(t(as.matrix(expr_all[target_genes, primary_ids,
                                             drop = FALSE])))
rm(expr_all)
manifest <- data.frame(sample_id = primary_ids, patient_id = patient_ids,
                       stringsAsFactors = FALSE)
manifest <- cbind(manifest, clinical[matched,
                                    c("project_id", "pathologic_stage",
                                      "pathologic_t", "pathologic_n")])
manifest$FAP <- expr$FAP
manifest$FAP13 <- score_zmean(expr, FAP13)
manifest$matrix4 <- score_zmean(expr, matrix4)
manifest$receptor2 <- score_zmean(expr, receptor2)
manifest$SDC4 <- expr$SDC4
manifest$CD44 <- expr$CD44
manifest$proxy6 <- score_zmean(expr, c(matrix4, receptor2))
write.csv(manifest, file.path(out_dir, "TCGA_full_primary_manifest.csv"),
          row.names = FALSE)

accounting <- data.frame(
  item = c("expression_profiles", "primary_tumor_code_01",
           "unique_primary_tumor_patients", "matched_GDC_cases",
           "TCGA_COAD", "TCGA_READ", "known_pathologic_stage",
           "known_pathologic_T", "known_pathologic_N"),
  n = c(length(sample_ids), sum(sample_type == "01"), nrow(manifest),
        sum(!is.na(manifest$project_id)), sum(manifest$project_id == "TCGA-COAD",
                                             na.rm = TRUE),
        sum(manifest$project_id == "TCGA-READ", na.rm = TRUE),
        sum(!is.na(manifest$pathologic_stage)),
        sum(!is.na(manifest$pathologic_t)),
        sum(!is.na(manifest$pathologic_n)))
)
write.csv(accounting,
          file.path(out_dir, "TCGA_full_primary_accounting.csv"),
          row.names = FALSE)

comparisons <- list(
  "FAP13 vs matrix4" = c("FAP13", "matrix4"),
  "FAP13 vs receptor2" = c("FAP13", "receptor2"),
  "FAP13 vs SDC4" = c("FAP13", "SDC4"),
  "FAP13 vs CD44" = c("FAP13", "CD44"),
  "FAP vs matrix4" = c("FAP", "matrix4")
)
cor_rows <- list()
counter <- 1L
for (cohort in c("TCGA-COADREAD", "TCGA-COAD", "TCGA-READ")) {
  subset <- if (cohort == "TCGA-COADREAD") manifest else
    manifest[manifest$project_id == cohort, ]
  for (comparison in names(comparisons)) {
    vars <- comparisons[[comparison]]
    values <- bootstrap_cor(subset[[vars[1]]], subset[[vars[2]]])
    cor_rows[[counter]] <- data.frame(cohort = cohort,
                                      comparison = comparison, t(values))
    counter <- counter + 1L
  }
}
cor_results <- do.call(rbind, cor_rows)
cor_results$fdr_bh_five_within_cohort <- ave(
  cor_results$p_value, cor_results$cohort,
  FUN = function(x) p.adjust(x, method = "BH"))
write.csv(cor_results,
          file.path(out_dir, "TCGA_full_primary_correlations.csv"),
          row.names = FALSE)

parse_t <- function(x) {
  out <- suppressWarnings(as.numeric(sub("^T([1-4]).*", "\\1", x)))
  out[!grepl("^T[1-4]", x)] <- NA_real_
  out
}
manifest$T_ordinal <- parse_t(manifest$pathologic_t)
trend_vars <- c("FAP", "FAP13", "matrix4", "receptor2")
trend <- do.call(rbind, lapply(trend_vars, function(variable) {
  values <- bootstrap_cor(manifest$T_ordinal, manifest[[variable]])
  data.frame(variable = variable, t(values))
}))
trend$fdr_bh_four <- p.adjust(trend$p_value, method = "BH")
write.csv(trend, file.path(out_dir, "TCGA_full_primary_T_trends.csv"),
          row.names = FALSE)

known_n <- grepl("^N[0-2]", manifest$pathologic_n)
nodal <- manifest[known_n, ]
nodal$node_positive <- grepl("^N[12]", nodal$pathologic_n)
threshold <- median(nodal$proxy6, na.rm = TRUE)
nodal$proxy_group <- ifelse(nodal$proxy6 >= threshold, "high", "low")
high_pos <- sum(nodal$proxy_group == "high" & nodal$node_positive)
high_neg <- sum(nodal$proxy_group == "high" & !nodal$node_positive)
low_pos <- sum(nodal$proxy_group == "low" & nodal$node_positive)
low_neg <- sum(nodal$proxy_group == "low" & !nodal$node_positive)
tab <- matrix(c(high_pos, high_neg, low_pos, low_neg), nrow = 2, byrow = TRUE,
              dimnames = list(group = c("high", "low"),
                              node = c("positive", "negative")))
fisher <- fisher.test(tab)
nodal$T_ordinal <- parse_t(nodal$pathologic_t)
model_data <- nodal[complete.cases(nodal[, c("node_positive", "proxy_group",
                                             "proxy6", "T_ordinal")]), ]
model_data$proxy_high <- as.integer(model_data$proxy_group == "high")
binary_model <- glm(node_positive ~ proxy_high + T_ordinal,
                    family = binomial(), data = model_data)
continuous_model <- glm(node_positive ~ proxy6 + T_ordinal,
                        family = binomial(), data = model_data)

model_row <- function(model, term, analysis) {
  beta <- coef(model)[term]
  se <- sqrt(diag(vcov(model)))[term]
  data.frame(
    analysis = analysis, n = nobs(model), threshold = threshold,
    high_positive = high_pos, high_negative = high_neg,
    low_positive = low_pos, low_negative = low_neg,
    odds_ratio = exp(beta), ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    p_value = summary(model)$coefficients[term, "Pr(>|z|)"]
  )
}

nodal_result <- rbind(
  data.frame(
    analysis = "unadjusted median split, Fisher exact", n = nrow(nodal),
    threshold = threshold, high_positive = high_pos, high_negative = high_neg,
    low_positive = low_pos, low_negative = low_neg,
    odds_ratio = unname(fisher$estimate), ci_low = fisher$conf.int[1],
    ci_high = fisher$conf.int[2], p_value = fisher$p.value
  ),
  model_row(binary_model, "proxy_high",
            "median split adjusted for ordinal T category"),
  model_row(continuous_model, "proxy6",
            "continuous proxy6 adjusted for ordinal T category")
)
write.csv(nodal_result,
          file.path(out_dir, "TCGA_full_primary_nodal_proxy.csv"),
          row.names = FALSE)

log_file <- file.path(out_dir, "TCGA_full_primary_audit.log")
sink(log_file)
cat("TCGA FULL PRIMARY-TUMOR AUDIT\n")
cat("Query UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
cat("GDC API cases returned:", length(hits), "\n")
cat("GDC data release endpoint checked separately from the API response.\n\n")
print(accounting, row.names = FALSE)
cat("\nCorrelations:\n")
print(cor_results, row.names = FALSE)
cat("\nT-category trends:\n")
print(trend, row.names = FALSE)
cat("\nNodal proxy:\n")
print(nodal_result, row.names = FALSE)
cat("\nSession information:\n")
print(sessionInfo())
sink()

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "TCGA_full_primary_sessionInfo.txt"))
