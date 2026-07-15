options(stringsAsFactors = FALSE, timeout = 180)

project_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
input_dir <- file.path(project_root, "work", "reproducibility", "inputs", "ualcan")
result_dir <- file.path(project_root, "work", "reproducibility", "results", "UALCAN_CPTAC")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(project_root, "work", "Rlib"), .libPaths()))
if (!requireNamespace("curl", quietly = TRUE)) {
  stop("The curl R package is required.")
}

genes <- c("FAP", "COL1A1", "COL1A2", "FN1", "CLDN1", "CLDN4")
base_url <- "https://ualcan.path.uab.edu/cgi-bin/CPTAC-Result.pl?genenam=%s&ctype=Colon"

fetch_page <- function(gene) {
  url <- sprintf(base_url, gene)
  raw_path <- file.path(input_dir, sprintf("%s.html", gene))
  if (file.exists(raw_path) && file.info(raw_path)$size > 1000L) {
    status_code <- 200L
    bytes <- file.info(raw_path)$size
  } else {
    handle <- curl::new_handle(
      connecttimeout = 30,
      timeout = 180,
      useragent = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) Codex-reproducibility-audit"
    )
    response <- curl::curl_fetch_memory(url, handle = handle)
    if (response$status_code != 200L) {
      stop(sprintf("UALCAN returned HTTP %s for %s", response$status_code, gene))
    }
    writeBin(response$content, raw_path)
    status_code <- response$status_code
    bytes <- length(response$content)
  }
  checksum <- if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(file = raw_path, algo = "sha256")
  } else {
    unname(tools::md5sum(raw_path))
  }
  list(
    gene = gene,
    url = url,
    status_code = status_code,
    bytes = bytes,
    raw_path = file.path("work", "reproducibility", "inputs", "ualcan", basename(raw_path)),
    checksum_algorithm = if (requireNamespace("digest", quietly = TRUE)) "SHA-256" else "MD5",
    checksum = checksum
  )
}

extract_number <- function(text, field) {
  pattern <- sprintf("%s\\s*:\\s*([-+]?[0-9]*\\.?[0-9]+(?:[Ee][-+]?[0-9]+)?)", field)
  hit <- regexec(pattern, text, perl = TRUE)
  match <- regmatches(text, hit)[[1]]
  if (length(match) < 2L) return(NA_real_)
  as.numeric(match[2])
}

parse_page <- function(gene, raw_path) {
  html <- paste(readLines(raw_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  category_hit <- regexpr("categories\\s*:\\s*\\[[^]]+\\]", html, perl = TRUE)
  if (category_hit[1] < 0L) {
    return(data.frame(gene = gene, available = FALSE, reason = "Primary tumour boxplot not found"))
  }
  category_text <- regmatches(html, category_hit)
  category_matches <- gregexpr("'([^']+)'", category_text, perl = TRUE)
  categories <- regmatches(category_text, category_matches)[[1]]
  categories <- gsub("^'|'$", "", categories)
  categories <- gsub("<br>", " ", categories, fixed = TRUE)
  category_n <- suppressWarnings(as.integer(sub(".*\\(n=([0-9]+)\\).*", "\\1", categories)))

  remainder <- substring(html, category_hit[1] + attr(category_hit, "match.length"))
  series_hit <- regexpr("series\\s*:\\s*\\[\\s*\\{data\\s*:\\s*\\[", remainder, perl = TRUE)
  if (series_hit[1] < 0L) {
    return(data.frame(gene = gene, available = FALSE, reason = "Primary tumour series not found"))
  }
  series_text <- substring(remainder, series_hit[1])
  series_end <- regexpr("\\]\\s*,\\s*tooltip", series_text, perl = TRUE)
  if (series_end[1] < 0L) {
    return(data.frame(gene = gene, available = FALSE, reason = "Primary tumour series end not found"))
  }
  series_text <- substring(series_text, 1L, series_end[1])
  object_hits <- gregexpr("\\{[^{}]+\\}", series_text, perl = TRUE)
  objects <- regmatches(series_text, object_hits)[[1]]
  if (length(objects) < 2L || length(categories) < 2L) {
    return(data.frame(gene = gene, available = FALSE, reason = "Fewer than two boxplot groups"))
  }

  p_pattern <- "Normal-vs-Primary</td><td[^>]*>\\s*([0-9.Ee+-]+)"
  p_hit <- regexec(p_pattern, html, perl = TRUE)
  p_match <- regmatches(html, p_hit)[[1]]
  p_value <- if (length(p_match) >= 2L) as.numeric(p_match[2]) else NA_real_

  normal <- objects[[1]]
  tumour <- objects[[2]]
  normal_median <- extract_number(normal, "median")
  tumour_median <- extract_number(tumour, "median")
  data.frame(
    gene = gene,
    available = TRUE,
    normal_n = category_n[1],
    tumour_n = category_n[2],
    normal_low = extract_number(normal, "low"),
    normal_q1 = extract_number(normal, "q1"),
    normal_median = normal_median,
    normal_q3 = extract_number(normal, "q3"),
    normal_high = extract_number(normal, "high"),
    tumour_low = extract_number(tumour, "low"),
    tumour_q1 = extract_number(tumour, "q1"),
    tumour_median = tumour_median,
    tumour_q3 = extract_number(tumour, "q3"),
    tumour_high = extract_number(tumour, "high"),
    median_difference = tumour_median - normal_median,
    direction = ifelse(tumour_median > normal_median, "Higher in primary tumour",
      ifelse(tumour_median < normal_median, "Lower in primary tumour", "No median difference")),
    p_value = p_value,
    test_reported_by_ualcan = "Unpaired two-sample t-test",
    source_note = "UALCAN display of CPTAC colon proteomics; not an independent cohort from CPTAC",
    reason = NA_character_
  )
}

manifest_rows <- lapply(genes, fetch_page)
manifest <- do.call(rbind, lapply(manifest_rows, as.data.frame))
write.csv(manifest, file.path(result_dir, "ualcan_download_manifest.csv"), row.names = FALSE)

summary_rows <- lapply(seq_len(nrow(manifest)), function(index) {
  parse_page(manifest$gene[index], manifest$raw_path[index])
})
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(result_dir, "ualcan_cptac_primary_vs_normal.csv"), row.names = FALSE)

cat("Downloaded", nrow(manifest), "UALCAN CPTAC pages to", input_dir, "\n")
cat("Manifest:", file.path(result_dir, "ualcan_download_manifest.csv"), "\n")
cat("Parsed results:", file.path(result_dir, "ualcan_cptac_primary_vs_normal.csv"), "\n")
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
