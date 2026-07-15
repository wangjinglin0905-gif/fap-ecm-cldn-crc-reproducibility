#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, timeout = 180)

stopifnot(requireNamespace("curl", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))
stopifnot(requireNamespace("data.table", quietly = TRUE))

input_file <- "work/revision_round3/extracted/prior_manuscript.txt"
out_dir <- "work/revision_round3/reports/reference_audit"
raw_dir <- file.path(out_dir, "crossref_raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

fetch_text <- function(url) {
  handle <- curl::new_handle(
    connecttimeout = 30,
    timeout = 180,
    useragent = paste0(
      "FAP-CLDN-JCMM-reference-audit/2.0 ",
      "(mailto:wangjinglin@gz5055.com)"
    )
  )
  response <- curl::curl_fetch_memory(url, handle = handle)
  if (response$status_code != 200L) {
    stop("HTTP ", response$status_code, " for ", url)
  }
  rawToChar(response$content)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

normalize_text <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^a-z0-9]+", " ", x))
  trimws(gsub("[[:space:]]+", " ", x))
}

lines <- readLines(input_file, encoding = "UTF-8", warn = FALSE)
reference_lines <- lines[grepl("^P01(2[0-9]|3[0-9]|4[0-9]|5[0-4])\\t", lines)]
reference_lines <- reference_lines[grepl("doi:", reference_lines, ignore.case = TRUE)]

extract_one <- function(line) {
  text <- sub("^[^\\t]+\\t\\[[^]]+\\]\\t", "", line)
  number <- as.integer(sub("^([0-9]+)\\..*$", "\\1", text))
  doi <- sub("^.*doi:([^[:space:]]+).*$", "\\1", text, ignore.case = TRUE)
  doi <- sub("[.;,]+$", "", doi)
  cited_first_author <- sub("^[0-9]+\\.\\s+([^[:space:]]+).*$", "\\1", text)
  data.frame(
    reference_number = number,
    cited_reference = text,
    cited_doi = tolower(doi),
    cited_first_author = cited_first_author
  )
}

cited <- data.table::rbindlist(lapply(reference_lines, extract_one))
cited <- rbind(
  cited,
  data.frame(
    reference_number = 36L,
    cited_reference = paste(
      "36. Boeckx B, De Roock W, De Craene B, De Hertogh G, Van den Eynde M,",
      "Demetter P, et al. Fibroblast activation protein identifies consensus molecular",
      "subtype 4 in colorectal cancer and allows its detection by 68Ga-FAPI-PET imaging.",
      "Br J Cancer. 2022;127:145-155. doi:10.1038/s41416-022-01748-z"
    ),
    cited_doi = "10.1038/s41416-022-01748-z",
    cited_first_author = "Boeckx"
  )
)
cited <- cited[order(cited$reference_number), ]

crossref_rows <- vector("list", nrow(cited))
for (i in seq_len(nrow(cited))) {
  doi <- cited$cited_doi[i]
  url <- paste0(
    "https://api.crossref.org/works/",
    utils::URLencode(doi, reserved = TRUE),
    "?mailto=wangjinglin%40gz5055.com"
  )
  cache_path <- file.path(raw_dir, sprintf("ref_%02d.json", cited$reference_number[i]))
  if (file.exists(cache_path)) {
    raw <- paste(readLines(cache_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  } else {
    raw <- fetch_text(url)
    writeLines(raw, cache_path, useBytes = TRUE)
  }
  parsed <- jsonlite::fromJSON(raw, simplifyVector = FALSE)$message
  authors <- vapply(parsed$author %||% list(), function(author) {
    trimws(paste(author$given %||% "", author$family %||% ""))
  }, character(1))
  date_parts <- parsed$published$`date-parts` %||% parsed$issued$`date-parts`
  year <- if (!is.null(date_parts) && length(date_parts)) date_parts[[1]][[1]] else NA_integer_
  authoritative_title <- paste(unlist(parsed$title %||% ""), collapse = " ")
  crossref_rows[[i]] <- data.frame(
    reference_number = cited$reference_number[i],
    cited_doi = cited$cited_doi[i],
    resolved_doi = tolower(parsed$DOI %||% ""),
    doi_exact_match = identical(tolower(parsed$DOI %||% ""), cited$cited_doi[i]),
    authoritative_title = authoritative_title,
    title_present_in_citation = grepl(
      normalize_text(authoritative_title),
      normalize_text(cited$cited_reference[i]),
      fixed = TRUE
    ),
    authoritative_first_author = if (length(parsed$author %||% list())) {
      parsed$author[[1]]$family %||% NA_character_
    } else {
      NA_character_
    },
    journal = paste(unlist(parsed$`container-title` %||% ""), collapse = " "),
    year = year,
    volume = parsed$volume %||% NA_character_,
    issue = parsed$issue %||% NA_character_,
    pages_or_article = parsed$page %||% parsed$`article-number` %||% NA_character_,
    authors = paste(authors, collapse = "; "),
    crossref_url = parsed$URL %||% NA_character_
  )
  if (!file.exists(cache_path)) Sys.sleep(0.12)
}

crossref <- data.table::rbindlist(crossref_rows, fill = TRUE)
audit <- merge(cited, crossref, by = c("reference_number", "cited_doi"), all.x = TRUE)
audit$first_author_match <- normalize_text(audit$cited_first_author) ==
  normalize_text(audit$authoritative_first_author)

# Resolve DOI-to-PMID mappings in a single NCBI request. Records not indexed in
# PubMed are retained as NA and remain verifiable through Crossref.
idconv_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=",
  utils::URLencode(paste(audit$cited_doi, collapse = ","), reserved = TRUE),
  "&format=json&tool=FAP_CLDN_JCMM_audit&email=wangjinglin%40gz5055.com"
)
idconv_cache <- file.path(out_dir, "ncbi_idconv_raw.json")
if (file.exists(idconv_cache)) {
  idconv_raw <- paste(readLines(idconv_cache, encoding = "UTF-8", warn = FALSE), collapse = "\n")
} else {
  idconv_raw <- fetch_text(idconv_url)
  writeLines(idconv_raw, idconv_cache, useBytes = TRUE)
}
idconv <- jsonlite::fromJSON(idconv_raw, simplifyVector = TRUE)$records
if (is.data.frame(idconv) && nrow(idconv)) {
  idconv$doi <- tolower(idconv$doi)
  audit <- merge(
    audit,
    idconv[, intersect(c("doi", "pmid", "pmcid", "status"), names(idconv)), drop = FALSE],
    by.x = "cited_doi",
    by.y = "doi",
    all.x = TRUE
  )
}

audit <- audit[order(audit$reference_number), ]
audit$verification_status <- ifelse(
  audit$doi_exact_match %in% TRUE & audit$title_present_in_citation %in% TRUE &
    audit$first_author_match %in% TRUE,
  "verified_doi_title_first_author",
  ifelse(audit$doi_exact_match %in% TRUE, "verified_doi_metadata_review_needed", "unresolved")
)

data.table::fwrite(audit, file.path(out_dir, "all_reference_verification.csv"), na = "")
summary <- data.frame(
  metric = c(
    "references_checked", "doi_exact_matches", "titles_present_in_citation",
    "verified_doi_title_first_author", "first_author_mismatches", "pubmed_indexed", "unresolved"
  ),
  value = c(
    nrow(audit),
    sum(audit$doi_exact_match %in% TRUE),
    sum(audit$title_present_in_citation %in% TRUE),
    sum(audit$verification_status == "verified_doi_title_first_author"),
    sum(!(audit$first_author_match %in% TRUE)),
    if ("pmid" %in% names(audit)) sum(!is.na(audit$pmid) & nzchar(audit$pmid)) else 0L,
    sum(audit$verification_status == "unresolved")
  )
)
data.table::fwrite(summary, file.path(out_dir, "reference_verification_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

print(summary, row.names = FALSE)
if (any(audit$verification_status != "verified_doi_title_first_author")) {
  print(audit[audit$verification_status != "verified_doi_title_first_author",
              c("reference_number", "cited_doi", "cited_first_author",
                "authoritative_first_author", "authoritative_title", "verification_status")])
}
