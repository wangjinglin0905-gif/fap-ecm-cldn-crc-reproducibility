options(stringsAsFactors = FALSE, timeout = 180)
user_library <- Sys.getenv("FAP_R_LIBRARY", unset = "")
if (nzchar(user_library)) .libPaths(c(user_library, .libPaths()))

if (!requireNamespace("curl", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("The curl and jsonlite R packages are required.")
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_root, "work", "citation_audit", "new_references")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fetch_text <- function(url) {
  handle <- curl::new_handle(
    connecttimeout = 30,
    timeout = 180,
    useragent = "FAP-CLDN-manuscript-reference-audit/1.0"
  )
  response <- curl::curl_fetch_memory(url, handle = handle)
  if (response$status_code != 200L) stop("HTTP ", response$status_code, " for ", url)
  rawToChar(response$content)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

geo_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166555&targ=self&form=text&view=full"
geo_text <- fetch_text(geo_url)
writeLines(geo_text, file.path(output_dir, "GSE166555_SOFT.txt"), useBytes = TRUE)

extract_soft <- function(field) {
  pattern <- paste0("(?m)^!Series_", field, " = (.+)$")
  hit <- regexec(pattern, geo_text, perl = TRUE)
  match <- regmatches(geo_text, hit)[[1]]
  if (length(match) >= 2L) trimws(match[2]) else NA_character_
}

geo_summary <- data.frame(
  accession = "GSE166555",
  title = extract_soft("title"),
  pubmed_id = extract_soft("pubmed_id"),
  submission_date = extract_soft("submission_date"),
  last_update_date = extract_soft("last_update_date"),
  sample_count = length(gregexpr("(?m)^\\^SAMPLE = ", geo_text, perl = TRUE)[[1]])
)
write.csv(geo_summary, file.path(output_dir, "GSE166555_metadata.csv"), row.names = FALSE)

pubmed_url <- paste0(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=",
  geo_summary$pubmed_id,
  "&retmode=json"
)
pubmed_json <- jsonlite::fromJSON(fetch_text(pubmed_url), simplifyVector = FALSE)
pubmed_record <- pubmed_json$result[[geo_summary$pubmed_id]]
pubmed_summary <- data.frame(
  pubmed_id = geo_summary$pubmed_id,
  title = pubmed_record$title %||% NA_character_,
  journal = pubmed_record$fulljournalname %||% pubmed_record$source %||% NA_character_,
  publication_date = pubmed_record$pubdate %||% NA_character_,
  volume = pubmed_record$volume %||% NA_character_,
  issue = pubmed_record$issue %||% NA_character_,
  pages = pubmed_record$pages %||% NA_character_,
  authors = paste(vapply(pubmed_record$authors, function(author) author$name %||% "", character(1)), collapse = "; "),
  doi = {
    ids <- pubmed_record$articleids
    doi_hits <- vapply(ids, function(identifier) identical(identifier$idtype, "doi"), logical(1))
    if (any(doi_hits)) ids[[which(doi_hits)[1]]]$value else NA_character_
  }
)
write.csv(pubmed_summary, file.path(output_dir, "GSE166555_pubmed_reference.csv"), row.names = FALSE)

dois <- c(
  Valdeolivas_spatial = "10.1038/s41698-023-00488-4",
  CPTAC_colon = "10.1016/j.cell.2019.03.030",
  UALCAN = "10.1016/j.neo.2017.05.002"
)

crossref_rows <- lapply(names(dois), function(label) {
  doi <- dois[[label]]
  url <- paste0("https://api.crossref.org/works/", utils::URLencode(doi, reserved = TRUE))
  parsed <- jsonlite::fromJSON(fetch_text(url), simplifyVector = FALSE)$message
  authors <- vapply(parsed$author, function(author) {
    paste(c(author$given %||% "", author$family %||% ""), collapse = " ")
  }, character(1))
  year <- parsed$published$`date-parts`[[1]][[1]]
  data.frame(
    label = label,
    doi = parsed$DOI,
    title = paste(unlist(parsed$title), collapse = " "),
    journal = paste(unlist(parsed$`container-title`), collapse = " "),
    year = year,
    volume = parsed$volume %||% NA_character_,
    issue = parsed$issue %||% NA_character_,
    pages = parsed$page %||% parsed$`article-number` %||% NA_character_,
    authors = paste(authors, collapse = "; "),
    url = parsed$URL
  )
})

crossref <- do.call(rbind, crossref_rows)
write.csv(crossref, file.path(output_dir, "crossref_verified_new_references.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

print(geo_summary)
print(pubmed_summary)
print(crossref[, c("label", "doi", "title", "journal", "year")])
