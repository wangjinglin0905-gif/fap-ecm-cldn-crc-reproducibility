#!/usr/bin/env Rscript

# Verify manuscript references against Crossref and PubMed official metadata.
# Crossref is queried first by DOI; PubMed is then queried by DOI [AID] or,
# for DOI-less references, by exact title. No identifier is inferred locally.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript verify_references_crossref_pubmed.R EXTRACTED_MD OUTPUT_DIR")
}
input_md <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(jsonlite))

user_agent <- "FAP-manuscript-reference-audit/1.0 (mailto:manuscript.audit@example.org)"
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

get_json <- function(url, tries = 3L) {
  last_error <- NULL
  for (attempt in seq_len(tries)) {
    tmp <- tempfile(fileext = ".json")
    ok <- tryCatch({
      suppressWarnings(download.file(url, tmp, quiet = TRUE, mode = "wb",
                                     method = "libcurl",
                                     headers = c("User-Agent" = user_agent)))
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    })
    if (ok) {
      ans <- tryCatch(jsonlite::fromJSON(tmp, simplifyVector = TRUE),
                      error = function(e) { last_error <<- conditionMessage(e); NULL })
      unlink(tmp)
      if (!is.null(ans)) return(ans)
    } else {
      unlink(tmp)
    }
    Sys.sleep(attempt)
  }
  structure(NULL, error = last_error)
}

normalise_text <- function(x) {
  x <- gsub("<[^>]+>", " ", x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(x)
  gsub("[^a-z0-9]+", " ", x)
}

extract_pubmed_doi <- function(articleids) {
  if (is.null(articleids)) return(NA_character_)
  if (is.data.frame(articleids)) {
    id_col <- if ("idtype" %in% names(articleids)) "idtype" else if ("idtype" %in% tolower(names(articleids))) names(articleids)[match("idtype", tolower(names(articleids)))] else NA
    value_col <- if ("value" %in% names(articleids)) "value" else NA
    if (!is.na(id_col) && !is.na(value_col)) {
      hit <- which(tolower(articleids[[id_col]]) == "doi")
      if (length(hit)) return(as.character(articleids[[value_col]][hit[[1]]]))
    }
  }
  NA_character_
}

lines <- readLines(input_md, warn = FALSE, encoding = "UTF-8")
ref_start <- grep("(style=Heading1) References", lines, fixed = TRUE)
if (length(ref_start) != 1L) stop("Could not uniquely locate References heading")
ref_lines <- lines[seq.int(ref_start + 1L, length(lines))]
ref_lines <- ref_lines[grepl("^\\[P[0-9]+\\] [0-9]+\\. ", ref_lines)]

citations <- sub("^\\[P[0-9]+\\] ", "", ref_lines)
numbers <- as.integer(sub("^([0-9]+)\\..*$", "\\1", citations))
citations <- citations[order(numbers)]
numbers <- sort(numbers)
if (!identical(numbers, seq_len(length(numbers)))) stop("Reference numbering is not consecutive")

doi_match <- regexpr("https://doi\\.org/[^ ]+$", citations, perl = TRUE)
dois <- rep(NA_character_, length(citations))
has_doi <- doi_match > 0
dois[has_doi] <- sub("^https://doi\\.org/", "", regmatches(citations, doi_match), ignore.case = TRUE)

# Exact title is needed only for the DOI-less Anticancer Research reference.
title_overrides <- rep(NA_character_, length(citations))
title_overrides[numbers == 8L] <- "Selective up-regulation of claudin-1 and claudin-2 in colorectal cancer"

rows <- vector("list", length(citations))
for (i in seq_along(citations)) {
  number <- numbers[[i]]
  citation <- citations[[i]]
  doi <- dois[[i]]
  message("Verifying reference ", number, "/", length(citations), ifelse(is.na(doi), " (title query)", paste0(" DOI ", doi)))

  cr_status <- "not_queried_no_doi"
  cr_title <- cr_journal <- NA_character_
  cr_year <- NA_integer_
  title_match <- NA
  if (!is.na(doi)) {
    cr_url <- paste0("https://api.crossref.org/works/", URLencode(doi, reserved = TRUE))
    cr <- get_json(cr_url)
    if (!is.null(cr) && identical(as.character(cr$status), "ok")) {
      cr_status <- "verified"
      msg <- cr$message
      if (!is.null(msg$title) && length(msg$title)) cr_title <- as.character(msg$title[[1]])
      if (!is.null(msg$container.title) && length(msg$container.title)) cr_journal <- as.character(msg$container.title[[1]])
      date_parts <- NULL
      if (!is.null(msg$published$`date-parts`)) date_parts <- msg$published$`date-parts`
      if (is.null(date_parts) && !is.null(msg$issued$`date-parts`)) date_parts <- msg$issued$`date-parts`
      if (!is.null(date_parts)) cr_year <- suppressWarnings(as.integer(unlist(date_parts)[[1]]))
      if (!is.na(cr_title)) {
        title_match <- grepl(trimws(normalise_text(cr_title)), normalise_text(citation), fixed = TRUE)
      }
    } else {
      cr_status <- "not_found_or_request_failed"
    }
  }

  term <- if (!is.na(doi)) paste0(doi, "[AID]") else paste0('"', title_overrides[[i]], '"[Title]')
  esearch_url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=5&term=",
    URLencode(term, reserved = TRUE)
  )
  esearch <- get_json(esearch_url)
  pmids <- character()
  if (!is.null(esearch) && !is.null(esearch$esearchresult$idlist)) {
    pmids <- as.character(unlist(esearch$esearchresult$idlist))
  }

  pmid <- pm_title <- pm_journal <- pm_doi <- NA_character_
  pm_year <- NA_integer_
  retraction_flag <- FALSE
  pubmed_status <- if (length(pmids)) "verified" else "not_indexed_or_not_found"
  if (length(pmids)) {
    pmid <- pmids[[1]]
    esummary_url <- paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=",
      URLencode(pmid, reserved = TRUE)
    )
    summary_json <- get_json(esummary_url)
    rec <- if (!is.null(summary_json$result[[pmid]])) summary_json$result[[pmid]] else NULL
    if (!is.null(rec)) {
      pm_title <- as.character(rec$title %||% NA_character_)
      pm_journal <- as.character(rec$fulljournalname %||% rec$source %||% NA_character_)
      pm_year <- suppressWarnings(as.integer(sub("^([0-9]{4}).*$", "\\1", as.character(rec$pubdate %||% ""))))
      pm_doi <- extract_pubmed_doi(rec$articleids)
      pubtypes <- paste(as.character(unlist(rec$pubtype)), collapse = "; ")
      retraction_flag <- grepl("Retracted Publication|Retraction of Publication", pubtypes, ignore.case = TRUE) ||
        grepl("retracted", pm_title, ignore.case = TRUE)
    }
  }

  doi_consistent <- if (!is.na(doi) && !is.na(pm_doi)) tolower(doi) == tolower(pm_doi) else NA
  status <- if (!is.na(doi) && cr_status != "verified") {
    "FAIL_CROSSREF"
  } else if (!is.na(title_match) && !title_match) {
    "REVIEW_TITLE"
  } else if (isTRUE(retraction_flag)) {
    "FAIL_RETRACTION_FLAG"
  } else if (!is.na(doi_consistent) && !doi_consistent) {
    "REVIEW_DOI_MISMATCH"
  } else if (pubmed_status == "verified") {
    "VERIFIED_CROSSREF_PUBMED"
  } else if (cr_status == "verified") {
    "VERIFIED_CROSSREF_NOT_IN_PUBMED"
  } else {
    "REVIEW_NO_IDENTIFIER_MATCH"
  }

  rows[[i]] <- data.frame(
    reference_number = number,
    manuscript_citation = citation,
    manuscript_doi = doi,
    crossref_status = cr_status,
    crossref_title = cr_title,
    crossref_journal = cr_journal,
    crossref_year = cr_year,
    crossref_title_in_citation = title_match,
    pubmed_status = pubmed_status,
    pmid = pmid,
    pubmed_title = pm_title,
    pubmed_journal = pm_journal,
    pubmed_year = pm_year,
    pubmed_doi = pm_doi,
    pubmed_doi_matches_manuscript = doi_consistent,
    retraction_flag = retraction_flag,
    verification_status = status,
    stringsAsFactors = FALSE
  )
  Sys.sleep(0.12)
}

results <- do.call(rbind, rows)
write.csv(results, file.path(output_dir, "reference_verification_crossref_pubmed.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

summary_lines <- c(
  "# Reference verification summary",
  "",
  paste0("- References checked: ", nrow(results)),
  paste0("- Verified in both Crossref and PubMed: ", sum(results$verification_status == "VERIFIED_CROSSREF_PUBMED")),
  paste0("- Verified in Crossref but not found/indexed in PubMed: ", sum(results$verification_status == "VERIFIED_CROSSREF_NOT_IN_PUBMED")),
  paste0("- DOI-less references verified in PubMed: ", sum(is.na(results$manuscript_doi) & results$pubmed_status == "verified")),
  paste0("- Retraction flags: ", sum(results$retraction_flag)),
  "",
  "## Records requiring manual review",
  ""
)
review_rows <- results[!grepl("^VERIFIED_", results$verification_status), ]
if (!nrow(review_rows)) {
  summary_lines <- c(summary_lines, "None.")
} else {
  for (j in seq_len(nrow(review_rows))) {
    summary_lines <- c(summary_lines, paste0("- Reference ", review_rows$reference_number[[j]], ": ", review_rows$verification_status[[j]]))
  }
}
summary_lines <- c(summary_lines, "", "## PubMed-unindexed or unmatched records", "")
not_pm <- results[results$pubmed_status != "verified", ]
if (!nrow(not_pm)) {
  summary_lines <- c(summary_lines, "None.")
} else {
  for (j in seq_len(nrow(not_pm))) {
    summary_lines <- c(summary_lines, paste0("- Reference ", not_pm$reference_number[[j]], ": ", not_pm$manuscript_citation[[j]]))
  }
}
writeLines(summary_lines, file.path(output_dir, "reference_verification_summary.md"), useBytes = TRUE)

print(table(results$verification_status, useNA = "ifany"))
