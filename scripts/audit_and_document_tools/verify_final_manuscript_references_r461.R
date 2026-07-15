#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
stopifnot(requireNamespace("data.table", quietly = TRUE))

manuscript_file <- "work/revision_round3/final/JCMM_Manuscript_Revised_Round3.md"
metadata_file <- "work/revision_round3/reports/reference_audit/all_reference_verification.csv"
out_dir <- "work/revision_round3/reports/reference_audit"

normalize_text <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^a-z0-9]+", " ", x))
  trimws(gsub("[[:space:]]+", " ", x))
}

lines <- readLines(manuscript_file, encoding = "UTF-8", warn = FALSE)
reference_heading <- which(lines == "## References")
stopifnot(length(reference_heading) == 1L)
body_lines <- lines[seq_len(reference_heading - 1L)]
reference_lines <- lines[(reference_heading + 1L):length(lines)]
reference_lines <- reference_lines[grepl("^[0-9]+\\. .+doi:10\\.", reference_lines)]

parse_reference <- function(text) {
  number <- as.integer(sub("^([0-9]+)\\..*$", "\\1", text))
  doi <- tolower(sub("^.*doi:([^[:space:]]+).*$", "\\1", text, ignore.case = TRUE))
  doi <- sub("[.;,]+$", "", doi)
  first_author <- sub("^[0-9]+\\.\\s+([^[:space:]]+).*$", "\\1", text)
  data.frame(
    final_reference_number = number,
    final_reference = text,
    cited_doi = doi,
    final_first_author = first_author
  )
}
final_refs <- data.table::rbindlist(lapply(reference_lines, parse_reference))
metadata <- data.table::fread(metadata_file, data.table = FALSE)
metadata <- metadata[, c(
  "cited_doi", "resolved_doi", "authoritative_title", "authoritative_first_author",
  "journal", "year", "volume", "issue", "pages_or_article", "pmid", "crossref_url"
)]
metadata <- metadata[!duplicated(metadata$cited_doi), ]
audit <- merge(final_refs, metadata, by = "cited_doi", all.x = TRUE)
audit <- audit[order(audit$final_reference_number), ]
audit$doi_resolved <- !is.na(audit$resolved_doi) & audit$resolved_doi == audit$cited_doi
audit$first_author_match <- normalize_text(audit$final_first_author) ==
  normalize_text(audit$authoritative_first_author)

# The TCGA consortium record does not expose a conventional Crossref first author.
audit$first_author_match[audit$final_reference_number == 2L &
                           grepl("Cancer Genome Atlas Network", audit$final_reference)] <- TRUE

data.table::fwrite(audit, file.path(out_dir, "final_manuscript_reference_verification.csv"), na = "")

expand_citation <- function(token) {
  token <- gsub("[", "", token, fixed = TRUE)
  token <- gsub("]", "", token, fixed = TRUE)
  token <- gsub(" ", "", token, fixed = TRUE)
  pieces <- strsplit(token, ",", fixed = TRUE)[[1]]
  values <- integer()
  for (piece in pieces) {
    if (grepl("-", piece, fixed = TRUE)) {
      bounds <- as.integer(strsplit(piece, "-", fixed = TRUE)[[1]])
      if (length(bounds) == 2L && all(is.finite(bounds))) {
        values <- c(values, seq.int(bounds[1], bounds[2]))
      }
    } else {
      value <- suppressWarnings(as.integer(piece))
      if (is.finite(value)) values <- c(values, value)
    }
  }
  values
}

matches <- regmatches(
  body_lines,
  gregexpr("\\[[0123456789,-]+\\]", body_lines, perl = TRUE)
)
citation_sequence <- unlist(lapply(matches, function(tokens) {
  if (length(tokens) == 0L) return(integer())
  unlist(lapply(tokens, expand_citation))
}))
first_appearance <- unique(citation_sequence)
expected <- seq_len(nrow(final_refs))
missing_citations <- setdiff(expected, citation_sequence)
extra_citations <- setdiff(citation_sequence, expected)
first_appearance_in_order <- identical(first_appearance, sort(first_appearance))

summary <- data.frame(
  metric = c(
    "references_in_list", "unique_dois", "doi_resolved", "first_author_matches",
    "references_cited_in_text", "missing_in_text", "out_of_range_citations",
    "first_appearance_in_numeric_order"
  ),
  value = c(
    nrow(final_refs), length(unique(final_refs$cited_doi)),
    sum(audit$doi_resolved %in% TRUE), sum(audit$first_author_match %in% TRUE),
    length(unique(citation_sequence)), length(missing_citations), length(extra_citations),
    first_appearance_in_order
  )
)
data.table::fwrite(summary, file.path(out_dir, "final_manuscript_reference_summary.csv"))
writeLines(
  c(
    paste0("First-appearance sequence: ", paste(first_appearance, collapse = ",")),
    paste0("Missing citations: ", paste(missing_citations, collapse = ",")),
    paste0("Out-of-range citations: ", paste(extra_citations, collapse = ","))
  ),
  file.path(out_dir, "final_manuscript_citation_sequence.txt")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "final_reference_check_sessionInfo.txt"))

print(summary, row.names = FALSE)
if (any(!(audit$doi_resolved %in% TRUE) | !(audit$first_author_match %in% TRUE))) {
  print(audit[!(audit$doi_resolved %in% TRUE) | !(audit$first_author_match %in% TRUE),
              c("final_reference_number", "cited_doi", "final_first_author",
                "authoritative_first_author", "doi_resolved", "first_author_match")])
}
