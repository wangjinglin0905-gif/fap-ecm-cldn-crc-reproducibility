#!/usr/bin/env Rscript

stopifnot(requireNamespace("curl", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

query <- utils::URLencode(
  "Fibroblast activation protein identifies consensus molecular subtype 4 in colorectal cancer",
  reserved = TRUE
)
url <- paste0(
  "https://api.crossref.org/works?query.title=", query,
  "&rows=10&select=DOI,title,author,published,container-title,volume,page"
)
handle <- curl::new_handle(
  connecttimeout = 30,
  timeout = 180,
  useragent = "FAP-CLDN-JCMM-reference-audit/2.0 (mailto:wangjinglin@gz5055.com)"
)
response <- curl::curl_fetch_memory(url, handle = handle)
stopifnot(response$status_code == 200L)
raw <- rawToChar(response$content)
dir.create("work/revision_round3/reports/reference_audit", recursive = TRUE, showWarnings = FALSE)
writeLines(raw, "work/revision_round3/reports/reference_audit/reference_36_title_query.json", useBytes = TRUE)
items <- jsonlite::fromJSON(raw, simplifyVector = FALSE)$message$items
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
rows <- lapply(seq_along(items), function(i) {
  item <- items[[i]]
  data.frame(
    rank = i,
    doi = item$DOI %||% NA_character_,
    title = paste(unlist(item$title %||% ""), collapse = " "),
    journal = paste(unlist(item$`container-title` %||% ""), collapse = " "),
    volume = item$volume %||% NA_character_,
    pages = item$page %||% NA_character_
  )
})
result <- do.call(rbind, rows)
write.csv(result, "work/revision_round3/reports/reference_audit/reference_36_title_query.csv", row.names = FALSE)
print(result, row.names = FALSE)
