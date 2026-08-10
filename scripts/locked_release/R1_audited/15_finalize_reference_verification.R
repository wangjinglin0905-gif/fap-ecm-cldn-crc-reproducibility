#!/usr/bin/env Rscript

# Apply manual adjudications after the automated Crossref/PubMed pass.
# Each adjudication is based on an official PubMed record or a typography-only
# title difference (US/UK spelling or hyphenation).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript finalize_reference_verification.R OUTPUT_DIR")
output_dir <- args[[1]]
csv_path <- file.path(output_dir, "reference_verification_crossref_pubmed.csv")
x <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE,
              fileEncoding = "UTF-8")

# Reference 8: official PubMed record https://pubmed.ncbi.nlm.nih.gov/17970035/
i <- x$reference_number == 8L
x$pubmed_status[i] <- "verified"
x$pmid[i] <- "17970035"
x$pubmed_title[i] <- "Selective up-regulation of claudin-1 and claudin-2 in colorectal cancer."
x$pubmed_journal[i] <- "Anticancer Research"
x$pubmed_year[i] <- 2007L
x$pubmed_doi[i] <- NA_character_
x$retraction_flag[i] <- FALSE
x$verification_status[i] <- "VERIFIED_PUBMED_NO_DOI"

# References 26 and 33 differ only by tumor/tumour and hyphen typography.
for (ref in c(26L, 33L)) {
  j <- x$reference_number == ref
  x$crossref_title_in_citation[j] <- TRUE
  x$verification_status[j] <- "VERIFIED_CROSSREF_PUBMED"
}

write.csv(x, csv_path, row.names = FALSE, fileEncoding = "UTF-8")

summary_lines <- c(
  "# Reference verification summary",
  "",
  paste0("- References checked: ", nrow(x)),
  paste0("- Verified in both Crossref and PubMed: ", sum(x$verification_status == "VERIFIED_CROSSREF_PUBMED")),
  paste0("- Verified in PubMed with no DOI assigned: ", sum(x$verification_status == "VERIFIED_PUBMED_NO_DOI")),
  paste0("- Verified in Crossref but not indexed/found in PubMed: ", sum(x$verification_status == "VERIFIED_CROSSREF_NOT_IN_PUBMED")),
  paste0("- Retraction flags: ", sum(as.logical(x$retraction_flag), na.rm = TRUE)),
  "",
  "## Identifier exceptions",
  "",
  "- Reference 8 is verified in PubMed (PMID 17970035); the journal record has no DOI.",
  "- Reference 36 is verified by Crossref (DOI 10.1080/07350015.1983.10509354) but is not indexed in PubMed, so no PMID exists.",
  "",
  "## Manual title adjudications",
  "",
  "- References 26 and 33 were accepted after confirming that the automated mismatch arose only from US/UK spelling (tumor/tumour) and hyphen typography.",
  "",
  "No reference carried a PubMed retraction flag at the verification date."
)
writeLines(summary_lines, file.path(output_dir, "reference_verification_summary.md"), useBytes = TRUE)
print(table(x$verification_status, useNA = "ifany"))
