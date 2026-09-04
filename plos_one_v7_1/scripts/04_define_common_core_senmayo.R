#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop(paste(
    "Usage: Rscript 04_define_common_core_senmayo.R",
    "<GSE132 gene_availability.csv> <GSE166 analysis_ledger.json>",
    "<GSE280 H5> <GSE334 features.tsv.gz> <Valdeolivas features.tsv.gz>",
    "<SenMayo source txt> <output directory>"
  ))
}

suppressPackageStartupMessages({
  library(BPCells)
  library(jsonlite)
})

paths <- normalizePath(args[1:6], mustWork = TRUE)
out_dir <- args[[7]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

availability <- read.csv(paths[[1]], check.names = FALSE)
row_132 <- availability[availability$set == "SenMayo_nonoverlap_represented", , drop = FALSE]
if (nrow(row_132) != 1L) stop("GSE132 represented-gene row not found uniquely")
g132 <- strsplit(row_132$genes[[1]], ";", fixed = TRUE)[[1]]

ledger_166 <- jsonlite::fromJSON(paths[[2]], simplifyVector = TRUE)
g166 <- unname(ledger_166$SenMayo119_available)

sen_source <- unique(trimws(readLines(paths[[6]], warn = FALSE)))
sen_source <- sen_source[nzchar(sen_source)]
fap13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2", "MMP9",
           "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
sen_nonoverlap <- setdiff(sen_source, fap13)

g280 <- trimws(BPCells:::read_hdf5_string_cpp(paths[[3]], "matrix/features/name", 16384L))

read_symbols <- function(path) {
  x <- read.delim(gzfile(path), header = FALSE, sep = "\t", quote = "", comment.char = "",
                  stringsAsFactors = FALSE)
  if (ncol(x) < 2L) stop("Feature file has fewer than two columns: ", path)
  trimws(x[[2]])
}

g334 <- read_symbols(paths[[4]])
gval <- read_symbols(paths[[5]])

represented <- list(
  GSE132465 = intersect(sen_nonoverlap, unique(g132)),
  GSE166555 = intersect(sen_nonoverlap, unique(g166)),
  GSE280315 = intersect(sen_nonoverlap, unique(g280)),
  GSE334323 = intersect(sen_nonoverlap, unique(g334)),
  Valdeolivas_Visium = intersect(sen_nonoverlap, unique(gval))
)
common_core <- Reduce(intersect, represented)

coverage <- do.call(rbind, lapply(names(represented), function(nm) {
  g <- represented[[nm]]
  data.frame(
    cohort = nm,
    source_nonoverlap_n = length(sen_nonoverlap),
    represented_n = length(g),
    common_core_n = length(common_core),
    missing_from_source = paste(setdiff(sen_nonoverlap, g), collapse = ";"),
    removed_to_common_core = paste(setdiff(g, common_core), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

writeLines(common_core, file.path(out_dir, "SenMayo_common_core_genes.txt"), useBytes = TRUE)
write.csv(coverage, file.path(out_dir, "SenMayo_common_core_coverage.csv"), row.names = FALSE)
write.csv(
  data.frame(order = seq_along(common_core), gene = common_core),
  file.path(out_dir, "SenMayo_common_core_gene_ledger.csv"), row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_common_core_definition.txt"))

cat("Common core:", length(common_core), "genes\n")
print(coverage, row.names = FALSE)
