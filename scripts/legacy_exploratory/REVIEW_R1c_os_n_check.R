# Quick check: GSE39582 OS model n discrepancy (codex n=557 vs my n=551)
suppressPackageStartupMessages({library(dplyr); library(survival)})
PROJ <- "."
DATA <- file.path(PROJ, "data")

expr <- read.csv(file.path(DATA, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(DATA, "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]; probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}
sm_file <- file.path(DATA, "GSE39582", "GSE39582_series_matrix.txt.gz")
con <- gzfile(sm_file, "rt"); field_rows <- list(); sample_names <- NULL
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_characteristics_ch1")) {
    cells <- strsplit(substr(line, nchar("!Sample_characteristics_ch1") + 1, nchar(line)), "\t")[[1]]
    cells <- gsub('"', '', trimws(cells)); field_rows[[length(field_rows) + 1]] <- cells
  } else if (startsWith(line, "!Sample_geo_accession")) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar("!Sample_geo_accession") + 1, nchar(line)), "\t")[[1]])
    sample_names <- sample_names[sample_names != ""]
  }
}
close(con)
field_names <- sapply(field_rows, function(r) { hit <- which(grepl(":", r))[1]; if (!is.na(hit)) sub(":.*", "", r[hit]) else "" })
get_field <- function(fname) {
  idx <- which(field_names == fname); if (length(idx) == 0) return(NULL)
  vals <- field_rows[[idx[1]]]; vals <- vals[vals != ""]
  vals <- sub(paste0(fname, ":"), "", vals, fixed = TRUE); trimws(vals)
}
n <- length(sample_names)
gse <- data.frame(
  sample = sample_names, dataset = get_field("dataset"),
  tnm_t = get_field("tnm.t"), tnm_n = get_field("tnm.n"), tnm_m = get_field("tnm.m"),
  stage = get_field("tnm.stage"),
  os_event = as.numeric(get_field("os.event")), os_delay = as.numeric(get_field("os.delay (months)")),
  age = as.numeric(get_field("age.at.diagnosis (year)")), Sex = get_field("Sex"),
  stringsAsFactors = FALSE)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)
gse_t <- gse[!is.na(gse$dataset) & gse$dataset != "Non Tumoral", ]
FAP13 <- c("FAP","POSTN","THY1","PDPN","TAGLN","ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
zmean <- function(df, genes) { a <- genes[genes %in% colnames(df)]; s <- as.data.frame(lapply(a, function(g) (df[[g]] - mean(df[[g]], na.rm=TRUE))/sd(df[[g]], na.rm=TRUE))); rowMeans(s, na.rm=TRUE) }
gse_t$FAP13 <- zmean(gse_t, FAP13)

cat("stage dist (tumor-only):\n"); print(table(gse_t$stage, useNA = "ifany"))
cat("\nos valid (delay>0, event not NA):",
    sum(!is.na(gse_t$os_delay) & gse_t$os_delay > 0 & !is.na(gse_t$os_event)), "\n")
cat("FAP13 non-NA:", sum(!is.na(gse_t$FAP13)), "\n")
cat("age non-NA:", sum(!is.na(gse_t$age)), "| Sex non-NA:", sum(!is.na(gse_t$Sex)), "\n")

# Model A: stage 1-4 only, require all covariates (my R1b)
m1 <- gse_t[!is.na(gse_t$os_delay) & gse_t$os_delay > 0 & !is.na(gse_t$os_event) &
            !is.na(gse_t$FAP13) & !is.na(gse_t$age) & !is.na(gse_t$Sex) &
            !is.na(gse_t$stage) & gse_t$stage %in% c("1","2","3","4"), ]
cat("\nModel A (stage 1-4, full covars): n =", nrow(m1), "events =", sum(m1$os_event), "\n")

# Model B: include stage 0? codex says exclude stage-0. Keep 1-4.
# Model C: relax age/sex completeness (drop rows only where model needs)
m2 <- gse_t[!is.na(gse_t$os_delay) & gse_t$os_delay > 0 & !is.na(gse_t$os_event) &
            !is.na(gse_t$FAP13) & !is.na(gse_t$age) & !is.na(gse_t$Sex) &
            !is.na(gse_t$tnm_t) & gse_t$tnm_t %in% c("T1","T2","T3","T4"), ]
cat("Model C (T-based, no stage field): n =", nrow(m2), "events =", sum(m2$os_event), "\n")

# codex claims n=557, events=190. Try: stage 0-4 all included as factor
m3 <- gse_t[!is.na(gse_t$os_delay) & gse_t$os_delay > 0 & !is.na(gse_t$os_event) &
            !is.na(gse_t$FAP13) & !is.na(gse_t$age) & !is.na(gse_t$Sex) & !is.na(gse_t$stage), ]
cat("Model D (stage all incl 0/N/A non-NA): n =", nrow(m3), "events =", sum(m3$os_event), "\n")
print(table(m3$stage, useNA = "ifany"))
