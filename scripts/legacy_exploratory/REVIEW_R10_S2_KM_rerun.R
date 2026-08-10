# =============================================================================
# REVIEW_R10: S2 KM re-run on the 566-tumor analysis set (GSE39582)
# v6.0 issue 3.2: S2a n=579 / S2b n=574 vs text "566 tumor samples"
# Fix: exclude 19 Non Tumoral arrays, re-run OS/RFS KM by median FAP-CAF score
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(survival)})
PROJ <- "."
DATA <- file.path(PROJ, "data")
OUT  <- file.path(PROJ, "output", "review_r10")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- 1. expression matrix (gene-level, same as R1b pipeline) ----
expr <- read.csv(file.path(DATA, "GSE39582", "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(DATA, "GSE39582", "gene2probe.csv"))
sample_ids <- as.character(expr[[1]])
gene_expr <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (g in unique(probe_map$gene)) {
  probes <- probe_map$probe[probe_map$gene == g]
  probes <- probes[probes %in% colnames(expr)]
  if (length(probes) > 0) gene_expr[[g]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}

# ---- 2. phenotype from series matrix (streaming parse) ----
sm_file <- file.path(DATA, "GSE39582", "GSE39582_series_matrix.txt.gz")
con <- gzfile(sm_file, "rt")
field_rows <- list(); sample_names <- NULL
while (TRUE) {
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) break
  if (startsWith(line, "!Sample_characteristics_ch1")) {
    cells <- strsplit(substr(line, nchar("!Sample_characteristics_ch1") + 1, nchar(line)), "\t")[[1]]
    cells <- gsub('"', '', trimws(cells))
    field_rows[[length(field_rows) + 1]] <- cells
  } else if (startsWith(line, "!Sample_geo_accession")) {
    sample_names <- gsub('"', '', strsplit(substr(line, nchar("!Sample_geo_accession") + 1, nchar(line)), "\t")[[1]])
    sample_names <- sample_names[sample_names != ""]
  }
}
close(con)
field_names <- sapply(field_rows, function(r) {
  hit <- which(grepl(":", r))[1]
  if (!is.na(hit)) sub(":.*", "", r[hit]) else ""
})
get_field <- function(fname) {
  idx <- which(field_names == fname)
  if (length(idx) == 0) return(NULL)
  vals <- field_rows[[idx[1]]]
  vals <- vals[vals != ""]
  vals <- sub(paste0(fname, ":"), "", vals, fixed = TRUE)
  trimws(vals)
}
gse <- data.frame(
  sample   = sample_names,
  dataset  = get_field("dataset"),
  os_event = suppressWarnings(as.numeric(get_field("os.event"))),
  os_delay = suppressWarnings(as.numeric(get_field("os.delay (months)"))),
  rfs_event = suppressWarnings(as.numeric(get_field("rfs.event"))),
  rfs_delay = suppressWarnings(as.numeric(get_field("rfs.delay"))),
  stringsAsFactors = FALSE
)
gse <- merge(gse, gene_expr, by = "sample", all.x = TRUE)

# ---- 3. FAP-CAF score (17-gene z-mean, as in v3.x KM figures) ----
FAP17 <- c("FAP","COL1A1","COL1A2","COL3A1","FN1","POSTN","THY1","PDPN","TAGLN",
           "ACTA2","MMP2","MMP9","CXCL12","TGFB1","INHBA","WNT2","WNT5A")
zmean <- function(df, genes) {
  a <- genes[genes %in% colnames(df)]
  s <- as.data.frame(lapply(a, function(g) (df[[g]] - mean(df[[g]], na.rm = TRUE)) / sd(df[[g]], na.rm = TRUE)))
  rowMeans(s, na.rm = TRUE)
}
gse$FAP_CAF <- zmean(gse, FAP17)

# ---- 4. TUMOR-ONLY analysis set (exclude Non Tumoral) ----
gse_t <- gse[!is.na(gse$dataset) & gse$dataset != "Non Tumoral", ]
cat("Total arrays:", nrow(gse), "| tumor samples:", nrow(gse_t), "\n")

run_km <- function(df, time, event, label) {
  d <- df[!is.na(df[[time]]) & df[[time]] > 0 & !is.na(df[[event]]) & !is.na(df$FAP_CAF), ]
  med <- median(d$FAP_CAF)
  d$grp <- factor(ifelse(d$FAP_CAF >= med, "FAP-CAF high", "FAP-CAF low"),
                  levels = c("FAP-CAF low", "FAP-CAF high"))
  fit <- survfit(Surv(d[[time]], d[[event]]) ~ grp, data = d)
  sd_ <- survdiff(Surv(d[[time]], d[[event]]) ~ grp, data = d)
  p <- 1 - pchisq(sd_$chisq, 1)
  cox <- coxph(Surv(d[[time]], d[[event]]) ~ grp, data = d)
  hr <- exp(coef(cox)[1]); ci <- exp(confint(cox)[1, ])
  cat(sprintf("[%s] n=%d (low=%d, high=%d) | events low=%d high=%d | log-rank P=%.3f | HR(high vs low)=%.2f (%.2f-%.2f)\n",
              label, nrow(d), sum(d$grp == "FAP-CAF low"), sum(d$grp == "FAP-CAF high"),
              sum(d[[event]][d$grp == "FAP-CAF low"]), sum(d[[event]][d$grp == "FAP-CAF high"]),
              p, hr, ci[1], ci[2]))
  list(data = d, fit = fit, p = p, hr = hr, ci = ci, n_low = sum(d$grp == "FAP-CAF low"),
       n_high = sum(d$grp == "FAP-CAF high"))
}

os <- run_km(gse_t, "os_delay", "os_event", "OS tumor-only")
rfs <- run_km(gse_t, "rfs_delay", "rfs_event", "RFS tumor-only")

# also run on FULL array set (old v3.x behavior) for comparison
os_full <- run_km(gse, "os_delay", "os_event", "OS full-array (old)")
rfs_full <- run_km(gse, "rfs_delay", "rfs_event", "RFS full-array (old)")

# ---- 5. redraw KM figures (PNG 300dpi + TIFF 600dpi) ----
draw_km <- function(km, time_unit, title, file_base, xmax = 150) {
  d <- km$data
  col <- c("FAP-CAF low" = "#2166AC", "FAP-CAF high" = "#B2182B")
  for (fmt in c("png", "tiff")) {
    if (fmt == "png") {
      png(file.path(OUT, paste0(file_base, ".png")), width = 2100, height = 1950, res = 300)
    } else {
      tiff(file.path(OUT, paste0(file_base, ".tiff")), width = 7, height = 6.5,
           units = "in", res = 600, compression = "lzw")
    }
    par(mar = c(8.5, 5.2, 3, 1.5), mgp = c(3.2, 0.8, 0), xpd = NA)
    plot(km$fit, col = col, lwd = 2.2, xlab = "", ylab = "", xaxt = "n",
         xlim = c(0, xmax), main = title, cex.main = 0.95)
    axis(1, at = seq(0, xmax, 25))
    title(xlab = paste0("Time (", time_unit, ")"), ylab = "Survival probability")
    legend("topright", legend = c(paste0("FAP-CAF low (n=", km$n_low, ")"),
                                  paste0("FAP-CAF high (n=", km$n_high, ")")),
           col = col, lwd = 2.2, bty = "n", cex = 0.9)
    txt <- sprintf("log-rank P = %.2f\nHR (high vs low) = %.2f (95%% CI %.2f\u2013%.2f)",
                   km$p, km$hr, km$ci[1], km$ci[2])
    text(xmax * 0.30, 0.16, txt, adj = c(0, 0.5), cex = 0.9)
    # risk table below x-axis
    times <- seq(0, xmax, 25)
    nr <- summary(km$fit, times = times, extend = TRUE)
    strata_levels <- levels(nr$strata)
    mtext("Number at risk", side = 1, line = 5.2, at = -12, cex = 0.8, adj = 1)
    for (si in seq_along(strata_levels)) {
      idx <- which(as.numeric(nr$strata) == si | nr$strata == strata_levels[si])
      lab <- gsub("grp=", "", strata_levels[si])
      mtext(lab, side = 1, line = 5.2 + si * 1.1, at = -12, cex = 0.75, adj = 1, col = col[lab])
      vals <- nr$n.risk[idx]
      mtext(vals, side = 1, line = 5.2 + si * 1.1, at = times, cex = 0.75, col = col[lab])
    }
    dev.off()
  }
  cat("saved:", file_base, "\n")
}

draw_km(os, "months", sprintf("GSE39582 overall survival by FAP-CAF score (tumor samples, n = %d)", nrow(os$data)), "S2a_GSE39582_OS_KM_v2")
draw_km(rfs, "months", sprintf("GSE39582 relapse-free survival by FAP-CAF score (tumor samples, n = %d)", nrow(rfs$data)), "S2b_GSE39582_RFS_KM_v2")

# summary csv
summ <- data.frame(
  analysis = c("OS tumor-only", "RFS tumor-only", "OS full-array", "RFS full-array"),
  n = c(nrow(os$data), nrow(rfs$data), nrow(os_full$data), nrow(rfs_full$data)),
  n_low = c(os$n_low, rfs$n_low, os_full$n_low, rfs_full$n_low),
  n_high = c(os$n_high, rfs$n_high, os_full$n_high, rfs_full$n_high),
  logrank_P = c(os$p, rfs$p, os_full$p, rfs_full$p),
  HR = c(os$hr, rfs$hr, os_full$hr, rfs_full$hr),
  HR_lo = c(os$ci[1], rfs$ci[1], os_full$ci[1], rfs_full$ci[1]),
  HR_hi = c(os$ci[2], rfs$ci[2], os_full$ci[2], rfs_full$ci[2])
)
write.csv(summ, file.path(OUT, "S2_KM_rerun_summary.csv"), row.names = FALSE)
print(summ)
cat("\n=== [DONE] REVIEW_R10 S2 KM re-run ===\n")
