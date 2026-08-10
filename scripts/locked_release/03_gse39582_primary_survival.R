options(stringsAsFactors = FALSE, warn = 1)
set.seed(2026080606)
suppressPackageStartupMessages(library(survival))

root <- normalizePath(getwd(), winslash = "/")
public_dir <- file.path(root, "data", "public", "GSE39582")
out_dir <- file.path(root, "results", "analysis")
expr <- read.csv(file.path(public_dir, "GSE39582_expr.csv"), check.names = FALSE)
probe_map <- read.csv(file.path(public_dir, "gene2probe.csv"), check.names = FALSE)
meta <- read.csv(file.path(out_dir, "GSE39582_parsed_metadata.csv"), check.names = FALSE)

FAP13 <- c("FAP", "POSTN", "THY1", "PDPN", "TAGLN", "ACTA2", "MMP2",
           "MMP9", "CXCL12", "TGFB1", "INHBA", "WNT2", "WNT5A")
matrix4 <- c("COL1A1", "COL1A2", "COL3A1", "FN1")
receptor2 <- c("SDC4", "CD44")
gene_expr <- data.frame(sample = as.character(expr[[1]]))
for (gene in unique(probe_map$gene)) {
  probes <- intersect(probe_map$probe[probe_map$gene == gene], names(expr))
  if (length(probes)) gene_expr[[gene]] <- rowMeans(expr[, probes, drop = FALSE], na.rm = TRUE)
}
dat <- merge(meta, gene_expr, by = "sample", all.x = TRUE, sort = FALSE)
dat <- dat[!is.na(dat$dataset) & dat$dataset != "Non Tumoral", ]
score <- function(x, genes) rowMeans(scale(x[, genes, drop = FALSE]), na.rm = TRUE)
dat$FAP13 <- score(dat, FAP13); dat$matrix4 <- score(dat, matrix4)
dat$receptor2 <- score(dat, receptor2); dat$FAP17 <- score(dat, c(FAP13, matrix4))
stopifnot(nrow(dat) == 566L)

cor_one <- function(y, seed) {
  set.seed(seed); x <- dat$FAP13
  boot <- replicate(5000L, { i <- sample.int(nrow(dat), replace = TRUE); suppressWarnings(cor(x[i], y[i], method = "spearman")) })
  test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(n = nrow(dat), rho = unname(test$estimate),
             ci_low = unname(quantile(boot, 0.025, type = 6, na.rm = TRUE)),
             ci_high = unname(quantile(boot, 0.975, type = 6, na.rm = TRUE)), p_value = test$p.value)
}
corr <- rbind(cbind(comparison = "FAP13 vs matrix4", cor_one(dat$matrix4, 2026080613L)),
              cbind(comparison = "FAP13 vs receptor2", cor_one(dat$receptor2, 2026080614L)))
corr$fdr_bh_two <- p.adjust(corr$p_value, method = "BH")
write.csv(corr, file.path(out_dir, "GSE39582_primary_correlations.csv"), row.names = FALSE)

run_survival <- function(time, event, endpoint) {
  keep <- is.finite(time) & is.finite(event) & is.finite(dat$FAP17) & time > 0
  d <- dat[keep, ]; threshold <- median(d$FAP17)
  d$group <- factor(ifelse(d$FAP17 >= threshold, "high", "low"), levels = c("low", "high"))
  fit <- coxph(Surv(time[keep], event[keep]) ~ group, data = d)
  lr <- survdiff(Surv(time[keep], event[keep]) ~ group, data = d)
  ci <- exp(confint(fit))[1, ]
  data.frame(endpoint = endpoint, n = nrow(d), n_low = sum(d$group == "low"),
             n_high = sum(d$group == "high"), events = sum(event[keep] == 1),
             median_cutpoint = threshold, HR_high_vs_low = exp(coef(fit))[1],
             ci_low = ci[1], ci_high = ci[2],
             logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE))
}
km <- rbind(run_survival(dat$os_time, dat$os_event, "Overall survival"),
            run_survival(dat$rfs_time, dat$rfs_event, "Relapse-free survival"))
write.csv(km, file.path(out_dir, "GSE39582_tumor_only_KM_final.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "GSE39582_sessionInfo.txt"))
