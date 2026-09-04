#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript 01_sc_cluster_and_correlation_lock.R",
    "<GSE132465 cell_scores.csv> <GSE132465 patient_scores.csv>",
    "<GSE166555 patient_compartment_scores.csv> <output_dir>"
  ))
}

cell_path <- normalizePath(args[[1]], mustWork = TRUE)
gse132_patient_path <- normalizePath(args[[2]], mustWork = TRUE)
gse166_patient_path <- normalizePath(args[[3]], mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

seed_cluster <- 2026090301L
seed_correlation <- 2026090302L
n_cluster_boot <- 5000L
n_correlation_boot <- 10000L

coef_table <- function(fit, model_id, term) {
  tab <- as.data.frame(summary(fit)$coefficients)
  if (!term %in% rownames(tab)) stop("Term not found in ", model_id, ": ", term)
  ci <- suppressMessages(confint(fit, parm = term, method = "Wald"))
  data.frame(
    model_id = model_id,
    estimator = "random-intercept lmer; Satterthwaite inference",
    term = term,
    estimate = tab[term, "Estimate"],
    std_error = tab[term, "Std. Error"],
    df = tab[term, "df"],
    statistic = tab[term, "t value"],
    p_value = tab[term, "Pr(>|t|)"],
    ci_low = ci[1, 1],
    ci_high = ci[1, 2],
    stringsAsFactors = FALSE
  )
}

# Fast patient-cluster bootstrap for a patient fixed-effects estimator. The
# within-patient model is represented by patient-specific cross-products, so a
# bootstrap draw only has to sum 23 small matrices rather than refit cell rows.
cluster_fe_bootstrap <- function(formula, data, patient, term, model_id,
                                 draws = n_cluster_boot, seed = seed_cluster) {
  mf <- model.frame(formula, data = data, na.action = na.omit)
  used <- as.integer(rownames(mf))
  pid <- droplevels(factor(patient[used]))
  tt <- terms(formula)
  y <- model.response(mf)
  x <- model.matrix(tt, mf)
  if ("(Intercept)" %in% colnames(x)) x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  if (!term %in% colnames(x)) stop("Term not present in model matrix: ", term)

  # Apply the within transform separately to outcome and every design column.
  y_w <- y - ave(y, pid, FUN = mean)
  x_w <- x
  for (j in seq_len(ncol(x_w))) x_w[, j] <- x_w[, j] - ave(x_w[, j], pid, FUN = mean)

  patients <- levels(pid)
  xtx <- vector("list", length(patients))
  xty <- vector("list", length(patients))
  for (i in seq_along(patients)) {
    keep <- pid == patients[[i]]
    xp <- x_w[keep, , drop = FALSE]
    yp <- y_w[keep]
    xtx[[i]] <- crossprod(xp)
    xty[[i]] <- crossprod(xp, yp)
  }

  solve_beta <- function(indices) {
    xx <- Reduce(`+`, xtx[indices])
    xy <- Reduce(`+`, xty[indices])
    beta <- tryCatch(qr.solve(xx, xy), error = function(e) rep(NA_real_, ncol(x_w)))
    names(beta) <- colnames(x_w)
    unname(beta[[term]])
  }

  original <- solve_beta(seq_along(patients))
  set.seed(seed)
  boot <- numeric(draws)
  for (b in seq_len(draws)) {
    boot[[b]] <- solve_beta(sample.int(length(patients), length(patients), replace = TRUE))
  }
  valid <- is.finite(boot)
  boot <- boot[valid]
  ci <- unname(quantile(boot, c(0.025, 0.975), type = 6, na.rm = TRUE))
  p_sign <- min(1, 2 * min(
    (1 + sum(boot <= 0)) / (length(boot) + 1),
    (1 + sum(boot >= 0)) / (length(boot) + 1)
  ))
  list(
    summary = data.frame(
      model_id = model_id,
      estimator = "patient fixed effects with patient-cluster percentile bootstrap",
      term = term,
      estimate = original,
      std_error = sd(boot),
      df = length(patients) - 1L,
      statistic = NA_real_,
      p_value = p_sign,
      ci_low = ci[[1]],
      ci_high = ci[[2]],
      clusters = length(patients),
      bootstrap_requested = draws,
      bootstrap_valid = length(boot),
      stringsAsFactors = FALSE
    ),
    draws = data.frame(model_id = model_id, draw = seq_along(boot), estimate = boot)
  )
}

cell <- read.csv(cell_path, check.names = FALSE)
fib <- subset(cell, compartment == "Fibroblast")
rownames(fib) <- seq_len(nrow(fib))
fib$Patient <- factor(fib$Patient)
fib$FAP_status <- factor(fib$FAP_status, levels = c("FAP-", "FAP+"))
fib$Cell_subtype <- factor(fib$Cell_subtype)

if (nrow(fib) != 1501L || nlevels(fib$Patient) != 23L) {
  warning("Expected 1,501 tumour fibroblasts from 23 patients; observed ",
          nrow(fib), " cells and ", nlevels(fib$Patient), " patients")
}

specs <- list(
  list(
    id = "GSE132_FAPdet_SenMayo_depth_subtype",
    formula = SenMayo_zmean ~ FAP_status + MKI67_z + log_nCount_RNA + Cell_subtype,
    term = "FAP_statusFAP+"
  ),
  list(
    id = "GSE132_FAPdet_SASP_depth_subtype",
    formula = SASP_zmean ~ FAP_status + MKI67_z + log_nCount_RNA + Cell_subtype,
    term = "FAP_statusFAP+"
  ),
  list(
    id = "GSE132_matrix4_on_FAP_depth_subtype",
    formula = matrix4_zmean ~ FAP_expr + log_nCount_RNA + Cell_subtype,
    term = "FAP_expr"
  )
)

lmer_rows <- list()
cluster_rows <- list()
cluster_draws <- list()
for (i in seq_along(specs)) {
  s <- specs[[i]]
  mixed_formula <- as.formula(paste(deparse(s$formula), "+ (1 | Patient)"))
  fit <- lmerTest::lmer(
    mixed_formula, data = fib, REML = FALSE,
    control = lmerControl(optimizer = "bobyqa")
  )
  lmer_row <- coef_table(fit, s$id, s$term)
  lmer_row$clusters <- nlevels(fib$Patient)
  lmer_row$bootstrap_requested <- NA_integer_
  lmer_row$bootstrap_valid <- NA_integer_
  lmer_rows[[i]] <- lmer_row

  cb <- cluster_fe_bootstrap(
    s$formula, fib, fib$Patient, s$term, s$id,
    draws = n_cluster_boot, seed = seed_cluster + i
  )
  cluster_rows[[i]] <- cb$summary
  cluster_draws[[i]] <- cb$draws
}

model_results <- rbind(do.call(rbind, lmer_rows), do.call(rbind, cluster_rows))
write.csv(model_results, file.path(out_dir, "gse132465_cluster_inference.csv"), row.names = FALSE)
gz <- gzfile(file.path(out_dir, "gse132465_cluster_bootstrap_draws.csv.gz"), "wt")
write.csv(do.call(rbind, cluster_draws), gz, row.names = FALSE)
close(gz)

# Build the aligned patient-level plotting source. SenMayo uses gene-wise zmean
# so the figure and headline estimand live on the same scale.
aggregate_patient_compartment <- function(df) {
  aggregate(
    cbind(SenMayo_zmean, SASP_zmean, MKI67) ~ Patient + compartment,
    data = df, FUN = mean
  )
}

aligned132 <- aggregate_patient_compartment(cell)
names(aligned132)[names(aligned132) == "Patient"] <- "patient"
aligned132$cohort <- "GSE132465"
aligned132$SenMayo_score <- aligned132$SenMayo_zmean
aligned132$SASP_score <- aligned132$SASP_zmean
aligned132$MKI67_score <- aligned132$MKI67
aligned132 <- aligned132[, c("cohort", "patient", "compartment", "SenMayo_score", "SASP_score", "MKI67_score")]

aligned166 <- read.csv(gse166_patient_path, check.names = FALSE)
aligned166$cohort <- "GSE166555"
aligned166$SenMayo_score <- aligned166$SenMayo119
aligned166$SASP_score <- aligned166$SASP25
aligned166$MKI67_score <- aligned166$MKI67
aligned166 <- aligned166[, c("cohort", "patient", "compartment", "SenMayo_score", "SASP_score", "MKI67_score")]
aligned <- rbind(aligned132, aligned166)
write.csv(aligned, file.path(out_dir, "scrna_patient_compartment_scores_aligned.csv"), row.names = FALSE)

paired_rows <- list()
k <- 1L
for (cc in unique(aligned$cohort)) {
  for (ep in c("SenMayo_score", "SASP_score", "MKI67_score")) {
    d <- aligned[aligned$cohort == cc, , drop = FALSE]
    wide <- reshape(d[, c("patient", "compartment", ep)], idvar = "patient",
                    timevar = "compartment", direction = "wide")
    fcol <- paste0(ep, ".Fibroblast")
    ecol <- paste0(ep, ".Epithelial")
    keep <- complete.cases(wide[, c(fcol, ecol)])
    delta <- wide[[fcol]][keep] - wide[[ecol]][keep]
    wt <- suppressWarnings(wilcox.test(delta, mu = 0, exact = length(delta) < 50L, correct = FALSE))
    paired_rows[[k]] <- data.frame(
      cohort = cc, endpoint = ep, n_patients = length(delta),
      median_fibroblast = median(wide[[fcol]][keep]),
      median_epithelial = median(wide[[ecol]][keep]),
      median_paired_difference = median(delta),
      positive_direction = sum(delta > 0),
      p_two_sided = wt$p.value,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}
write.csv(do.call(rbind, paired_rows), file.path(out_dir, "scrna_paired_compartment_summary_aligned.csv"), row.names = FALSE)

safe_rho <- function(x, y) {
  keep <- complete.cases(x, y)
  if (sum(keep) < 4L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

correlation_contrasts <- function(df, cohort, x_name, m_name, s_name,
                                  draws = n_correlation_boot, seed = seed_correlation) {
  x <- df[[x_name]]
  m <- df[[m_name]]
  s <- df[[s_name]]
  keep <- complete.cases(x, m, s)
  x <- x[keep]; m <- m[keep]; s <- s[keep]
  n <- length(x)
  observed <- c(
    rho_FAP_matrix = safe_rho(x, m),
    rho_FAP_SenMayo = safe_rho(x, s),
    rho_SenMayo_matrix = safe_rho(s, m)
  )
  contrasts <- c(
    FAP_matrix_minus_FAP_SenMayo = observed[[1]] - observed[[2]],
    FAP_matrix_minus_SenMayo_matrix = observed[[1]] - observed[[3]]
  )

  set.seed(seed)
  boots <- matrix(NA_real_, nrow = draws, ncol = 5L)
  colnames(boots) <- c(names(observed), names(contrasts))
  for (b in seq_len(draws)) {
    idx <- sample.int(n, n, replace = TRUE)
    rb <- c(safe_rho(x[idx], m[idx]), safe_rho(x[idx], s[idx]), safe_rho(s[idx], m[idx]))
    boots[b, ] <- c(rb, rb[[1]] - rb[[2]], rb[[1]] - rb[[3]])
  }

  rows <- list()
  j <- 1L
  for (nm in colnames(boots)) {
    vals <- boots[, nm]
    vals <- vals[is.finite(vals)]
    est <- if (nm %in% names(observed)) observed[[nm]] else contrasts[[nm]]
    ci <- unname(quantile(vals, c(0.025, 0.975), type = 6, na.rm = TRUE))
    p_sign <- min(1, 2 * min(
      (1 + sum(vals <= 0)) / (length(vals) + 1),
      (1 + sum(vals >= 0)) / (length(vals) + 1)
    ))
    rows[[j]] <- data.frame(
      cohort = cohort,
      estimand = nm,
      n_patients = n,
      estimate = est,
      ci_low = ci[[1]],
      ci_high = ci[[2]],
      bootstrap_sign_p = p_sign,
      bootstrap_requested = draws,
      bootstrap_valid = length(vals),
      stringsAsFactors = FALSE
    )
    j <- j + 1L
  }
  list(
    summary = do.call(rbind, rows),
    draws = data.frame(cohort = cohort, draw = seq_len(draws), boots, check.names = FALSE)
  )
}

g132 <- read.csv(gse132_patient_path, check.names = FALSE)
c132 <- correlation_contrasts(
  g132, "GSE132465", "FAP_expr", "matrix4_mean", "SenMayo_mean",
  seed = seed_correlation + 1L
)
g166all <- read.csv(gse166_patient_path, check.names = FALSE)
g166 <- subset(g166all, compartment == "Fibroblast")
c166 <- correlation_contrasts(
  g166, "GSE166555", "FAP13", "matrix4", "SenMayo119",
  seed = seed_correlation + 2L
)

write.csv(rbind(c132$summary, c166$summary),
          file.path(out_dir, "dependent_correlation_contrasts.csv"), row.names = FALSE)
gz <- gzfile(file.path(out_dir, "dependent_correlation_bootstrap_draws.csv.gz"), "wt")
write.csv(rbind(c132$draws, c166$draws), gz, row.names = FALSE)
close(gz)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo_sc_cluster_correlation.txt"))
cat("Completed single-cell cluster and dependent-correlation lock analyses.\n")
print(model_results, digits = 6, row.names = FALSE)
print(rbind(c132$summary, c166$summary), digits = 6, row.names = FALSE)
