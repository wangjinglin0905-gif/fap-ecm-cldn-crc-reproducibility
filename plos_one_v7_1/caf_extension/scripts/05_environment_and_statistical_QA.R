# Scoped environment check: do not enumerate/load the user's R package library.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1], winslash="/", mustWork=TRUE)
extra_library <- Sys.getenv("FAP_R_LIBRARY", "")
if (nzchar(extra_library)) .libPaths(c(extra_library, .libPaths()))
pkgs <- c("ggplot2", "patchwork", "ragg", "svglite")
for (p in pkgs) {
  cat(p, as.character(packageVersion(p)), "at", find.package(p), "\n")
  loadNamespace(p)
}
capture.output(sessionInfo(), file=file.path(root,"qa/R_figure_environment.txt"))
res <- file.path(root,"results/adapted_v02")
allp <- read.csv(file.path(res,"patient_associations.csv"),check.names=FALSE)
s <- read.csv(file.path(res,"primary_single_cell_results.csv"),check.names=FALSE)
checks <- list()
for (i in seq_len(nrow(s))) {
  r <- s[i,]
  cells <- read.csv(file.path(res,paste0(r$cohort,"_cell_program_scores.csv")),check.names=FALSE)
  ps <- allp[allp$cohort==r$cohort & allp$program==r$program & allp$variant=="purged" & allp$model=="adjusted",]
  recomputed <- numeric(nrow(ps))
  for (j in seq_len(nrow(ps))) {
    d <- cells[cells$patient==ps$patient[j],]
    # Reconstruct MKI67 normalized counts from the exported metadata when present.
    # GSE166555 gets exact normalized MKI67 in the compact numerical QA input.
    gene <- read.csv(file.path(root,"derived",paste0(r$cohort,"MKI67_QA.csv")))
    d$mki <- gene$MKI67[match(d$cell,gene$cell)]
    xx <- rank(d[[paste0(r$program,"_purged")]])
    yy <- rank(d$SenMayo)
    cov <- cbind(1,rank(log1p(d$nCount_RNA)),rank(d$mki))
    recomputed[j] <- cor(qr.resid(qr(cov),xx),qr.resid(qr(cov),yy))
  }
  z <- atanh(recomputed)
  est <- tanh(mean(z)); tt <- t.test(z)
  checks[[i]] <- data.frame(cohort=r$cohort,program=r$program,
      max_patient_r_difference=max(abs(recomputed-ps$r)),summary_r_difference=abs(est-r$r),
      p_difference=abs(tt$p.value-r$p),t_ci_low=tanh(tt$conf.int[1]),t_ci_high=tanh(tt$conf.int[2]))
}
q <- do.call(rbind,checks)
write.csv(q,file.path(root,"qa/independent_R_statistical_check.csv"),row.names=FALSE)
stopifnot(all(q$max_patient_r_difference < 1e-8), all(q$summary_r_difference<1e-8),all(q$p_difference<1e-8))
print(q)
cat("INDEPENDENT R CHECK PASSED\n")
