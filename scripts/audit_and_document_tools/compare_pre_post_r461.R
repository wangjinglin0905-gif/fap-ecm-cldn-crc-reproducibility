#!/usr/bin/env Rscript

# Compare all pre-R-4.6.1 and post-R-4.6.1 result files.
# This script performs only read-only reconciliation and writes a machine-readable
# audit table. It intentionally does not modify any scientific result.

old_root <- normalizePath(
  "work/revision_round3/results/pre_r461/results",
  winslash = "/",
  mustWork = TRUE
)
new_root <- normalizePath(
  "work/reproducibility/results",
  winslash = "/",
  mustWork = TRUE
)
out_dir <- "work/revision_round3/reports"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(requireNamespace("data.table", quietly = TRUE))
stopifnot(requireNamespace("digest", quietly = TRUE))

relative_files <- function(root) {
  full <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  full <- full[!file.info(full)$isdir]
  substring(normalizePath(full, winslash = "/"), nchar(root) + 2L)
}

old_files <- relative_files(old_root)
new_files <- relative_files(new_root)
all_files <- sort(union(old_files, new_files))

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

is_table <- function(path) {
  grepl("\\.(csv|tsv|txt)(\\.gz)?$", path, ignore.case = TRUE) &&
    !grepl("(warnings|sessionInfo)\\.txt$", path, ignore.case = TRUE)
}

read_table_safe <- function(path) {
  tryCatch(
    data.table::fread(path, data.table = FALSE, showProgress = FALSE),
    error = function(e) structure(list(error = conditionMessage(e)), class = "read_error")
  )
}

compare_tables <- function(old_path, new_path) {
  old <- read_table_safe(old_path)
  new <- read_table_safe(new_path)
  if (inherits(old, "read_error") || inherits(new, "read_error")) {
    return(list(
      table_status = "read_error",
      old_nrow = NA_integer_, old_ncol = NA_integer_,
      new_nrow = NA_integer_, new_ncol = NA_integer_,
      numeric_diff_cells = NA_integer_, character_diff_cells = NA_integer_,
      max_abs_numeric_diff = NA_real_, differing_columns = NA_character_,
      read_error = paste(
        if (inherits(old, "read_error")) old$error else "",
        if (inherits(new, "read_error")) new$error else "",
        sep = " | "
      )
    ))
  }

  ans <- list(
    table_status = "comparable",
    old_nrow = nrow(old), old_ncol = ncol(old),
    new_nrow = nrow(new), new_ncol = ncol(new),
    numeric_diff_cells = 0L, character_diff_cells = 0L,
    max_abs_numeric_diff = 0, differing_columns = "", read_error = ""
  )

  if (!identical(dim(old), dim(new)) || !identical(names(old), names(new))) {
    ans$table_status <- "schema_or_dimension_changed"
    ans$differing_columns <- paste(
      setdiff(union(names(old), names(new)), intersect(names(old), names(new))),
      collapse = ";"
    )
    return(ans)
  }

  differing <- character()
  for (nm in names(old)) {
    x <- old[[nm]]
    y <- new[[nm]]
    if (is.numeric(x) && is.numeric(y)) {
      both_na <- is.na(x) & is.na(y)
      one_na <- xor(is.na(x), is.na(y))
      delta <- abs(x - y)
      changed <- one_na | (!both_na & !is.na(delta) & delta > 1e-12)
      n_changed <- sum(changed)
      ans$numeric_diff_cells <- ans$numeric_diff_cells + n_changed
      if (any(is.finite(delta))) {
        ans$max_abs_numeric_diff <- max(ans$max_abs_numeric_diff, max(delta, na.rm = TRUE))
      }
    } else {
      xs <- as.character(x)
      ys <- as.character(y)
      both_na <- is.na(xs) & is.na(ys)
      changed <- !both_na & (xor(is.na(xs), is.na(ys)) | (!is.na(xs) & !is.na(ys) & xs != ys))
      n_changed <- sum(changed)
      ans$character_diff_cells <- ans$character_diff_cells + n_changed
    }
    if (n_changed > 0L) differing <- c(differing, nm)
  }
  ans$differing_columns <- paste(differing, collapse = ";")
  if (ans$numeric_diff_cells == 0L && ans$character_diff_cells == 0L) {
    ans$table_status <- "identical_values"
  } else {
    ans$table_status <- "values_changed"
  }
  ans
}

rows <- vector("list", length(all_files))
for (i in seq_along(all_files)) {
  rel <- all_files[[i]]
  old_path <- file.path(old_root, rel)
  new_path <- file.path(new_root, rel)
  old_exists <- file.exists(old_path)
  new_exists <- file.exists(new_path)
  row <- list(
    file = gsub("\\\\", "/", rel),
    old_exists = old_exists,
    new_exists = new_exists,
    old_bytes = if (old_exists) file.info(old_path)$size else NA_real_,
    new_bytes = if (new_exists) file.info(new_path)$size else NA_real_,
    old_sha256 = if (old_exists) sha256_file(old_path) else NA_character_,
    new_sha256 = if (new_exists) sha256_file(new_path) else NA_character_,
    byte_identical = if (old_exists && new_exists) {
      identical(sha256_file(old_path), sha256_file(new_path))
    } else NA
  )
  tab <- list(
    table_status = if (is_table(rel)) "not_comparable_missing_file" else "not_tabular",
    old_nrow = NA_integer_, old_ncol = NA_integer_,
    new_nrow = NA_integer_, new_ncol = NA_integer_,
    numeric_diff_cells = NA_integer_, character_diff_cells = NA_integer_,
    max_abs_numeric_diff = NA_real_, differing_columns = NA_character_,
    read_error = NA_character_
  )
  if (old_exists && new_exists && is_table(rel)) {
    tab <- compare_tables(old_path, new_path)
  }
  rows[[i]] <- c(row, tab)
}

comparison <- data.table::rbindlist(rows, fill = TRUE)
data.table::fwrite(
  comparison,
  file.path(out_dir, "pre_post_r461_result_comparison.csv"),
  na = ""
)

summary <- data.frame(
  metric = c(
    "files_total", "files_byte_identical", "files_byte_changed",
    "tables_identical_values", "tables_values_changed",
    "tables_schema_or_dimension_changed", "tables_read_error"
  ),
  value = c(
    nrow(comparison),
    sum(comparison$byte_identical %in% TRUE),
    sum(comparison$byte_identical %in% FALSE),
    sum(comparison$table_status == "identical_values", na.rm = TRUE),
    sum(comparison$table_status == "values_changed", na.rm = TRUE),
    sum(comparison$table_status == "schema_or_dimension_changed", na.rm = TRUE),
    sum(comparison$table_status == "read_error", na.rm = TRUE)
  )
)
data.table::fwrite(summary, file.path(out_dir, "pre_post_r461_result_summary.csv"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "pre_post_compare_sessionInfo.txt"))
print(summary, row.names = FALSE)
