# =============================================================================
# Phase 1 生信分析：FAP-ECM-SDC4/CD44 轴在早期 CRC 间质转变中的作用
# 配置文件 & 全局参数
# =============================================================================

# --- Paths ---
PROJ_ROOT  <- "."
OUT_DIR    <- file.path(PROJ_ROOT, "output")
FIG_DIR    <- file.path(OUT_DIR, "figures")
TAB_DIR    <- file.path(OUT_DIR, "tables")
DATA_DIR   <- file.path(PROJ_ROOT, "data")

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# --- R binary (R 4.6.1, use x64 version to avoid segfault) ---
R_BIN <- "Rscript"

# --- Global thresholds ---
FDR_THRESHOLD    <- 0.05
LOGFC_THRESHOLD  <- 0.5
COR_PVAL         <- 0.05

# --- Gene sets for scoring ---
FAP_CAF_GENES <- c("FAP", "COL1A1", "COL1A2", "FN1", "ACTA2")
ECM_SDC4_CD44_GENES <- c("COL1A1", "COL1A2", "COL3A1", "FN1", "SDC4", "CD44")
PROLIFERATION_CONFOUNDER <- "MKI67"  # for correcting proliferation confounding

# --- Claudin gene list (≥10 members) ---
CLAUDIN_GENES <- c("CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN7",
                   "CLDN8", "CLDN12", "CLDN15", "CLDN18", "CLDN23")

# --- Staging groups (AJCC 8th) ---
STAGE_EARLY <- "I"
STAGE_LATE  <- c("II", "III", "IV")

# --- Required packages ---
REQUIRED_PKGS <- c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2", "edgeR",
  "GSVA", "GSEABase", "survival", "survminer", "rstatix",
  "ggplot2", "ggpubr", "pheatmap", "ComplexHeatmap",
  "dplyr", "tidyr", "tibble", "stringr", "corrplot"
)

for (pkg in REQUIRED_PKGS) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

message("[CONFIG] All paths and packages loaded. Ready for analysis.")
