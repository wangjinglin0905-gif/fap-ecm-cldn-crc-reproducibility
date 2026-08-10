# =============================================================================
# Master orchestrator: Phase 1 Analysis Pipeline
# Run: Rscript master_run.R
# =============================================================================
SCRIPT_DIR <- "./scripts"
PYTHON_BIN <- "python"
R_BIN      <- "Rscript"

message("============================================")
message("Phase 1: FAP-ECM-SDC4/CD44 in Early CRC")
message("Analysis Pipeline")
message("============================================")

# --- A1: TCGA bulk (R) ---
message("\n[1/5] Running A1: TCGA-COAD Staging-Stratified Bulk Analysis...")
system2(R_BIN, file.path(SCRIPT_DIR, "01_A1_TCGA_bulk.R"))

# --- A2: Single-cell (Python) ---
message("\n[2/5] Running A2: Single-Cell FAP+ CAF Analysis...")
system2(PYTHON_BIN, file.path(SCRIPT_DIR, "02_A2_single_cell.py"))

# --- A2b: CellChat (R, conditional on A2 output) ---
cellchat_r <- file.path(SCRIPT_DIR, "02b_A2_CellChat.R")
if (file.exists(cellchat_r)) {
  message("\n[2b/5] Running A2b: CellChat inference...")
  system2(R_BIN, cellchat_r)
}

# --- A3: Spatial (Python, conditional on data) ---
message("\n[3/5] Running A3: Spatial proximity analysis...")
system2(PYTHON_BIN, file.path(SCRIPT_DIR, "03_A3_spatial.py"))

# --- A4: Claudin individual (R) ---
message("\n[4/5] Running A4: Claudin family individual analysis...")
system2(R_BIN, file.path(SCRIPT_DIR, "04_A4_claudin_individual.R"))

# --- A5: Pseudotime (R) ---
message("\n[5/5] Running A5: FAP+ CAF pseudotime trajectory...")
system2(R_BIN, file.path(SCRIPT_DIR, "05_A5_pseudotime.R"))

message("\n============================================")
message("Phase 1 pipeline complete!")
message("Output: ./output/")
message("============================================")
