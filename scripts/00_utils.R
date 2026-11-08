# ==============================================================================
# Script: 00_utils.R
# Description: Shared utility functions, logging, paths, and clinical lab palette
# Project: BioMarker-R (Bulk RNA-Seq & Biomarker Discovery)
# ==============================================================================

# Suppress noisy package startup messages globally when sourced
suppressPackageStartupMessages({
  # Standard libraries checked when available
})

# ------------------------------------------------------------------------------
# 1. Deterministic Reproducibility
# ------------------------------------------------------------------------------
set_reproducible_seed <- function(seed = 42) {
  set.seed(seed)
  if (exists("RNGkind")) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  }
}

# ------------------------------------------------------------------------------
# 2. Precision Logging (Console & File)
# ------------------------------------------------------------------------------
log_info <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[INFO]  [%s] %s\n", timestamp, msg))
}

log_success <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[SUCCESS] [%s] %s\n", timestamp, msg))
}

log_warn <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[WARN]  [%s] %s\n", timestamp, msg))
}

log_error <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[ERROR] [%s] %s\n", timestamp, msg))
}

# ------------------------------------------------------------------------------
# 3. Path Management & Directory Assertions
# ------------------------------------------------------------------------------
get_project_paths <- function() {
  # Root is assumed to be project root
  root_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  
  paths <- list(
    root      = root_dir,
    scripts   = file.path(root_dir, "scripts"),
    data_raw  = file.path(root_dir, "data", "raw"),
    data_proc = file.path(root_dir, "data", "processed"),
    figures   = file.path(root_dir, "figures"),
    reports   = file.path(root_dir, "reports"),
    docs      = file.path(root_dir, "docs")
  )
  
  # Ensure target directories exist
  for (p in paths) {
    if (!dir.exists(p)) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(paths)
}

# ------------------------------------------------------------------------------
# 4. Clinical Laboratory Chromatic Standards (Conforms strictly to design.md)
#    No purple-to-blue gradients; high-contrast editorial lab palette.
# ------------------------------------------------------------------------------
LAB_PALETTE <- list(
  upregulated   = "#B91C1C", # Deep Crimson Red (high effect size / overexpression)
  downregulated = "#1D4ED8", # Cobalt Navy (downregulation)
  neutral       = "#94A3B8", # Cool Slate Gray (non-significant genes)
  normal_cohort = "#0D9488", # Teal Diagnostic (Solid Tissue Normal Control)
  tumor_cohort  = "#BE185D", # Deep Rose Plum (TNBC Cohort)
  biomarker     = "#047857", # Emerald (Triangulated prognostic candidate)
  border        = "#E2E8F0", # Clean 1px solid border
  bg_panel      = "#F8F9FA", # Lab off-white container
  text_primary  = "#0F172A"  # High-contrast Slate 900
)

# ------------------------------------------------------------------------------
# 5. Cohort & Count Matrix Verification Checkpoints
# ------------------------------------------------------------------------------
assert_sample_parity <- function(count_matrix, metadata, sample_col = "sample_id") {
  count_samples <- colnames(count_matrix)
  meta_samples  <- metadata[[sample_col]]
  
  if (is.null(count_samples) || is.null(meta_samples)) {
    stop("Sample names missing in count matrix columns or metadata table.")
  }
  
  # Check identical ordering
  if (identical(count_samples, meta_samples)) {
    log_success(sprintf("Sample parity verified: %d matched samples in exact order.", length(count_samples)))
    return(TRUE)
  }
  
  # Check matching set
  common_samples <- intersect(count_samples, meta_samples)
  if (length(common_samples) == 0) {
    stop("Zero overlapping samples between count matrix and metadata!")
  }
  
  log_warn(sprintf("Counts (%d) and Metadata (%d) samples misaligned. Synchronizing on %d shared samples.",
                   length(count_samples), length(meta_samples), length(common_samples)))
  return(FALSE)
}
