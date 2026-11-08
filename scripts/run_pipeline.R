# ==============================================================================
# Script: run_pipeline.R
# Description: Master Pipeline Orchestrator for BioMarker-R
# Project: BioMarker-R (Bulk RNA-Seq & Biomarker Discovery)
# ==============================================================================

cat("\n======================================================================\n")
cat("      BIOMARKER-R: END-TO-END TRANSCRIPTOMICS & BIOMARKER DISCOVERY    \n")
cat("======================================================================\n\n")

start_time <- Sys.time()

# 1. Source utility helpers
source("scripts/00_utils.R")
set_reproducible_seed(42)
paths <- get_project_paths()

# 2. Sequential execution of pipeline modules
pipeline_steps <- c(
  "scripts/01_data_download.R",
  "scripts/02_preprocessing_qc.R",
  "scripts/03_deseq2_diffexp.R",
  "scripts/04_pathway_enrichment.R",
  "scripts/05_survival_analysis.R"
)

for (step_script in pipeline_steps) {
  step_start <- Sys.time()
  log_info(sprintf("Executing: %s", step_script))
  
  if (!file.exists(step_script)) {
    log_error(sprintf("Script not found: %s", step_script))
    stop("Aborting pipeline execution.")
  }
  
  source(step_script)
  step_elapsed <- round(as.numeric(difftime(Sys.time(), step_start, units = "secs")), 2)
  log_success(sprintf("Finished %s in %.2f seconds.\n", basename(step_script), step_elapsed))
}

total_elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)

cat("\n======================================================================\n")
cat(sprintf("  PIPELINE EXECUTION COMPLETE IN %.2f SECONDS\n", total_elapsed))
cat("======================================================================\n")
cat("Artifacts generated:\n")
cat(sprintf("  - Processed Data: %s\n", paths$data_proc))
cat(sprintf("  - Figures:        %s\n", paths$figures))
cat("Next Step: Compile interactive Quarto/Rmd report in reports/\n\n")
