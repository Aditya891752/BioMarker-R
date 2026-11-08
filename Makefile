# ==============================================================================
# Makefile for BioMarker-R Computational Transcriptomics Pipeline
# ==============================================================================

.PHONY: all pipeline report clean help

all: pipeline report
	@echo "======================================================================"
	@echo " BioMarker-R pipeline and interactive report built successfully!"
	@echo " View the report at docs/index.html"
	@echo "======================================================================"

pipeline:
	@echo "==> Running core transcriptomics pipeline (scripts 01 to 05)..."
	Rscript scripts/run_pipeline.R

report:
	@echo "==> Compiling interactive Quarto / R Markdown report to docs/index.html..."
	quarto render reports/BioMarker_Analysis_Report.qmd --output-dir ../docs --output index.html || \
	Rscript -e "rmarkdown::render('reports/BioMarker_Analysis_Report.qmd', output_dir = 'docs', output_file = 'index.html')"

clean:
	@echo "==> Cleaning temporary R build caches and session logs..."
	rm -rf .Rhistory .RData reports/.quarto reports/*.html
	@echo "Clean completed."

help:
	@echo "Available make commands:"
	@echo "  make all       - Run complete pipeline and render interactive report"
	@echo "  make pipeline  - Run R scripts 01 to 05 (data -> PCA -> DEA -> ORA -> Survival)"
	@echo "  make report    - Render Quarto report to docs/index.html"
	@echo "  make clean     - Remove temporary files"
