# System Architecture (`architecture.md`)

## 1. System Overview
BioMarker-R follows a modular, reproducible computational genomics pipeline architecture. The system executes deterministically from raw public data to validated publication figures and an interactive HTML report hosted on GitHub Pages.

```
+-----------------------------------------------------------------------------------+
|                         DATA ACQUISITION & HARMONIZATION                          |
|                                                                                   |
|  [TCGA GDC Portal via TCGAbiolinks]        [NCBI GEO via GEOquery]                |
|                     |                                     |                       |
|                     +------------------+------------------+                       |
|                                        v                                          |
|                           [scripts/01_data_download.R]                            |
|                                        |                                          |
|                                        v                                          |
|                          data/raw/ (Ignored by Git)                               |
+----------------------------------------|------------------------------------------+
                                         v
+-----------------------------------------------------------------------------------+
|                           CORE BIOINFORMATICS PIPELINE                            |
|                                                                                   |
|  [02_preprocessing_qc.R]       - Filter low counts, map Ensembl->Symbol/Entrez,   |
|                                  perform VST, run PCA & sample distance matrix.   |
|            |                                                                      |
|            v                                                                      |
|  [03_deseq2_diffexp.R]         - Negative binomial GLM fitting, apeglm shrinkage, |
|                                  significance filtering, edgeR cross-check.       |
|            |                                                                      |
|            v                                                                      |
|  [04_pathway_enrichment.R]     - clusterProfiler GO (BP/MF/CC) and KEGG ORA/GSEA. |
|            |                                                                      |
|            v                                                                      |
|  [05_survival_analysis.R]      - Kaplan-Meier curves, log-rank tests, Cox PH,     |
|                                  triangulating candidate prognostic biomarkers.   |
+----------------------------------------|------------------------------------------+
                                         v
+-----------------------------------------------------------------------------------+
|                        REPORT GENERATION & PUBLICATION                            |
|                                                                                   |
|  [reports/BioMarker_Analysis_Report.qmd / .Rmd]                                   |
|   ├── Compiles narrative, methods, interactive DT tables, and Plotly graphics    |
|   ├── Generates static publication vector graphics (PDF/PNG in figures/)         |
|   └── Renders standalone docs/index.html                                          |
+----------------------------------------|------------------------------------------+
                                         v
+-----------------------------------------------------------------------------------+
|                          CONTINUOUS INTEGRATION & HOSTING                         |
|                                                                                   |
|  GitHub Actions Workflow (.github/workflows/deploy-report.yml)                    |
|   ├── R environment setup via renv cache                                          |
|   ├── Pipeline execution and report compilation                                  |
|   └── Deploys docs/ to GitHub Pages (https://<username>.github.io/BioMarker-R)     |
+-----------------------------------------------------------------------------------+
```

---

## 2. Directory Layout

```
BioMarker-R/
├── .github/
│   └── workflows/
│       ├── test-pipeline.yml         # CI workflow testing R code and scripts
│       └── deploy-report.yml         # Compiles Quarto report & deploys to GH Pages
├── .gitignore                        # Strict shielding of raw data, caches, secrets
├── Dockerfile                        # Hermetic container environment for pipeline
├── Makefile                          # 1-command execution: `make all`, `make report`
├── README.md                         # Comprehensive portfolio & methodology readme
├── prd.md                            # Product Requirements Document
├── architecture.md                   # System Architecture
├── architecture_diagram.md           # Pipeline & GitHub Pages Deployment Guide
├── techstack.md                      # Technology Stack
├── aistack.md                        # AI Assistance & Guardrails
├── rules.md                          # Visual & Security Prohibitions
├── design.md                         # Publication & Report Visual Design
├── task.md                           # GSD Phase Execution Checklist
├── memory.md                         # Project Context & Progress Log
│
├── renv.lock                         # Exact package version lockfile
├── renv/                             # renv infrastructure scripts
│
├── scripts/                          # Modular analytical R scripts
│   ├── 00_utils.R                    # Helper logging & directory checks
│   ├── 01_data_download.R            # Ingests TCGA & GEO cohorts
│   ├── 02_preprocessing_qc.R         # Filtering, VST, PCA, sample clustering
│   ├── 03_deseq2_diffexp.R           # DESeq2 GLM & edgeR benchmark
│   ├── 04_pathway_enrichment.R       # clusterProfiler GO/KEGG analysis
│   ├── 05_survival_analysis.R        # KM survival, Cox PH & biomarker ranking
│   └── run_pipeline.R                # Master controller executing scripts 01 to 05
│
├── reports/
│   ├── BioMarker_Analysis_Report.qmd # Quarto document (or .Rmd)
│   ├── _quarto.yml                   # Quarto styling & configuration
│   └── references.bib                # Literature citations
│
├── figures/                          # 300 DPI publication plots (PDF + PNG)
│   ├── fig1_pca_qc.pdf
│   ├── fig2_volcano_plot.pdf
│   ├── fig3_deg_heatmap.pdf
│   ├── fig4_pathway_enrichment.pdf
│   └── fig5_kaplan_meier_biomarkers.pdf
│
├── data/
│   ├── raw/                          # Raw count matrices (ignored by Git)
│   └── processed/                    # Clean analytical tables (CSV & RDS)
│       ├── metadata_harmonized.csv
│       ├── deg_all_genes.csv
│       ├── deg_significant.csv
│       ├── pathway_enrichment_results.csv
│       └── candidate_biomarkers.csv
│
└── docs/                             # Built HTML site for GitHub Pages (index.html)
```

---

## 3. Data Processing & Artifact Flow

1. **Raw Acquisition**: `scripts/01_data_download.R` queries GDC via `TCGAbiolinks` and caches the expression matrix and clinical metadata in `data/raw/`.
2. **Quality Control & VST**: `scripts/02_preprocessing_qc.R` removes low-count genes, applies variance stabilizing transformation (`vst()`), and outputs `figures/fig1_pca_qc.pdf`.
3. **Statistical Modeling**: `scripts/03_deseq2_diffexp.R` constructs `DESeqDataSet`, performs dispersion estimation and Wald testing, applies `lfcShrink()`, cross-validates against `edgeR`, and outputs `figures/fig2_volcano_plot.pdf` and `figures/fig3_deg_heatmap.pdf`.
4. **Pathway Enrichment**: `scripts/04_pathway_enrichment.R` translates genes to Entrez IDs, runs GO and KEGG tests, simplifies terms, and generates `figures/fig4_pathway_enrichment.pdf`.
5. **Survival Analysis**: `scripts/05_survival_analysis.R` links patient survival metrics, stratifies cohorts by gene expression, computes log-rank tests and hazard ratios, isolates final biomarkers, and outputs `figures/fig5_kaplan_meier_biomarkers.pdf`.
6. **Report Compilation**: `reports/BioMarker_Analysis_Report.qmd` imports the processed tables and figures, embedding interactive searchable data tables and tooltips, outputting `docs/index.html`.

---

## 4. Reproducibility & Environment Isolation
- **R Environment**: Controlled by `renv`. Running `renv::restore()` on any machine reproduces the exact library environment.
- **Container Environment**: The repository includes a `Dockerfile` based on `rocker/r-ver:4.3.2` with all Bioconductor and CRAN packages pre-installed.
- **GitHub Actions**: Automated pipeline on every push to `main` verifies script execution, builds the Quarto report, and publishes the static HTML artifact directly to the `gh-pages` branch.
