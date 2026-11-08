# Project Tasks & Execution Roadmap (`task.md`)

This roadmap follows the **GSD (Get Stuff Done)** structured execution model. Each phase is broken down into atomic, testable tasks with concrete verification steps and deliverables.

---

## Phase 1: Environment Setup & Data Acquisition
**Goal**: Establish hermetic reproducibility with `renv` and fetch matched TNBC/Normal cohorts from TCGA and GEO.

- [x] **Task 1.1: Initialize Project Directory & R Environment**
  - Create directory structure (`scripts/`, `reports/`, `figures/`, `data/raw/`, `data/processed/`, `docs/`).
  - Initialize R dependency lockfile and folder structure.
  - Configure `.gitignore` to prevent any data cache or secret leaks.
  - *Verification*: Git initialized with clean initial commit; `.gitignore` covers `.env*`, `data/raw/`.

- [ ] **Task 1.2: Programmatic Data Ingestion (TCGA-BRCA)**
  - Write `scripts/01_data_download.R`.
  - Use `TCGAbiolinks` to query GDC for TCGA-BRCA RNA-seq raw count matrices and clinical data.
  - Filter for Triple-Negative Breast Cancer samples (ER-, PR-, HER2- by IHC/FISH) and solid tissue normal samples.
  - *Verification*: Script outputs count matrix of $>100$ samples and clinical metadata table with matched barcodes.

- [ ] **Task 1.3: Secondary GEO Cohort Ingestion (Cross-Validation)**
  - Ingest validation dataset (e.g., GSE58135 / GSE76275) using `GEOquery`.
  - Parse expression series and extract survival/subtype phenotypes.
  - *Verification*: Saved RDS object with synchronized sample IDs.

---

## Phase 2: Preprocessing, Harmonization & Quality Control
**Goal**: Filter unexpressed genes, harmonize Ensembl to HGNC symbols, and detect batch effects or outliers via PCA.

- [ ] **Task 2.1: Gene Identifier Mapping & Filtering**
  - Map Ensembl IDs to HGNC symbols and Entrez IDs using `org.Hs.eg.db`.
  - Resolve duplicates by keeping the highest mean count transcript.
  - Filter low-count genes (count $\ge 10$ in at least $N$ samples, where $N$ = smallest group size).
  - *Verification*: Summary printout showing pre-filtering vs. post-filtering gene counts (~14,000–16,000 retained).

- [ ] **Task 2.2: Variance Stabilization & Sample Distance Matrix**
  - Construct `DESeqDataSet` object.
  - Run `vst(blind = TRUE)` for unsupervised quality assessment.
  - Compute sample-to-sample Euclidean distance matrix and hierarchical clustering.
  - *Verification*: Generated sample distance heatmap showing clear clustering of normal vs. tumor.

- [ ] **Task 2.3: Principal Component Analysis (PCA)**
  - Calculate PCA on the top 500 most variable genes.
  - Plot PC1 vs. PC2 colored by condition (`TNBC` vs. `Normal`) with percentage variance explained.
  - Save publication plot: `figures/fig1_pca_qc.pdf` and `figures/fig1_pca_qc.png`.
  - *Verification*: PC1 separates Tumor from Normal ($> 40\%$ variance explained); no extreme outliers detected.

---

## Phase 3: Differential Expression Analysis (DEA)
**Goal**: Fit negative binomial generalized linear models, apply empirical Bayes shrinkage, and benchmark against edgeR.

- [ ] **Task 3.1: Fit DESeq2 Model**
  - Define model design formula: `design = ~ condition`.
  - Execute `DESeq()` estimation (size factors, dispersions, Wald tests).
  - Extract contrast: `results(dds, contrast = c("condition", "TNBC", "Normal"))`.
  - Apply `lfcShrink()` with `apeglm` or `normal` prior.
  - *Verification*: Result table contains `baseMean`, `log2FoldChange`, `lfcSE`, `pvalue`, `padj`.

- [ ] **Task 3.2: Secondary Benchmark with edgeR**
  - Run parallel DGEList workflow with TMM normalization and Quasi-Likelihood F-test (`glmQLFit`).
  - Calculate Spearman correlation between DESeq2 and edgeR $\log_2\text{FC}$ values.
  - *Verification*: Correlation coefficient $r > 0.95$, confirming statistical robustness.

- [ ] **Task 3.3: Significance Filtering & Visualization**
  - Apply strict cutoffs: $\text{FDR} < 0.01$ and $|\log_2\text{FC}| > 2.0$.
  - Render publication-grade static Volcano Plot using `ggplot2`: `figures/fig2_volcano_plot.pdf`.
  - Generate top 50 DEG expression heatmap using `ComplexHeatmap`: `figures/fig3_deg_heatmap.pdf`.
  - Save clean CSV tables: `data/processed/deg_all_genes.csv` and `data/processed/deg_significant.csv`.
  - *Verification*: Identified ~45–150 highly significant DEGs including known TNBC markers.

---

## Phase 4: Functional & Pathway Enrichment Analysis
**Goal**: Identify dysregulated biological processes (GO) and canonical molecular pathways (KEGG).

- [ ] **Task 4.1: Gene Ontology (GO) Over-Representation Analysis**
  - Run `enrichGO()` via `clusterProfiler` for BP (Biological Process), MF (Molecular Function), and CC (Cellular Component).
  - Use Benjamini-Hochberg adjustment ($p_{\text{adj}} < 0.05$).
  - Simplify redundant terms with `simplify()`.
  - *Verification*: Enriched mitotic cell cycle, cell division, and extracellular matrix organization terms detected.

- [ ] **Task 4.2: KEGG Pathway Enrichment**
  - Convert candidate gene symbols to Entrez IDs.
  - Run `enrichKEGG()` and extract pathways (Cell Cycle, p53 signaling, PI3K-Akt signaling).
  - *Verification*: Enriched pathways match established triple-negative breast oncogenic biology.

- [ ] **Task 4.3: Export Enrichment Visualizations & Tables**
  - Generate publication dot plots and network maps: `figures/fig4_pathway_enrichment.pdf`.
  - Export serialized pathway records to `data/processed/pathway_enrichment_results.csv`.
  - *Verification*: CSV contains term descriptions, gene ratios, adjusted p-values, and associated gene symbols.

---

## Phase 5: Clinical Survival Modeling & Biomarker Shortlist
**Goal**: Stratify patients by candidate gene expression, compute Kaplan-Meier curves, and isolate prognostic biomarkers.

- [ ] **Task 5.1: Clinical Cohort Alignment**
  - Merge patient survival metrics (overall survival days/months, vital status) with normalized expression matrix.
  - Filter out cases with missing follow-up or survival times $< 30$ days.
  - *Verification*: Matched cohort with $> 100$ complete clinical records.

- [ ] **Task 5.2: Kaplan-Meier Estimation & Log-Rank Testing**
  - Stratify patients into High vs. Low expression using median or `surv_cutpoint()` (optimal cutoff).
  - Fit `survfit(Surv(time, status) ~ strata)` using `survival` and `survminer`.
  - Perform Log-Rank test; calculate hazard ratios via Cox Proportional Hazards (`coxph`).
  - Render publication plots: `figures/fig5_kaplan_meier_biomarkers.pdf`.
  - *Verification*: Generated KM step plots with risk tables and log-rank p-values.

- [ ] **Task 5.3: Biomarker Triangulation & Shortlist Export**
  - Intersect criteria:
    1. Highly significant DEG ($|\log_2\text{FC}| > 2.0$, $\text{FDR} < 0.01$).
    2. Over-represented in key oncogenic pathways.
    3. Significant overall survival association ($p < 0.05$).
  - Export candidate profiles to `data/processed/candidate_biomarkers.csv`.
  - *Verification*: Shortlist of 5–15 high-confidence candidate prognostic biomarkers.

---

## Phase 6: Interactive Quarto / R Markdown Report Authoring
**Goal**: Build a self-contained, publication-grade interactive report.

- [ ] **Task 6.1: Report Structure & Narrative Drafting**
  - Author `reports/BioMarker_Analysis_Report.qmd` with structured sections: Introduction & Biological Context, Methods, Quality Control, Differential Expression, Functional Enrichment, Survival Analysis, Discussion & Limitations.
  - Adhere strictly to `design.md` and `rules.md` (no emojis in headings, no purple-blue gradients, high-contrast typography).
  - *Verification*: Clean document compiling locally with zero rendering errors.

- [ ] **Task 6.2: Embed Interactive DataTables & Plotly Graphics**
  - Embed interactive `DT::datatable()` for `deg_significant.csv` and `pathway_enrichment_results.csv` with instant search and download buttons.
  - Embed interactive `plotly` volcano plot with gene hover tooltips.
  - *Verification*: Table searches and plot tooltips operate smoothly without external network calls.

- [ ] **Task 6.3: Render Standalone Static Site**
  - Render report to `docs/index.html`.
  - *Verification*: Opening `docs/index.html` in a web browser works completely offline with all assets self-contained.

---

## Phase 7: Containerization & GitHub Pages CI/CD
**Goal**: Automate testing and deployment to GitHub Pages.

- [ ] **Task 7.1: Containerization with Docker**
  - Write `Dockerfile` based on `rocker/r-ver:4.3.2`.
  - Create `Makefile` with targets: `make setup`, `make pipeline`, `make report`, `make all`.
  - *Verification*: `docker build -t biomarker-r .` builds successfully and reproduces outputs.

- [ ] **Task 7.2: GitHub Actions Workflow for GitHub Pages**
  - Author `.github/workflows/deploy-report.yml`.
  - Configure automatic compilation on push to `main` and deployment to GitHub Pages (`docs/` folder or `gh-pages` branch).
  - *Verification*: GitHub Action passes and publishes live report to `https://<username>.github.io/BioMarker-R/`.

---

## Phase 8: Final Review & Portfolio Packaging
**Goal**: Prepare project for technical interviews and public demonstration.

- [ ] **Task 8.1: Publication Figures Review**
  - Audit all figures in `figures/` for font legibility, label clarity, 300 DPI resolution, and color contrast.
  - *Verification*: Figures meet Nature / BMC Bioinformatics submission standards.

- [ ] **Task 8.2: Comprehensive README & Technical Documentation**
  - Write executive `README.md` detailing: Project Overview, Research Questions, Pipeline Architecture, Key Biological Discoveries, How to Reproduce (Local / Docker), Interview Talking Points.
  - *Verification*: External reviewer can clone repo and run pipeline with 1 command.
