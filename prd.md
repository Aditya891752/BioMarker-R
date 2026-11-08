# Product Requirements Document (PRD)

## Project Title
**End-to-End Bulk RNA-Seq Analysis & Biomarker Discovery Platform (BioMarker-R)**

---

## 1. Context, Purpose & Delivery Model
BioMarker-R is a computational biology transcriptomics project engineered to demonstrate an end-to-end, fully reproducible analysis workflow from raw expression data to clinical interpretation. Designed as an industry-grade portfolio piece for Bioinformatics and Data Analyst roles, the project adheres to strict statistical rigor, publication standards, and transparent documentation.

Instead of introducing unnecessary cloud microservices or external hosting overhead, the project delivers:
1. **A Hermetically Reproducible Pipeline**: Master R orchestrator backed by `renv` and Docker for 1-command replication.
2. **Interactive HTML Report (GitHub Pages)**: A standalone, publication-grade Quarto / R Markdown report containing interactive data tables (DT), zoomable volcano/PCA plots (plotly), publication heatmaps (ComplexHeatmap), and survival curves.
3. **Publication-Ready Figures**: High-resolution vector (PDF) and raster (PNG, 300 DPI) figures formatted for peer-reviewed journal submission.

The chosen disease model is **Triple-Negative Breast Cancer (TNBC)** contrasted against **adjacent normal breast tissue**.

---

## 2. Biological Questions & Objectives
1. **Differential Expression**: Which transcripts exhibit statistically and biologically significant expression alterations in TNBC compared to matched normal breast tissue?
2. **Pathway & Gene Ontology Enrichment**: Which biological processes (GO BP/MF/CC) and canonical molecular pathways (KEGG / Reactome) are perturbed in TNBC?
3. **Prognostic Biomarker Identification**: Does stratified expression (high vs. low) of top candidate biomarkers correlate significantly with patient overall survival (OS) outcomes?

---

## 3. Data Sources & Cohort Requirements
- **Primary Cohort**: TCGA-BRCA (The Cancer Genome Atlas - Breast Invasive Carcinoma) downloaded programmatically via `TCGAbiolinks`, filtered for TNBC samples (ER negative, PR negative, HER2 negative) and matched solid tissue normal controls.
- **Secondary Validation Cohort**: NCBI GEO accession (e.g., GSE58135 or GSE76275) ingested via `GEOquery`.
- **Target Size**: 100+ patient samples across disease and control groups.
- **Metadata Fields**: Sample Barcode, Clinical Condition (`Tumor` vs. `Normal`), Subtype, Age at Diagnosis, Overall Survival Time (days/months), Vital Status (0 = Censored/Alive, 1 = Event/Deceased), Pathological Stage.

---

## 4. Analytical Workflow & Methodology

### 4.1 Data Ingestion & Identifier Harmonization
- Programmatic acquisition via `TCGAbiolinks` and `GEOquery`.
- Gene identifier translation: Harmonize Ensembl Gene IDs (`ENSG...`) to HGNC approved gene symbols and Entrez Gene IDs via `org.Hs.eg.db`.
- Verify sample matching: Ensure strict 1-to-1 correspondence between count matrix column names and clinical metadata rows.

### 4.2 Quality Control & Exploratory Analysis
- Filter unexpressed and low-count genes (retain genes with count $\ge 10$ in at least $N$ samples, where $N$ equals the smallest group size).
- Apply Variance Stabilizing Transformation (`vst`) or regularized logarithm (`rlog`) via `DESeq2`.
- Principal Component Analysis (PCA) on top 500–1000 variable genes to evaluate cluster separation and rule out batch effects or outliers.
- Sample-to-sample distance matrix and unsupervised hierarchical clustering (Euclidean distance, Ward.D2 linkage).

### 4.3 Differential Expression Analysis (DEA)
- Primary engine: `DESeq2` fitting negative binomial generalized linear models:
  $$\text{design} = \sim \text{condition}$$
- Empirical Bayes shrinkage using `lfcShrink()` with `apeglm` or `normal` prior.
- Secondary cross-check: `edgeR` quasi-likelihood F-test (`glmQLFit`) to verify gene ranking and fold-change consistency.
- Significance Thresholds:
  - Adjusted p-value (Benjamini-Hochberg FDR) $< 0.01$
  - Absolute effect size: $|\log_2 \text{Fold Change}| > 2.0$ (4-fold expression change)

### 4.4 Functional Enrichment Analysis
- Over-representation analysis (ORA) with `clusterProfiler`:
  - Gene Ontology (GO): Biological Process (BP), Molecular Function (MF), Cellular Component (CC).
  - KEGG pathway analysis using Entrez IDs.
- Significance cutoff: $p_{\text{adj}} < 0.05$.
- Gene Set Enrichment Analysis (GSEA) ranking genes by $\text{sign}(\log_2\text{FC}) \times (-\log_{10} p_{\text{adj}})$.

### 4.5 Survival Modeling & Clinical Stratification
- Stratify patients into High vs. Low expression groups using median or maximally selected rank statistics (`surv_cutpoint`).
- Fit Kaplan-Meier curves and test for statistical divergence via Log-Rank test ($p < 0.05$) using `survival` and `survminer`.
- Fit Univariate and Multivariate Cox Proportional Hazards regression (`coxph`) adjusting for patient age and stage.

### 4.6 Biomarker Triangulation
- Final candidate biomarkers must satisfy the triple requirement:
  1. High effect size and statistical significance ($|\log_2\text{FC}| > 2.0$, $p_{\text{adj}} < 0.01$).
  2. Direct participation in enriched oncogenic/mitotic pathways.
  3. Statistically significant correlation with overall patient survival ($p_{\text{log-rank}} < 0.05$).

---

## 5. Deliverables & Presentation
1. **Fully Reproducible R Pipeline**: Executable via a single command (`Rscript run_pipeline.R` or `make all`).
2. **Self-Contained Quarto / R Markdown Interactive Report**:
   - Deployed automatically to **GitHub Pages** via GitHub Actions.
   - Embeds interactive data tables (`DT`) with client-side sorting/filtering.
   - Features interactive plotly scatter/volcano charts alongside static publication figures.
3. **Publication-Quality Figures Folder**:
   - `figures/pca_plot.{pdf,png}`
   - `figures/volcano_plot.{pdf,png}`
   - `figures/deg_heatmap.{pdf,png}`
   - `figures/go_kegg_enrichment.{pdf,png}`
   - `figures/kaplan_meier_survival.{pdf,png}`
4. **Structured Clean Artifacts**:
   - `data/processed/deg_results_all.csv`
   - `data/processed/deg_significant.csv`
   - `data/processed/pathway_enrichment.csv`
   - `data/processed/biomarker_candidates.csv`
5. **Comprehensive Project README**: Documenting scientific background, methods, key findings, interview talking points, and instructions to reproduce.
