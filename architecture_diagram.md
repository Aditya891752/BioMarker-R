# Architecture & Deployment Diagram (`architecture_diagram.md`)

This document details the end-to-end execution workflow and provides step-by-step instructions for running the pipeline locally, compiling the interactive report, and publishing it freely to **GitHub Pages** with zero cloud hosting dependencies.

---

## 1. End-to-End Pipeline & Publishing Workflow

```mermaid
flowchart TD
    subgraph Data_Sources["1. Public Biological Repositories"]
        A1["TCGA GDC Portal (TCGA-BRCA)"]
        A2["NCBI GEO (Validation Cohort)"]
    end

    subgraph R_Pipeline["2. Computational Biology Pipeline (R / Bioconductor)"]
        B1["scripts/01_data_download.R<br/>(TCGAbiolinks / GEOquery)"]
        B2["scripts/02_preprocessing_qc.R<br/>(org.Hs.eg.db, VST, PCA, Clustering)"]
        B3["scripts/03_deseq2_diffexp.R<br/>(DESeq2 NB-GLM, apeglm, edgeR check)"]
        B4["scripts/04_pathway_enrichment.R<br/>(clusterProfiler: GO BP/MF/CC & KEGG)"]
        B5["scripts/05_survival_analysis.R<br/>(survival, survminer: KM & Cox PH)"]

        B1 --> B2 --> B3 --> B4 --> B5
    end

    subgraph Artifacts["3. Processed Data & Publication Outputs"]
        C1["data/processed/<br/>(deg_significant.csv, biomarker_candidates.csv)"]
        C2["figures/<br/>(fig1_pca, fig2_volcano, fig3_heatmap, fig5_km.pdf/png)"]
        B5 --> C1
        B5 --> C2
    end

    subgraph Report_Engine["4. Quarto / R Markdown Publishing"]
        D1["reports/BioMarker_Analysis_Report.qmd"]
        C1 --> D1
        C2 --> D1
        D1 -->|"quarto render / rmarkdown::render()"| D2["docs/index.html<br/>(Interactive DT Tables, Plotly, Publication Figures)"]
    end

    subgraph GitHub_Hosting["5. Version Control & GitHub Pages"]
        E1["GitHub Repository (BioMarker-R)"]
        E2[".github/workflows/deploy-report.yml"]
        E3["GitHub Pages Branch (gh-pages / docs)"]
        U["Public Viewers / Recruiters / Hiring Managers"]

        D2 --> E1
        E1 --> E2 --> E3 --> U
    end

    A1 --> B1
    A2 --> B1
```

---

## 2. Step-by-Step Instructions: Local Run to GitHub Pages Deployment

### Step A: Initialize the Reproducible Environment
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/BioMarker-R.git
   cd BioMarker-R
   ```
2. Restore package dependencies in R:
   ```R
   # Inside an R terminal or RStudio
   renv::restore()
   ```

---

### Step B: Run the Complete Bioinformatics Pipeline
Execute the master pipeline script, which runs all sequential modules from data download to biomarker triangulation:
```bash
Rscript scripts/run_pipeline.R
```
Or use the Makefile:
```bash
make all
```
**Outputs generated**:
- Clean tables saved in `data/processed/` (`deg_significant.csv`, `pathway_enrichment_results.csv`, `candidate_biomarkers.csv`).
- High-resolution publication plots saved in `figures/` (PDF and 300 DPI PNG).

---

### Step C: Compile the Interactive Report
Render the comprehensive Quarto / R Markdown report into a standalone, searchable HTML page located in `docs/`:
```bash
quarto render reports/BioMarker_Analysis_Report.qmd --output-dir ../docs --output index.html
```
Or via R console:
```R
rmarkdown::render("reports/BioMarker_Analysis_Report.Rmd", output_dir = "docs", output_file = "index.html")
```
*Tip*: Open `docs/index.html` in your local browser to inspect the interactive tables and figures before publishing.

---

### Step D: Publish Automatically on GitHub Pages

1. **Commit and Push to GitHub**:
   ```bash
   git add .
   git commit -m "feat: complete reproducible transcriptomics pipeline and report"
   git push origin main
   ```
2. **Enable GitHub Pages**:
   - Go to your repository on GitHub: `https://github.com/<your-username>/BioMarker-R`.
   - Click **Settings** $\rightarrow$ **Pages** (in the left sidebar).
   - Under **Build and deployment** $\rightarrow$ **Source**:
     - Choose **Deploy from a branch**.
     - **Branch**: Select `main`, folder: `/docs`.
     - Click **Save**.
     *(Alternatively, use the included GitHub Actions workflow `.github/workflows/deploy-report.yml` to build and deploy automatically).*
3. **Verify the Live Deployment**:
   - Within 1–2 minutes, your report will be live at:
     `https://<your-username>.github.io/BioMarker-R/`
   - Test that interactive tables (search, sort, pagination) and hover tooltips work seamlessly without any hosting costs or backend maintenance.
