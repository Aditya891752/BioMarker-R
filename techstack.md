# Technology Stack (`techstack.md`)

## 1. Stack Matrix Overview

| Domain / Layer | Technology | Version / Standard | Role & Justification |
| :--- | :--- | :--- | :--- |
| **Statistical Computing** | **R** | `>= 4.3.0` | Gold standard programming language for statistical genetics and transcriptomics. |
| **Package Hermeticity** | **renv** | `>= 1.0.0` | Creates exact dependency lockfiles (`renv.lock`) to eliminate dependency drift across operating systems. |
| **Data Ingestion** | **TCGAbiolinks** / **GEOquery** | Bioconductor | Direct programmatic API extraction of GDC Harmonized BRCA datasets and NCBI GEO accessions. |
| **Differential Expression** | **DESeq2** | Bioconductor | Empirical Bayes dispersion shrinkage and negative binomial generalized linear modeling. |
| **Secondary Verification** | **edgeR** | Bioconductor | Quasi-likelihood F-tests for cross-validating DEG stability, effect sizes, and ranking. |
| **Enrichment Analysis** | **clusterProfiler** | Bioconductor | Over-representation analysis (ORA) and GSEA for Gene Ontology (BP/MF/CC) and KEGG pathways. |
| **Genome Annotation** | **org.Hs.eg.db** | Bioconductor | Authoritative mapping between Ensembl Gene IDs (`ENSG...`), HGNC symbols, and Entrez IDs. |
| **Survival Modeling** | **survival** / **survminer** | CRAN | Kaplan-Meier estimation, Log-Rank statistical testing, and Cox Proportional Hazards regression. |
| **Publication Graphics** | **ComplexHeatmap** / **ggplot2** | Bioconductor / CRAN | Generation of vector (PDF) and 300 DPI (PNG) publication-grade heatmaps, volcano, and PCA plots. |
| **Interactive Visuals** | **plotly** | CRAN | Interactive zoom, pan, and hover inspection of volcano plots and PCA embeddings in the report. |
| **Interactive Data Tables**| **DT (DataTables)** | CRAN | Searchable, paginated, client-side interactive tables for browsing thousands of DEGs and pathways. |
| **Report Generation** | **Quarto / R Markdown** | `>= 1.4` | Literate programming environment compiling code, narrative, math formulas, and plots into standalone HTML. |
| **Containerization** | **Docker** (`rocker/r-ver`) | `>= 4.3.2` | Containerizes the entire R runtime, system libraries, and Bioconductor dependencies for platform independence. |
| **Continuous Delivery** | **GitHub Actions** | Latest | Automatically runs linting, tests pipeline execution, and publishes compiled HTML to GitHub Pages. |
| **Report Hosting** | **GitHub Pages** | Free / Static Edge | Serves the interactive report with zero maintenance, zero server cold starts, and zero hosting costs. |

---

## 2. Bioinformatics & Statistical Stack Justification

### R & Bioconductor over Python
- **Industry & Academic Authority**: While Python has gained ground in single-cell analysis (`Scanpy`), bulk RNA-seq differential expression and survival modeling are overwhelmingly benchmarked against `DESeq2`, `edgeR`, and `clusterProfiler`. Utilizing native Bioconductor packages demonstrates industry-standard knowledge directly applicable to clinical trials and biotech research roles.
- **Empirical Bayes Shrinkage**: `DESeq2`'s `lfcShrink()` with `apeglm` prevents low-count genes with high variance from dominating top DEG lists, producing biologically meaningful effect sizes.

---

## 3. Reporting & Publishing Architecture Justification

### Quarto / R Markdown over Heavy Web Frameworks
- **Direct Scientific Storytelling**: Allows code, biological interpretation, mathematical notation, and visualizations to coexist in a single reproducible narrative document.
- **Interactive DT Tables**: Instead of spinning up a separate database and API server to browse 15,000 genes, `DT` bundles searchable, sortable JavaScript tables directly into the static HTML file, delivering sub-millisecond client-side filtering with zero server dependencies.
- **Permanent Availability**: Hosting on GitHub Pages guarantees that hiring managers and recruiters can access your live interactive analysis without encountering sleeping cloud instances or expired hosting trials.

---

## 4. Reproducibility Guarantee
- **Two Tiers of Reproducibility**:
  1. **Native R**: Replicable with `renv::restore()`.
  2. **Containerized**: Fully encapsulated Docker container executing via `docker run -v $(pwd):/workspace biomarker-r make all`.
