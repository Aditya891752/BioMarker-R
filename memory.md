# Project Memory & Context (`memory.md`)

## 1. Project Identity & Purpose
- **Project Name**: BioMarker-R (End-to-End Bulk RNA-Seq Analysis & Biomarker Discovery)
- **Author Persona**: Undergraduate B.Tech Biotechnology Engineering student.
- **Strategic Goal**: Serve as an industry-grade computational biology portfolio and resume piece for Data Analyst and Bioinformatics roles, demonstrating mastery of statistical genomics, reproducible code, and scientific communication.
- **Disease Focus**: **Triple-Negative Breast Cancer (TNBC)** vs. matched solid tissue normal controls.

---

## 2. Key Biological Context & Prior Knowledge
- **Why TNBC?**:
  - TNBC lacks Estrogen Receptor (ER), Progesterone Receptor (PR), and HER2 expression.
  - Lacks targeted endocrine therapies (like Tamoxifen) or HER2-targeted therapies (like Trastuzumab/Herceptin).
  - High risk of recurrence, aggressive clinical progression, and significant heterogeneity.
  - Finding prognostic and predictive biomarkers is an urgent, active area of oncological translational research.
- **Expected Biological Signatures**:
  - Down-regulated: *FOXA1*, *GATA3*, *ESR1*, *PGR*, *ERBB2* (luminal and HER2 lineage markers).
  - Up-regulated: *MKI67* (proliferation), *FOXM1*, *SOX10*, *VGLL1*, basal cytokeratins (*KRT5*, *KRT14*, *KRT17*).
  - Pathways: Mitotic Cell Cycle, DNA Replication, p53 Signaling, Extracellular Matrix Degradation.

---

## 3. Data Lineage & Pipeline Parameters
- **Primary Source**: TCGA-BRCA RNA-Seq count matrix downloaded programmatically via `TCGAbiolinks`.
- **Secondary Validation**: NCBI GEO series (e.g., GSE58135 or GSE76275) ingested via `GEOquery`.
- **Target Cohort Size**: $> 100$ samples (TNBC tumor and matched normal breast tissue).
- **Differential Expression Cutoffs**:
  - Adjusted $p$-value (Benjamini-Hochberg FDR) $< 0.01$.
  - $|\log_2\text{Fold Change}| > 2.0$ (representing a minimum 4-fold biological change).
- **Statistical Model**:
  - Negative binomial Generalized Linear Model with empirical Bayes dispersion shrinkage via `DESeq2`.
  - Independent cross-check using `edgeR` quasi-likelihood F-test.
- **Survival Stratification**:
  - Kaplan-Meier estimator with Log-Rank test ($p < 0.05$).
  - Stratification by median or maximally selected rank statistics (`surv_cutpoint`).
  - Univariate and multivariate Cox proportional hazards modeling (`coxph`).

---

## 4. Architectural & Delivery Strategy (Pure R & GitHub Pages)
- **Delivery Model Decision**:
  - Cloud microservices (Render/Vercel) were intentionally bypassed to eliminate external hosting maintenance, sleeping backend servers, and unnecessary complexity.
  - The project is architected as an **end-to-end, hermetic R computational biology pipeline** with automated compilation of an interactive **Quarto / R Markdown report deployed to GitHub Pages**.
- **Interactive Capabilities on Static Edge**:
  - Leverages client-side JavaScript libraries directly embedded by R:
    - `DT (DataTables)` for sub-millisecond searching, filtering, and CSV export of 15,000+ gene records.
    - `plotly` for zoomable, hover-inspectable volcano and PCA plots.
  - Accompanied by static, 300 DPI vector publication figures (`figures/`).

---

## 5. Security & Aesthetic Constraints (Non-Negotiable)
- **Aesthetic**:
  - No purple-to-blue gradients, no gradient hero text, no emojis in headings, no Inter font everywhere, no colored border cards, no glassmorphism, no low-contrast dark mode, no 3 icon boxes in a row, no floating badges above headlines, no Lucide icons plastered everywhere, no untouched Shadcn defaults, no fade-in on scroll, no cursor beams, no button opacity fades, no inconsistent spacing, no em dashes everywhere, no generic passwords, no serif italic accents, no Space Grotesk + Instrument Serif, no grain over gradients.
- **Security & Data Safety**:
  - Zero API keys or secrets in repository.
  - `.gitignore` strictly protects `.env*`, `.Renviron`, and raw data caches (`data/raw/`).
  - Hermetic reproducibility guaranteed through `renv.lock` and Dockerfile.

---

## 6. Current Project State
- **Stage**: Project documentation, governance, and planning complete.
- **Active Documents**:
  - `prd.md`: Product Requirements Document (Pure R + Quarto + GitHub Pages).
  - `architecture.md`: Pipeline modular structure, artifact flow, and CI/CD.
  - `architecture_diagram.md`: Visual workflow & step-by-step GitHub Pages deployment guide.
  - `techstack.md`: Streamlined R, Bioconductor, Quarto, DT, and Docker stack.
  - `aistack.md`: AI governance, prompt grounding, and safeguards.
  - `rules.md`: Strict aesthetic, bioinformatics, and security rules.
  - `design.md`: Clinical laboratory UI specification for publication plots & report.
  - `task.md`: GSD-style phase execution plan (Phases 1–8).
  - `memory.md`: Project memory and context reference.
- **Next Step**: Phase 1 execution (directory creation, R `renv` setup, and data acquisition script authoring).
