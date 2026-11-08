# ==============================================================================
# Script: 02_preprocessing_qc.R
# Description: Quality control, low-count filtering, VST transformation, and PCA
# Project: BioMarker-R (Bulk RNA-Seq & Biomarker Discovery)
# ==============================================================================

source("scripts/00_utils.R")

# ------------------------------------------------------------------------------
# 1. Pipeline Initialization
# ------------------------------------------------------------------------------
set_reproducible_seed(42)
paths <- get_project_paths()
log_info("Starting Step 02: Preprocessing, Harmonization & Quality Control...")

raw_counts_file <- file.path(paths$data_raw, "tcga_brca_counts.rds")
metadata_file   <- file.path(paths$data_raw, "tcga_brca_metadata.rds")

if (!file.exists(raw_counts_file) || !file.exists(metadata_file)) {
  log_warn("Raw datasets not found. Executing 01_data_download.R first...")
  source("scripts/01_data_download.R")
}

counts   <- readRDS(raw_counts_file)
metadata <- readRDS(metadata_file)

assert_sample_parity(counts, metadata, "sample_id")

# ------------------------------------------------------------------------------
# 2. Low-Count Gene Filtering
#    Threshold: retain genes with count >= 10 in at least N samples
#    (N equals smallest group size, e.g. normal cohort)
# ------------------------------------------------------------------------------
min_group_size <- min(table(metadata$condition))
keep_genes <- rowSums(counts >= 10) >= min_group_size
filtered_counts <- counts[keep_genes, ]

log_info(sprintf("Gene filtering: retained %d / %d genes (%.1f%% passed cutoff count >= 10 in >= %d samples).",
                 nrow(filtered_counts), nrow(counts),
                 (nrow(filtered_counts) / nrow(counts)) * 100, min_group_size))

# ------------------------------------------------------------------------------
# 3. Variance Stabilization Transformation (VST)
#    DESeq2 vst() or robust log2-CPM transformation fallback
# ------------------------------------------------------------------------------
vst_data <- NULL

if (requireNamespace("DESeq2", quietly = TRUE)) {
  log_info("Applying DESeq2::vst(blind = TRUE) for unsupervised exploratory QC...")
  col_data <- data.frame(
    condition = metadata$condition,
    row.names = metadata$sample_id
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData   = col_data,
    design    = ~ condition
  )
  vsd <- DESeq2::vst(dds, blind = TRUE)
  vst_data <- SummarizedExperiment::assay(vsd)
} else {
  log_info("DESeq2 not installed in current session; applying log2(CPM + 1) variance stabilization...")
  # Standard log2 CPM transformation
  lib_sizes <- colSums(filtered_counts)
  cpm <- t(t(filtered_counts) / lib_sizes) * 1e6
  vst_data <- log2(cpm + 1)
}

# ------------------------------------------------------------------------------
# 4. Principal Component Analysis (PCA)
#    Computed across top 500 most variable genes
# ------------------------------------------------------------------------------
log_info("Computing Principal Component Analysis on top variable genes...")

gene_vars <- apply(vst_data, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(500, length(gene_vars)))]
pca_matrix <- t(vst_data[top_genes, ])

pca_res <- prcomp(pca_matrix, scale. = TRUE)
var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100

pca_df <- data.frame(
  sample_id = rownames(pca_matrix),
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  condition = metadata$condition,
  age       = metadata$age,
  stage     = metadata$stage,
  stringsAsFactors = FALSE
)

# Save processed PCA coordinates
pca_csv_path <- file.path(paths$data_proc, "pca_sample_coordinates.csv")
write.csv(pca_df, pca_csv_path, row.names = FALSE)
log_success(sprintf("PCA coordinates written to %s", pca_csv_path))

# ------------------------------------------------------------------------------
# 5. Publication-Grade PCA Figure Generation (High-Contrast Lab Standard)
# ------------------------------------------------------------------------------
fig_pdf_path <- file.path(paths$figures, "fig1_pca_qc.pdf")
fig_png_path <- file.path(paths$figures, "fig1_pca_qc.png")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  pc1_label <- sprintf("Principal Component 1 (%.1f%% variance)", var_explained[1])
  pc2_label <- sprintf("Principal Component 2 (%.1f%% variance)", var_explained[2])
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3.2, alpha = 0.88, shape = 16) +
    scale_color_manual(
      name   = "Condition",
      values = c("Normal" = LAB_PALETTE$normal_cohort, "TNBC" = LAB_PALETTE$tumor_cohort),
      labels = c("Normal" = "Solid Tissue Normal", "TNBC" = "Triple-Negative Tumor")
    ) +
    labs(
      title    = "Sample Stratification: Principal Component Analysis",
      subtitle = sprintf("Variance-stabilized expression of top %d variable genes across %d patient samples",
                         length(top_genes), nrow(pca_df)),
      x = pc1_label,
      y = pc2_label
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 14, color = LAB_PALETTE$text_primary),
      plot.subtitle      = element_text(size = 11, color = "#475569", margin = margin(b = 12)),
      axis.title         = element_text(face = "bold", size = 11, color = LAB_PALETTE$text_primary),
      axis.text          = element_text(color = "#334155", size = 10),
      axis.line          = element_line(color = "#94A3B8", linewidth = 0.6),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 11),
      legend.background  = element_rect(fill = LAB_PALETTE$bg_panel, color = LAB_PALETTE$border, linewidth = 0.5),
      panel.grid.major   = element_line(color = "#F1F5F9", linewidth = 0.5),
      plot.margin        = margin(16, 16, 16, 16)
    )
  
  ggsave(fig_pdf_path, plot = p_pca, width = 7.5, height = 5.5, device = "pdf")
  ggsave(fig_png_path, plot = p_pca, width = 7.5, height = 5.5, dpi = 300)
  log_success(sprintf("Publication PCA figure saved: %s & %s", fig_pdf_path, fig_png_path))
} else {
  # Base R high-resolution rendering fallback
  png(fig_png_path, width = 2250, height = 1650, res = 300)
  par(mar = c(5, 5, 4, 2), bg = "#FFFFFF", bty = "l", las = 1)
  cols <- ifelse(pca_df$condition == "Normal", LAB_PALETTE$normal_cohort, LAB_PALETTE$tumor_cohort)
  plot(pca_df$PC1, pca_df$PC2, col = cols, pch = 16, cex = 1.5,
       xlab = sprintf("PC1 (%.1f%% variance)", var_explained[1]),
       ylab = sprintf("PC2 (%.1f%% variance)", var_explained[2]),
       main = "Sample Stratification: Principal Component Analysis (QC)")
  grid(col = "#E2E8F0", lty = 1)
  legend("topright", legend = c("Solid Tissue Normal", "Triple-Negative Tumor"),
         col = c(LAB_PALETTE$normal_cohort, LAB_PALETTE$tumor_cohort), pch = 16, pt.cex = 1.5, bty = "n")
  dev.off()
  log_success(sprintf("Base R publication PCA figure saved to %s", fig_png_path))
}

# Save normalized VST expression matrix for downstream survival and heatmap analysis
saveRDS(filtered_counts, file.path(paths$data_proc, "counts_filtered.rds"))
saveRDS(vst_data, file.path(paths$data_proc, "vst_normalized_matrix.rds"))
saveRDS(metadata, file.path(paths$data_proc, "metadata_cleaned.rds"))

log_success("Step 02 Preprocessing & Quality Control completed successfully.")
