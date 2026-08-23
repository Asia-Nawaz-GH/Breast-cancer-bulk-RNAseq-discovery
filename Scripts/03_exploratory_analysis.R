# ==============================================================================
# Script: 03_pca_analysis.R
# Purpose: Exploratory Data Analysis & Principal Component Analysis (PCA)
# Author: Asia Nawaz
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Loading Cleaned Data Objects & Required Libraries
# ------------------------------------------------------------------------------
library(ggplot2)
library(ggrepel)

# Load binary objects from Results directory
expr_matrix <- readRDS("Results/clean_expr_matrix.rds")
pheno_data  <- readRDS("Results/pheno_metadata.rds")

# ------------------------------------------------------------------------------
# 2. Computing Principal Component Analysis (PCA)
# ------------------------------------------------------------------------------
# Transpose matrix so samples are rows and probes are columns
pca_res <- prcomp(t(expr_matrix), scale. = TRUE)

# Calculate percentage of total variance explained by PC1 and PC2
var_explained <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 1)

# ------------------------------------------------------------------------------
# 3. Structuring PCA Coordinate Data Frame
# ------------------------------------------------------------------------------
pca_df <- as.data.frame(pca_res$x)
pca_df$Sample_ID <- rownames(pca_df)

# Merge metadata directly from pheno_data
pca_df$Group     <- pheno_data$Group[match(pca_df$Sample_ID, pheno_data$Sample_ID)]
pca_df$Replicate <- pheno_data$Replicate[match(pca_df$Sample_ID, pheno_data$Sample_ID)]
pca_df$Label     <- paste(pca_df$Group, pca_df$Replicate)

# ------------------------------------------------------------------------------
# 4. Rendering ggplot2 PCA Scatter Plot
# ------------------------------------------------------------------------------
group_palette <- c("Control" = "#4E79A7", "DPBQ" = "#E15759")

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  # Add grid origin lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", linewidth = 0.5) +
  
  # Plot sample points
  geom_point(size = 5, alpha = 0.9) +
  
  # Repel text labels dynamically to eliminate text overlapping
  geom_text_repel(
    aes(label = Label),
    fontface      = "bold",
    size          = 4,
    box.padding   = 0.6,
    point.padding = 0.5,
    max.overlaps  = Inf,
    show.legend   = FALSE
  ) +
  
  # Color scaling and limits
  scale_color_manual(values = group_palette) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
  
  # Labels and themes
  labs(
    title    = "Exploratory Data Analysis: Principal Component Analysis (PCA)",
    subtitle = "Dataset: GSE73710 (Breast Cancer Cell Lines)",
    x        = paste0("PC1: ", var_explained[1], "% Variance Explained"),
    y        = paste0("PC2: ", var_explained[2], "% Variance Explained"),
    color    = "Treatment Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, margin = margin(b = 4)),
    plot.subtitle    = element_text(face = "italic", color = "gray30", size = 10, margin = margin(b = 10)),
    axis.title       = element_text(face = "bold", size = 10),
    legend.title     = element_text(face = "bold", size = 10),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------------------
# 5. Saving High-Resolution Plot
# ------------------------------------------------------------------------------
ggsave(
  filename = "Figures/02_pca_plot.png",
  plot     = pca_plot,
  width    = 7.5,
  height   = 5.5,
  dpi      = 300
)

cat("\nScript 03 execution complete. Saved non-overlapping PCA plot to Figures/02_pca_plot.png\n")