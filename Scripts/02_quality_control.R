# ==============================================================================
# Script: 02_quality_control.R
# Purpose: Quality Control & Intensity Distribution Analysis
# Author: Asia Nawaz
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Loading Cleaned Data Objects from Binary Storage
# ------------------------------------------------------------------------------
# 1. Load cleaned data objects from the Results folder
expr_matrix <- readRDS("Results/clean_expr_matrix.rds")
pheno_data  <- readRDS("Results/pheno_metadata.rds")

# ------------------------------------------------------------------------------
# 2. Summary Statistics & Data Range Inspection
# ------------------------------------------------------------------------------
# 2. Inspect overall summary statistics across all matrix values
summary(as.vector(expr_matrix))

# 3. Inspect the exact minimum and maximum values directly
range(expr_matrix)

# 4. Re-check range by explicitly ignoring missing (NA) values
range(expr_matrix, na.rm = TRUE)

# ------------------------------------------------------------------------------
# 3. Cleaning NA Values & Resaving Matrix
# ------------------------------------------------------------------------------
# 5. Remove probe rows containing NA values to ensure clean matrix
expr_matrix <- na.omit(expr_matrix)

# 6. Verify that NAs are completely gone
sum(is.na(expr_matrix))

# Update clean matrix on disk to guarantee zero NAs in downstream analyses
saveRDS(expr_matrix, "Results/clean_expr_matrix.rds")
# ------------------------------------------------------------------------------
# 7. Defining Plotting Parameters & Reorder Samples Group-Wise
# ------------------------------------------------------------------------------
# Set colorblind-friendly publication palette
group_palette <- c("Control" = "#4E79A7", "DPBQ" = "#E15159")

# Construct full sample labels
raw_sample_labels <- paste(pheno_data$Group, pheno_data$Replicate)

# Create desired sample order: Controls first, then DPBQ treatments
desired_order <- c("Control 1", "Control 2", "Control 3", "DPBQ 1", "DPBQ 2", "DPBQ 3")
order_indices <- match(desired_order, raw_sample_labels)

# Reorder matrix columns, metadata groups, and colors to group controls together
expr_matrix   <- expr_matrix[, order_indices]
sample_labels <- raw_sample_labels[order_indices]
sample_groups <- as.character(pheno_data$Group)[order_indices]
sample_colors <- group_palette[sample_groups]

# ------------------------------------------------------------------------------
# 8. Initializing High-Resolution Graphic Device (300 DPI PNG)
# ------------------------------------------------------------------------------
# Ensure Figures directory exists
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

png(
  filename = "Figures/01_qc_intensity_boxplot.png",
  width    = 7,
  height   = 5.5,
  units    = "in",
  res      = 300
)

# Set publication margins: c(bottom, left, top, right)
par(
  mar      = c(5.5, 4.5, 4.8, 1.5),
  mgp      = c(2.5, 0.7, 0),
  font.lab = 2,
  cex.lab  = 1.0,
  cex.axis = 0.95
)

# ------------------------------------------------------------------------------
# 9. Rendering Intensity Boxplot
# ------------------------------------------------------------------------------
# Draw primary boxplot structure with expanded Y-limits (1 to 19 for headroom)
boxplot(
  expr_matrix,
  names   = sample_labels,
  col     = sample_colors,
  border  = "gray20",
  outline = FALSE,
  las     = 2,
  ylim    = c(1, 19),
  main    = "",
  xlab    = "",
  ylab    = expression(Log[2] ~ "Signal Intensity")
)

# Add subtle background gridlines for precise median comparison
grid(nx = NA, ny = NULL, col = "gray85", lty = "dotted", lwd = 1)

# Redraw boxplot outlines over gridlines for crisp rendering
boxplot(
  expr_matrix,
  names   = sample_labels,
  col     = sample_colors,
  border  = "gray20",
  outline = FALSE,
  las     = 2,
  ylim    = c(1, 19),
  add     = TRUE
)

# Add main title and top subtitle with custom vertical spacing
mtext("Quality Control: Log2 Signal Intensity Distribution", side = 3, line = 2.4, font = 2, cex = 1.1)
mtext("Dataset: GSE73710 (Breast Cancer Cell Lines)", side = 3, line = 0.9, font = 3, cex = 0.85, col = "gray30")

# Add custom X-axis label centered below sample names
mtext("Experimental Samples", side = 1, line = 4.2, font = 2, cex = 1.0)

# Add publication legend with white background fill
legend(
  "topright",
  legend     = names(group_palette),
  fill       = group_palette,
  border     = "gray20",
  title      = "Treatment Group",
  title.font = 2,
  bg         = "white",
  box.col    = "gray70",
  cex        = 0.80
)

# ------------------------------------------------------------------------------
# 10. Finalizing Image File Connection
# ------------------------------------------------------------------------------
# Close PNG graphics device and flush image output to disk
dev.off()

cat("\nScript 02 execution complete. Saved publication boxplot to Figures/01_qc_intensity_boxplot.png\n")