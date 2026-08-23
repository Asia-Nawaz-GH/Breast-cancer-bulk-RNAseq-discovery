# ==============================================================================
# Script: 05_visualization.R
# Purpose: Volcano Plot & Heatmap generation
# Author: Asia Nawaz
# ==============================================================================

# ------------------------------------------------------------------------------
# Loading required Packages, Libraries & Adjusting values
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) BiocManager::install("hgu133plus2.db")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(hgu133plus2.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)

# Load data objects
full_results <- readRDS("Results/deg_full_results.rds")
expr_matrix  <- readRDS("Results/clean_expr_matrix.rds")
pheno_data   <- readRDS("Results/pheno_metadata.rds")

# Remove missing values
full_results <- full_results[!is.na(full_results$logFC) & !is.na(full_results$P.Value), ]

# Map Probe IDs to official HUGO Gene Symbols
full_results$Gene_Symbol <- mapIds(
  hgu133plus2.db,
  keys      = full_results$Probe_ID,
  column    = "SYMBOL",
  keytype   = "PROBEID",
  multiVals = "first"
)

# Manual fallback override for specific unannotated probes
full_results$Gene_Symbol[full_results$Probe_ID == "210524_x_at"] <- "CYTOR"
full_results$Gene_Symbol[full_results$Probe_ID == "217165_x_at"] <- "TNFRSF10D"
full_results$Gene_Symbol[full_results$Probe_ID == "213629_x_at"] <- "AKR1C3"

# Fill missing gene symbols
missing_genes <- is.na(full_results$Gene_Symbol) | full_results$Gene_Symbol == ""
full_results$Gene_Symbol[missing_genes] <- full_results$Probe_ID[missing_genes]

# Define screening thresholds
p_raw_cutoff <- 0.01
lfc_cutoff   <- 0.585

# Calculate -log10(P.Value)
full_results$neg_log10_p <- -log10(full_results$P.Value)

# Assign categorical expression status based on Adjusted P-value (FDR)
full_results$Expression <- "Not Significant"
full_results$Expression[full_results$adj.P.Val < 0.05 & full_results$logFC >= lfc_cutoff]  <- "Upregulated"
full_results$Expression[full_results$adj.P.Val < 0.05 & full_results$logFC <= -lfc_cutoff] <- "Downregulated"

full_results$Expression <- factor(
  full_results$Expression, 
  levels = c("Upregulated", "Downregulated", "Not Significant")
)

# Calculate -log10(adj.P.Val) for proper y-axis scaling
full_results$neg_log10_fdr <- -log10(full_results$adj.P.Val)

# Deduplicate by Gene_Symbol for top labeling (take lowest FDR per gene)
full_results_sorted <- full_results[order(full_results$adj.P.Val), ]
unique_genes_df    <- full_results_sorted[!duplicated(full_results_sorted$Gene_Symbol), ]
unique_genes_df    <- unique_genes_df[!grepl("_at$", unique_genes_df$Gene_Symbol), ]

# Select Top 8 Upregulated and Top 8 Downregulated genes (16 total)
top_up_labels <- unique_genes_df %>% 
  filter(Expression == "Upregulated") %>% 
  arrange(adj.P.Val) %>% 
  head(8)

top_down_labels <- unique_genes_df %>% 
  filter(Expression == "Downregulated") %>% 
  arrange(adj.P.Val) %>% 
  head(8)

top_genes <- rbind(top_up_labels, top_down_labels)

# ------------------------------------------------------------------------------
# 2. Creating Volcano Plot
# ------------------------------------------------------------------------------

y_max <- max(full_results$neg_log10_fdr[!is.infinite(full_results$neg_log10_fdr)], na.rm = TRUE) * 1.15

volcano_plot <- ggplot(full_results, aes(x = logFC, y = neg_log10_fdr, color = Expression)) +
  geom_point(alpha = 0.60, size = 1.8) +
  scale_color_manual(values = c(
    "Upregulated"     = "#B2182C", 
    "Downregulated"   = "#2B5C8F", 
    "Not Significant" = "gray80"
  )) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.5) +
  # Separate left labeling to push labels outward to the left
  geom_text_repel(
    data               = filter(top_genes, Expression == "Downregulated"),
    aes(label          = Gene_Symbol),
    color              = "black",
    size               = 3.3,
    fontface           = "bold",
    nudge_x            = -0.8,
    direction          = "y",
    hjust              = 1,
    box.padding        = 0.35,
    point.padding      = 0.3,
    segment.color      = "gray30",
    segment.size       = 0.35,
    max.overlaps       = Inf,
    min.segment.length = 0
  ) +
  # Separate right labeling to push labels outward to the right
  geom_text_repel(
    data               = filter(top_genes, Expression == "Upregulated"),
    aes(label          = Gene_Symbol),
    color              = "black",
    size               = 3.3,
    fontface           = "bold",
    nudge_x            = 0.8,
    direction          = "y",
    hjust              = 0,
    box.padding        = 0.35,
    point.padding      = 0.3,
    segment.color      = "gray30",
    segment.size       = 0.35,
    max.overlaps       = Inf,
    min.segment.length = 0
  ) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0.01, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(
    title    = "Differential Expression Landscape: DPBQ vs. Control",
    subtitle = "Dataset: GSE73710 | Thresholds: Adj. P < 0.05 & |Log2FC| >= 0.585",
    x        = expression(bold("Log"[2]*" Fold Change")),
    y        = expression(bold("-Log"[10]*"(FDR)")),
    color    = "Status"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle     = element_text(face = "italic", size = 9.5, hjust = 0.5, color = "gray30"),
    axis.title        = element_text(face = "bold", size = 11),
    axis.text         = element_text(color = "black", size = 10),
    legend.position   = "top",
    legend.title      = element_text(face = "bold", size = 10),
    legend.background = element_rect(fill = "gray98", color = "gray80"),
    legend.key        = element_blank(),
    plot.margin       = margin(t = 15, r = 15, b = 10, l = 10)
  )

ggsave(
  filename = "Figures/03_volcano_plot.png",
  plot     = volcano_plot,
  width    = 8.5,
  height   = 6.5,
  dpi      = 300
)

cat("Volcano plot successfully exported to Figures/03_volcano_plot.png\n")

# ------------------------------------------------------------------------------
# 3. Creating Top DEGs Heatmap
# ------------------------------------------------------------------------------

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap", update = FALSE)
}
library(ComplexHeatmap)
library(circlize)

# Select top 20 Upregulated and top 20 Downregulated genes
top_up <- unique_genes_df %>% 
  filter(Expression == "Upregulated") %>% 
  arrange(desc(logFC)) %>% 
  head(20)

top_down <- unique_genes_df %>% 
  filter(Expression == "Downregulated") %>% 
  arrange(logFC) %>% 
  head(20)

top_40_degs <- rbind(top_up, top_down)

# Extract expression matrix rows using Probe IDs
heatmap_matrix <- expr_matrix[top_40_degs$Probe_ID, ]
rownames(heatmap_matrix) <- top_40_degs$Gene_Symbol

# Z-score scaling across rows
heatmap_matrix_z <- t(scale(t(heatmap_matrix)))

# Order samples explicitly
sample_order <- c("GSM1901163", "GSM1901165", "GSM1901167", "GSM1901164", "GSM1901166", "GSM1901168")
heatmap_matrix_z <- heatmap_matrix_z[, sample_order]

# Square-root height scaling for dendrogram branches
row_dist <- dist(heatmap_matrix_z, method = "euclidean")
row_hclust <- hclust(row_dist, method = "ward.D2")
row_hclust$height <- sqrt(row_hclust$height)

# Set color mapping for Z-Score gradient (-1.5 to +1.5)
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2B5C8F", "#F7F7F7", "#B2182C"))

# Construct Top Annotation (Treatment Bar above expression grid)
top_ann <- HeatmapAnnotation(
  Treatment = pheno_data[sample_order, "Group"],
  col = list(Treatment = c("Control" = "#2B5C8F", "DPBQ" = "#B2182C")),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Treatment = list(
      title = "Treatment",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  )
)

# Define clean sample labels mapping directly to sample_order
clean_sample_labels <- c(
  "GSM1901163" = "Control 1",
  "GSM1901165" = "Control 2",
  "GSM1901167" = "Control 3",
  "GSM1901164" = "DPBQ 1",
  "GSM1901166" = "DPBQ 2",
  "GSM1901168" = "DPBQ 3"
)

# Render publication-ready heatmap to file
png("Figures/04_expression_heatmap.png", width = 7.0, height = 10.0, units = "in", res = 300)

Heatmap(
  heatmap_matrix_z,
  name = "Z-Score",  # Continuous color bar title
  col = col_fun,
  top_annotation = top_ann,
  cluster_rows = row_hclust,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_labels = clean_sample_labels[colnames(heatmap_matrix_z)], # Displays Control 1-3, DPBQ 1-3
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,
  row_dend_width = unit(2.2, "cm"),
  column_title = "Hierarchical Clustering of Top 40 DEGs",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

dev.off()

cat("Publication-ready heatmap generated with precise legend separation.\n")