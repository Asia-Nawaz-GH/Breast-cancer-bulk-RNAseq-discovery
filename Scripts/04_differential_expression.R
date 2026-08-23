# ==============================================================================
# Script: 04_differential_expression.R
# Purpose: Differential Gene Expression Analysis via Empirical Bayes Moderation
# Author: Asia Nawaz
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Loading Cleaned Data Objects & Microarray Differential Analysis Suite
# ------------------------------------------------------------------------------
library(limma)

# Load binary objects from Results directory
expr_matrix <- readRDS("Results/clean_expr_matrix.rds")
pheno_data  <- readRDS("Results/pheno_metadata.rds")

# ------------------------------------------------------------------------------
# 2. Constructing Model Matrix & Parameterize Experimental Design
# ------------------------------------------------------------------------------
# Build full-rank design matrix specifying group assignments without intercept
design <- model.matrix(~ 0 + Group, data = pheno_data)

# Fix coefficient column names to clean group labels
colnames(design) <- levels(pheno_data$Group)

cat("=== DESIGN MATRIX STRUCTURE ===\n")
print(design)

# ------------------------------------------------------------------------------
# 3. Formulating Contrast Matrix for Target Comparison (DPBQ vs. Control)
# ------------------------------------------------------------------------------
# Define explicit contrast: DPBQ (treatment) minus Control (baseline)
contrast_matrix <- makeContrasts(
  DPBQ_vs_Control = DPBQ - Control,
  levels = design
)

cat("\n=== CONTRAST MATRIX ===\n")
print(contrast_matrix)

# ------------------------------------------------------------------------------
# 4. Fitting Linear Models & Executing Empirical Bayes Moderation
# ------------------------------------------------------------------------------
# Fit gene-wise linear model across all probes
fit <- lmFit(expr_matrix, design)

# Estimate contrast coefficients from fitted model
fit_contrasts <- contrasts.fit(fit, contrast_matrix)

# Compute Empirical Bayes statistics for moderated t-tests
ebayes_fit <- eBayes(fit_contrasts)

# ------------------------------------------------------------------------------
# 5. Extracting Differential Expression Table & Saving Results
# ------------------------------------------------------------------------------
# Extract full genome results table sorted by adjusted p-value
deg_results <- topTable(
  ebayes_fit,
  coef = "DPBQ_vs_Control",
  number = Inf,
  adjust.method = "BH"
)

# Preserve Probe IDs as explicit column
deg_results$Probe_ID <- rownames(deg_results)

# Print top 10 differentially expressed genes to console
cat("\n=== TOP 10 DIFFERENTIALLY EXPRESSED PROBES ===\n")
print(head(deg_results[, c("Probe_ID", "logFC", "P.Value", "adj.P.Val")], 10))

# Export complete results table to disk as binary .rds and tabular text
saveRDS(deg_results, "Results/deg_full_results.rds")
write.table(
  deg_results,
  file = "Results/deg_full_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nScript 04 execution complete. Saved differential expression results to Results/\n")