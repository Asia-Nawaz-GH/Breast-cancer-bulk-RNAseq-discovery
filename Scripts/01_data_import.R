# ==============================================================================
# Script: 01_data_import_setup.R
# Purpose: Data Import, Metadata Structuring & Binary Serialization
# Author: Asia Nawaz
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Parsing Raw File Header & Extracting Sample Barcodes
# ------------------------------------------------------------------------------
# Read raw metadata headers from the expression dataset
header_lines <- readLines("Data/BC_Dataset.txt", n = 2)

# Extract GEO accession barcodes (GSM IDs) from header line
raw_header <- unlist(strsplit(header_lines[1], "\t"))
sample_ids <- raw_header[grepl("GSM", raw_header)]

# Print extracted sample IDs to console for verification
cat("Extracted Sample IDs:\n")
print(sample_ids)

# ------------------------------------------------------------------------------
# 2. Building Phenotype Metadata Data Frame (pheno_data)
# ------------------------------------------------------------------------------
# Map experimental groups alternating per GSM accession order (3 Control, 3 DPBQ)
# GSM1901163=Ctrl1, GSM1901164=DPBQ1, GSM1901165=Ctrl2, GSM1901166=DPBQ2...
pheno_data <- data.frame(
  Sample_ID = sample_ids,
  Group     = factor(rep(c("Control", "DPBQ"), times = 3), levels = c("Control", "DPBQ")),
  Replicate = c(1, 1, 2, 2, 3, 3),
  row.names = sample_ids,
  stringsAsFactors = FALSE
)

# Print pheno_data structure to console for verification
cat("\nPhenotype Metadata Structure:\n")
print(pheno_data)

# ------------------------------------------------------------------------------
# 3. Loading Expression Intensity Matrix (expr_matrix)
# ------------------------------------------------------------------------------
# Read full tabular expression dataset into memory
raw_data <- read.delim(
  "Data/BC_Dataset.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Set probe ID column as row names and isolate numerical signal matrix
rownames(raw_data) <- raw_data[[1]]
expr_matrix <- as.matrix(raw_data[, sample_ids])

# Remove missing probe values (NAs) to guarantee numerical matrix clean state
expr_matrix <- na.omit(expr_matrix)

# Ensure numerical matrix dimensions and data type integrity
storage.mode(expr_matrix) <- "numeric"

# Print matrix dimensions to console (Probes x Samples)
cat("\nExpression Matrix Dimensions (Probes x Samples):\n")
print(dim(expr_matrix))

# ------------------------------------------------------------------------------
# 4. Serializing Cleaned Objects to Disk (.rds)
# ------------------------------------------------------------------------------
# Create Results directory if it does not already exist
if (!dir.exists("Results")) {
  dir.create("Results")
}

# Save phenotype metadata and expression matrix as binary .rds files
saveRDS(pheno_data, file = "Results/pheno_metadata.rds")
saveRDS(expr_matrix, file = "Results/clean_expr_matrix.rds")

cat("\nScript 01 execution complete. Saved binary files to Results/ directory.\n")