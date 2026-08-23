# ==============================================================================
# Script: 06_pathway_enrichment.R
# Purpose: Functional Enrichment Analysis (GO, KEGG, and GSEA)
# Author: Asia Nawaz
# ==============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(hgu133plus2.db)
  library(msigdbr)
  library(ggplot2)
  library(enrichplot)
  library(dplyr)
  library(stringr)
})
# ------------------------------------------------------------------------------
# 1. Loading Annotated Differential Expression Results
# ------------------------------------------------------------------------------

full_results <- readRDS("Results/deg_full_results.rds")
mapped_ids <- bitr(full_results$Probe_ID, fromType = "PROBEID", toType = "ENTREZID", OrgDb = hgu133plus2.db)
full_annotated <- merge(full_results, mapped_ids, by.x = "Probe_ID", by.y = "PROBEID")
full_annotated <- full_annotated[!is.na(full_annotated$ENTREZID) & !duplicated(full_annotated$ENTREZID), ]

# Filter DEGs for ORA
fdr_cutoff <- 0.05
lfc_cutoff <- 0.585
deg_genes  <- full_annotated[full_annotated$adj.P.Val < fdr_cutoff & abs(full_annotated$logFC) >= lfc_cutoff, ]
valid_entrez <- unique(deg_genes$ENTREZID)

# Helper function to parse ORA data frames
parse_ora_df <- function(ora_obj) {
  df <- as.data.frame(ora_obj)
  if (nrow(df) == 0) return(df)
  df$GeneRatio_Numeric <- sapply(df$GeneRatio, function(x) {
    parts <- as.numeric(strsplit(x, "/")[[1]])
    return(parts[1] / parts[2])
  })
  df$neg_log10_p <- -log10(df$p.adjust)
  return(df)
}
# ------------------------------------------------------------------------------
# 2. GO Biological Process (ORA) Analysis
# ------------------------------------------------------------------------------

go_bp_raw <- enrichGO(gene = valid_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

# Prune redundant parent/child GO terms (cutoff 0.7 removes overlapping child terms)
go_bp <- clusterProfiler::simplify(go_bp_raw, cutoff = 0.7, by = "p.adjust", select_fun = min)
top_go <- parse_ora_df(go_bp)[1:min(15, nrow(as.data.frame(go_bp))), ]

p_go <- ggplot(top_go, aes(x = GeneRatio_Numeric, y = reorder(Description, GeneRatio_Numeric))) +
  geom_point(aes(size = Count, color = neg_log10_p)) +
  scale_color_gradient(low = "#2B5C8F", high = "#B2182C", name = "-log10(Adj P)") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38)) +
  theme_classic(base_size = 11) +
  labs(title = "Top Enriched GO Biological Processes", x = "Gene Ratio", y = NULL, size = "Gene Count") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title  = element_text(face = "bold", size = 11))

ggsave("Figures/05_go_bp_enrichment.png", plot = p_go, width = 8.5, height = 7.0, dpi = 300)

# ------------------------------------------------------------------------------
# 3. KEGG Pathway Analysis (ORA) Analysis
# ------------------------------------------------------------------------------

kegg_res <- enrichKEGG(gene = valid_entrez, organism = "hsa", pvalueCutoff = 0.05)
top_kegg <- parse_ora_df(kegg_res)[1:min(15, nrow(as.data.frame(kegg_res))), ]

p_kegg <- ggplot(top_kegg, aes(x = GeneRatio_Numeric, y = reorder(Description, GeneRatio_Numeric))) +
  geom_point(aes(size = Count, color = neg_log10_p)) +
  scale_color_gradient(low = "#2B5C8F", high = "#B2182C", name = "-log10(Adj P)") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38)) +
  theme_classic(base_size = 11) +
  labs(title = "Top Enriched KEGG Pathways", x = "Gene Ratio", y = NULL, size = "Gene Count") +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title  = element_text(face = "bold", size = 11))

ggsave("Figures/06_kegg_enrichment.png", plot = p_kegg, width = 8.5, height = 7.0, dpi = 300)

# ------------------------------------------------------------------------------
# 4. Gene Set Enrichment Analysis (GSEA)
# ------------------------------------------------------------------------------

full_annotated$rank_metric <- sign(full_annotated$logFC) * (-log10(full_annotated$adj.P.Val))
gene_list <- sort(setNames(full_annotated$rank_metric, full_annotated$ENTREZID), decreasing = TRUE)
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

gsea_res <- GSEA(geneList = gene_list, TERM2GENE = hallmark_sets, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)

# Native Title Annotation Exporter for GSEA
export_gsea_native <- function(gsea_obj, gene_set_id, title_text, output_path) {
  res_df <- as.data.frame(gsea_obj)
  
  target_row <- res_df[res_df$ID == gene_set_id | res_df$Description == gene_set_id, ]
  
  if (nrow(target_row) == 0) {
    stop(paste("Gene set ID not found in GSEA object:", gene_set_id))
  }
  
  # Format NES explicitly to 2 decimal places (ensures trailing zeros aren't dropped, e.g., 2.40)
  nes_val <- sprintf("%.2f", target_row$NES[1])
  raw_fdr <- target_row$p.adjust[1]
  
  # Handle permutation floor threshold
  if (is.na(raw_fdr) || raw_fdr <= 1e-08) {
    fdr_str <- "< 0.001"
  } else {
    fdr_str <- paste0("= ", formatC(raw_fdr, format = "e", digits = 2))
  }
  
  full_title <- paste0(title_text, " (NES = ", nes_val, ", FDR ", fdr_str, ")")
  
  p <- gseaplot2(
    gsea_obj,
    geneSetID = target_row$ID[1],
    title     = full_title,
    color     = "#B2182C",
    base_size = 11,
    subplots  = 1:3
  )
  
  ggsave(filename = output_path, plot = p, width = 8.5, height = 6.5, dpi = 300)
}
export_gsea_native(gsea_res, "HALLMARK_HYPOXIA", "GSEA: Hallmark Hypoxia", "Figures/07_gsea_hypoxia.png")
export_gsea_native(gsea_res, "HALLMARK_P53_PATHWAY", "GSEA: Hallmark p53 Pathway", "Figures/08_gsea_p53_pathway.png")