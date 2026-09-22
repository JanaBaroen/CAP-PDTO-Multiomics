# ============================================================
# Install and load packages
# ============================================================

# Run this block once to install all required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "org.Hs.eg.db", "fgsea"))
install.packages(c("tidyverse", "ggplot2"))

# Load packages
library(tidyverse)
library(ggplot2)
library(fgsea)
library(org.Hs.eg.db)

# =============================================================================
# SETUP — same BASE_PATH as your main pipeline script
# =============================================================================

BASE_PATH <- "/your/path/to/RNAseq_output"   # CHANGE THIS to your output directory

# !! Make sure this .gmt file is in your BASE_PATH folder !!
# Download from https://www.gsea-msigdb.org/gsea/msigdb (free account needed)
#   Hallmark:  h.all.vX.X.Hs.symbols.gmt          (H collection)

GMT_HALLMARK <- file.path(BASE_PATH, "h.all.vX.X.Hs.symbols.gmt")             # H collection

# =============================================================================
# NOTE: This script assumes you have already run the main pipeline and have
# res.treatment available in your R environment. If not, load it from the TSV:
# =============================================================================

# Option A: if res.treatment is already in your environment, skip this block
# Option B: load from saved TSV if starting fresh
res.treatment.df <- read_tsv(file.path(BASE_PATH, "de_results_treatment.tsv"))

# =============================================================================
# Build GSEA ranking vector
# log2FC * -log10(pvalue) — captures both effect size and significance
# =============================================================================

res.for.gsea <- res.treatment.df %>%
  filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
  mutate(gene_clean = sub("\\..*", "", gene))  # strip Ensembl version numbers

# Map Ensembl IDs to gene symbols
id.map.gsea <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = res.for.gsea$gene_clean,
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)

res.for.gsea <- res.for.gsea %>%
  left_join(id.map.gsea, by = c("gene_clean" = "ENSEMBL")) %>%
  filter(!is.na(SYMBOL), SYMBOL != "")

ranking <- res.for.gsea$log2FoldChange * -log10(res.for.gsea$pvalue)
names(ranking) <- res.for.gsea$SYMBOL
ranking <- sort(ranking, decreasing = TRUE)

cat("Genes in ranking:", length(ranking), "\n")
cat("Top 5:\n"); print(head(ranking, 5))

# =============================================================================
# Load gene sets
# =============================================================================

pathways.hallmark <- gmtPathways(GMT_HALLMARK)
cat("Hallmark gene sets loaded:", length(pathways.hallmark), "\n")

# =============================================================================
# Run GSEA - Hallmark collection
# =============================================================================

set.seed(42)
res.fgsea.hallmark <- fgsea(pathways.hallmark, ranking, minSize = 15, maxSize = 500)
cat("Hallmark GSEA done.\n")

# =============================================================================
# Helper function to make and save GSEA barplot
# =============================================================================

make_gsea_plot <- function(res.fgsea, collection_name, filename_prefix,
                           strip_prefix = NULL, pval_threshold = 0.10,
                           top_n = NULL, width = 8, height = 6) {
  
  df <- res.fgsea %>%
    filter(pval < pval_threshold) %>%
    arrange(NES)
  
  if (!is.null(top_n) && nrow(df) > top_n) {
    df <- df %>% slice_max(abs(NES), n = top_n)
  }
  
  if (nrow(df) == 0) {
    cat("No pathways below pval <", pval_threshold, "for", collection_name, "\n")
    return(NULL)
  }
  
  df <- df %>%
    arrange(NES) %>%
    mutate(
      pathway     = if (!is.null(strip_prefix))
        str_replace_all(pathway, strip_prefix, "") else pathway,
      pathway     = str_replace_all(pathway, "_", " ") %>% str_to_title(),
      direction   = ifelse(NES > 0, "Up in Treated", "Up in Untreated"),
      significant = ifelse(padj < 0.05, "padj < 0.05", "padj >= 0.05")
    )
  
  p <- ggplot(df, aes(x = NES, y = reorder(pathway, NES),
                      fill = direction, alpha = significant)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Up in Treated"   = "firebrick",
                                 "Up in Untreated" = "steelblue")) +
    scale_alpha_manual(values = c("padj < 0.05" = 1.0,
                                  "padj >= 0.05" = 0.4)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(
      title    = paste("GSEA - Treatment effect (", collection_name, ")"),
      subtitle = "Solid = padj < 0.05; transparent = pval < 0.10 (uncorrected)",
      x        = "Normalised Enrichment Score (NES)",
      y        = NULL, fill = NULL, alpha = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  
  print(p)
  ggsave(file.path(BASE_PATH, paste0(filename_prefix, ".pdf")),
         p, width = width, height = height)
  ggsave(file.path(BASE_PATH, paste0(filename_prefix, ".png")),
         p, width = width, height = height, dpi = 300)
  
  cat(collection_name, "barplot saved.\n")
  return(p)
}

# =============================================================================
# Hallmark — results and plot
# =============================================================================

cat("\n--- HALLMARK ---\n")
cat("Significant pathways (padj < 0.05):\n")
res.fgsea.hallmark %>%
  filter(padj < 0.05) %>% arrange(padj) %>%
  dplyr::select(pathway, NES, pval, padj, size) %>%
  as.data.frame() %>% print()

make_gsea_plot(res.fgsea.hallmark, "Hallmark", "gsea_barplot_hallmark",
               strip_prefix = "HALLMARK_", width = 8, height = 6)

res.fgsea.hallmark %>% arrange(pval) %>%
  write_tsv(file.path(BASE_PATH, "gsea_hallmark_full.tsv"))
res.fgsea.hallmark %>% filter(pval < 0.10) %>% arrange(pval) %>%
  write_tsv(file.path(BASE_PATH, "gsea_hallmark_trending.tsv"))

cat("\nGSEA analysis complete! Output files saved to:", BASE_PATH, "\n")
cat("Note: res.fgsea.hallmark and pathways.hallmark must remain in the R session\n")
cat("      if you plan to run rnaseq_analysis_pipeline_heatmap.R next.\n")
cat("Next step: run rnaseq_analysis_pipeline_heatmap.R to generate manuscript heatmaps.\n")
cat("Next step: run rnaseq_treatment_per_cancer_type.R for per-cancer-type analysis.\n")
