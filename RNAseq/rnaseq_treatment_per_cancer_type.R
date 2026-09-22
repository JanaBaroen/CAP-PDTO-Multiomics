# ============================================================
# Treatment effect analysis per cancer type (HNSCC and PDAC)
#
# Runs DESeq2 and Hallmark GSEA separately within each
# cancer type, producing annotated DE results and GSEA output files.
#
# Output files are named with the cancer type suffix
# e.g. gsea_treatment_HNSCC_hallmark_full.tsv
# ============================================================

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
library(DESeq2)
library(ggplot2)
library(fgsea)
library(org.Hs.eg.db)

# =============================================================================
# SETUP
# =============================================================================

BASE_PATH <- "/your/path/to/RNAseq_output"   # CHANGE THIS to your output directory

COL_PATIENT   <- "Patient"
COL_CANCER    <- "Type"
COL_TREATMENT <- "Condition"

LEVEL_TREATED   <- "Treated"
LEVEL_UNTREATED <- "Untreated"
LEVEL_HNSCC     <- "HNSCC"
LEVEL_PDAC      <- "PDAC"

# Download from https://www.gsea-msigdb.org/gsea/msigdb (free account needed)
GMT_HALLMARK <- file.path(BASE_PATH, "h.all.vX.X.Hs.symbols.gmt")

# =============================================================================
# 0. Load data
# =============================================================================

exp.mat  <- readRDS(file.path(BASE_PATH, "exp.mat.rds"))
pdat.mat <- readRDS(file.path(BASE_PATH, "pdat.mat.rds"))

pdat.mat[[COL_TREATMENT]] <- factor(pdat.mat[[COL_TREATMENT]],
                                    levels = c(LEVEL_UNTREATED, LEVEL_TREATED))
pdat.mat[[COL_CANCER]]    <- factor(pdat.mat[[COL_CANCER]],
                                    levels = c(LEVEL_HNSCC, LEVEL_PDAC))
pdat.mat[[COL_PATIENT]]   <- factor(pdat.mat[[COL_PATIENT]])

# =============================================================================
# Helper functions
# =============================================================================

# Annotate DE results using org.Hs.eg.db
annotate_results <- function(res.df) {
  res.df <- res.df %>%
    mutate(gene_clean = sub("\\..*", "", gene))

  annotations <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = res.df$gene_clean,
    columns = c("ENSEMBL", "SYMBOL", "GENENAME"),
    keytype = "ENSEMBL"
  ) %>%
    filter(!is.na(SYMBOL)) %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  res.df %>%
    left_join(annotations, by = c("gene_clean" = "ENSEMBL")) %>%
    dplyr::select(gene, SYMBOL, GENENAME,
                  log2FoldChange, pvalue, padj, baseMean, lfcSE) %>%
    arrange(padj)
}

# Make GSEA ranking vector
make_ranking <- function(res) {
  res.df <- res %>%
    data.frame() %>%
    rownames_to_column("gene") %>%
    filter(!is.na(pvalue), !is.na(log2FoldChange)) %>%
    mutate(gene_clean = sub("\\..*", "", gene))

  id.map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = res.df$gene_clean,
    columns = c("ENSEMBL", "SYMBOL"),
    keytype = "ENSEMBL"
  ) %>%
    filter(!is.na(SYMBOL), SYMBOL != "") %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  res.df <- res.df %>%
    left_join(id.map, by = c("gene_clean" = "ENSEMBL")) %>%
    filter(!is.na(SYMBOL), SYMBOL != "")

  ranking <- res.df$log2FoldChange * -log10(res.df$pvalue)
  names(ranking) <- res.df$SYMBOL
  sort(ranking, decreasing = TRUE)
}

# Make and save GSEA barplot
make_gsea_plot <- function(res.fgsea, collection_name, filename_prefix,
                           cancer_type, strip_prefix = NULL,
                           pval_threshold = 0.10, top_n = 20,
                           width = 8, height = 6) {
  df <- res.fgsea %>%
    filter(pval < pval_threshold)

  if (nrow(df) == 0) {
    cat("No pathways below pval <", pval_threshold, "for", collection_name, "\n")
    return(NULL)
  }

  if (!is.null(top_n) && nrow(df) > top_n) {
    df <- df %>% slice_max(abs(NES), n = top_n)
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
    scale_alpha_manual(values = c("padj < 0.05" = 1.0, "padj >= 0.05" = 0.4)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    labs(
      title    = paste("GSEA - Treatment effect in", cancer_type,
                       "(", collection_name, ")"),
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
  cat(collection_name, cancer_type, "GSEA plot saved.\n")
  return(p)
}

# Load gene sets once
pathways.hallmark <- gmtPathways(GMT_HALLMARK)

# =============================================================================
# Run analysis for each cancer type
# =============================================================================

for (cancer in c(LEVEL_HNSCC, LEVEL_PDAC)) {

  cat("\n\n=================================================\n")
  cat("Analysing treatment effect in:", cancer, "\n")
  cat("=================================================\n")

  # --- Subset to this cancer type -------------------------------------------

  samples.cancer <- rownames(pdat.mat)[pdat.mat[[COL_CANCER]] == cancer]
  exp.sub        <- exp.mat[, samples.cancer]
  pdat.sub       <- pdat.mat[samples.cancer, ]

  cat("Samples:", ncol(exp.sub), "\n")

  # --- DESeq2 object --------------------------------------------------------

  dds.sub <- DESeqDataSetFromMatrix(
    countData = exp.sub,
    colData   = pdat.sub,
    design    = as.formula(paste("~", COL_PATIENT, "+", COL_TREATMENT))
  )

  # Collapse technical replicates
  pdat.sub$group <- paste(pdat.sub[[COL_PATIENT]], pdat.sub[[COL_TREATMENT]], sep = "_")
  dds.sub <- collapseReplicates(dds.sub, groupby = pdat.sub$group, renameCols = TRUE)
  cat("Samples after collapsing:", ncol(dds.sub), "\n")

  # --- DESeq2 DE ------------------------------------------------------------

  dds.sub <- DESeq(dds.sub)

  res.sub <- results(
    dds.sub,
    contrast = c(COL_TREATMENT, LEVEL_TREATED, LEVEL_UNTREATED),
    alpha    = 0.05
  )

  cat("\nDE summary for", cancer, ":\n")
  print(summary(res.sub))

  # --- Annotate -------------------------------------------------------------

  res.sub.df <- res.sub %>%
    data.frame() %>%
    rownames_to_column("gene") %>%
    arrange(padj)

  res.sub.annotated <- annotate_results(res.sub.df)

  sig.sub.annotated <- res.sub.annotated %>% filter(padj < 0.05)

  cat("Significant genes (padj<0.05) in", cancer, ":", nrow(sig.sub.annotated), "\n")
  cat("  Up in Treated:  ", sum(sig.sub.annotated$log2FoldChange > 0, na.rm = TRUE), "\n")
  cat("  Up in Untreated:", sum(sig.sub.annotated$log2FoldChange < 0, na.rm = TRUE), "\n")

  # Save results
  res.sub.annotated %>%
    write_tsv(file.path(BASE_PATH, paste0("de_results_treatment_", cancer, "_annotated.tsv")))
  sig.sub.annotated %>%
    write_tsv(file.path(BASE_PATH, paste0("de_results_treatment_", cancer, "_significant.tsv")))

  # --- GSEA -----------------------------------------------------------------

  ranking.sub <- make_ranking(res.sub)
  cat("Genes in ranking for", cancer, ":", length(ranking.sub), "\n")

  set.seed(42)
  fgsea.hallmark <- fgsea(pathways.hallmark, ranking.sub, minSize = 15, maxSize = 500)

  cat("\nGSEA significant pathways (padj<0.05) for", cancer, ":\n")
  cat("Hallmark:", sum(fgsea.hallmark$padj < 0.05, na.rm = TRUE), "\n")

  # Print top trending
  cat("\nHallmark trending (pval<0.10):\n")
  fgsea.hallmark %>% filter(pval < 0.10) %>% arrange(pval) %>%
    dplyr::select(pathway, NES, pval, padj, size) %>%
    as.data.frame() %>% print()

  # Make plots
  make_gsea_plot(fgsea.hallmark, "Hallmark",
                 paste0("gsea_treatment_", cancer, "_hallmark"),
                 cancer_type  = cancer,
                 strip_prefix = "HALLMARK_",
                 pval_threshold = 0.10, width = 8, height = 6)

  # Save GSEA results
  fgsea.hallmark %>% arrange(pval) %>%
    write_tsv(file.path(BASE_PATH, paste0("gsea_treatment_", cancer, "_hallmark_full.tsv")))
  
  fgsea.hallmark %>% filter(pval < 0.10) %>% arrange(pval) %>%
    write_tsv(file.path(BASE_PATH, paste0("gsea_treatment_", cancer, "_hallmark_trending.tsv")))
  
  cat("\nCompleted analysis for:", cancer, "\n")
}

cat("\nAll analyses complete!\n")
cat("Output files saved to:", BASE_PATH, "\n")
