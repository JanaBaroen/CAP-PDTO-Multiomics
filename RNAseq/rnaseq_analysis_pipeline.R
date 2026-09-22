# ============================================================
# Install and load packages
# ============================================================

# Run this block once to install all required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "org.Hs.eg.db"))
install.packages(c("tidyverse", "ggplot2", "fgsea"))

# Load packages
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(fgsea)
library(org.Hs.eg.db)

# =============================================================================
# SETUP — adjust these variable names to match your pdat.mat column names
# =============================================================================

# !! SET YOUR BASE PATH HERE — only need to change this one line !!
BASE_PATH <- "/your/path/to/RNAseq_output"   # CHANGE THIS to your output directory

# Column names from your pdat.mat
COL_PATIENT    <- "Patient"        # Column name for patient ID in pdat.mat
COL_CANCER     <- "Type"           # Cancer type (HNSCC / PDAC)
COL_TREATMENT  <- "Condition"      # Treatment condition

# Level labels as they appear in your data
LEVEL_TREATED   <- "Treated"       # Treated level in Condition column
LEVEL_UNTREATED <- "Untreated"     # Reference level in Condition column
LEVEL_HNSCC     <- "HNSCC"         # Cancer type 1 (reference)
LEVEL_PDAC      <- "PDAC"          # Cancer type 2

# =============================================================================
# 0. Load data
# =============================================================================

exp.mat  <- readRDS(file.path(BASE_PATH, "exp.mat.rds"))   # genes x samples count matrix
pdat.mat <- readRDS(file.path(BASE_PATH, "pdat.mat.rds"))  # samples x variables metadata data.frame

dim(exp.mat)   # rows = genes, cols = samples
dim(pdat.mat)  # rows = samples

# Critical: sample order must match between expression matrix and metadata
all(colnames(exp.mat) == rownames(pdat.mat))  # Must be TRUE before continuing!

# =============================================================================
# 1. Factor setup
#    Set reference levels so contrasts are interpretable:
#      - Treatment: Untreated is reference  → positive LFC = up in Treated
#      - Cancer type: HNSCC is reference    → positive LFC = up in PDAC
#      - Patient: blocking factor (no reference needed for interpretation)
# =============================================================================

pdat.mat[[COL_TREATMENT]] <- factor(pdat.mat[[COL_TREATMENT]],
                                    levels = c(LEVEL_UNTREATED, LEVEL_TREATED))

pdat.mat[[COL_CANCER]]    <- factor(pdat.mat[[COL_CANCER]],
                                    levels = c(LEVEL_HNSCC, LEVEL_PDAC))

pdat.mat[[COL_PATIENT]]   <- factor(pdat.mat[[COL_PATIENT]])

str(pdat.mat)  # Confirm all key columns are factors

# =============================================================================
# 2. DESeq2 object + VST normalisation (for PCA / visualisation)
#
#    Design: ~ patient + treatment
#      - patient:   blocks for within-patient pairing (absorbs inter-patient
#                   variation AND cancer type, since patient and cancer type
#                   are perfectly confounded — each patient has only one type)
#      - treatment: primary variable of interest
#
#    NOTE: cancer_type cannot be included separately in this model because
#    it is fully determined by patient (all samples from a patient are the
#    same cancer type). To test cancer type effect, run a separate model
#    using only one sample per patient (e.g. untreated only).
# =============================================================================

design_formula <- as.formula(
  paste("~", COL_PATIENT, "+", COL_TREATMENT)
)

dds <- DESeqDataSetFromMatrix(
  countData = exp.mat,
  colData   = pdat.mat,
  design    = design_formula
)

# =============================================================================
# 3. COLLAPSE TECHNICAL REPLICATES
# Each patient has 3 technical replicates per condition → sum into one sample
# This reduces 36 samples to 12 (6 patients × 2 conditions)
# =============================================================================

# Create a grouping variable combining Patient and Condition
pdat.mat$group <- paste(pdat.mat[[COL_PATIENT]], pdat.mat[[COL_TREATMENT]], sep = "_")

# Collapse replicates using DESeq2's built-in function
# This sums counts across technical replicates
dds <- collapseReplicates(
  dds,
  groupby    = pdat.mat$group,
  renameCols = TRUE
)

cat("Samples after collapsing:", ncol(dds), "\n")  # Should be 12

# Extract collapsed metadata — use this instead of pdat.mat for all plots
pdat.collapsed <- as.data.frame(colData(dds)) %>%
  dplyr::select(all_of(c(COL_PATIENT, COL_CANCER, COL_TREATMENT))) %>%
  mutate(
    across(all_of(COL_TREATMENT), ~ factor(., levels = c(LEVEL_UNTREATED, LEVEL_TREATED))),
    across(all_of(COL_CANCER),    ~ factor(., levels = c(LEVEL_HNSCC, LEVEL_PDAC))),
    across(all_of(COL_PATIENT),   ~ factor(.))
  )

# =============================================================================
# 4. VST normalisation (for PCA and heatmaps - NOT for DE testing)
# =============================================================================

dds.vst <- vst(dds, blind = FALSE)
vst.mat <- assay(dds.vst)

# =============================================================================
# 5. PCA
# =============================================================================

# Remove zero-variance genes (would break scaling)
gene_vars        <- apply(vst.mat, 1, var)
vst.mat.filtered <- vst.mat[gene_vars > 0, ]
cat("Genes retained for PCA:", nrow(vst.mat.filtered), "of", nrow(vst.mat), "\n")

pc            <- prcomp(t(vst.mat.filtered), center = TRUE, scale. = TRUE)
var_explained <- pc$sdev^2
prop_var      <- var_explained / sum(var_explained)
cum_var       <- cumsum(prop_var)

cat("PCs needed to explain 70% variance:", which(cum_var >= 0.70)[1], "\n")
cat("PC1:", round(prop_var[1]*100, 1), "%  |  PC2:", round(prop_var[2]*100, 1), "%\n")

pdat.collapsed$PC1 <- pc$x[, 1]
pdat.collapsed$PC2 <- pc$x[, 2]

pca_plot <- pdat.collapsed %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_line(aes(group = .data[[COL_PATIENT]]),
            color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = .data[[COL_TREATMENT]],
                 shape = .data[[COL_CANCER]]),
             size = 4, alpha = 0.9) +
  scale_color_manual(values = c("steelblue", "firebrick"), name = "Treatment") +
  scale_shape_manual(values = c(16, 17), name = "Cancer type") +
  labs(
    title = "PCA - VST normalised expression",
    x = paste0("PC1 (", round(prop_var[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(prop_var[2]*100, 1), "%)")
  ) +
  theme_bw(base_size = 13)

print(pca_plot)
ggsave(file.path(BASE_PATH, "pca_plot.pdf"), pca_plot, width = 7, height = 6)
ggsave(file.path(BASE_PATH, "pca_plot.png"), pca_plot, width = 7, height = 6, dpi = 300)

# =============================================================================
# 6. Differential expression - DESeq2
# =============================================================================

dds <- DESeq(dds)
resultsNames(dds)

res.treatment <- results(
  dds,
  contrast = c(COL_TREATMENT, LEVEL_TREATED, LEVEL_UNTREATED),
  alpha    = 0.05
)

summary(res.treatment)

# Volcano plot — all significant genes (padj < 0.05) coloured
volcano_plot <- res.treatment %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(size = 2, alpha = 0.6,
             aes(color = case_when(
               padj < 0.05 & log2FoldChange > 0 ~ "Up in Treated",
               padj < 0.05 & log2FoldChange < 0 ~ "Up in Untreated",
               TRUE                              ~ "NS"
             ))) +
  scale_color_manual(
    name   = NULL,
    values = c("Up in Treated" = "firebrick", "Up in Untreated" = "steelblue", "NS" = "grey60")
  ) +
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey40") +
  labs(title = "Treatment effect: Treated vs Untreated",
       subtitle = "Adjusted for patient pairing — coloured = padj < 0.05") +
  theme_bw(base_size = 13)

print(volcano_plot)
ggsave(file.path(BASE_PATH, "volcano_treatment_all_significant_genes.pdf"),
       volcano_plot, width = 7, height = 6)
ggsave(file.path(BASE_PATH, "volcano_treatment_all_significant_genes.png"),
       volcano_plot, width = 7, height = 6, dpi = 300)

# Save full DE results
res.treatment %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  write_tsv(file.path(BASE_PATH, "de_results_treatment.tsv"))

# =============================================================================
# 7. Annotate DE results using org.Hs.eg.db (no internet needed)
# =============================================================================

res.treatment.df <- res.treatment %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj) %>%
  mutate(gene_clean = sub("\\..*", "", gene))  # strip Ensembl version numbers

# Map Ensembl IDs to gene symbols, names, Entrez IDs
annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = res.treatment.df$gene_clean,
  columns = c("ENSEMBL", "SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cat("Genes successfully annotated:", nrow(annotations), "\n")

res.annotated <- res.treatment.df %>%
  left_join(annotations, by = c("gene_clean" = "ENSEMBL")) %>%
  dplyr::select(gene, SYMBOL, GENENAME, ENTREZID,
                log2FoldChange, pvalue, padj, baseMean, lfcSE) %>%
  arrange(padj)

# --- Significant genes (strict): padj < 0.05 ---------------

sig.genes <- res.treatment.df %>%
  filter(padj < 0.05)    # all 58 significant genes, no LFC filter

sig.annotated <- res.annotated %>%
  filter(padj < 0.05)    # same for annotated version

cat("Total significant DE genes (padj<0.05):", nrow(sig.annotated), "\n")
cat("  Up in Treated:  ", sum(sig.annotated$log2FoldChange > 0), "\n")
cat("  Up in Untreated:", sum(sig.annotated$log2FoldChange < 0), "\n")

sig.annotated %>%
  dplyr::select(SYMBOL, GENENAME, log2FoldChange, padj) %>%
  print()

# Save annotated results
res.annotated %>%
  write_tsv(file.path(BASE_PATH, "de_results_treatment_annotated.tsv"))

sig.annotated %>%
  write_tsv(file.path(BASE_PATH, "de_results_significant_annotated.tsv"))

cat("Annotated results saved.\n")

# =============================================================================
# 8. Heatmaps
# =============================================================================

# Shared column annotation using collapsed metadata
col.annotation <- data.frame(
  "Cancer type" = pdat.collapsed[[COL_CANCER]],
  "Treatment"   = pdat.collapsed[[COL_TREATMENT]],
  row.names     = rownames(pdat.collapsed),
  check.names   = FALSE
)

annotation.colors <- list(
  "Cancer type" = setNames(c("#7B68EE", "#FFA07A"), c(LEVEL_HNSCC, LEVEL_PDAC)),
  "Treatment"   = setNames(c("steelblue", "firebrick"),
                           c(LEVEL_UNTREATED, LEVEL_TREATED))
)

# Order columns: cancer type then treatment
col.order <- pdat.collapsed %>%
  rownames_to_column("sample") %>%
  arrange(.data[[COL_CANCER]], .data[[COL_TREATMENT]]) %>%
  pull(sample)

col.annotation.ordered <- col.annotation[col.order, , drop = FALSE]

# --- 8a. Heatmap: significant genes (strict) ---------------------------------

# Manual overrides for genes missing a symbol in org.Hs.eg.db,
# resolved via biomaRt / Ensembl / HGNC lookup
manual_symbols <- c(
  "ENSG00000203804" = "ADAMTSL4-AS1",
  "ENSG00000254708" = "MMADHCP2",
  "ENSG00000262074" = "SNORD3B-2",
  "ENSG00000265185" = "SNORD3B-1"
  # ENSG00000270022 = "RNA, U12 small nuclear" -- no unique HGNC symbol confirmed; left as ENSG
  # remaining IDs are genuine "novel transcript" loci with no assigned symbol anywhere -- kept as ENSG
)

heatmap.mat <- vst.mat[sig.genes$gene, ]

rownames(heatmap.mat) <- sig.annotated %>%
  mutate(
    gene_clean = sub("\\..*", "", gene),                  # strip version suffix for matching/display
    label = case_when(
      !is.na(SYMBOL)                          ~ SYMBOL,
      gene_clean %in% names(manual_symbols)   ~ manual_symbols[gene_clean],
      TRUE                                    ~ gene_clean  # fall back to version-stripped ENSG ID
    )
  ) %>%
  pull(label)

heatmap.mat.scaled <- t(scale(t(heatmap.mat)))[, col.order]

cat("\nMain pipeline complete! Output files saved to:", BASE_PATH, "\n")
cat("Next steps:\n")
cat("  1. Run rnaseq_analysis_pipeline_GSEA.R to generate GSEA results\n")
cat("  2. Run rnaseq_analysis_pipeline_heatmap.R to generate manuscript heatmaps\n")
cat("  3. Run rnaseq_treatment_per_cancer_type.R for per-cancer-type analysis\n")
