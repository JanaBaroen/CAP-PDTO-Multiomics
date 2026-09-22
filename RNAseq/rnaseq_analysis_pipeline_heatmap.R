# ============================================================
# Manuscript heatmap figures — treatment effect
#
# Produces two heatmaps used in the manuscript:
#   1. heatmap_curated_main_figure     — curated gene selection with
#      pathway annotation per gene (main figure)
#   2. supplementary_heatmap_all_significant_DE_genes — all significant
#      DE genes (padj < 0.05) for the supplementary figures
#
# Prerequisites — run these scripts first in the same R session:
#   1. rnaseq_analysis_pipeline.R  (provides: vst.mat, res.annotated,
#      pdat.collapsed, col.annotation.ordered, annotation.colors,
#      col.order, manual_symbols)
#   2. rnaseq_analysis_pipeline_GSEA.R  (provides: res.fgsea.hallmark,
#      pathways.hallmark)
# ============================================================

# ============================================================
# Install and load packages
# ============================================================

# Run this block once to install all required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap"))
install.packages(c("tidyverse", "circlize"))

# Load packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# -----------------------------------------------------------------------
# 1. Define significance filter and base path
# -----------------------------------------------------------------------
# !! SET YOUR BASE PATH HERE — must match rnaseq_analysis_pipeline.R !!
BASE_PATH <- "/your/path/to/RNAseq_output"   # CHANGE THIS to your output directory

LFC_THRESHOLD <- 0   # set to 0 for padj<0.05 only; set to 1 for |log2FC|>1 additionally

sig.annotated <- res.annotated %>%
  filter(padj < 0.05, abs(log2FoldChange) > LFC_THRESHOLD)

filter_label <- if (LFC_THRESHOLD > 0) {
  paste0("padj<0.05, |log2FC|>", LFC_THRESHOLD)
} else {
  "padj<0.05"
}

cat("Significant genes (", filter_label, "):", nrow(sig.annotated), "\n")
cat("  Up in Treated:  ", sum(sig.annotated$log2FoldChange > 0), "\n")
cat("  Down in Treated:", sum(sig.annotated$log2FoldChange < 0), "\n")

# Resolve display labels the same way the original script does
sig.annotated <- sig.annotated %>%
  mutate(
    gene_clean = sub("\\..*", "", gene),
    label = case_when(
      !is.na(SYMBOL)                          ~ SYMBOL,
      gene_clean %in% names(manual_symbols)   ~ manual_symbols[gene_clean],
      TRUE                                    ~ gene_clean
    ),
    direction = ifelse(log2FoldChange > 0, "Induced in treated", "Reduced in treated")
  )

# -----------------------------------------------------------------------
# 2. Build Hallmark pathway-membership lookup
# -----------------------------------------------------------------------
# Maps each significant gene to all trending Hallmark pathways it belongs to.

if (!exists("res.fgsea.hallmark") || !exists("pathways.hallmark")) {
  stop("res.fgsea.hallmark and/or pathways.hallmark not found in this session. ",
       "Run rnaseq_analysis_pipeline_GSEA.R through the Hallmark fgsea() call ",
       "AND the gmtPathways(GMT_HALLMARK) call before this script.")
}

res.fgsea.hallmark <- as_tibble(res.fgsea.hallmark)  # in case it's still a data.table

# Display labels for Hallmark pathways — shortened for figure legibility.
# Must match the exact HALLMARK_* pathway names in your .gmt file as keys.
# If a pathway name doesn't match, the lookup will produce NA and a warning
# will be shown below.
hallmark_display_names <- c(
  "HALLMARK_UV_RESPONSE_UP"               = "UV response up",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"      = "TNFa via NFkB",
  "HALLMARK_HYPOXIA"                      = "Hypoxia",
  "HALLMARK_MYOGENESIS"                   = "Myogenesis",
  "HALLMARK_P53_PATHWAY"                  = "p53 pathway",
  "HALLMARK_UV_RESPONSE_DN"               = "UV response dn",
  "HALLMARK_BILE_ACID_METABOLISM"         = "Bile acid metab.",
  "HALLMARK_ANDROGEN_RESPONSE"            = "Androgen resp.",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY"      = "Estrogen early",
  "HALLMARK_E2F_TARGETS"                  = "E2F targets",
  "HALLMARK_PROTEIN_SECRETION"            = "Protein secr.",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"    = "IFN-\u03b3 response",  # IFN-γ response
  "HALLMARK_G2M_CHECKPOINT"               = "G2M checkpoint",
  "HALLMARK_MITOTIC_SPINDLE"              = "Mitotic spindle",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"    = "IFN-\u03b1 response"   # IFN-α response
)

trending_pathways <- res.fgsea.hallmark %>%
  dplyr::filter(pval < 0.10) %>%               # trending cutoff
  dplyr::arrange(pval) %>%
  dplyr::pull(pathway)

cat("Trending Hallmark pathways (pval < 0.10):", length(trending_pathways), "\n")

unmapped <- setdiff(trending_pathways, names(hallmark_display_names))
if (length(unmapped) > 0) {
  warning("These trending pathways have no exact display label in ",
          "hallmark_display_names and will show their raw HALLMARK_ name instead: ",
          paste(unmapped, collapse = ", "),
          ". Add them to hallmark_display_names above with the wording you want.")
}

gene_pathway_map <- tibble(pathway = trending_pathways, rank = seq_along(trending_pathways)) %>%
  mutate(member_genes = pathways.hallmark[trending_pathways]) %>%
  tidyr::unnest(member_genes) %>%
  dplyr::rename(SYMBOL = member_genes) %>%
  mutate(pathway_short = dplyr::coalesce(hallmark_display_names[pathway], pathway)) %>%
  arrange(SYMBOL, rank) %>%
  group_by(SYMBOL) %>%
  summarise(gsea_pathways = list(unique(pathway_short)), .groups = "drop")

# Tag each significant gene with its GSEA-confirmed pathway memberships
sig.annotated <- sig.annotated %>%
  left_join(gene_pathway_map, by = c("label" = "SYMBOL")) %>%
  mutate(gsea_pathways = ifelse(is.na(gsea_pathways), list(character(0)), gsea_pathways))

sig.annotated$is_pathway_member <- lengths(sig.annotated$gsea_pathways) > 0

cat("Significant genes related to a trending Hallmark pathway:",
    sum(sig.annotated$is_pathway_member), "of", nrow(sig.annotated), "\n")
cat("  Up + pathway member:  ",
    sum(sig.annotated$log2FoldChange > 0 & sig.annotated$is_pathway_member), "\n")
cat("  Down + pathway member:",
    sum(sig.annotated$log2FoldChange < 0 & sig.annotated$is_pathway_member), "\n")

# -----------------------------------------------------------------------
# 3. Curate genes for the main-text figure
# -----------------------------------------------------------------------
# Manually selected genes of biological interest from the significant DE gene list.
# Edit curated_up and curated_down to change the gene selection.

curated_up   <- c("DUSP1", "VEGFA", "HBEGF", "SOCS1", "ARID5A",
                  "MYADM", "LSMEM1", "NEU1", "PTCH2")
curated_down <- c("H1-3", "H2AC4", "H3C1", "H4C9", "H2BC3",
                  "SELE", "PBX1", "ENC1", "TET3", "ARFGEF2")

cat("\nCurated gene selection:\n")
cat("  Up genes requested:  ", length(curated_up), "\n")
cat("  Down genes requested:", length(curated_down), "\n")

curated_genes <- sig.annotated %>%
  filter(label %in% c(curated_up, curated_down))

requested <- c(curated_up, curated_down)
missing   <- setdiff(requested, curated_genes$label)
cat("\nCurated genes matched:", nrow(curated_genes), "of", length(requested), "requested\n")
if (length(missing) > 0) {
  cat("Not found in sig.annotated (check spelling or significance):",
      paste(missing, collapse = ", "), "\n")
}
cat("Of these, GSEA-confirmed pathway members:", sum(curated_genes$is_pathway_member),
    "of", nrow(curated_genes), "\n")

# -----------------------------------------------------------------------
# 4. Hand-curated pathway display per gene
# -----------------------------------------------------------------------
# Confidence notes:
#   HIGH:     DUSP1, ARID5A, SELE (TNFa/NFkB); SOCS1 (IFN response);
#             VEGFA (Hypoxia); ARFGEF2 (Protein secretion)
#   MODERATE: DUSP1, HBEGF (TNFa/NFkB)
#   LOW:      PBX1 (Androgen response — documented mainly in prostate cancer)
#   NOTE:     Histone genes assigned to E2F/G2M/Mitotic Spindle reflect
#             general biology; these specific symbols may not be formal
#             members of those Hallmark gene sets

pathway_display <- list(
  "DUSP1"   = c("TNFa via NFkB", "UV response up"),
  "HBEGF"   = c("TNFa via NFkB"),
  "SOCS1"   = c("IFN-\u03b3 response", "IFN-\u03b1 response"),
  "VEGFA"   = c("Hypoxia"),
  "ARID5A"  = c("TNFa via NFkB"),
  "SELE"    = c("TNFa via NFkB"),
  "ARFGEF2" = c("Protein secr."),
  "PBX1"    = c("Androgen resp."),
  "H1-3"    = c("E2F targets", "G2M checkpoint"),
  "H2AC4"   = c("E2F targets", "G2M checkpoint"),
  "H3C1"    = c("E2F targets", "G2M checkpoint"),
  "H4C9"    = c("E2F targets", "Mitotic spindle"),
  "H2BC3"   = c("E2F targets", "G2M checkpoint")
  # Add/edit/remove entries freely. Any curated gene not listed here (or
  # listed with an empty vector) gets zero pathways and is dropped in
  # Section 4b below. Use the exact label strings from hallmark_display_names
  # (Section 2) so pathway_colors matches them up correctly.
)

curated_genes <- curated_genes %>%
  rowwise() %>%
  mutate(
    combined_pathways = list(if (label %in% names(pathway_display))
      pathway_display[[label]] else character(0))
  ) %>%
  ungroup()

curated_genes$n_pathways <- lengths(curated_genes$combined_pathways)
cat("\nHand-curated pathway assignment applied.\n")
cat("Genes with at least one displayed pathway:",
    sum(curated_genes$n_pathways > 0), "of", nrow(curated_genes), "\n")
cat("Genes with 2 displayed pathways:", sum(curated_genes$n_pathways == 2), "\n")
cat("Genes with 3+ displayed pathways:", sum(curated_genes$n_pathways >= 3), "\n")

# -----------------------------------------------------------------------
# 4b. Keep ONLY genes that have at least one related pathway
# -----------------------------------------------------------------------
# Drops curated genes with no pathway assignment in pathway_display above.
# To include a gene, add it to pathway_display with its relevant pathway(s).

dropped_genes <- curated_genes$label[curated_genes$n_pathways == 0]
curated_genes <- curated_genes %>% filter(n_pathways > 0)

cat("\nDropped (no related pathway):", paste(dropped_genes, collapse = ", "), "\n")
cat("Remaining curated genes:", nrow(curated_genes), "\n")
cat("  Induced in treated:", sum(curated_genes$direction == "Induced in treated"), "\n")
cat("  Reduced in treated:", sum(curated_genes$direction == "Reduced in treated"), "\n")

# -----------------------------------------------------------------------
# 5. Recolor the annotations: Cancer type, Treatment
# -----------------------------------------------------------------------
# Cancer type colours: HNSCC = orange, PDAC = purple
# Must match the level names used in annotation.colors from the main pipeline

if (!all(c("HNSCC", "PDAC") %in% names(annotation.colors$`Cancer type`))) {
  stop("Expected names 'HNSCC' and 'PDAC' in annotation.colors$`Cancer type`, found: ",
       paste(names(annotation.colors$`Cancer type`), collapse = ", "),
       ". Update the names below to match.")
}
annotation.colors$`Cancer type`["HNSCC"] <- "#FFA07A"   # orange
annotation.colors$`Cancer type`["PDAC"]  <- "#7B68EE"   # purple

# Treatment colours: teal = Untreated, mustard/gold = Treated
annotation.colors$Treatment <- setNames(
  c("#2A9D8F", "#E9C46A"),      # teal = Untreated, mustard/gold = Treated
  c(LEVEL_UNTREATED, LEVEL_TREATED)
)

# -----------------------------------------------------------------------
# 6. Build the multi-pathway display columns and colors
# -----------------------------------------------------------------------
# Each gene shows up to MAX_PATHWAY_SLOTS pathway columns.
# Increase MAX_PATHWAY_SLOTS if any gene has more pathways than the current limit.

MAX_PATHWAY_SLOTS <- 2   # increase if any gene has more than 2 pathways

pathway_slot_names <- c("Pathway 1", "Pathway 2")

overflow <- curated_genes$n_pathways > MAX_PATHWAY_SLOTS
if (any(overflow)) {
  warning("These genes relate to more than ", MAX_PATHWAY_SLOTS, " pathways — only the ",
          "first ", MAX_PATHWAY_SLOTS, " are shown: ",
          paste(curated_genes$label[overflow], collapse = ", "),
          ". Add more explicit Pathway N columns below if you want all of them displayed.")
}

# Each gene's pathways fill slots in order; unused slots are NA (shown as white).
pathway_wide <- curated_genes %>%
  rowwise() %>%
  mutate(
    `Pathway 1` = if (length(combined_pathways) >= 1) combined_pathways[[1]] else NA_character_,
    `Pathway 2` = if (length(combined_pathways) >= 2) combined_pathways[[2]] else NA_character_
  ) %>%
  ungroup() %>%
  dplyr::select(label, `Pathway 1`, `Pathway 2`)

all_pathway_labels <- sort(unique(unlist(curated_genes$combined_pathways)))
cat("\nDistinct related-pathway labels across curated genes:", length(all_pathway_labels), "\n")

# Colour rule: avoid white, blue, and red; avoid exact hues used for
# Cancer type (#FFA07A, #7B68EE), Treatment (#2A9D8F, #E9C46A),
# and Z-score (steelblue - white - firebrick)

pathway_palette <- c(
  "#6A1B9A",  # deep purple      — distinct from Cancer type's lighter lavender #7B68EE
  "#2E7D32",  # forest green
  "#BF5B04",  # burnt orange/rust — distinct from Cancer type's salmon #FFA07A
  "#C2185B",  # magenta/raspberry
  "#827717",  # dark olive/mustard — distinct from Treatment's gold #E9C46A
  "#00695C",  # dark teal        — distinct from Treatment's teal #2A9D8F
  "#6D4C41",  # brown
  "#9E9E9E",  # gray
  "#9ACD32",  # yellow-green/lime
  "#880E4F"   # dark wine/magenta — distinct from the raspberry above (darker, more toward wine)
)

# Completeness check: stops if more pathway labels than palette colours
if (length(all_pathway_labels) > length(pathway_palette)) {
  stop("More distinct pathway labels (", length(all_pathway_labels), ") than colors in ",
       "pathway_palette (", length(pathway_palette), "): ",
       paste(all_pathway_labels, collapse = ", "),
       ". Add more hex codes to pathway_palette above before proceeding — ",
       "continuing would silently leave some labels uncolored (rendered as white).")
}
pathway_colors <- setNames(pathway_palette[seq_along(all_pathway_labels)], all_pathway_labels)
stopifnot(all(all_pathway_labels %in% names(pathway_colors)))  # belt-and-suspenders

# -----------------------------------------------------------------------
# 7. Build the curated matrix, split by direction, columns grouped
# -----------------------------------------------------------------------
# Rows ordered by significance, split by direction (induced/reduced).
# Gene names on the left, no clustering.

curated.mat <- vst.mat[curated_genes$gene, , drop = FALSE]
rownames(curated.mat) <- curated_genes$label
curated.mat.scaled <- t(scale(t(curated.mat)))[, col.order, drop = FALSE]

row_order_match <- match(rownames(curated.mat.scaled), curated_genes$label)
pathway_wide_ordered <- pathway_wide[match(rownames(curated.mat.scaled), pathway_wide$label), ]

row_split <- factor(
  curated_genes$direction[row_order_match],
  levels = c("Induced in treated", "Reduced in treated")
)

# Column grouping: Cancer type x Treatment, in the same order as col.order
col_group_labels <- paste(col.annotation.ordered[["Cancer type"]],
                          col.annotation.ordered[["Treatment"]])
col_split <- factor(col_group_labels, levels = unique(col_group_labels))

col_fun <- colorRamp2(seq(-2, 2, length.out = 11),
                      colorRampPalette(c("steelblue", "white", "firebrick"))(11))

ha_top <- HeatmapAnnotation(
  df  = col.annotation.ordered,
  col = annotation.colors,
  annotation_legend_param = list(
    Treatment     = list(title = "Treatment"),
    `Cancer type` = list(title = "Cancer type")
  ),
  show_annotation_name = FALSE
)

# Build pathway annotation columns. Legend is built manually to ensure all
# pathways appear regardless of which slot they occupy.
pathway_annotation_args <- setNames(
  lapply(pathway_slot_names, function(nm) pathway_wide_ordered[[nm]]),
  pathway_slot_names
)
pathway_col_list <- setNames(
  rep(list(pathway_colors), MAX_PATHWAY_SLOTS),
  pathway_slot_names
)

ha_right <- do.call(rowAnnotation, c(
  pathway_annotation_args,
  list(
    col = pathway_col_list,
    na_col = "white",    # empty pathway slots shown as white (intentional)
    show_legend = rep(FALSE, MAX_PATHWAY_SLOTS)   # suppressed - see explicit Legend() below
  )
))

pathway_legend <- Legend(
  labels    = names(pathway_colors),
  legend_gp = gpar(fill = unname(pathway_colors)),
  title     = "Likely related pathway(s)"
)

ht_curated <- Heatmap(
  curated.mat.scaled,
  name                 = "Z-score",
  col                  = col_fun,
  top_annotation       = ha_top,
  right_annotation     = ha_right,
  row_split            = row_split,
  row_title_gp         = gpar(fontsize = 11, fontface = "bold"),
  row_title_rot        = 90,   # vertical row title
  row_gap              = unit(3, "mm"),
  cluster_row_slices   = FALSE,
  cluster_rows         = FALSE,    # rows stay in significance-selection order
  row_names_side       = "left",   # gene names on the left
  column_split         = col_split,
  column_title         = NULL,     # cancer type and treatment shown by annotation bars above
  column_gap           = unit(2, "mm"),
  cluster_columns      = FALSE,
  show_column_names    = FALSE,
  show_row_names       = TRUE,
  row_names_gp         = gpar(fontsize = 10, fontface = "italic"),
  heatmap_legend_param = list(title = "Z-score")
)

for (ext in c("pdf", "png")) {
  if (ext == "pdf") {
    pdf(file.path(BASE_PATH, "heatmap_curated_main_figure.pdf"), width = 8.5, height = 7)
  } else {
    png(file.path(BASE_PATH, "heatmap_curated_main_figure.png"),
        width = 8.5, height = 7, units = "in", res = 300)
  }
  draw(ht_curated,
       heatmap_legend_side     = "bottom",
       annotation_legend_side  = "bottom",
       annotation_legend_list  = list(pathway_legend),
       merge_legend            = TRUE,
       column_title            = paste0("Curated DE genes - Treatment effect\n(",
                                        filter_label, ", Z-scored VST, n=",
                                        nrow(curated_genes), " of ", nrow(sig.annotated), ")"))
  dev.off()
}
cat("Curated main-figure heatmap saved.\n")

# -----------------------------------------------------------------------
# 8. Full significant-gene heatmap -> supplementary figure
# -----------------------------------------------------------------------
# Full heatmap of all significant DE genes for supplementary figure

full.mat <- vst.mat[sig.annotated$gene, , drop = FALSE]
rownames(full.mat) <- sig.annotated$label
full.mat.scaled <- t(scale(t(full.mat)))[, col.order, drop = FALSE]

ht_full <- Heatmap(
  full.mat.scaled,
  name                 = "Z-score",
  col                  = col_fun,
  top_annotation       = ha_top,
  show_column_names    = FALSE,
  show_row_names       = TRUE,
  cluster_columns      = FALSE,
  cluster_rows         = TRUE,
  row_names_gp         = gpar(fontsize = 7),
  heatmap_legend_param = list(title = "Z-score")
)

for (ext in c("pdf", "png")) {
  if (ext == "pdf") {
    pdf(file.path(BASE_PATH, "supplementary_heatmap_all_significant_DE_genes.pdf"),
        width = 10, height = 14)
  } else {
    png(file.path(BASE_PATH, "supplementary_heatmap_all_significant_DE_genes.png"),
        width = 10, height = 14, units = "in", res = 300)
  }
  draw(ht_full,
       heatmap_legend_side    = "bottom",
       annotation_legend_side = "bottom",
       merge_legend           = TRUE,
       column_title           = paste0("All significant DE genes - Treatment effect\n(",
                                       filter_label, ", Z-scored VST, n=",
                                       nrow(sig.annotated), ")"))
  dev.off()
}
cat("\nHeatmap figures complete! Output files saved to:", BASE_PATH, "\n")
cat("Next step: run rnaseq_treatment_per_cancer_type.R for per-cancer-type analysis.\n")

