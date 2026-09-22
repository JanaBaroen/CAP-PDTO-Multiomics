# =============================================================================
# Step 4 — Gene Set Enrichment Analysis (GSEA)
#
# Method:         fgsea (fast preranked GSEA)
# Ranking metric: limma moderated t-statistic (from step 3)
# Gene sets:      MSigDB Hallmarks, KEGG Medicus, GO Biological Process
#
# Input:  step3_DE_<contrast>.tsv files
# Output: step4_GSEA_Treat_PDAC.tsv   — primary output used in manuscript
#         step4_GSEA_Treat_HNSCC.tsv  — primary output used in manuscript
#         step4_GSEA_HNSCC_vs_PDAC.tsv
#         step4_GSEA_Interaction.tsv
#         step4_GSEA_summary.tsv
#
# BEFORE RUNNING: download .gmt files from https://www.gsea-msigdb.org/gsea/msigdb/
#   1. Hallmarks:    h.all.vX.X.Hs.symbols.gmt
#   2. KEGG Medicus: c2.cp.kegg_medicus.vX.X.Hs.symbols.gmt
#   3. GO BP:        c5.go.bp.vX.X.Hs.symbols.gmt  (optional but recommended)
# Save all .gmt files to the folder set in gsea_data_dir below.
# =============================================================================

# ---- Packages ----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_cran <- c("dplyr", "stringr", "readr", "tibble", "forcats")
pkgs_bioc <- c("fgsea")

new_cran <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
new_bioc <- pkgs_bioc[!sapply(pkgs_bioc, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) install.packages(new_cran)
if (length(new_bioc) > 0) BiocManager::install(new_bioc)

suppressPackageStartupMessages(lapply(c(pkgs_cran, pkgs_bioc), library,
                                       character.only = TRUE))

# Prevent namespace conflicts
select <- dplyr::select
filter <- dplyr::filter


# ==============================================================================
# SECTION A — Configuration  (EDIT THIS BLOCK)
# ==============================================================================

output_dir    <- "/path/to/output"          # <-- CHANGE THIS
gsea_data_dir <- "/path/to/gsea_genesets"  # <-- CHANGE THIS (folder with .gmt files)

dir.create(gsea_data_dir, showWarnings = FALSE, recursive = TRUE)


# ==============================================================================
# SECTION B — Load gene sets from .gmt files
# ==============================================================================

gmt_files <- list.files(gsea_data_dir, pattern = "\\.gmt$", full.names = TRUE)

if (length(gmt_files) == 0) {
  stop("No .gmt files found in: ", gsea_data_dir,
       "\nDownload from https://www.gsea-msigdb.org/gsea/msigdb/ and save to: ",
       gsea_data_dir)
}

message("Loading gene sets from: ", gsea_data_dir)
gene_sets <- list()
for (f in gmt_files) {
  gs <- gmtPathways(f)
  bn <- basename(f)
  nm <- if      (grepl("^h\\.",        bn, ignore.case = TRUE)) "Hallmarks"
        else if (grepl("kegg",         bn, ignore.case = TRUE)) "KEGG"
        else if (grepl("c5\\.go\\.bp", bn, ignore.case = TRUE)) "GO_BP"
        else bn %>% str_remove("\\.v[0-9]+.*$") %>% str_replace_all("\\.", "_")
  gene_sets[[nm]] <- gs
  message("  ", nm, ": ", length(gs), " gene sets loaded")
}

# Combine all collections into one list, prefixing names to avoid collisions
all_gene_sets <- list()
for (coll in names(gene_sets)) {
  gs_named  <- gene_sets[[coll]]
  new_names <- ifelse(
    str_detect(names(gs_named), "^HALLMARK_|^KEGG_|^GOBP_|^REACTOME_"),
    names(gs_named),
    paste0(coll, "_", names(gs_named))
  )
  names(gs_named) <- new_names
  all_gene_sets   <- c(all_gene_sets, gs_named)
}
message("Total gene sets: ", length(all_gene_sets))


# ==============================================================================
# SECTION C — Load DE results and build ranked protein lists
# ==============================================================================
# Ranking metric: moderated t-statistic from limma.
# Preferred over logFC because it accounts for measurement precision.

contrast_names <- c("Treat_PDAC", "Treat_HNSCC", "HNSCC_vs_PDAC", "Interaction")

de_results <- lapply(contrast_names, function(ct) {
  f <- file.path(output_dir, paste0("step3_DE_", ct, ".tsv"))
  if (!file.exists(f)) stop("File not found: ", f, "\nRun step3 first.")
  read_tsv(f, show_col_types = FALSE)
})
names(de_results) <- contrast_names

ranked_lists <- lapply(de_results, function(res) {
  ranks <- setNames(res$t_stat, res$SYMBOL)
  sort(ranks[!is.na(ranks)], decreasing = TRUE)
})

cat("\n=== Ranked list sizes ===\n"); print(sapply(ranked_lists, length))


# ==============================================================================
# SECTION D — Run fgsea for all contrasts
# ==============================================================================

set.seed(42)
gsea_results <- list()

for (ct in contrast_names) {
  message("\nRunning fgsea: ", ct, " (", length(ranked_lists[[ct]]),
          " proteins, ", length(all_gene_sets), " gene sets)")

  res <- fgsea(
    pathways    = all_gene_sets,
    stats       = ranked_lists[[ct]],
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000
  ) %>%
    as_tibble() %>%
    arrange(pval) %>%
    mutate(
      Contrast      = ct,
      Significant   = padj < 0.05,
      Direction     = case_when(
        Significant & NES > 0 ~ "Up",
        Significant & NES < 0 ~ "Down",
        TRUE                  ~ "NS"
      ),
      # Clean pathway names for readability
      pathway_clean = pathway %>%
        str_remove("^HALLMARK_") %>%
        str_remove("^KEGG_MEDICUS_REFERENCE_") %>%
        str_remove("^KEGG_MEDICUS_PARENT_") %>%
        str_remove("^KEGG_") %>%
        str_remove("^GOBP_") %>%
        str_replace_all("_", " ") %>%
        str_to_title() %>%
        str_trunc(60)
    )

  message("  Significant pathways (padj<0.05): ", sum(res$Significant, na.rm = TRUE))
  gsea_results[[ct]] <- res
}


# ==============================================================================
# SECTION E — Save results
# ==============================================================================

for (ct in contrast_names) {
  gsea_results[[ct]] %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) %>%
    write_tsv(file.path(output_dir, paste0("step4_GSEA_", ct, ".tsv")))
  message("Saved: step4_GSEA_", ct, ".tsv")
}

gsea_summary <- bind_rows(lapply(gsea_results, function(res) {
  tibble(
    Contrast     = res$Contrast[1],
    n_tested     = nrow(res),
    n_sig_padj05 = sum(res$padj < 0.05, na.rm = TRUE),
    n_up         = sum(res$Direction == "Up",   na.rm = TRUE),
    n_down       = sum(res$Direction == "Down", na.rm = TRUE),
    n_pval01     = sum(res$pval < 0.01, na.rm = TRUE),
    n_pval05     = sum(res$pval < 0.05, na.rm = TRUE)
  )
}))

cat("\n=== GSEA summary ===\n"); print(gsea_summary)
write_tsv(gsea_summary, file.path(output_dir, "step4_GSEA_summary.tsv"))
message("Saved: step4_GSEA_summary.tsv")

cat("\n=== Step 4 complete ===\n")
cat("Primary outputs used in manuscript:\n")
cat("  - step4_GSEA_Treat_PDAC.tsv\n")
cat("  - step4_GSEA_Treat_HNSCC.tsv\n")
