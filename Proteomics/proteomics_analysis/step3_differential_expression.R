# =============================================================================
# Step 3 — Differential expression
# Model: limma + duplicateCorrelation (patient as blocking factor)
#
# Design: group-means parameterisation with 4 groups
#   (PDAC_UT, PDAC_T, HNSCC_UT, HNSCC_T)
# Contrasts:
#   1. Treatment effect in PDAC   (PDAC_T vs PDAC_UT)   [exploratory]
#   2. Treatment effect in HNSCC  (HNSCC_T vs HNSCC_UT) [exploratory]
#   3. Cancer type at baseline    (HNSCC_UT vs PDAC_UT)  [confirmed]
#   4. Interaction: treatment difference between cancer types
#
# Input:  step2_matrix_normalized.tsv
#         step2_sample_metadata.tsv
# Output: step3_DE_<contrast>.tsv  — one file per contrast (used in step 4)
#         step3_DE_summary.tsv     — counts of significant proteins
# =============================================================================

# ---- Packages ----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("dplyr", "readr", "tibble", "limma")
new  <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new) > 0) {
  if ("limma" %in% new) BiocManager::install("limma")
  install.packages(new[new != "limma"])
}
suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))


# ==============================================================================
# SECTION A — Configuration  (EDIT THIS BLOCK)
# ==============================================================================

output_dir <- "/path/to/output"   # <-- CHANGE THIS

mat <- read_tsv(file.path(output_dir, "step2_matrix_normalized.tsv"),
                show_col_types = FALSE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

meta <- read_tsv(file.path(output_dir, "step2_sample_metadata.tsv"),
                 show_col_types = FALSE) %>%
  mutate(
    Group      = paste0(CancerType, "_", Condition),
    Condition  = factor(Condition,  levels = c("UT", "T")),
    CancerType = factor(CancerType, levels = c("PDAC", "HNSCC"))
  )

mat <- mat[, meta$Sample]


# ==============================================================================
# SECTION B — Design matrix and contrasts
# ==============================================================================

Group  <- factor(meta$Group, levels = c("PDAC_UT", "PDAC_T", "HNSCC_UT", "HNSCC_T"))
design <- model.matrix(~ 0 + Group)
colnames(design) <- levels(Group)

contrast_mat <- makeContrasts(
  Treat_PDAC    = PDAC_T   - PDAC_UT,
  Treat_HNSCC   = HNSCC_T  - HNSCC_UT,
  HNSCC_vs_PDAC = HNSCC_UT - PDAC_UT,
  Interaction   = (HNSCC_T - HNSCC_UT) - (PDAC_T - PDAC_UT),
  levels = design
)

cat("=== Design matrix (first 6 rows) ===\n"); print(head(design))
cat("\n=== Contrasts ===\n");                  print(contrast_mat)


# ==============================================================================
# SECTION C — Estimate within-patient correlation (duplicateCorrelation)
# ==============================================================================
# Accounts for 5 biological replicates nested within each patient.
# Typical intra-patient correlation range: 0.1-0.6.

message("\nEstimating within-patient correlation...")
Patient <- factor(meta$Patient)
corfit  <- duplicateCorrelation(mat, design, block = Patient)
cat(sprintf("Estimated intra-patient correlation: %.3f\n",
            corfit$consensus.correlation))


# ==============================================================================
# SECTION D — Fit linear model
# ==============================================================================

message("Fitting linear model...")
fit1 <- lmFit(mat, design, block = Patient,
              correlation = corfit$consensus.correlation)
fit2 <- contrasts.fit(fit1, contrast_mat)
fit3 <- eBayes(fit2, trend = TRUE, robust = TRUE)
# trend = TRUE: accounts for mean-variance trend (important for proteomics)
# robust = TRUE: downweights outlier proteins


# ==============================================================================
# SECTION E — Extract results per contrast
# ==============================================================================
# STRICT  (HNSCC_vs_PDAC, Interaction): FDR < 0.05, |logFC| > 1
# RELAXED (Treat_PDAC, Treat_HNSCC):    nominal p < 0.05, |logFC| > 0.5
#   Rationale: n=3 patients per cancer type gives low power for treatment
#   contrasts. Relaxed threshold is exploratory and clearly labelled.

strict_contrasts <- c("HNSCC_vs_PDAC", "Interaction")
contrast_names   <- colnames(contrast_mat)
de_results       <- list()
de_summary       <- tibble()

for (ct in contrast_names) {

  if (ct %in% strict_contrasts) {
    fc_thresh    <- 1.0
    thresh_label <- "FDR<0.05, |logFC|>1"
  } else {
    fc_thresh    <- 0.5
    thresh_label <- "p<0.05 (nominal), |logFC|>0.5 [exploratory]"
  }

  res <- topTable(fit3, coef = ct, number = Inf, sort.by = "P") %>%
    as.data.frame() %>%
    rownames_to_column("SYMBOL") %>%
    as_tibble() %>%
    rename(t_stat = t, p_value = P.Value, FDR = adj.P.Val) %>%
    mutate(
      Contrast    = ct,
      ThreshLabel = thresh_label,
      Significant = if (ct %in% strict_contrasts) {
        FDR < 0.05 & abs(logFC) > fc_thresh
      } else {
        p_value < 0.05 & abs(logFC) > fc_thresh
      },
      Direction = case_when(
        Significant & logFC > 0 ~ "Up",
        Significant & logFC < 0 ~ "Down",
        TRUE                    ~ "NS"
      )
    )

  de_results[[ct]] <- res

  de_summary <- bind_rows(de_summary, tibble(
    Contrast   = ct, Threshold = thresh_label,
    n_tested   = nrow(res),
    n_sig      = sum(res$Significant),
    n_up       = sum(res$Direction == "Up"),
    n_down     = sum(res$Direction == "Down"),
    n_FDR05    = sum(res$FDR < 0.05),
    n_p05_FC05 = sum(res$p_value < 0.05 & abs(res$logFC) > 0.5)
  ))

  write_tsv(res, file.path(output_dir, paste0("step3_DE_", ct, ".tsv")))
  message("  Saved: step3_DE_", ct, ".tsv  |  sig=", sum(res$Significant))
}

cat("\n=== Differential expression summary ===\n"); print(de_summary)
write_tsv(de_summary, file.path(output_dir, "step3_DE_summary.tsv"))

cat("\n=== Step 3 complete. Proceed to step4_GSEA.R ===\n")
