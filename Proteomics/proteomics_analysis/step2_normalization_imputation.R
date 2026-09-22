# =============================================================================
# Step 2 — Normalisation and imputation
#
# Input:  step1_protein_matrix_filtered.tsv
#         step1_sample_metadata.tsv
# Output: step2_matrix_normalized.tsv  — log2, quantile-normalised,
#                                        MinProb-imputed matrix
#         step2_sample_metadata.tsv
# =============================================================================

# ---- Packages ----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs_cran <- c("dplyr", "tidyr", "readr", "tibble")
pkgs_bioc <- c("limma")

new_cran <- pkgs_cran[!sapply(pkgs_cran, requireNamespace, quietly = TRUE)]
new_bioc <- pkgs_bioc[!sapply(pkgs_bioc, requireNamespace, quietly = TRUE)]
if (length(new_cran) > 0) install.packages(new_cran)
if (length(new_bioc) > 0) BiocManager::install(new_bioc)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(limma)
})


# ==============================================================================
# SECTION A — Configuration  (EDIT THIS BLOCK)
# ==============================================================================

# Directory where step 1 outputs were saved and where step 2 outputs will go
output_dir <- "/path/to/output"   # <-- CHANGE THIS

# Read step 1 outputs
mat_raw <- read_tsv(file.path(output_dir, "step1_protein_matrix_filtered.tsv"),
                    show_col_types = FALSE) %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

meta <- read_tsv(file.path(output_dir, "step1_sample_metadata.tsv"),
                 show_col_types = FALSE) %>%
  mutate(Group = paste0(CancerType, "_", Condition))

mat_raw <- mat_raw[, meta$Sample]


# ==============================================================================
# SECTION B — Log2 transformation
# ==============================================================================

mat_raw[mat_raw == 0] <- NA
mat_log2 <- log2(mat_raw)
message("Log2 transformation done.")


# ==============================================================================
# SECTION C — Quantile normalisation
# ==============================================================================

mat_norm <- normalizeBetweenArrays(mat_log2, method = "quantile")
message("Quantile normalisation done.")


# ==============================================================================
# SECTION D — MinProb imputation
# ==============================================================================
# Imputes missing values by drawing from a low-intensity Gaussian
# (1st percentile, width = 0.3 * SD of observed values per sample).
# Simulates proteins near the detection limit.

set.seed(42)

minprob_impute <- function(mat, q = 0.01, width = 0.3) {
  mat_imp <- mat
  for (j in seq_len(ncol(mat))) {
    col       <- mat[, j]
    obs       <- col[!is.na(col)]
    if (length(obs) < 10) next
    mu_imp    <- quantile(obs, probs = q, na.rm = TRUE)
    sd_imp    <- sd(obs, na.rm = TRUE) * width
    n_missing <- sum(is.na(col))
    mat_imp[is.na(col), j] <- rnorm(n_missing, mean = mu_imp, sd = sd_imp)
  }
  mat_imp
}

mat_imputed <- minprob_impute(mat_norm, q = 0.01, width = 0.3)
message("MinProb imputation done. Missing values remaining: ", sum(is.na(mat_imputed)))


# ==============================================================================
# SECTION E — Save outputs
# ==============================================================================

mat_imputed %>%
  as.data.frame() %>%
  rownames_to_column("SYMBOL") %>%
  write_tsv(file.path(output_dir, "step2_matrix_normalized.tsv"))
message("Saved: step2_matrix_normalized.tsv")

write_tsv(meta, file.path(output_dir, "step2_sample_metadata.tsv"))
message("Saved: step2_sample_metadata.tsv")

cat("\n=== Step 2 complete. Proceed to step3_differential_expression.R ===\n")
