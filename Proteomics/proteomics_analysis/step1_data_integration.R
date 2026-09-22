# =============================================================================
# Step 1 — Data integration (PEAKS lfq.proteins input)
#
# Merges 6 separate PEAKS LFQ runs (one per patient) into one protein matrix.
# Proteins detected in fewer than min_patients_detected out of 3 patients
# within a cancer type are excluded.
#
# Input:  one lfq.proteins.csv per patient (PEAKS Studio output)
# Output: step1_protein_matrix_filtered.tsv  — merged, filtered intensity matrix
#         step1_sample_metadata.tsv          — sample annotation table
#         step1_summary.tsv                  — run statistics
#
# Patient naming convention used throughout:
#   PDTO1–PDTO3 = PDAC patients
#   PDTO4–PDTO6 = HNSCC patients
# =============================================================================

# ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(tibble)
})


# ==============================================================================
# SECTION A — Configuration  (EDIT THIS BLOCK)
# ==============================================================================

# Base path to the folder containing one subfolder per patient.
# Each subfolder must contain an lfq.proteins.csv file from PEAKS Studio.
# Example structure:
#   /path/to/data/PDTO1/lfq.proteins.csv
#   /path/to/data/PDTO2/lfq.proteins.csv
#   ...
base_path <- "/path/to/PEAKS_output"   # <-- CHANGE THIS

# List all 6 patient files.
# Names must match the patient IDs used in the PEAKS sample columns
# (i.e. whatever prefix appears before the underscore in "Sample X Area" columns).
peaks_files <- list(
  PDTO1 = file.path(base_path, "PDTO1", "lfq.proteins.csv"),
  PDTO2 = file.path(base_path, "PDTO2", "lfq.proteins.csv"),
  PDTO3 = file.path(base_path, "PDTO3", "lfq.proteins.csv"),
  PDTO4 = file.path(base_path, "PDTO4", "lfq.proteins.csv"),
  PDTO5 = file.path(base_path, "PDTO5", "lfq.proteins.csv"),
  PDTO6 = file.path(base_path, "PDTO6", "lfq.proteins.csv")
)

# Patient -> cancer type mapping.
# Update values to match your cancer types.
patient_cancer_type <- c(
  PDTO1 = "PDAC",
  PDTO2 = "PDAC",
  PDTO3 = "PDAC",
  PDTO4 = "HNSCC",
  PDTO5 = "HNSCC",
  PDTO6 = "HNSCC"
)

# PEAKS assigns global sample numbers across all patients in one project.
# Pattern used here: within each patient block of 10 samples,
#   odd numbers  = Treated (T)
#   even numbers = Untreated (UT)
# Replicate = ceiling(position_within_block / 2)
#
# Update peaks_num values to match the "Sample X Area" column numbers
# that appear in YOUR lfq.proteins.csv files.
sample_map <- tribble(
  ~peaks_num, ~Patient, ~Condition, ~Replicate,
  # PDTO1: samples 1-10
   1, "PDTO1", "T",  1,   2, "PDTO1", "UT", 1,
   3, "PDTO1", "T",  2,   4, "PDTO1", "UT", 2,
   5, "PDTO1", "T",  3,   6, "PDTO1", "UT", 3,
   7, "PDTO1", "T",  4,   8, "PDTO1", "UT", 4,
   9, "PDTO1", "T",  5,  10, "PDTO1", "UT", 5,
  # PDTO2: samples 21-30
  21, "PDTO2", "T",  1,  22, "PDTO2", "UT", 1,
  23, "PDTO2", "T",  2,  24, "PDTO2", "UT", 2,
  25, "PDTO2", "T",  3,  26, "PDTO2", "UT", 3,
  27, "PDTO2", "T",  4,  28, "PDTO2", "UT", 4,
  29, "PDTO2", "T",  5,  30, "PDTO2", "UT", 5,
  # PDTO3: samples 41-50
  41, "PDTO3", "T",  1,  42, "PDTO3", "UT", 1,
  43, "PDTO3", "T",  2,  44, "PDTO3", "UT", 2,
  45, "PDTO3", "T",  3,  46, "PDTO3", "UT", 3,
  47, "PDTO3", "T",  4,  48, "PDTO3", "UT", 4,
  49, "PDTO3", "T",  5,  50, "PDTO3", "UT", 5,
  # PDTO4: samples 11-20
  11, "PDTO4", "T",  1,  12, "PDTO4", "UT", 1,
  13, "PDTO4", "T",  2,  14, "PDTO4", "UT", 2,
  15, "PDTO4", "T",  3,  16, "PDTO4", "UT", 3,
  17, "PDTO4", "T",  4,  18, "PDTO4", "UT", 4,
  19, "PDTO4", "T",  5,  20, "PDTO4", "UT", 5,
  # PDTO5: samples 31-40
  31, "PDTO5", "T",  1,  32, "PDTO5", "UT", 1,
  33, "PDTO5", "T",  2,  34, "PDTO5", "UT", 2,
  35, "PDTO5", "T",  3,  36, "PDTO5", "UT", 3,
  37, "PDTO5", "T",  4,  38, "PDTO5", "UT", 4,
  39, "PDTO5", "T",  5,  40, "PDTO5", "UT", 5,
  # PDTO6: samples 51-60
  51, "PDTO6", "T",  1,  52, "PDTO6", "UT", 1,
  53, "PDTO6", "T",  2,  54, "PDTO6", "UT", 2,
  55, "PDTO6", "T",  3,  56, "PDTO6", "UT", 3,
  57, "PDTO6", "T",  4,  58, "PDTO6", "UT", 4,
  59, "PDTO6", "T",  5,  60, "PDTO6", "UT", 5
) %>%
  mutate(
    Sample     = paste0(Patient, "_", Condition, "_", Replicate),
    CancerType = patient_cancer_type[Patient]
  )

# Minimum number of patients per cancer type in which a protein must be
# detected (non-NA, non-zero in >= 1 replicate) to be retained.
# Recommended: 2 out of 3.
min_patients_detected <- 2

# Output directory for all downstream steps
output_dir <- "/path/to/output"   # <-- CHANGE THIS


# ==============================================================================
# SECTION B — Load and tidy each patient file
# ==============================================================================

load_patient_peaks <- function(file_path, patient_id) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  raw <- read_csv(file_path, show_col_types = FALSE)
  colnames(raw) <- str_trim(colnames(raw))

  # --- Gene symbol -----------------------------------------------------------
  # Primary: "Gene" column (e.g. KLHL4)
  # Fallback: parse from "Accession" column (e.g. sp|Q9C0H6|KLHL4_HUMAN)
  if ("Gene" %in% colnames(raw)) {
    raw <- raw %>%
      mutate(
        SYMBOL = str_trim(str_extract(Gene, "^[^;|]+")),
        SYMBOL = if_else(
          is.na(SYMBOL) | SYMBOL == "",
          str_extract(Accession, "(?<=\\|)[^|]+(?=_HUMAN)"),
          SYMBOL
        )
      )
  } else if ("Accession" %in% colnames(raw)) {
    raw <- raw %>%
      mutate(SYMBOL = str_extract(Accession, "(?<=\\|)[^|]+(?=_HUMAN)"))
  } else {
    stop("Cannot find Gene or Accession column in: ", file_path,
         "\nActual columns: ", paste(colnames(raw)[1:10], collapse = ", "))
  }

  # --- Find all "Sample X Area" columns in this file ------------------------
  area_cols <- grep("^Sample [0-9]+ Area$", colnames(raw), value = TRUE)

  if (length(area_cols) == 0) {
    stop("No 'Sample X Area' columns found in: ", file_path)
  }

  # Extract sample numbers and look up in the global sample_map
  sample_nums <- as.integer(str_extract(area_cols, "[0-9]+"))
  col_info    <- sample_map %>%
    dplyr::filter(peaks_num %in% sample_nums, Patient == patient_id)

  if (nrow(col_info) == 0) {
    stop("No matching entries in sample_map for patient ", patient_id,
         " with sample numbers: ", paste(sample_nums, collapse = ", "))
  }

  matched_cols <- paste0("Sample ", col_info$peaks_num, " Area")

  dat <- raw %>%
    dplyr::select(SYMBOL, all_of(matched_cols)) %>%
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "")

  # --- Zero -> NA (PEAKS uses 0 for undetected proteins) --------------------
  dat <- dat %>%
    mutate(across(all_of(matched_cols), ~ na_if(as.numeric(.x), 0)))

  # --- Rename columns to clean sample names ---------------------------------
  # "Sample 1 Area" -> "PDTO1_T_1"
  rename_vec        <- setNames(col_info$Sample, matched_cols)
  colnames(dat)[-1] <- rename_vec[matched_cols]

  # --- Collapse duplicate gene symbols --------------------------------------
  # Keep the row with the most non-NA values per symbol
  dat <- dat %>%
    group_by(SYMBOL) %>%
    slice_max(
      order_by = rowSums(!is.na(across(where(is.numeric)))),
      n = 1, with_ties = FALSE
    ) %>%
    ungroup()

  message("  Loaded ", patient_id, ": ", nrow(dat), " proteins, ",
          ncol(dat) - 1, " samples")
  dat
}

message("Loading patient files...")
patient_data <- mapply(
  load_patient_peaks,
  file_path  = peaks_files,
  patient_id = names(peaks_files),
  SIMPLIFY   = FALSE
)


# ==============================================================================
# SECTION C — Merge all patients into one matrix
# ==============================================================================
# Full outer join on SYMBOL.
# Proteins absent from a patient get NA for all of that patient's columns.

message("\nMerging patients...")
merged <- Reduce(
  function(a, b) full_join(a, b, by = "SYMBOL"),
  patient_data
)

message("After merging: ", nrow(merged), " proteins x ",
        ncol(merged) - 1, " samples")


# ==============================================================================
# SECTION D — Verify sample metadata
# ==============================================================================

metadata <- sample_map %>%
  dplyr::select(Sample, Patient, Condition, Replicate, CancerType)

sample_cols_in_matrix <- colnames(merged)[-1]
missing_from_matrix   <- setdiff(metadata$Sample, sample_cols_in_matrix)
extra_in_matrix       <- setdiff(sample_cols_in_matrix, metadata$Sample)

if (length(missing_from_matrix) > 0) {
  warning("In metadata but missing from matrix: ",
          paste(missing_from_matrix, collapse = ", "))
}
if (length(extra_in_matrix) > 0) {
  warning("In matrix but not in metadata: ",
          paste(extra_in_matrix, collapse = ", "))
}

cat("\n=== Sample metadata ===\n")
print(metadata, n = nrow(metadata))


# ==============================================================================
# SECTION E — Filter: keep proteins in >= 2/3 patients per cancer type
# ==============================================================================

is_detected <- function(x) !is.na(x) & x > 0

mat <- merged %>%
  column_to_rownames("SYMBOL") %>%
  as.matrix()

mat <- mat[, metadata$Sample[metadata$Sample %in% colnames(mat)]]

message("\nApplying coverage filter (>=", min_patients_detected,
        "/3 patients per cancer type)...")

filter_coverage <- function(mat, meta, min_patients = 2) {

  cancer_types <- unique(na.omit(meta$CancerType))
  passes <- matrix(
    FALSE,
    nrow = nrow(mat), ncol = length(cancer_types),
    dimnames = list(rownames(mat), cancer_types)
  )

  for (ct in cancer_types) {
    patients_in_type <- unique(meta$Patient[meta$CancerType == ct])

    detected_per_patient <- sapply(patients_in_type, function(p) {
      pt_samples <- intersect(meta$Sample[meta$Patient == p], colnames(mat))
      if (length(pt_samples) == 0) return(rep(FALSE, nrow(mat)))
      pt_mat <- mat[, pt_samples, drop = FALSE]
      rowSums(apply(pt_mat, 2, is_detected)) >= 1
    })

    if (is.vector(detected_per_patient)) {
      detected_per_patient <- matrix(detected_per_patient, ncol = 1)
    }

    passes[, ct] <- rowSums(detected_per_patient) >= min_patients
    message("  ", ct, ": ", sum(passes[, ct]), " proteins pass")
  }

  rowSums(passes) >= 1
}

keep         <- filter_coverage(mat, metadata, min_patients_detected)
mat_filtered <- mat[keep, ]

message("Proteins retained: ", nrow(mat_filtered), " / ", nrow(mat))


# ==============================================================================
# SECTION F — Missing value report
# ==============================================================================

pct_missing_total <- mean(is.na(mat_filtered)) * 100
cat(sprintf("\nOverall missing rate: %.1f%%\n", pct_missing_total))

pct_per_sample <- sort(colMeans(is.na(mat_filtered)) * 100, decreasing = TRUE)
cat("\nPer-sample missing rate (%):\n")
print(round(pct_per_sample, 1))

high_missing <- names(pct_per_sample)[pct_per_sample > 50]
if (length(high_missing) > 0) {
  warning("Samples >50% missing — inspect before imputation:\n  ",
          paste(high_missing, collapse = "\n  "))
} else {
  message("All samples <50% missing — good to proceed.")
}

missing_summary <- metadata %>%
  mutate(pct_missing = round(pct_per_sample[Sample], 1)) %>%
  group_by(CancerType, Condition) %>%
  summarise(mean_pct_missing = mean(pct_missing, na.rm = TRUE),
            .groups = "drop")

cat("\nMean missing rate by group:\n")
print(missing_summary)


# ==============================================================================
# SECTION G — Save outputs
# ==============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mat_filtered %>%
  as.data.frame() %>%
  rownames_to_column("SYMBOL") %>%
  write_tsv(file.path(output_dir, "step1_protein_matrix_filtered.tsv"))
message("\nSaved: step1_protein_matrix_filtered.tsv")

write_tsv(metadata, file.path(output_dir, "step1_sample_metadata.tsv"))
message("Saved: step1_sample_metadata.tsv")

tibble(
  total_proteins_before_filter = nrow(mat),
  proteins_after_filter        = nrow(mat_filtered),
  pct_missing_overall          = round(pct_missing_total, 1),
  n_samples_total              = ncol(mat_filtered),
  n_patients                   = length(peaks_files),
  filter_rule                  = paste0(">=", min_patients_detected,
                                        "/3 patients per cancer type")
) %>%
  write_tsv(file.path(output_dir, "step1_summary.tsv"))
message("Saved: step1_summary.tsv")

cat("\n=== Step 1 complete. Proceed to step2_qc_normalization_pca.R ===\n")
