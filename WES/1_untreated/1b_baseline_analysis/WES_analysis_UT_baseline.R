# ============================================================
# WES Somatic Variant Analysis — HNSCC & PDAC
# Untreated (UT) tumor samples — baseline mutation landscape
# Tools: maftools, clusterProfiler, ggplot2
# ============================================================

# ============================================================
# SECTION 1: Load packages
# ============================================================

library(maftools)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpubr)

# ============================================================
# SECTION 2: Load and filter MAF files
# ============================================================

# Set directories
# CHANGE THIS: directory containing your final filtered MAF files
MAF_DIR_UT <- "/your/path/to/MAF_files/final"

# CHANGE THIS: directory where output figures and tables will be saved
output_dir <- "/your/path/to/WES_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Define patient metadata for UT samples
# Tumor_Sample_Barcode must match what is inside the MAF files
patient_info_ut <- data.frame(
  Tumor_Sample_Barcode = c(
    "PDTO4", "PDTO5", "PDTO6",  # HNSCC untreated
    "PDTO1",  "PDTO2",  "PDTO3"    # PDAC untreated
  ),
  patient = c("PDTO4", "PDTO5", "PDTO6",
              "PDTO1",  "PDTO2",  "PDTO3"),
  cancer_type = c("HNSCC", "HNSCC", "HNSCC",
                  "PDAC",  "PDAC",  "PDAC"),
  status = rep("Untreated", 6),
  stringsAsFactors = FALSE
)

# MAF file paths
maf_files_ut <- list(
  PDTO4_UT = file.path(MAF_DIR_UT, "PDTO4.final.maf"),
  PDTO5_UT = file.path(MAF_DIR_UT, "PDTO5.final.maf"),
  PDTO6_UT = file.path(MAF_DIR_UT, "PDTO6.final.maf"),
  PDTO1_UT  = file.path(MAF_DIR_UT, "PDTO1.final.maf"),
  PDTO2_UT  = file.path(MAF_DIR_UT, "PDTO2.final.maf"),
  PDTO3_UT  = file.path(MAF_DIR_UT, "PDTO3.final.maf")
)

# Lookup table: list name -> Tumor_Sample_Barcode (must match MAF internals)
barcode_lookup <- c(
  PDTO4_UT = "PDTO4",
  PDTO5_UT = "PDTO5",
  PDTO6_UT = "PDTO6",
  PDTO1_UT  = "PDTO1",
  PDTO2_UT  = "PDTO2",
  PDTO3_UT  = "PDTO3"
)

# Read MAF files
# IMPORTANT: These are tumor-only MAFs so we use PASS variants only
# to remove germline, contamination, and other non-somatic variants
mafs_ut <- lapply(names(maf_files_ut), function(patient) {
  cat("Loading", patient, "...\n")
  read.maf(
    maf = maf_files_ut[[patient]],
    clinicalData = patient_info_ut[
      patient_info_ut$Tumor_Sample_Barcode == barcode_lookup[patient], ],
    vc_nonSyn = c(
      "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
      "Translation_Start_Site", "Nonsense_Mutation",
      "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
      "Missense_Mutation"
    ),
    verbose = FALSE
  )
})
names(mafs_ut) <- names(maf_files_ut)

# Merge all UT MAFs
maf_all_ut <- merge_mafs(mafs_ut)

# Separate by cancer type
maf_hnscc_ut <- merge_mafs(mafs_ut[c("PDTO4_UT", "PDTO5_UT", "PDTO6_UT")])
maf_pdac_ut  <- merge_mafs(mafs_ut[c("PDTO1_UT", "PDTO2_UT", "PDTO3_UT")])

# Check variant counts after loading
cat("\n=== Variant counts after loading ===\n")
cat("All UT samples:", nrow(maf_all_ut@data), "variants\n")
cat("HNSCC UT:", nrow(maf_hnscc_ut@data), "variants\n")
cat("PDAC UT:", nrow(maf_pdac_ut@data), "variants\n")

# ============================================================
# SECTION 3: Tumor Mutational Burden (TMB)
# ============================================================
tmb_all_ut <- tmb(
  maf = maf_all_ut,
  captureSize = 38,  # CHANGE THIS: update if using a different exome capture kit
  logScale = FALSE
)
print("TMB per UT sample:")
print(tmb_all_ut)

# Add cancer type info
tmb_all_ut$cancer_type <- ifelse(
  tmb_all_ut$Tumor_Sample_Barcode %in% c("PDTO4", "PDTO5", "PDTO6"),
  "HNSCC", "PDAC"
)
tmb_all_ut$patient <- tmb_all_ut$Tumor_Sample_Barcode

# Save TMB results
write.csv(tmb_all_ut,
          file.path(output_dir, "tmb_per_sample_UT_filtered.csv"), row.names = FALSE)


# ============================================================
# SECTION 4: Extract gene lists for enrichment analysis
# ============================================================

genes_all_ut   <- getGeneSummary(maf_all_ut)$Hugo_Symbol
genes_hnscc_ut <- getGeneSummary(maf_hnscc_ut)$Hugo_Symbol
genes_pdac_ut  <- getGeneSummary(maf_pdac_ut)$Hugo_Symbol

# Get recurrently mutated genes (>=2 samples)
gene_summary_hnscc_ut <- getGeneSummary(maf_hnscc_ut)
gene_summary_pdac_ut  <- getGeneSummary(maf_pdac_ut)

genes_hnscc_ut_recurrent <- gene_summary_hnscc_ut %>%
  filter(MutatedSamples >= 2) %>%
  pull(Hugo_Symbol)

genes_pdac_ut_recurrent <- gene_summary_pdac_ut %>%
  filter(MutatedSamples >= 2) %>%
  pull(Hugo_Symbol)

cat("Recurrent HNSCC UT genes (>=2 samples):", length(genes_hnscc_ut_recurrent), "\n")
cat("Recurrent PDAC UT genes (>=2 samples):", length(genes_pdac_ut_recurrent), "\n")
print(genes_hnscc_ut_recurrent)
print(genes_pdac_ut_recurrent)

# Save gene lists
write.table(genes_hnscc_ut,
            file.path(output_dir, "mutated_genes_HNSCC_UT_filtered.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genes_pdac_ut,
            file.path(output_dir, "mutated_genes_PDAC_UT_filtered.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
