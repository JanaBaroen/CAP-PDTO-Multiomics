# ============================================================
# WES Somatic Variant Analysis — HNSCC & PDAC
# Treated (T) vs Untreated (UT) tumor pairs
# Tools: maftools, clusterProfiler, ggplot2
# ============================================================

# ============================================================
# SECTION 1: Install and load packages
# ============================================================

# Run this block once to install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "maftools",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ReactomePA"
))

install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "ggpubr"))

# Load packages
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
# SECTION 2: Load MAF files
# ============================================================

# CHANGE THIS: directory containing your MAF files from vcf2maf
MAF_DIR <- "/your/path/to/MAF_files"

# Define patient metadata
patient_info <- data.frame(
  Tumor_Sample_Barcode = c(
    "PDTO4_treated", "PDTO5_treated", "PDTO6_treated",  # HNSCC treated
    "PDTO1_treated", "PDTO2_treated", "PDTO3_treated"   # PDAC treated
  ),
  patient = c("PDTO4", "PDTO5", "PDTO6",
              "PDTO1", "PDTO2", "PDTO3"),
  cancer_type = c("HNSCC", "HNSCC", "HNSCC",
                  "PDAC",  "PDAC",  "PDAC"),
  status = rep("Treated", 6),
  stringsAsFactors = FALSE
)

# Load individual MAF files
maf_files <- list(
  PDTO4 = file.path(MAF_DIR, "PDTO4_treated_vs_PDTO4_untreated.maf"),
  PDTO5 = file.path(MAF_DIR, "PDTO5_treated_vs_PDTO5_untreated.maf"),
  PDTO6 = file.path(MAF_DIR, "PDTO6_treated_vs_PDTO6_untreated.maf"),
  PDTO1 = file.path(MAF_DIR, "PDTO1_treated_vs_PDTO1_untreated.maf"),
  PDTO2 = file.path(MAF_DIR, "PDTO2_treated_vs_PDTO2_untreated.maf"),
  PDTO3 = file.path(MAF_DIR, "PDTO3_treated_vs_PDTO3_untreated.maf")
)

# Read all MAF files
mafs <- lapply(names(maf_files), function(patient) {
  read.maf(
    maf = maf_files[[patient]],
    clinicalData = patient_info[grep(
      substr(patient, 1, 4),
      patient_info$patient
    ), ],
    verbose = FALSE
  )
})
names(mafs) <- names(maf_files)

# Merge all MAFs into one combined object
maf_all <- merge_mafs(mafs)

# Separate by cancer type
maf_hnscc <- merge_mafs(mafs[c("PDTO4", "PDTO5", "PDTO6")])
maf_pdac  <- merge_mafs(mafs[c("PDTO1", "PDTO2", "PDTO3")])

# CHANGE THIS: directory where output files will be saved
output_dir <- "/your/path/to/WES_output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# SECTION 3: Tumor Mutational Burden (TMB)
# ============================================================

# Calculate TMB for each sample
# WES captures approximately 38 Mb of coding sequence
tmb_all <- tmb(
  maf = maf_all,
  captureSize = 38,
  logScale = FALSE
)
dev.copy(png, file.path(output_dir, "tmb_plot.png"), width = 1200, height = 800, res = 120)
dev.off()

print("TMB per sample:")
print(tmb_all)

# Add cancer type and patient info
tmb_all$cancer_type <- ifelse(
  tmb_all$Tumor_Sample_Barcode %in% c("PDTO4_treated", "PDTO5_treated", "PDTO6_treated"),
  "HNSCC", "PDAC"
)
tmb_all$patient <- c(
  "PDTO4", "PDTO5", "PDTO6",
  "PDTO1", "PDTO2", "PDTO3"
)[match(
  tmb_all$Tumor_Sample_Barcode,
  c("PDTO4_treated", "PDTO5_treated", "PDTO6_treated",
    "PDTO1_treated", "PDTO2_treated", "PDTO3_treated")
)]

# Plot TMB comparison
ggplot(tmb_all, aes(x = patient, y = total_perMB, fill = cancer_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("HNSCC" = "#4E79A7", "PDAC" = "#F28E2B")) +
  labs(
    title = "Tumor Mutational Burden per Patient",
    x = "Patient",
    y = "Mutations per Mb",
    fill = "Cancer type"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  scale_y_continuous(limits = c(0, 0.35))
ggsave(file.path(output_dir, "tmb_per_patient.png"), width = 8, height = 6, dpi = 120)

# ============================================================
# SECTION 4: Extract gene lists for enrichment analysis
# ============================================================

# Get all mutated genes per cohort
genes_all   <- getGeneSummary(maf_all)$Hugo_Symbol
genes_hnscc <- getGeneSummary(maf_hnscc)$Hugo_Symbol
genes_pdac  <- getGeneSummary(maf_pdac)$Hugo_Symbol

# Get high confidence genes (mutated in >1 sample per cohort)
gene_summary_hnscc <- getGeneSummary(maf_hnscc)
gene_summary_pdac  <- getGeneSummary(maf_pdac)

genes_hnscc_recurrent <- gene_summary_hnscc %>%
  filter(MutatedSamples >= 2) %>%
  pull(Hugo_Symbol)

genes_pdac_recurrent <- gene_summary_pdac %>%
  filter(MutatedSamples >= 2) %>%
  pull(Hugo_Symbol)

cat("Recurrent HNSCC genes (>=2 samples):", length(genes_hnscc_recurrent), "\n")
cat("Recurrent PDAC genes (>=2 samples):", length(genes_pdac_recurrent), "\n")
print(genes_hnscc_recurrent)
print(genes_pdac_recurrent)

# ============================================================
# SECTION 5: Pathway enrichment analysis — clusterProfiler
# ============================================================

# Convert gene symbols to Entrez IDs
convert_to_entrez <- function(genes) {
  result <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  return(result$ENTREZID)
}

entrez_all   <- convert_to_entrez(genes_all)
entrez_hnscc <- convert_to_entrez(genes_hnscc)
entrez_pdac  <- convert_to_entrez(genes_pdac)

# Reactome pathway enrichment
run_reactome_enrichment <- function(entrez_ids, label) {
  ereact <- enrichPathway(
    gene          = entrez_ids,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    readable      = TRUE
  )
  
  if (!is.null(ereact) && nrow(as.data.frame(ereact)) > 0) {
    print(paste("Reactome enrichment results for", label))
    print(head(as.data.frame(ereact), 20))
    print(dotplot(ereact, showCategory = 15,
                  title = paste("Reactome pathways —", label)))
    dev.copy(png, file.path(output_dir, paste0("Reactome_dotplot_", label, ".png")),
             width = 1200, height = 1000, res = 120)
    dev.off()
  } else {
    cat("No significant Reactome pathways for", label, "\n")
  }
  return(ereact)
}

react_all   <- run_reactome_enrichment(entrez_all,   "All_samples")

react_hnscc <- run_reactome_enrichment(entrez_hnscc, "HNSCC")
if (!is.null(react_hnscc) && nrow(as.data.frame(react_hnscc)) > 0) {
  dotplot(react_hnscc, showCategory = 15, title = "Reactome pathways — HNSCC")
  dev.copy(png, file.path(output_dir, "Reactome_dotplot_HNSCC.png"), 
           width = 1200, height = 1000, res = 120)
  dev.off()
}

react_pdac  <- run_reactome_enrichment(entrez_pdac,  "PDAC")
if (!is.null(react_pdac) && nrow(as.data.frame(react_pdac)) > 0) {
  dotplot(react_pdac, showCategory = 15, title = "Reactome pathways — PDAC")
  dev.copy(png, file.path(output_dir, "Reactome_dotplot_PDAC.png"), 
           width = 1200, height = 1000, res = 120)
  dev.off()
}

write.csv(as.data.frame(react_hnscc), 
          file.path(output_dir, "Reactome_enrichment_HNSCC.csv"), 
          row.names = FALSE)
write.csv(as.data.frame(react_pdac), 
          file.path(output_dir, "Reactome_enrichment_PDAC.csv"), 
          row.names = FALSE)

# ============================================================
# SECTION 6: Save results to files
# ============================================================

# Save gene summaries
write.csv(
  getGeneSummary(maf_all),
  file.path(output_dir, "gene_summary_all.csv"),
  row.names = FALSE
)
write.csv(
  getGeneSummary(maf_hnscc),
  file.path(output_dir, "gene_summary_hnscc.csv"),
  row.names = FALSE
)
write.csv(
  getGeneSummary(maf_pdac),
  file.path(output_dir, "gene_summary_pdac.csv"),
  row.names = FALSE
)

# Save TMB results
write.csv(
  tmb_all,
  file.path(output_dir, "tmb_per_sample.csv"),
  row.names = FALSE
)


# Save mutated gene lists
write.table(
  genes_hnscc,
  file.path(output_dir, "mutated_genes_HNSCC.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
write.table(
  genes_pdac,
  file.path(output_dir, "mutated_genes_PDAC.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

cat("\n=== Analysis complete! ===\n")
cat("Results saved to:", output_dir, "\n")
