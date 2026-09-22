# ============================================================
# Step 5: Hotspot whitelisting + consequence filtering in R
# ============================================================
library(maftools)
library(data.table)

# CHANGE THIS: your root scratch/working directory
SCRATCH <- "/your/scratch/directory"

# CHANGE THIS: directory containing filtered MAF files from vcf2maf
MAF_DIR <- "/your/path/to/MAF_files/filtered"

# CHANGE THIS: output directory for final filtered MAFs
OUT_DIR <- "/your/path/to/MAF_files/final"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- Load all filtered MAFs ----
# Sample IDs (corresponding to associated manuscript)
patients <- c("PDTO1","PDTO2","PDTO3",
              "PDTO4","PDTO5","PDTO6")

maf_list <- lapply(patients, function(p) {
    f <- file.path(MAF_DIR, paste0(p, ".filtered.maf"))
    if (!file.exists(f)) { 
        warning(paste("Missing:", f))
        return(NULL) 
    }
    fread(f, skip = "Hugo_Symbol")
})
combined <- rbindlist(maf_list, fill = TRUE)
cat("Total variants loaded:", nrow(combined), "\n")

# ---- Consequence filtering ----
keep_variants <- c(
    "Missense_Mutation", "Nonsense_Mutation",
    "Frame_Shift_Del",   "Frame_Shift_Ins",
    "Splice_Site",       "In_Frame_Del",
    "In_Frame_Ins",      "Nonstop_Mutation",
    "Translation_Start_Site"
)
maf_func <- combined[Variant_Classification %in% keep_variants]
cat("After consequence filter:", nrow(maf_func), "\n")

# ---- Cancer hotspot whitelisting ----
# Key HNSCC + PDAC hotspots to always protect
# Gene + amino acid position (NA = whole gene protected)
# See README for details
hotspots <- data.table(
    Hugo_Symbol = c(
        "KRAS","KRAS","KRAS","KRAS",
        "TP53","TP53","TP53","TP53","TP53","TP53",
        "PIK3CA","PIK3CA","PIK3CA",
        "EGFR","EGFR","EGFR",
        "CDKN2A","SMAD4","RNF43",
        "GNAS","BRCA2","BRCA1",
        "FAT1","NOTCH1","HRAS",
        "ARID1A","ATM",
        "KMT2D"
    ),
    AA_pos = c(
        12, 13, 61, 146,
        175, 245, 248, 249, 273, 282,
        542, 545, 1047,
        858, 719, 746,
        58, 361, 117,
        201, 1813, 1775,
        NA, NA, 61,
        NA, NA,
        NA
    )
)

# Mark hotspot variants (by gene + amino acid position)
maf_func[, AA_num := as.numeric(gsub("[^0-9]", "", HGVSp_Short))]
maf_func[, is_hotspot := FALSE]

for (i in seq_len(nrow(hotspots))) {
    gene <- hotspots$Hugo_Symbol[i]
    pos  <- hotspots$AA_pos[i]
    if (is.na(pos)) {
        # Gene-level hotspot — protect any mutation in this gene
        maf_func[Hugo_Symbol == gene, is_hotspot := TRUE]
    } else {
        # Position-level hotspot
        maf_func[Hugo_Symbol == gene & AA_num == pos, is_hotspot := TRUE]
    }
}
cat("Hotspot-protected variants:", sum(maf_func$is_hotspot), "\n")

# ---- Final gnomAD safety check ----
if ("vep_gnomADe_AF" %in% colnames(maf_func)) {
    maf_func[, gnomad_num := suppressWarnings(as.numeric(vep_gnomADe_AF))]
    maf_final <- maf_func[
        is_hotspot == TRUE |
        is.na(gnomad_num) |
        gnomad_num <= 0.001
    ]
    maf_func[, gnomad_num := NULL]
} else {
    warning("gnomAD AF column not found - skipping safety check")
    maf_final <- maf_func
}

# ---- Clean up helper columns ----
maf_final[, c("AA_num", "is_hotspot") := NULL]
cat("Final variants after all filtering:", nrow(maf_final), "\n")

# ---- Save per-sample MAFs ----
for (p in patients) {
    sample_maf <- maf_final[Tumor_Sample_Barcode == p]
    fwrite(sample_maf,
           file.path(OUT_DIR, paste0(p, ".final.maf")),
           sep = "\t")
    cat(p, "->", nrow(sample_maf), "variants\n")
}

# ---- Save combined MAF ----
fwrite(maf_final,
       file.path(OUT_DIR, "all_samples.final.maf"),
       sep = "\t")

# ---- Load into maftools and verify ----
clinical <- data.frame(
    Tumor_Sample_Barcode = patients,
    cancer_type = c("PDAC","PDAC","PDAC","HNSCC","HNSCC","HNSCC")
)

maf_obj <- read.maf(
    maf          = file.path(OUT_DIR, "all_samples.final.maf"),
    clinicalData = clinical
)

# Print summary
cat("\n=== FINAL SUMMARY ===\n")
print(getSampleSummary(maf_obj))
tmb_result <- tmb(maf_obj)
print(tmb_result)

# Save summary plot
pdf(file.path(OUT_DIR, "maf_summary_filtered.pdf"), width=12, height=8)
plotmafSummary(maf_obj, rmOutlier=TRUE, addStat="median")
dev.off()
