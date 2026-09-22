###############################################################################
# RNA-seq pipeline from FASTQ → gene-level counts
# Paired-end; merges technical lanes; runs FastQC, TrimGalore (optional),
# Salmon quantification, tximport → exp.mat, pdat.mat.
###############################################################################

############################### 0. USER SETTINGS ##############################

options(repos = c(CRAN = "https://cloud.r-project.org"))
external_drive <- "/your/path/to/scratch"   # CHANGE THIS
Sys.setenv(PATH = paste(Sys.getenv("PATH"), sep = ":"))

# Top-level RNAseq folder containing patient subfolders
# Structure expected:
#   RNAseq/
#     PDTO1/
#       UT/   <- untreated FASTQs
#       T/    <- treated FASTQs
#     PDTO2/
#       UT/
#       T/
#     ... etc.

rnaseq_root <- file.path(external_drive, "RNAseq")

# Project directory for pipeline outputs
proj_dir <- file.path(external_drive, "RNAseq_project")

trim_threads <- 8
salmon_threads <- 16

# Output directories
trim_dir   <- file.path(proj_dir, "trimmed")
merged_dir <- file.path(proj_dir, "merged")
quant_dir  <- file.path(proj_dir, "quant")
qc_dir     <- file.path(proj_dir, "qc")
res_dir    <- file.path(proj_dir, "results")

dirs <- c(trim_dir, merged_dir, quant_dir, qc_dir, res_dir)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Sample sheet
sample_sheet_csv <- file.path(proj_dir, "samples.csv")

# Reference files — GENCODE v49
ref_dir        <- file.path(external_drive, "RNAseq", "RNAseq_ref")
transcripts_fa <- file.path(ref_dir, "gencode.v49.transcripts.fa")
annotation_gtf <- file.path(ref_dir, "gencode.v49.primary_assembly.annotation.gtf")
salmon_index   <- file.path(ref_dir, "salmon_index")

# Options
do_trimming    <- TRUE
salmon_libtype <- "A"

###############################################################################
########################### 1. PACKAGE SETUP ##################################
###############################################################################

need <- c("tximport", "DESeq2", "readr", "dplyr", "stringr",
          "tibble", "rtracklayer", "purrr")

to_install <- need[!need %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(setdiff(to_install, c("DESeq2","tximport","rtracklayer")))
  if (!"BiocManager" %in% rownames(installed.packages())) {
    install.packages("BiocManager")
  }
  BiocManager::install(intersect(to_install, c("DESeq2","tximport","rtracklayer")))
}

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(rtracklayer)
  library(purrr)
})

###############################################################################
############################ 2. CHECK BINARIES ################################
###############################################################################

require_bin <- function(bin) {
  if (Sys.which(bin) == "") stop(sprintf("ERROR: '%s' not found in PATH.", bin))
}

require_bin("salmon")
require_bin("fastqc")
require_bin("multiqc")
if (do_trimming) require_bin("trim_galore")

###############################################################################
############################ 3. LOAD SAMPLE SHEET #############################
###############################################################################

if (!file.exists(sample_sheet_csv)) {
  stop("Sample sheet missing: ", sample_sheet_csv)
}

pdat.mat <- read_csv(sample_sheet_csv, show_col_types = FALSE) |>
  mutate(
    sample_id = as.character(sample_id),
    Type      = factor(Type,      levels = c("HNSCC", "PDAC")),
    Condition = factor(Condition, levels = c("Untreated", "Treated"))
  ) |>
  distinct(sample_id, .keep_all = TRUE) |>
  column_to_rownames("sample_id")

###############################################################################
############### 4. MERGE TECHNICAL LANES PER SAMPLE ##########################
#
# Folder structure on drive:
#   RNAseq/<Patient>/<UT or T>/<sample_id>_L00#_R1_001.fastq.gz
#
# The function looks up each sample's Patient and Condition from pdat.mat,
# maps Condition → subfolder name (Untreated → UT, Treated → T),
# then searches that subfolder for all lane files for the sample.
###############################################################################

merge_fastqs <- function(sample) {
  
  patient   <- pdat.mat[sample, "Patient"]
  condition <- pdat.mat[sample, "Condition"]
  
  # Map condition to folder name
  cond_folder <- ifelse(as.character(condition) == "Untreated", "UT", "T")
  
  # Build path to the correct subfolder
  sample_dir <- file.path(rnaseq_root, patient, cond_folder)
  
  if (!dir.exists(sample_dir)) {
    stop("Directory not found: ", sample_dir)
  }
  
  # Find all lane files for this sample (matched by GC number prefix)
  r1 <- list.files(sample_dir, pattern = paste0("^", sample, ".*R1.*\\.fastq(\\.gz)?$"),
                   full.names = TRUE)
  r2 <- list.files(sample_dir, pattern = paste0("^", sample, ".*R2.*\\.fastq(\\.gz)?$"),
                   full.names = TRUE)
  
  if (length(r1) == 0) stop("No R1 reads found for sample ", sample, " in ", sample_dir)
  if (length(r2) == 0) stop("No R2 reads found for sample ", sample, " in ", sample_dir)
  
  # Sort to ensure lanes are concatenated in consistent order
  r1 <- sort(r1)
  r2 <- sort(r2)
  
  message("Merging ", length(r1), " lane file(s) for sample: ", sample,
          " [", patient, " / ", cond_folder, "]")
  
  out_r1 <- file.path(merged_dir, paste0(sample, "_R1.fastq.gz"))
  out_r2 <- file.path(merged_dir, paste0(sample, "_R2.fastq.gz"))
  
  if (!file.exists(out_r1))
    system(paste("cat", paste(shQuote(r1), collapse = " "), ">", shQuote(out_r1)))
  
  if (!file.exists(out_r2))
    system(paste("cat", paste(shQuote(r2), collapse = " "), ">", shQuote(out_r2)))
  
  return(c(R1 = out_r1, R2 = out_r2))
}

merged_pairs <- lapply(rownames(pdat.mat), merge_fastqs)
names(merged_pairs) <- rownames(pdat.mat)

###############################################################################
########################### 5. FASTQC + MULTIQC ###############################
###############################################################################

for (s in names(merged_pairs)) {
  fqc_out <- file.path(qc_dir, paste0(s, "_R1_fastqc.html"))
  if (!file.exists(fqc_out)) {
    system2("fastqc", c(merged_pairs[[s]][["R1"]], merged_pairs[[s]][["R2"]], "-o", qc_dir))
  } else {
    message("Skipping FastQC for ", s, " (already done)")
  }
}

if (!file.exists(file.path(qc_dir, "multiqc_report.html"))) {
  system2("multiqc", c(qc_dir, "-o", qc_dir))
}

###############################################################################
########################## 6. OPTIONAL: TRIMMING ##############################
###############################################################################

get_input_files <- function(sample) {
  if (!do_trimming) return(merged_pairs[[sample]])
  
  r1 <- merged_pairs[[sample]][["R1"]]
  r2 <- merged_pairs[[sample]][["R2"]]
  
  tr1_check <- file.path(trim_dir, paste0(sample, "_R1_val_1.fq.gz"))
  if (!file.exists(tr1_check)) {
    system2("trim_galore",
            c("--paired", "--cores", trim_threads, "--output_dir", trim_dir, r1, r2))
  } else {
    message("Skipping trimming for ", sample, " (already done)")
  }
  
  tr1 <- list.files(trim_dir, pattern = paste0(sample, ".*val_1.fq.gz"), full.names = TRUE)
  tr2 <- list.files(trim_dir, pattern = paste0(sample, ".*val_2.fq.gz"), full.names = TRUE)
  
  return(c(R1 = tr1, R2 = tr2))
}

inputs_for_quant <- lapply(rownames(pdat.mat), get_input_files)
names(inputs_for_quant) <- rownames(pdat.mat)

###############################################################################
############################ 7. BUILD SALMON INDEX ############################
###############################################################################

if (!dir.exists(salmon_index) || length(list.files(salmon_index)) == 0) {
  system2("salmon", c("index", "-t", transcripts_fa, "-i", salmon_index, "-k", 31))
}

###############################################################################
############################ 8. SALMON QUANT ##################################
###############################################################################

run_salmon <- function(sample) {
  out <- file.path(quant_dir, sample)
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  
  r1 <- inputs_for_quant[[sample]][["R1"]]
  r2 <- inputs_for_quant[[sample]][["R2"]]
  
  system2("salmon",
          c("quant", "-i", salmon_index, "-l", salmon_libtype,
            "-1", r1, "-2", r2,
            "-p", salmon_threads, "--gcBias", "--validateMappings",
            "-o", out))
  return(file.path(out, "quant.sf"))
}

quant_files <- sapply(rownames(pdat.mat), run_salmon)

###############################################################################
############################ 9. BUILD TX2GENE MAP #############################
###############################################################################

gtf <- rtracklayer::import(annotation_gtf)
df  <- as.data.frame(mcols(gtf))

tx2gene <- df |>
  dplyr::select(transcript_id, gene_id) |>
  filter(!is.na(transcript_id), !is.na(gene_id)) |>
  distinct()

colnames(tx2gene) <- c("TXNAME", "GENEID")

###############################################################################
############################ 10. TXIMPORT #####################################
###############################################################################

names(quant_files) <- rownames(pdat.mat)

txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreAfterBar = TRUE)

exp.mat <- round(txi$counts)

pdat.mat <- pdat.mat[colnames(exp.mat), , drop = FALSE]

###############################################################################
############################ 11. SAVE OUTPUTS #################################
###############################################################################

saveRDS(exp.mat,  file = file.path(res_dir, "exp.mat.rds"))
saveRDS(pdat.mat, file = file.path(res_dir, "pdat.mat.rds"))

write_tsv(as.data.frame(exp.mat) |> rownames_to_column("gene_id"),
          file.path(res_dir, "exp.mat.tsv"))
write_tsv(pdat.mat |> rownames_to_column("sample_id"),
          file.path(res_dir, "pdat.mat.tsv"))

message("✅ Pipeline complete! Files saved to: ", res_dir)