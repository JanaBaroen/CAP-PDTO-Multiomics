# Proteomics Analysis Pipeline — LFQ Mass Spectrometry

R analysis pipeline for label-free quantification (LFQ) proteomics data from PEAKS Studio, applied to patient-derived tumour organoid cultures from two cancer types (PDAC and HNSCC).

## Study design

- **6 patients**: PDTO1–PDTO3 (PDAC), PDTO4–PDTO6 (HNSCC)
- **2 conditions** per patient: Treated (T) and Untreated (UT)
- **5 biological replicates** per patient per condition
- **60 samples** total
- Protein identification and LFQ quantification: **PEAKS Studio**
- Patient-specific FASTA databases (canonical human proteome + patient-specific variant peptides)

## Scripts

Run scripts in order. Each script reads the output of the previous step.

| Script | Description |
|--------|-------------|
| `step1_data_integration.R` | Merge 6 separate PEAKS runs into one protein matrix; coverage filtering |
| `step2_normalization_imputation.R` | Log2 transform, quantile normalisation, MinProb imputation |
| `step3_differential_expression.R` | limma + duplicateCorrelation; 4 contrasts; DE results saved as TSV |
| `step4_GSEA.R` | Pre-ranked GSEA (fgsea); Hallmarks, KEGG, GO BP gene sets |

## Primary outputs used in manuscript

- `step4_GSEA_Treat_PDAC.tsv` — GSEA results, treatment effect in PDAC
- `step4_GSEA_Treat_HNSCC.tsv` — GSEA results, treatment effect in HNSCC

## Getting started

### 1. Set paths

Every script has a **Section A — Configuration** block at the top. Edit the lines marked `# <-- CHANGE THIS` before running:

```r
# All scripts (steps 1-4)
output_dir <- "/path/to/output"       # where all outputs will be saved

# Step 1 only
base_path  <- "/path/to/PEAKS_output" # folder containing one subfolder per patient

# Step 4 only
gsea_data_dir <- "/path/to/gsea_genesets"  # folder containing .gmt files
```

### 2. Update sample number mapping (step 1 only)

The `sample_map` table in step 1 maps PEAKS internal sample numbers to patient IDs and conditions. Update the `peaks_num` values to match the `Sample X Area` column numbers in your `lfq.proteins.csv` files. The pattern used here:

- **Odd** sample numbers = Treated (T)
- **Even** sample numbers = Untreated (UT)

### 3. Download gene set files (step 4 only)

Download `.gmt` files from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) (free account required) and save to `gsea_data_dir`:

- Hallmarks: `h.all.vX.X.Hs.symbols.gmt`
- KEGG Medicus: `c2.cp.kegg_medicus.vX.X.Hs.symbols.gmt`
- GO Biological Process: `c5.go.bp.vX.X.Hs.symbols.gmt`

The script detects `.gmt` files automatically — no filename changes needed.

## Statistical model

Differential expression was tested using **limma** with `duplicateCorrelation()` to account for 5 biological replicates nested within each patient.

Four contrasts were tested:

| Contrast | Threshold | Rationale |
|----------|-----------|-----------|
| Treatment effect — PDAC | nominal p < 0.05, \|logFC\| > 0.5 | exploratory; n=3 patients |
| Treatment effect — HNSCC | nominal p < 0.05, \|logFC\| > 0.5 | exploratory; n=3 patients |
| HNSCC vs PDAC (baseline) | FDR < 0.05, \|logFC\| > 1 | confirmed |
| Interaction | FDR < 0.05, \|logFC\| > 1 | confirmed |

GSEA was run on all four contrasts using the full ranked protein list (t-statistic). The treatment contrasts are the primary output.

## Required R packages

```r
# CRAN
install.packages(c("dplyr", "tidyr", "stringr", "readr", "tibble", "forcats"))

# Bioconductor
BiocManager::install(c("limma", "fgsea"))
```
