# RNAseq Analysis

This directory contains all scripts for the RNAseq analysis of untreated and CAP-treated patient-derived tumour organoids (PDTOs) from pancreatic ductal adenocarcinoma (PDAC; PDTO1–3) and head and neck squamous cell carcinoma (HNSCC; PDTO4–6), as described in the associated manuscript.

## Study design

Six patients (3 HNSCC + 3 PDAC) each provided paired untreated and CAP-treated tumour samples with 3 technical replicates per condition (36 samples total). Technical replicates were collapsed by summing counts prior to analysis, yielding 12 samples (6 patients × 2 conditions).

| | Details |
|---|---|
| **Cancer types** | HNSCC (PDTO4–6) and PDAC (PDTO1–3) |
| **Design** | Paired: each patient contributes both untreated and treated samples |
| **Technical replicates** | 3 per patient per condition (collapsed before analysis) |
| **Sequencing** | Bulk RNAseq, paired-end |
| **Reference genome** | GRCh38 / GENCODE v49 |

## Directory structure

```
RNAseq/
├── run_rnaseq_preprocessing.sh         — SLURM job script for preprocessing
├── rnaseq_pipeline.R                   — Stage 1: preprocessing FASTQ → count matrix
├── rnaseq_analysis_pipeline.R          — Stage 2: DESeq2 treatment effect analysis
├── rnaseq_analysis_pipeline_GSEA.R     — Stage 3: Hallmark GSEA
├── rnaseq_analysis_pipeline_heatmap.R  — Stage 4: manuscript heatmap figures
└── rnaseq_treatment_per_cancer_type.R  — Stage 5: per-cancer-type analysis
```

---

## Stage 1: Preprocessing — FASTQ to count matrix


### SLURM job (`run_rnaseq_preprocessing.sh`)

**Before running:** update the two paths marked `# CHANGE THIS` at the bottom of the script.

```bash
sbatch run_rnaseq_preprocessing.sh
```

**Modules loaded:**
| Module | Version |
|---|---|
| Miniconda3 | 23.5.2-0 |
| R-bundle-Bioconductor | 3.22-foss-2025a-R-4.5.1 |
| Salmon | 1.10.3-GCC-13.3.0 |
| MultiQC | 1.28-gfbf-2024a |
| trim-galore (conda) | via trimgalore environment |

### R preprocessing pipeline (`rnaseq_pipeline.R`)

Runs all preprocessing steps in sequence:

| Step | Description | Tool |
|---|---|---|
| 1 | Create output directories | R |
| 2 | Load sample sheet (`samples.csv`) | R |
| 3 | Merge technical sequencing lanes per sample | bash via R |
| 4 | Quality control on merged FASTQ files | FastQC + MultiQC |
| 5 | Adapter and quality trimming | Trim Galore |
| 6 | Build Salmon index (skipped if already exists) | Salmon |
| 7 | Quantify transcript expression | Salmon |
| 8 | Build transcript-to-gene map from GTF | rtracklayer |
| 9 | Summarise transcript counts to gene level | tximport |
| 10 | Save output matrices | R |

**Before running:** update `external_drive` at line 10 marked `# CHANGE THIS` to your working directory.

**Input folder structure required:**
```
RNAseq/
├── PDTO1/
│   ├── UT/    ← untreated FASTQ files (one or more lanes per sample)
│   └── T/     ← treated FASTQ files
├── PDTO2/
│   ├── UT/
│   └── T/
└── ... etc.
```

**Sample sheet format** (`samples.csv`) — one row per sample (not per lane):
```
sample_id,Patient,Type,Condition
SAMPLE01,PDTO1,PDAC,Untreated
SAMPLE02,PDTO1,PDAC,Untreated
SAMPLE03,PDTO1,PDAC,Untreated
SAMPLE04,PDTO1,PDAC,Treated
SAMPLE05,PDTO1,PDAC,Treated
SAMPLE06,PDTO1,PDAC,Treated
SAMPLE07,PDTO2,PDAC,Untreated
...
```

- `sample_id` — unique sample identifier, must match the FASTQ filename prefix
- `Patient` — patient identifier, must match the folder name under `RNAseq/`
- `Type` — cancer type: `HNSCC` or `PDAC`
- `Condition` — treatment status: `Untreated` (reference level) or `Treated`
- 3 rows per patient per condition (one per technical replicate)

**Reference files required:**
- GENCODE v49 transcript FASTA: `gencode.v49.transcripts.fa` — available from [GENCODE](https://www.gencodegenes.org/human/release_49.html)
- GENCODE v49 annotation GTF: `gencode.v49.primary_assembly.annotation.gtf` — available from [GENCODE](https://www.gencodegenes.org/human/release_49.html)
- Place both files in a `RNAseq_ref/` subfolder under your working directory

**Output files** (saved to `results/`):
| File | Description |
|---|---|
| `exp.mat.rds` | Raw integer count matrix, genes × samples (R format) |
| `exp.mat.tsv` | Raw integer count matrix (tab-separated) |
| `pdat.mat.rds` | Sample metadata, Patient/Type/Condition columns (R format) |
| `pdat.mat.tsv` | Sample metadata (tab-separated) |


---

## Stage 2–5: Downstream analysis

Stages 2–4 must be run **in order within the same R session** as objects are shared between scripts. Stage 5 can be run independently.

### Stage 2: DESeq2 treatment effect (`rnaseq_analysis_pipeline.R`)

Loads `exp.mat.rds` and `pdat.mat.rds`, collapses technical replicates, fits DESeq2 model (`~Patient + Condition`), and annotates DE results.

**Before running:** update `BASE_PATH` at the top of the script to your working directory.

**Output files used in the manuscript:**
- `pca_plot.pdf/.png` — PCA of all 12 collapsed samples
- `volcano_treatment_all_significant_genes.pdf/.png` — volcano plot, treatment effect
- `de_results_treatment_annotated.tsv` — full annotated DE results

### Stage 3: GSEA (`rnaseq_analysis_pipeline_GSEA.R`)

Runs Hallmark GSEA on the treatment effect. Loads DE results from Stage 2. **Must be run before Stage 4** as it produces objects (`res.fgsea.hallmark`, `pathways.hallmark`) required by the heatmap script.

**Gene set files required:**
Download the Hallmark collection from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb) (free account required). Select **H collection → Gene Symbols**.

```
h.all.vX.X.Hs.symbols.gmt
```

Place the file in your working directory and update `GMT_HALLMARK` at the top of the script to match the exact filename.

**Output files used in the manuscript:**
- `gsea_hallmark_full.tsv` — full Hallmark GSEA results

### Stage 4: Manuscript heatmaps (`rnaseq_analysis_pipeline_heatmap.R`)

Generates the two manuscript heatmap figures using ComplexHeatmap. Requires objects in memory from Stages 2 and 3.

**Output files used in the manuscript:**
- `heatmap_curated_main_figure.pdf/.png` — curated gene heatmap with pathway annotation (main figure)
- `supplementary_heatmap_all_significant_DE_genes.pdf/.png` — all significant DE genes (supplementary figure)

### Stage 5: Per-cancer-type analysis (`rnaseq_treatment_per_cancer_type.R`)

Runs DESeq2 and Hallmark GSEA separately within HNSCC and PDAC. Can be run independently from Stages 2–4.

**Before running:** update `BASE_PATH` and `GMT_HALLMARK` at the top of the script.

**Output files used in the manuscript:**
- `gsea_treatment_HNSCC_hallmark_full.tsv` — full Hallmark GSEA results for HNSCC
- `gsea_treatment_PDAC_hallmark_full.tsv` — full Hallmark GSEA results for PDAC

---

## How to run

1. Download raw data from EGA and organise into the folder structure described above
2. Download GENCODE v49 reference files and place in `RNAseq_ref/`
3. Set up the `trimgalore` conda environment: `conda create -n trimgalore -c bioconda trim-galore fastqc`
4. Update paths in `run_rnaseq_preprocessing.sh` and `rnaseq_pipeline.R`
5. Submit preprocessing: `sbatch run_rnaseq_preprocessing.sh`
6. Once `exp.mat.rds` and `pdat.mat.rds` are produced, open R or RStudio
7. Set `BASE_PATH` at the top of each analysis script to your working directory
8. Run Stages 2–4 in order within the same R session
9. Run Stage 5 independently if per-cancer-type results are needed

---

## Statistical approach

| Analysis | Model | Contrast |
|---|---|---|
| Treatment effect (all patients) | `~Patient + Condition` | Treated vs Untreated |
| Treatment effect per cancer type | `~Patient + Condition` | Treated vs Untreated (within HNSCC or PDAC) |

- Patient is included as a blocking factor to account for the paired design and the confounding of patient with cancer type
- Technical replicates are collapsed by summing counts using `DESeq2::collapseReplicates()` before fitting any model
- Significance threshold: padj < 0.05 (Benjamini-Hochberg correction)
- GSEA ranking metric: log2FoldChange × −log10(p-value)
- GSEA significance: padj < 0.05 reported as significant; pval < 0.10 reported as trending

---

## Software versions

| Software | Version |
|---|---|
| R | 4.5.1 |
| DESeq2 | via R-bundle-Bioconductor 3.22 |
| fgsea | via R-bundle-Bioconductor 3.22 |
| ComplexHeatmap | via R-bundle-Bioconductor 3.22 |
| tximport | via R-bundle-Bioconductor 3.22 |
| rtracklayer | via R-bundle-Bioconductor 3.22 |
| org.Hs.eg.db | via R-bundle-Bioconductor 3.22 |
| Salmon | 1.10.3 |
| Trim Galore | via conda |
| FastQC | via conda |
| MultiQC | 1.28 |

## Reference files

| Resource | Version |
|---|---|
| Reference genome | GRCh38 |
| Gene annotation | GENCODE v49 |
| Hallmark gene sets | MSigDB v2026.1 |
