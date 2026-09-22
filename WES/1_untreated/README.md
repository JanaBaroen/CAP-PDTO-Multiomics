# Untreated Sample WES Analysis

This directory contains all scripts for the untreated PDTO baseline mutation analysis, covering three stages: variant calling with nf-core/Sarek, a custom germline filtering pipeline, and downstream mutation analysis.

## Directory structure

```
1_untreated/
├── 1a_germline_filtering/
│   ├── filter_vcfs.sh          — Steps 1–3: BCFtools VAF/gnomAD/PON filtering
│   ├── vcf2maf.sh              — Step 4: VCF to MAF conversion
│   ├── run_step5_filter.sh     — SLURM wrapper for Step 5
│   └── germline_filter.R       — Step 5: consequence and hotspot filtering
└── 1b_baseline_analysis/
    └── WES_analysis_UT_baseline.R  — downstream TMB and gene list analysis
```

---

## Stage 1: Variant calling with nf-core/Sarek

Variant calling was performed using nf-core/Sarek v3.4.4 with Mutect2 in tumor-only mode (no matched normal available).

**Command:**
```bash
nextflow run nf-core/sarek -r 3.4.4 \
    -profile docker \
    -c ~/nextflow.config \
    --input /path/to/samplesheet.csv \
    --genome GATK.GRCh38 \
    --outdir /path/to/results \
    -w /path/to/sarek_work \
    --tools mutect2,strelka
```

**Resource parameters** (in `nextflow.config`):
```
params {
    max_cpus   = 12
    max_memory = '80.GB'
}
```

**Sample sheet format** (`samplesheet.csv`):
```
patient,sample,lane,fastq_1,fastq_2,sex,status
PDTO1,PDTO1_sample,L001,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz,NA,1
```

- `status = 1` indicates tumour sample (tumor-only mode)
- `lane` should match the lane identifier from your sequencing run (e.g. L006)
- `sex` set to NA if unknown
- Reference genome: GATK.GRCh38

The raw Mutect2 VCF output files (`.mutect2.filtered.vep.vcf`) are the input to the germline filtering pipeline below.

---

## Stage 2: Germline filtering pipeline

Scripts are in `1a_germline_filtering/`. Run in order: Steps 1–3 → Step 4 → Step 5.

### Steps 1–3: BCFtools filtering (`filter_vcfs.sh`)

Applies three sequential filters to remove germline and artifactual variants:

| Step | Filter | Tool |
|---|---|---|
| 1 | PASS filter + read depth (DP ≥ 10) + VAF filter | BCFtools filter |
| 2 | gnomAD population allele frequency (AF > 0.001) | BCFtools annotate |
| 3 | Panel of Normals (1000 Genomes PON) | BCFtools isec |

**VAF filter note:** variants with AF < 0.40 or AF > 0.65 are retained. The germline heterozygous range (0.40–0.65) is excluded to reduce germline contamination in tumor-only data. An initial upper exclusion band (AF ≥ 0.90) was removed after systematic investigation showed it incorrectly excluded loss-of-heterozygosity-consistent somatic mutations at tumour suppressor loci including TP53 and CDKN2A.

**Reference files required:**
- gnomAD af-only resource: `af-only-gnomad.hg38.vcf.gz` — available from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
- Panel of Normals: `1000g_pon.hg38.vcf.gz` — available from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
- GRCh38 reference FASTA: `GRCh38_no_alt.fa`

**Software required:**
- BCFtools v1.22

**Before running:** update all paths and SLURM settings marked `# CHANGE THIS` at the top of the script.

```bash
sbatch filter_vcfs.sh
```

---

### Step 4: VCF to MAF conversion (`vcf2maf.sh`)

Converts filtered VCFs to MAF format using vcf2maf with VEP re-annotation.

**Software required:**
- vcf2maf v1.6.22
- VEP v115.2
- BCFtools v1.22
- SAMtools v1.22.1

**Before running:** update all paths and SLURM settings marked `# CHANGE THIS` at the top of the script.

```bash
sbatch vcf2maf.sh
```

---

### Step 5: Consequence and hotspot filtering (`germline_filter.R`)

R-based filtering that:
- Retains only protein-coding consequences (missense, nonsense, frameshift, splice site, in-frame indels, nonstop, translation start site)
- Applies a curated driver gene hotspot whitelist to protect known cancer driver positions from downstream gnomAD filtering
- Applies a final gnomAD safety check (AF ≤ 0.001)

Run via SLURM:
```bash
sbatch run_step5_filter.sh
```

Or directly in R:
```r
source("germline_filter.R")
```

**Software required:**
- R v4.5.1
- R packages: data.table, maftools

**Known limitation:** KMT2D was identified as recurrently mutated across all three HNSCC samples (VAF 0.48–0.563, nonsense and missense variants with COSMIC annotation) but is excluded by the mid-band VAF filter (0.40–0.65). This is an inherent limitation of tumor-only WES analysis where clonal somatic mutations at moderate tumour purity are indistinguishable from germline heterozygous variants by VAF alone in the absence of a matched normal. KMT2D is included in the hotspot whitelist in `germline_filter.R` but cannot be rescued at Step 5 since it is removed upstream in Step 1.

**Output:** per-sample and combined final MAF files, used as input to the downstream analysis script.

---

## Stage 3: Downstream baseline analysis (`WES_analysis_UT_baseline.R`)

Script is in `1b_baseline_analysis/`.

Performs downstream analysis on the final filtered MAF files:
- Tumor mutational burden (TMB; 38 Mb exome capture size)
- Recurrently mutated gene lists per cancer type (HNSCC and PDAC)

**Software required:**
- R v4.5.1
- R packages: maftools, ggplot2, dplyr

**Before running:** update the two paths at lines 27–28 marked `# CHANGE THIS`:
- `MAF_DIR_UT` — directory containing final filtered MAF files from Step 5
- `output_dir` — directory where output files will be saved

**Output files used in the study:**
- `tmb_per_sample_UT_filtered.csv` — TMB values per sample
- `mutated_genes_HNSCC_UT_filtered.txt` — recurrently mutated genes in HNSCC PDTOs
- `mutated_genes_PDAC_UT_filtered.txt` — recurrently mutated genes in PDAC PDTOs

---

## Software versions

| Software | Version |
|---|---|
| nf-core/Sarek | 3.4.4 |
| BCFtools | 1.22 |
| SAMtools | 1.22.1 |
| VEP | 115.2 |
| vcf2maf | 1.6.22 |
| R | 4.5.1 |
| maftools | via R-bundle-Bioconductor 3.22 |

## Reference files

| Resource | Build/Version |
|---|---|
| Reference genome | GRCh38 (GATK.GRCh38) |
| gnomAD | af-only-gnomad.hg38.vcf.gz |
| Panel of Normals | 1000g_pon.hg38.vcf.gz |
| VEP cache | GRCh38 v115 |
