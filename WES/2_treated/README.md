# Treated Sample WES Analysis

This directory contains all scripts for the CAP-treated PDTO mutation analysis, covering two stages: paired variant calling with nf-core/Sarek (treated tumour vs. untreated tumour as reference) and downstream mutation analysis.

## Directory structure

```
2_treated/
├── run_vep_annotation.sh       — VCF to MAF conversion with VEP annotation
└── WES_analysis.R              — downstream TMB, gene list and Reactome analysis
```

---

## Stage 1: Variant calling with nf-core/Sarek

Variant calling was performed using nf-core/Sarek v3.4.4 with Mutect2 in paired mode, using the untreated tumour sample as the reference (normal). This design means only variants present in the treated sample but absent in the untreated sample are called, removing the need for a separate germline filtering pipeline.

**Command:**
```bash
nextflow run nf-core/sarek \
    -profile docker \
    -c /path/to/nextflow.config \
    --input /path/to/samplesheet.csv \
    --genome GATK.GRCh38 \
    --outdir /path/to/results \
    -w /path/to/sarek_work \
    --tools mutect2,strelka \
```

**Resource parameters** (in `nextflow.config`):
```
params {
    max_cpus   = 12
    max_memory = '80.GB'
}
```

**Sample sheet format** (`samplesheet.csv`) — one row per sample per patient, with treated (status=1) and untreated (status=0) listed together:
```
patient,sample,sex,status,fastq_1,fastq_2
PDTO1,PDTO1_treated,NA,1,/path/to/treated_R1.fastq.gz,/path/to/treated_R2.fastq.gz
PDTO1,PDTO1_untreated,NA,0,/path/to/untreated_R1.fastq.gz,/path/to/untreated_R2.fastq.gz
```

- `status = 1` indicates the treated tumour sample
- `status = 0` indicates the untreated tumour sample (used as reference/normal)
- One samplesheet per patient, run separately
- Reference genome: GATK.GRCh38

The raw Mutect2 VCF output files (`.mutect2.pass.vcf.gz`) are the input to the VCF to MAF conversion step below.

---

## Stage 2: VCF to MAF conversion (`run_vep_annotation.sh`)

Converts paired Mutect2 VCFs to MAF format using vcf2maf with VEP annotation.

**Software required:**
- vcf2maf v1.6.22
- VEP v115.2
- SAMtools v1.22.1

**Before running:** update all paths and SLURM settings marked `# CHANGE THIS` at the top of the script.

```bash
sbatch run_vep_annotation.sh
```

**Note on sample pairs:** the script processes treated/untreated pairs together. The `--tumor-id` corresponds to the treated sample and `--normal-id` corresponds to the untreated sample, matching the paired Sarek calling design.

**Output:** per-pair MAF files (one per patient) used as input to the downstream analysis script.

---

## Stage 3: Downstream analysis (`WES_analysis.R`)

Script performs downstream analysis on the paired MAF files:
- Tumor mutational burden (TMB; 38 Mb exome capture size)
- Recurrently mutated gene lists per cancer type (HNSCC and PDAC)
- Reactome pathway enrichment analysis
- Gene summaries per cancer type and combined

**Software required:**
- R v4.5.1
- R packages: maftools, clusterProfiler, ReactomePA, org.Hs.eg.db, enrichplot, ggplot2, dplyr, tidyr, RColorBrewer, ggpubr

**Before running:** update the two paths at the top of the script marked `# CHANGE THIS`:
- `MAF_DIR` — directory containing MAF files from `run_vep_annotation.sh`
- `output_dir` — directory where output files will be saved

**Output files used in the study:**
- `tmb_per_sample.csv` — TMB values per sample
- `gene_summary_all.csv` — gene-level mutation summary across all samples
- `gene_summary_hnscc.csv` — gene-level mutation summary for HNSCC samples
- `gene_summary_pdac.csv` — gene-level mutation summary for PDAC samples
- `mutated_genes_HNSCC.txt` — list of mutated genes in HNSCC PDTOs
- `mutated_genes_PDAC.txt` — list of mutated genes in PDAC PDTOs
- `Reactome_enrichment_HNSCC.csv` — Reactome pathway enrichment results for HNSCC
- `Reactome_enrichment_PDAC.csv` — Reactome pathway enrichment results for PDAC

---

## Software versions

| Software | Version |
|---|---|
| nf-core/Sarek | 3.4.4 |
| VEP | 115.2 |
| vcf2maf | 1.6.22 |
| SAMtools | 1.22.1 |
| R | 4.5.1 |
| maftools | via R-bundle-Bioconductor 3.22 |

## Reference files

| Resource | Build/Version |
|---|---|
| Reference genome | GRCh38 (GATK.GRCh38) |
| VEP cache | GRCh38 v115 |
