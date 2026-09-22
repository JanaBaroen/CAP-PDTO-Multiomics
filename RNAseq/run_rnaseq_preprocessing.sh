#!/bin/bash
#SBATCH --job-name=RNAseq_preprocessing
#SBATCH --account=your_HPC_account				#CHANGE THIS: your HPC account name
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1 --cpus-per-task=16
#SBATCH --mem=40gb
#SBATCH --output=/your/scratch/directory/rnaseq-%j.out	#CHANGE THIS: your scratch directory
#SBATCH --error=/your/scratch/directory/rnaseq-%j.err	#CHANGE THIS: your scratch directory

# NOTE: Module names below are specific to the CalcUA HPC cluster (University of Antwerp).
# Replace with the equivalent modules available on your HPC system.
module purge
module load calcua/all
module load Miniconda3/23.5.2-0
module load R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1
module load Salmon/1.10.3-GCC-13.3.0
module load MultiQC/1.28-gfbf-2024a

source activate trimgalore

cd /your/path/to/RNAseq_project         # CHANGE THIS: your project directory
Rscript /your/path/to/rnaseq_pipeline.R  # CHANGE THIS: path to the R script