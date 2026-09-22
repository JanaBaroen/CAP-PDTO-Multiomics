#!/bin/bash
#SBATCH --job-name=step5_filter
#SBATCH --output=/your/scratch/directory/step5_filter_%j.log	#CHANGE THIS: your scratch directory
#SBATCH --error=/your/scratch/directory/step5_filter_%j.err	#CHANGE THIS: your scratch directory
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=your_cluster_partition			#CHANGE THIS: your cluster partition name
#SBATCH --account=your_HPC_account				#CHANGE THIS: your HPC account name

module load R-bundle-Bioconductor/3.22-foss-2025a-R-4.5.1	#CHANGE THIS: your R module version

#CHANGE THIS: replace with the full path to germline_filter.R on your HPC system
Rscript /your/path/to/germline_filter.R
