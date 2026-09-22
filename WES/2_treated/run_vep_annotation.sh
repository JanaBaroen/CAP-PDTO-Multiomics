#!/bin/bash
#SBATCH --job-name=vep_annotation
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=/your/scratch/directory/vep_annotation_%j.log		# CHANGE THIS: your scratch directory
#SBATCH --error=/your/scratch/directory/vep_annotation_%j.err		# CHANGE THIS: your scratch directory
#SBATCH --partition=your_cluster_partition				# CHANGE THIS: your cluster partition name
#SBATCH --account=your_HPC_account					# CHANGE THIS: your HPC account name

module load calcua/2024a						# CHANGE THIS: your cluster module environment
module load VEP/115.2-GCC-13.3.0					# CHANGE THIS: your VEP module version

VCF2MAF=/your/path/to/vcf2maf.pl                    			# CHANGE THIS: path to vcf2maf.pl
SAMTOOLS=/your/path/to/samtools                      			# CHANGE THIS: path to samtools
VCF_DIR=/your/path/to/WES_data                       			# CHANGE THIS: directory containing paired Mutect2 VCFs
OUT_DIR=/your/path/to/MAF_files                      			# CHANGE THIS: output directory for MAF files
CACHE_DIR=/your/path/to/vep_cache                    			# CHANGE THIS: path to VEP cache directory
REF=/your/path/to/GRCh38_no_alt.fa                  			# CHANGE THIS: path to GRCh38 reference FASTA
VEP_PATH=/your/path/to/VEP                           			# CHANGE THIS: path to VEP installation

mkdir -p $OUT_DIR

# Sample pairs: treated tumour vs untreated tumour as reference
# Key: treated_vs_untreated, Value: "tumor_ID normal_ID"
# Replace PDTO names with your sample identifiers
declare -A PAIRS
PAIRS["PDTO1_treated_vs_PDTO1_untreated"]="PDTO1_treated PDTO1_untreated"
PAIRS["PDTO2_treated_vs_PDTO2_untreated"]="PDTO2_treated PDTO2_untreated"
PAIRS["PDTO3_treated_vs_PDTO3_untreated"]="PDTO3_treated PDTO3_untreated"
PAIRS["PDTO4_treated_vs_PDTO4_untreated"]="PDTO4_treated PDTO4_untreated"
PAIRS["PDTO5_treated_vs_PDTO5_untreated"]="PDTO5_treated PDTO5_untreated"
PAIRS["PDTO6_treated_vs_PDTO6_untreated"]="PDTO6_treated PDTO6_untreated"

for PATIENT in "${!PAIRS[@]}"; do
  read TUMOR NORMAL <<< "${PAIRS[$PATIENT]}"

  echo "Processing $PATIENT (tumor: $TUMOR, normal: $NORMAL)..."

  # Decompress VCF
  gunzip -c $VCF_DIR/${PATIENT}.mutect2.pass.vcf.gz \
    > $VCF_DIR/${PATIENT}.mutect2.pass.vcf

  # Run vcf2maf using full path
  perl $VCF2MAF \
    --input-vcf $VCF_DIR/${PATIENT}.mutect2.pass.vcf \
    --output-maf $OUT_DIR/${PATIENT}.maf \
    --tumor-id $TUMOR \
    --normal-id $NORMAL \
    --vep-path $VEP_PATH \
    --vep-data $CACHE_DIR \
    --vep-forks 4 \
    --ref-fasta $REF \
    --ncbi-build GRCh38 \
    --samtools-exec $SAMTOOLS \
    --verbose

  # Remove temporary uncompressed VCF
  rm $VCF_DIR/${PATIENT}.mutect2.pass.vcf

  echo "Done $PATIENT"
done

echo "All samples annotated!"
