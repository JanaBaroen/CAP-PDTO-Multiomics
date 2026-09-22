#!/bin/bash
#SBATCH --job-name=vcf2maf
#SBATCH --output=/your/scratch/directory/vcf2maf_%j.log		# CHANGE THIS: your scratch directory
#SBATCH --error=/your/scratch/directory/vcf2maf_%j.err		# CHANGE THIS: your scratch directory
#SBATCH --time=24:00:00
#SBATCH --mem=72G
#SBATCH --cpus-per-task=8
#SBATCH --partition=your_cluster_partition			# CHANGE THIS: your cluster partition name
#SBATCH --account=your_HPC_account				# CHANGE THIS: your HPC account name

module load VEP/115.2-GCC-13.3.0    # CHANGE THIS: your VEP module version
module load SAMtools/1.22.1-GCC-14.2.0    # CHANGE THIS: your SAMtools module version
module load BCFtools/1.22-GCC-14.2.0    # CHANGE THIS: your BCFtools module version

# CHANGE THIS: path to vcf2maf.pl in your software environment
VCF2MAF="/your/path/to/vcf2maf.pl"

# CHANGE THIS: path to your VEP installation
VEP_PATH="/your/path/to/VEP"

# CHANGE THIS: path to your VEP cache directory
VEP_DATA="/your/path/to/VEP_cache"

# CHANGE THIS: your root scratch/working directory
SCRATCH="/your/scratch/directory"

# CHANGE THIS: directory containing filtered VCFs from filter_vcfs.sh
VCF_DIR="$SCRATCH/WES_vcfs/filtered"

#CHANGE THIS: output directory for MAF files
MAF_DIR="$SCRATCH/MAF_files/filtered"

# CHANGE THIS: path to GRCh38 reference FASTA
REF="$SCRATCH/GRCh38_no_alt.fa"

mkdir -p "$MAF_DIR"

#sample IDs (corresponding to associated manuscript)
PATIENTS=("PDTO1" "PDTO2" "PDTO3" "PDTO4" "PDTO5" "PDTO6")

for PATIENT in "${PATIENTS[@]}"; do
    VCF_GZ="$VCF_DIR/${PATIENT}.final_filtered.vcf.gz"
    VCF="$VCF_DIR/${PATIENT}.final_filtered.vcf"
    MAF="$MAF_DIR/${PATIENT}.filtered.maf"

    echo "Converting $PATIENT to MAF..."

    if [ ! -f "$VCF_GZ" ]; then
        echo "ERROR: VCF not found for $PATIENT"
        continue
    fi

    # Decompress VCF for vcf2maf
    echo "Decompressing $PATIENT..."
    bcftools view "$VCF_GZ" -O v -o "$VCF"

    # Convert to MAF
    perl "$VCF2MAF" \
        --input-vcf "$VCF" \
        --output-maf "$MAF" \
        --tumor-id "$PATIENT" \
        --vep-path "$VEP_PATH" \
        --vep-data "$VEP_DATA" \
        --ref-fasta "$REF" \
        --ncbi-build GRCh38 \
        --inhibit-vep

    # Remove decompressed VCF to save space
    rm "$VCF"

    echo "Done: $PATIENT -> $MAF"
done

echo "All MAF conversions complete."
