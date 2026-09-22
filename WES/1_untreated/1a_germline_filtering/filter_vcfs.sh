#!/bin/bash
#SBATCH --job-name=germline_filter
#SBATCH --output=/your/scratch/directory/germline_filter_%j.log		# CHANGE THIS: your scratch directory
#SBATCH --error=/your/scratch/directory/germline_filter_%j.err		# CHANGE THIS: your scratch directory
#SBATCH --time=24:00:00
#SBATCH --mem=72G
#SBATCH --cpus-per-task=8
#SBATCH --partition=your_cluster_partition				# CHANGE THIS: your cluster partition name
#SBATCH --account=your_HPC_account					# CHANGE THIS: your HPC account name

module load BCFtools/1.22-GCC-14.2.0

# CHANGE THIS: your root scratch/working directory
SCRATCH="/your/scratch/directory"

# CHANGE THIS: directory containing raw Mutext2 VEP-annotated VCFs
VCF_DIR="$SCRATCH/WES_UT"

# CHANGE THIS: output directory for filtered VCFs
OUT_DIR="$SCRATCH/WES_vcfs/filtered"

# CHANGE THIS: path to gnomAD af-only resource file (GATK hg38 version)
GNOMAD="$SCRATCH/reference_files/af-only-gnomad.hg38.vcf.gz"

# CHANGE THIS: path to GATK 1000 Genomes Panel of Normals (hg38)
PON="$SCRATCH/reference_files/1000g_pon.hg38.vcf.gz"

mkdir -p "$OUT_DIR"

declare -A SAMPLES=(
    ["YOUR_RUN_ID_1"]="PDTO4"
    ["YOUR_RUN_ID_2"]="PDTO5"
    ["YOUR_RUN_ID_3"]="PDTO6"
    ["YOUR_RUN_ID_4"]="PDTO1"
    ["YOUR_RUN_ID_5"]="PDTO2"
    ["YOUR_RUN_ID_6"]="PDTO3"
)

for GC_ID in "${!SAMPLES[@]}"; do
    PATIENT="${SAMPLES[$GC_ID]}"
    VCF="$VCF_DIR/${GC_ID}.mutect2.filtered.vep.vcf"

    echo "=============================="
    echo "Processing: $GC_ID ($PATIENT)"
    echo "=============================="

    if [ ! -f "$VCF" ]; then
        echo "ERROR: VCF not found: $VCF"
        continue
    fi

    # Compress and index if needed
    if [ ! -f "${VCF}.gz" ]; then
        echo "Compressing $GC_ID..."
        bgzip -c "$VCF" > "${VCF}.gz"
        bcftools index "${VCF}.gz"
    fi

    BEFORE=$(bcftools view -H "${VCF}.gz" | wc -l)
    echo "  Variants before filtering: $BEFORE"

    # --- Step 1: PASS filter + depth + VAF filter ---
    echo "Applying PASS + depth + VAF filter..."
    bcftools view -f PASS "${VCF}.gz" | \
    
# VAF filter: excludes germline heterozygous range (0.40-0.65)
# NOTE: upper bound (>=0.90) removed to retain LOH consistent somatic mutations
# See README for pipeline details
    bcftools filter \
        --include 'FORMAT/DP>=10 && (FORMAT/AF<0.40 || FORMAT/AF>0.65)' \
        -Oz \
        -o "$OUT_DIR/${GC_ID}.pass_vaf_filtered.vcf.gz"

    bcftools index "$OUT_DIR/${GC_ID}.pass_vaf_filtered.vcf.gz"

    AFTER_PASS=$(bcftools view -H "$OUT_DIR/${GC_ID}.pass_vaf_filtered.vcf.gz" | wc -l)
    echo "  Variants after PASS+DP+VAF filter: $AFTER_PASS"

    # --- Step 2: Annotate with gnomAD genome AF ---
    echo "Annotating with gnomAD AF..."
    bcftools annotate \
        --annotations "$GNOMAD" \
        --columns "INFO/AF" \
        --output-type z \
        --output "$OUT_DIR/${GC_ID}.gnomad_annotated.vcf.gz" \
        "$OUT_DIR/${GC_ID}.pass_vaf_filtered.vcf.gz"

    bcftools index "$OUT_DIR/${GC_ID}.gnomad_annotated.vcf.gz"

    # --- Step 3: gnomAD AF filter ---
    echo "Applying gnomAD AF filter..."
    bcftools filter \
        --include 'INFO/AF[0] <= 0.001 || INFO/AF[0] = "."' \
        --output-type z \
        --output "$OUT_DIR/${GC_ID}.gnomad_filtered.vcf.gz" \
        "$OUT_DIR/${GC_ID}.gnomad_annotated.vcf.gz"

    bcftools index "$OUT_DIR/${GC_ID}.gnomad_filtered.vcf.gz"

    AFTER_GNOMAD=$(bcftools view -H "$OUT_DIR/${GC_ID}.gnomad_filtered.vcf.gz" | wc -l)
    echo "  Variants after gnomAD filter: $AFTER_GNOMAD"

    # --- Step 4: PON filter ---
    echo "Applying PON filter..."
    bcftools isec \
        "$OUT_DIR/${GC_ID}.gnomad_filtered.vcf.gz" \
        "$PON" \
        --complement \
        --write 1 \
        --output-type z \
        --output "$OUT_DIR/${GC_ID}.pon_filtered.vcf.gz"

    bcftools index "$OUT_DIR/${GC_ID}.pon_filtered.vcf.gz"

    AFTER_PON=$(bcftools view -H "$OUT_DIR/${GC_ID}.pon_filtered.vcf.gz" | wc -l)
    echo "  Variants after PON filter: $AFTER_PON"

    # --- Save with patient ID name ---
    cp "$OUT_DIR/${GC_ID}.pon_filtered.vcf.gz" \
       "$OUT_DIR/${PATIENT}.final_filtered.vcf.gz"
    bcftools index "$OUT_DIR/${PATIENT}.final_filtered.vcf.gz"

    # --- Clean up intermediate files ---
    rm "$OUT_DIR/${GC_ID}.pass_vaf_filtered.vcf.gz" \
       "$OUT_DIR/${GC_ID}.pass_vaf_filtered.vcf.gz.csi" \
       "$OUT_DIR/${GC_ID}.gnomad_annotated.vcf.gz" \
       "$OUT_DIR/${GC_ID}.gnomad_annotated.vcf.gz.csi" \
       "$OUT_DIR/${GC_ID}.gnomad_filtered.vcf.gz" \
       "$OUT_DIR/${GC_ID}.gnomad_filtered.vcf.gz.csi"

    echo "Done: $PATIENT"
    echo "  Total reduction: $BEFORE -> $AFTER_PON variants"
    echo ""
done

echo "=============================="
echo "All samples complete"
ls -lh "$OUT_DIR"/*.final_filtered.vcf.gz
echo "=============================="
