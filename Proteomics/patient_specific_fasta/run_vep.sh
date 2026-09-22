#!/bin/bash
# run_vep.sh
# ==========
# Submits one SLURM job per patient to run VEP annotation in parallel.
# Adjust module names, paths, and SLURM parameters for your HPC environment.
#
# Usage:
#   bash run_vep.sh

# Update these paths to match your environment
VCF_IN=/path/to/vcfs/filtered			# CHANGE THIS: directory containing filtered VCFs
VCF_OUT=/path/to/vep_annotated			# CHANGE THIS: output directory for VEP-annotated VCFs
CACHE=/path/to/vep_cache				# CHANGE THIS: path to VEP cache directory
SCRATCH=/path/to/scratch				# CHANGE THIS: your scratch/working directory

mkdir -p $VCF_OUT

for PATIENT in PDTO1 PDTO2 PDTO3 PDTO4 PDTO5 PDTO6; do	# CHANGE THIS: list your patient IDs

cat > ${SCRATCH}/vep_${PATIENT}.sh << EOF
#!/bin/bash
#SBATCH --job-name=vep_${PATIENT}
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --output=${SCRATCH}/vep_${PATIENT}_%j.out

# Adjust module load command for your HPC environment
module load VEP/115					# CHANGE THIS: adjust VEP module version for your HPC

echo "Annotating ${PATIENT}..."
vep \\
    --input_file  ${VCF_IN}/${PATIENT}.final_filtered.vcf.gz \\
    --output_file ${VCF_OUT}/${PATIENT}.vep.vcf \\
    --vcf \\
    --cache \\
    --dir_cache ${CACHE} \\
    --assembly GRCh38 \\
    --species homo_sapiens \\
    --symbol \\
    --canonical \\
    --protein \\
    --hgvs \\
    --hgvsp \\
    --sift b \\
    --polyphen b \\
    --offline \\
    --force_overwrite \\
    --fork 4
echo "Done: ${PATIENT}"
EOF

sbatch ${SCRATCH}/vep_${PATIENT}.sh
echo "Submitted job for $PATIENT"

done
