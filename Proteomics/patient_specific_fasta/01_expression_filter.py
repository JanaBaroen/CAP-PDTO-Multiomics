#!/usr/bin/env python3
"""
01_expression_filter.py
=======================
Filters expressed genes per patient from RNA-seq count data.

Reads exp.mat.tsv and pdat.mat.tsv, averages counts across untreated
replicates per patient, and outputs a per-patient list of expressed genes.

Usage:
    python 01_expression_filter.py \
        --exp        /path/to/exp.mat.tsv \
        --pdat       /path/to/pdat.mat.tsv \
        --outdir     /path/to/expression_filtered \
        --min_counts 10

Input:
    --exp        Gene x sample count matrix (from tximport/Salmon)
                 Rows = ENSEMBL gene IDs, columns = sample IDs
    --pdat       Sample metadata TSV with columns:
                 sample_id, Patient, Type, Condition
    --outdir     Output directory for per-patient gene lists
    --min_counts Minimum mean count across replicates to consider
                 a gene expressed (default: 10)

Output (per patient, in outdir):
    <PATIENT>_expressed_genes.txt   One ENSEMBL gene ID per line
    <PATIENT>_mean_counts.tsv       Full mean count table
    expression_summary.tsv          Summary across all patients

Notes:
    - PATIENT_MAP below maps anonymized patient IDs to the Patient
      column values in your pdat.mat file. Update this mapping to
      match your own sample sheet.
    - ENSEMBL version numbers are stripped automatically
      (e.g. ENSG00000000003.17 -> ENSG00000000003)
"""

import argparse
import os
import pandas as pd

###############################################################################
# PATIENT ID MAPPING
# Maps anonymized patient IDs to Patient column values in pdat.mat
# Update this mapping to match your own sample sheet
###############################################################################
PATIENT_MAP = {  # CHANGE THIS: update to match your patient IDs and sample sheet
    "PDTO1": "PDTO1",
    "PDTO2": "PDTO2",
    "PDTO3": "PDTO3",
    "PDTO4": "PDTO4",
    "PDTO5": "PDTO5",
    "PDTO6": "PDTO6",
}

###############################################################################
# ARGUMENT PARSING
###############################################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter expressed genes per patient from RNA-seq counts."
    )
    parser.add_argument("--exp",        required=True,
                        help="Path to exp.mat.tsv (gene x sample count matrix)")
    parser.add_argument("--pdat",       required=True,
                        help="Path to pdat.mat.tsv (sample metadata)")
    parser.add_argument("--outdir",     required=True,
                        help="Output directory for per-patient gene lists")
    parser.add_argument("--min_counts", type=float, default=10,
                        help="Minimum mean count to consider a gene expressed (default: 10)")
    return parser.parse_args()

###############################################################################
# MAIN
###############################################################################
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading expression matrix: {args.exp}")
    exp = pd.read_csv(args.exp, sep="\t", index_col=0)
    print(f"  -> {exp.shape[0]} genes x {exp.shape[1]} samples")

    print(f"Loading sample metadata: {args.pdat}")
    pdat = pd.read_csv(args.pdat, sep="\t")
    print(f"  -> {pdat.shape[0]} samples")

    # Keep only Untreated samples
    pdat_ut = pdat[pdat["Condition"] == "Untreated"].copy()
    print(f"  -> {pdat_ut.shape[0]} untreated samples retained")

    # Strip ENSEMBL version numbers (e.g. ENSG00000000003.17 -> ENSG00000000003)
    exp.index = exp.index.str.split(".").str[0]

    summary_rows = []

    for short_id, long_id in PATIENT_MAP.items():
        # Find untreated samples for this patient
        samples = pdat_ut.loc[pdat_ut["Patient"] == long_id, "sample_id"].tolist()

        if not samples:
            print(f"  WARNING: No untreated samples found for {short_id}, skipping.")
            continue

        # Check all samples exist in expression matrix
        missing = [s for s in samples if s not in exp.columns]
        if missing:
            print(f"  WARNING: Samples not found in exp.mat for {short_id}: {missing}")
            samples = [s for s in samples if s in exp.columns]

        if not samples:
            print(f"  ERROR: No valid samples for {short_id}, skipping.")
            continue

        print(f"\nPatient {short_id}: {len(samples)} untreated replicates")
        print(f"  Samples: {samples}")

        # Average counts across replicates
        mean_counts = exp[samples].mean(axis=1)

        # Filter by minimum mean count threshold
        expressed = mean_counts[mean_counts >= args.min_counts]
        print(f"  -> {len(expressed)} / {len(mean_counts)} genes expressed "
              f"(mean counts >= {args.min_counts})")

        # Save expressed gene list (ENSEMBL IDs, one per line)
        out_genes = os.path.join(args.outdir, f"{short_id}_expressed_genes.txt")
        expressed.index.to_series().to_csv(out_genes, index=False, header=False)
        print(f"  -> Saved: {out_genes}")

        # Save mean counts table
        out_counts = os.path.join(args.outdir, f"{short_id}_mean_counts.tsv")
        mean_counts.reset_index().rename(
            columns={"gene_id": "gene_id", 0: "mean_count"}
        ).to_csv(out_counts, sep="\t", index=False)

        summary_rows.append({
            "patient":         short_id,
            "n_replicates":    len(samples),
            "total_genes":     len(mean_counts),
            "expressed_genes": len(expressed),
            "threshold":       args.min_counts,
        })

    # Save summary
    summary_df = pd.DataFrame(summary_rows)
    summary_out = os.path.join(args.outdir, "expression_summary.tsv")
    summary_df.to_csv(summary_out, sep="\t", index=False)

    print("\n" + "="*60)
    print("EXPRESSION FILTERING SUMMARY")
    print("="*60)
    print(summary_df.to_string(index=False))
    print(f"\nAll outputs saved to: {args.outdir}")
    print("Done!")

if __name__ == "__main__":
    main()
