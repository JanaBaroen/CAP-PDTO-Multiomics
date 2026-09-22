#!/usr/bin/env python3
"""
02_generate_genemap.py
======================
Downloads a mapping table linking ENSEMBL gene IDs to UniProt accessions
via the Ensembl BioMart REST API. Required by 03_build_fasta.py.

This script only needs to be run once.

Usage:
    python 02_generate_genemap.py --outdir /path/to/proteomics_ref

Output:
    <outdir>/uniprot_genemap.tsv
    Columns: ensembl_gene_id | gene_name | uniprot_accession

Requirements:
    pip install requests pandas


"""

import argparse
import os
import sys
import requests
import pandas as pd
from io import StringIO

###############################################################################
# ARGUMENT PARSING
###############################################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate UniProt -> ENSEMBL gene ID mapping via BioMart."
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    return parser.parse_args()

###############################################################################
# MAIN
###############################################################################
def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    out_file = os.path.join(args.outdir, "uniprot_genemap.tsv")

    print("Downloading ENSEMBL -> UniProt mapping from Ensembl BioMart...")
    print("(This may take 1-2 minutes)")

    url = "https://www.ensembl.org/biomart/martservice"
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1"
       uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="uniprotswissprot"/>
  </Dataset>
</Query>"""

    response = requests.post(url, data={"query": xml}, timeout=300)
    if response.status_code != 200:
        print(f"ERROR: BioMart returned status {response.status_code}")
        print(response.text[:500])
        sys.exit(1)

    df = pd.read_csv(StringIO(response.text), sep="\t")
    df.columns = ["ensembl_gene_id", "gene_name", "uniprot_accession"]

    # Drop rows with no UniProt accession
    df = df[df["uniprot_accession"].notna() & (df["uniprot_accession"] != "")]
    df = df.drop_duplicates()

    df.to_csv(out_file, sep="\t", index=False)

    print(f"\nGene map saved: {out_file}")
    print(f"  {len(df)} rows")
    print(f"  {df['uniprot_accession'].nunique()} unique UniProt accessions")
    print(f"  {df['ensembl_gene_id'].nunique()} unique ENSEMBL gene IDs")
    print("\nPreview:")
    print(df.head(5).to_string(index=False))
    print("Done!")

if __name__ == "__main__":
    main()
