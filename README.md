# CAP-PDTO-Multiomics

Analysis scripts accompanying the manuscript:

> Multi-omics profiling reveals pronounced interpatient heterogeneity and low
> mutagenic impact of cold atmospheric plasma in patient-derived tumor
> organoids
>
> Jana Baroen, Axelle Van der Voort, Emma Peeters, Lisa Van der heyden,
> Sara Esbati, Louize Brants, Jasper Ott, Sofie Seghers,
> Felicia Rodrigues Fortes, Andreia Lopes, Daniel Flender, Hanne Verswyvel,
> Hannah Zaryouh, Geert Roeyen, Vera Hartman, Gilles Van Haesendonck,
> Hans Prenen, Senada Koljenovic, An Wouters, Christophe Deben,
> Jonas Van Audenaerde, Angela Privat-Maldonado, Inge Mertens,
> Evelien Smits, Annemie Bogaerts
>
> *(manuscript in preparation)*

This repository contains the code used to generate the multi-omics results
in the paper: RNA sequencing, whole-exome sequencing, proteomics, and
immunofluorescence imaging, across a panel of patient-derived tumor
organoids (PDTOs) treated with cold atmospheric plasma (CAP). Each
subfolder has its own README with setup instructions and usage examples
specific to that analysis.

## Contents

- **[`CRT_image_analysis/`](./CRT_image_analysis)** — organoid immunofluorescence
  imaging and perinuclear calreticulin (CRT) ring intensity quantification.
- **[`Proteomics/`](./Proteomics)** — label-free quantification (LFQ)
  proteomics analysis, including the patient-specific FASTA files used as
  the search database for each organoid.
- **[`RNAseq/`](./RNAseq)** — RNA sequencing analysis.
- **[`WES/`](./WES)** — whole-exome sequencing analysis.

## Data Availability

This repository contains analysis code only. The underlying sequencing and
proteomics data are deposited in public repositories (see manuscript for more details)

## Requirements

Each subfolder specifies its own dependencies (R or Python, as applicable)
in its README.

## Citation

A citation will be added here once the manuscript is published.
