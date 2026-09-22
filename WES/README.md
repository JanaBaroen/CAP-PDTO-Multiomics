# WES Analysis

This directory contains the scripts used for whole-exome sequencing (WES) analysis of untreated and CAP-treated patient-derived tumour organoids (PDTOs) from pancreatic ductal adenocarcinoma (PDAC; PDTO1–3) and head and neck squamous cell carcinoma (HNSCC; PDTO4–6), as described in the associated manuscript.

## Study design

Two separate variant calling and analysis pipelines were used depending on sample type:

| | Untreated samples | Treated samples |
|---|---|---|
| **Mutect2 mode** | Tumor-only (no matched normal) | Paired (treated tumour vs. untreated tumour as reference) |
| **Germline filtering** | Custom multi-step pipeline required | Not required — germline signal removed by paired calling |
| **Purpose** | Baseline mutation landscape | Treatment-induced mutations |

## Directory structure

```
WES/
├── 1_untreated/     — variant calling, germline filtering, and baseline analysis
│                      of untreated PDTO samples
└── 2_treated/       — paired variant calling and analysis of CAP-treated
                       PDTO samples
```

See the README in each subfolder for full details on the pipeline, software requirements, and how to run the scripts.
