# Organoid imaging & expression analysis pipeline

Scripts for perinuclear-ring intensity analysis of organoid immunofluorescence
imaging (CZI → TIFF → segmentation → per-nucleus/per-image quantification →
aggregation), plus a separate RNA-seq expressed-gene filtering step.

## Pipeline overview

**Imaging analysis** (run in this order):

1. **`extract_tiffs.py`** — reads raw `.czi` files, keeps the large tissue
   scenes, and writes each channel out as a 16-bit TIFF.
2. **`nuclei_ring_analysis.py`** — segments nuclei (Cellpose-SAM), builds a
   perinuclear "ring" mask around each nucleus, and measures signal-channel
   intensity in the ring and in the nucleus. Run with `--one <path>` for a
   single image or `--all` to process a whole `output_tiffs/` tree.
3. **`aggregate_results.py`** — rolls the per-image outputs of step 2 into
   condition-level summary tables and QC figures for a single dataset.
4. **`combine_results.py`** — combines two already-aggregated datasets
   (e.g. two separate imaging batches) into one combined table + figure set.

**Expression analysis** (independent of the imaging pipeline):

5. **`01_expression_filter.py`** — filters expressed genes per
   patient/organoid from an RNA-seq count matrix, averaging counts across
   untreated replicates.

Sample/organoid IDs throughout these scripts and examples use the format
`PDTO1`, `PDTO2`, etc.

## Setup

Before running any script on your own data, check the `# CHANGE THIS`
comments in each file — they mark values that need to be adapted to your
project:

| Script | What to configure |
|---|---|
| `extract_tiffs.py` | `ROOT` — your project directory (or just pass `--source` / `--output` on the command line) |
| `nuclei_ring_analysis.py` | `PIXEL_SIZE_UM` — your microscope's µm/pixel calibration; `ORG_RE` — the regex used to pull the sample ID out of each filename (defaults to matching IDs like `PDTO1`, `PDTO2`, ...) |
| `combine_results.py` | `DATASETS` — the labels and subfolder names of the result sets you want to combine |
| `01_expression_filter.py` | `PATIENT_MAP` — mapping of your sample IDs to the `Patient` column values in your metadata file |

## Directory layout expected by the imaging scripts

```
your_project/
├── source_file/                  # raw .czi files (input to extract_tiffs.py)
├── output_tiffs/                 # extracted TIFFs (output of extract_tiffs.py,
│                                  #   input to nuclei_ring_analysis.py)
│   ├── PDTO1/
│   ├── Controls/
│   └── ...
├── results/                      # per-image/per-nucleus outputs of
│                                  #   nuclei_ring_analysis.py + aggregate_results.py
│   ├── masks/
│   ├── qc/
│   └── per_image_summary.csv
└── results_final/                # staged datasets for combine_results.py
    ├── dataset1_source_file/
    └── dataset2_source_file_2/
```

## Dependencies

```
numpy
pandas
scipy
matplotlib
tifffile
scikit-image
aicspylibczi      # extract_tiffs.py only
cellpose          # nuclei_ring_analysis.py only (needs a GPU for practical runtimes)
```

Install with:

```bash
pip install numpy pandas scipy matplotlib tifffile scikit-image aicspylibczi cellpose
```

## Usage

```bash
# 1. Extract TIFFs from raw CZI files
python extract_tiffs.py --source ./source_file --output ./output_tiffs

# 2. Run nuclei/ring analysis on a single image
python nuclei_ring_analysis.py --one "output_tiffs/PDTO1/PDTO1 REP1 UT-01_S00_C1.tif" --qc

# 2b. Or run it on everything
python nuclei_ring_analysis.py --all --qc

# 3. Aggregate one dataset's per-image results
python aggregate_results.py --results-dir ./results

# 4. Combine two datasets that have already been aggregated
python combine_results.py

# 5. Filter expressed genes per patient from RNA-seq counts
python 01_expression_filter.py \
    --exp exp.mat.tsv \
    --pdat pdat.mat.tsv \
    --outdir expression_filtered \
    --min_counts 10
```

## Notes

- `nuclei_ring_analysis.py` expects paired TIFFs per image, named
  `<stem>_C1.tif` (nuclei channel) and `<stem>_C0.tif` (signal channel).
- The `condition` for each image is inferred from the filename: images in a
  `Controls/` folder are labeled `control`, filenames containing `ut` as a
  standalone token are labeled `untreated`, and everything else is labeled
  `treated`. Adjust the logic in `parse_meta()` if your naming convention
  differs.
