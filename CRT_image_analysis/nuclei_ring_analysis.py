"""
Nuclei perinuclear-ring intensity analysis.

Channel roles (set by NUCLEI_SUFFIX / SIGNAL_SUFFIX below):
  C1 = nuclei (segmented)        C0 = staining signal (measured)

For each image pair in output_tiffs/:
  mask1 : nuclei segmentation (Cellpose-SAM 'cpsam' model on the nuclei channel)
  mask2 : each nucleus grown so its area is +10% (Voronoi-limited, no overlap)
  mask3 : ring = mask2 - mask1  (perinuclear band, area ~= 10% of nucleus)
  measure signal-channel intensity inside mask3 (and, for reference, inside mask1)

Outputs (under results/):
  <stem>_mask1.tif / _mask2.tif / _mask3.tif   label images (uint16/uint32)
  <stem>_nuclei.csv                            one row per nucleus
  per_image_summary.csv                        one row per image (appended)
  qc/<stem>_qc.png                             QC overlay (optional)

Usage (pass the NUCLEI-channel tiff, i.e. *_C1.tif):
  python nuclei_ring_analysis.py --one "output_tiffs/PDTO4/PDTO4 REP1 UT-01_S00_C1.tif"
  python nuclei_ring_analysis.py --all
"""

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parent
TIFF_DIR = ROOT / "output_tiffs"
RESULTS_DIR = ROOT / "results"
MASK_DIR = RESULTS_DIR / "masks"
QC_DIR = RESULTS_DIR / "qc"

RING_AREA_FRAC = 0.10  # mask2 is +10% area -> ring (mask3) area = 10% of nucleus
CELLPOSE_DIAMETER = 30  # cpsam native scale; nuclei measured ~30 px here

# CHANGE THIS: microns per pixel, from your CZI/microscope metadata
PIXEL_SIZE_UM = 0.2287
SCALEBAR_UM = 50        # QC scale-bar length in microns

# Channel roles. Nuclei are segmented; signal is measured in the ring.
NUCLEI_SUFFIX = "_C1.tif"   # nuclei channel
SIGNAL_SUFFIX = "_C0.tif"   # staining / signal channel

_MODEL = None


def get_model():
    global _MODEL
    if _MODEL is None:
        from cellpose import models
        _MODEL = models.CellposeModel(gpu=True)
        print(f"  cellpose device: {_MODEL.device}", flush=True)
    return _MODEL


def tissue_bbox(c0: np.ndarray, pad: int = 200):
    """Foreground bounding box (y0, y1, x0, x1) with padding, to skip empty margins."""
    from skimage import filters, morphology

    ds = c0[::8, ::8].astype(np.float32)
    try:
        th = filters.threshold_otsu(ds)
    except ValueError:
        return 0, c0.shape[0], 0, c0.shape[1]
    fg = morphology.remove_small_objects(ds > th, 32)
    if not fg.any():
        return 0, c0.shape[0], 0, c0.shape[1]
    ys, xs = np.nonzero(fg)
    h, w = c0.shape
    y0 = max(0, ys.min() * 8 - pad)
    x0 = max(0, xs.min() * 8 - pad)
    y1 = min(h, ys.max() * 8 + pad)
    x1 = min(w, xs.max() * 8 + pad)
    return y0, y1, x0, x1


def segment_nuclei(c0: np.ndarray) -> np.ndarray:
    """Return uint32 label image of nuclei (mask1)."""
    model = get_model()
    masks, _flows, _styles = model.eval(c0, diameter=CELLPOSE_DIAMETER, normalize=True)
    return masks.astype(np.uint32)


def build_rings(labels: np.ndarray, area_frac: float = RING_AREA_FRAC):
    """
    Grow each nucleus so its area increases by `area_frac`, limited to its
    own Voronoi territory so rings never overlap and never sit on another
    nucleus. Returns (mask2_labels, mask3_ring_labels), both same dtype as labels.

    The ring for object i is the `ceil(area_frac * area_i)` background pixels
    closest to object i (and closer to i than to any other nucleus). This makes
    the ring area exactly area_frac of the nucleus area (per object).
    """
    labels = labels.astype(np.uint32, copy=False)
    bg = labels == 0
    if not bg.any() or labels.max() == 0:
        return labels.copy(), np.zeros_like(labels)

    # nearest-nucleus label + distance for every background pixel
    dist, (iy, ix) = ndi.distance_transform_edt(bg, return_indices=True)
    nearest = labels[iy, ix]  # nucleus label of nearest nucleus, for all pixels

    bg_ys, bg_xs = np.nonzero(bg)
    bg_lbl = nearest[bg_ys, bg_xs].astype(np.uint32)
    bg_dist = dist[bg_ys, bg_xs]

    # per-object pixel budget = ceil(frac * area)
    area = np.bincount(labels.ravel())  # area[label]; area[0] = background
    budget = np.ceil(area_frac * area).astype(np.int64)

    # order background pixels by (label, distance) and keep the closest `budget`
    order = np.lexsort((bg_dist, bg_lbl))
    sl = bg_lbl[order]
    change = np.empty(sl.shape, dtype=bool)
    change[0] = True
    change[1:] = sl[1:] != sl[:-1]
    grp_start = np.flatnonzero(change)
    counts = np.diff(np.append(grp_start, sl.size))
    rank = np.arange(sl.size) - np.repeat(grp_start, counts)
    keep = rank < budget[sl]

    keep_idx = order[keep]
    ring = np.zeros_like(labels)
    ring[bg_ys[keep_idx], bg_xs[keep_idx]] = bg_lbl[keep_idx]

    mask2 = labels.copy()
    rmask = ring > 0
    mask2[rmask] = ring[rmask]
    return mask2, ring


def _stats(values: np.ndarray):
    if values.size == 0:
        return dict(mean=np.nan, median=np.nan, sum=0.0, std=np.nan, n_px=0)
    return dict(
        mean=float(values.mean()),
        median=float(np.median(values)),
        sum=float(values.sum()),
        std=float(values.std()),
        n_px=int(values.size),
    )


def measure(sig: np.ndarray, nuc: np.ndarray, ring: np.ndarray) -> pd.DataFrame:
    """Per-nucleus signal stats in the ring (mask3) and in the nucleus (mask1)."""
    c1f = sig.astype(np.float64)
    n = int(nuc.max())
    if n == 0:
        return pd.DataFrame()

    # ring stats per label
    ring_sum = ndi.sum_labels(c1f, ring, index=np.arange(1, n + 1))
    ring_cnt = ndi.sum_labels(np.ones_like(c1f), ring, index=np.arange(1, n + 1))
    ring_mean = np.divide(ring_sum, ring_cnt, out=np.full(n, np.nan), where=ring_cnt > 0)
    ring_med = ndi.labeled_comprehension(
        c1f, ring, np.arange(1, n + 1), np.median, float, np.nan
    )
    # nucleus stats per label
    nuc_sum = ndi.sum_labels(c1f, nuc, index=np.arange(1, n + 1))
    nuc_cnt = ndi.sum_labels(np.ones_like(c1f), nuc, index=np.arange(1, n + 1))
    nuc_mean = np.divide(nuc_sum, nuc_cnt, out=np.full(n, np.nan), where=nuc_cnt > 0)

    # centroids
    centroids = ndi.center_of_mass(np.ones_like(c1f), nuc, index=np.arange(1, n + 1))
    cy = [c[0] for c in centroids]
    cx = [c[1] for c in centroids]

    return pd.DataFrame(
        {
            "label": np.arange(1, n + 1),
            "centroid_y": cy,
            "centroid_x": cx,
            "nucleus_area_px": nuc_cnt.astype(int),
            "ring_area_px": ring_cnt.astype(int),
            "ring_signal_mean": ring_mean,
            "ring_signal_median": ring_med,
            "ring_signal_sum": ring_sum,
            "nucleus_signal_mean": nuc_mean,
            "nucleus_signal_sum": nuc_sum,
        }
    )


# CHANGE THIS: pattern matching your organoid/sample ID naming scheme.
# Filenames are expected to contain an ID like "PDTO1", "PDTO2", etc.
ORG_RE = re.compile(r"(PDTO\d+)", re.I)


def parse_meta(nuc_path: Path):
    """Derive (group, image_name, condition) robustly from folder + filename.

    Works for both the sub-foldered first dataset and the flat second dataset,
    and for names using spaces or underscores as separators.
    """
    rel = nuc_path.relative_to(TIFF_DIR)
    folder = rel.parts[0] if len(rel.parts) > 1 else None
    name = nuc_path.name[: -len(NUCLEI_SUFFIX)]  # strip _C1.tif
    toks = re.split(r"[ _\-]+", name.lower())

    # controls live in a Controls/ folder in dataset 1
    if folder == "Controls":
        return "Controls", name, "control"

    m = ORG_RE.search(name)
    group = re.sub(r"\s+", "", m.group(1)).upper() if m else (folder or "unknown")

    condition = "untreated" if "ut" in toks else "treated"
    return group, name, condition


def process_pair(nuc_path: Path, save_masks=True, save_qc=False):
    nuc_path = Path(nuc_path).resolve()
    sig_path = nuc_path.with_name(nuc_path.name.replace(NUCLEI_SUFFIX, SIGNAL_SUFFIX))
    if not sig_path.exists():
        print(f"  [SKIP] no signal channel for {nuc_path.name}", file=sys.stderr)
        return None

    group, name, condition = parse_meta(nuc_path)
    nuc_full = tifffile.imread(nuc_path)   # nuclei channel (C1)
    sig_full = tifffile.imread(sig_path)   # signal channel (C0)
    H, W = nuc_full.shape

    t = time.time()
    y0, y1, x0, x1 = tissue_bbox(nuc_full)
    nuc_img = nuc_full[y0:y1, x0:x1]
    sig_img = sig_full[y0:y1, x0:x1]

    nuc = segment_nuclei(nuc_img)      # labels on the crop
    n = int(nuc.max())
    mask2, ring = build_rings(nuc)
    df = measure(sig_img, nuc, ring)   # signal measured in ring/nucleus
    dt = time.time() - t

    # centroids -> full-image coordinates
    if n:
        df["centroid_y"] = df["centroid_y"] + y0
        df["centroid_x"] = df["centroid_x"] + x0
    df.insert(0, "image", name)
    df.insert(1, "group", group)
    df.insert(2, "condition", condition)

    if save_masks:
        MASK_DIR.mkdir(parents=True, exist_ok=True)
        dt_save = np.uint16 if n < 65535 else np.uint32

        def place_full(crop_lbl):
            full = np.zeros((H, W), dtype=dt_save)
            full[y0:y1, x0:x1] = crop_lbl.astype(dt_save)
            return full

        kw = dict(compression="zlib")
        tifffile.imwrite(MASK_DIR / f"{name}_mask1.tif", place_full(nuc), **kw)
        tifffile.imwrite(MASK_DIR / f"{name}_mask2.tif", place_full(mask2), **kw)
        tifffile.imwrite(MASK_DIR / f"{name}_mask3.tif", place_full(ring), **kw)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS_DIR / f"{name}_nuclei.csv", index=False)

    if save_qc:
        save_qc_overlay(nuc_img, sig_img, nuc, ring, name)

    summary = dict(
        image=name,
        group=group,
        condition=condition,
        n_nuclei=n,
        total_nucleus_area_px=int(df["nucleus_area_px"].sum()) if n else 0,
        total_ring_area_px=int(df["ring_area_px"].sum()) if n else 0,
        ring_signal_mean=float(np.nansum(df["ring_signal_sum"]) / max(df["ring_area_px"].sum(), 1)) if n else np.nan,
        ring_signal_median_of_nuclei=float(np.nanmedian(df["ring_signal_mean"])) if n else np.nan,
        ring_signal_total=float(np.nansum(df["ring_signal_sum"])) if n else 0.0,
        nucleus_signal_mean=float(np.nansum(df["nucleus_signal_sum"]) / max(df["nucleus_area_px"].sum(), 1)) if n else np.nan,
        seconds=round(dt, 1),
    )
    print(f"  {name}: {n} nuclei, ring mean signal={summary['ring_signal_mean']:.1f} ({dt:.0f}s)", flush=True)
    return summary


def save_qc_overlay(nuc_img, sig_img, nuc, ring, name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries

    QC_DIR.mkdir(parents=True, exist_ok=True)
    # crop to a dense 1024 window for legibility
    ys, xs = np.nonzero(nuc)
    if ys.size:
        cy, cx = int(np.median(ys)), int(np.median(xs))
    else:
        cy, cx = nuc_img.shape[0] // 2, nuc_img.shape[1] // 2
    h = 512
    y0, x0 = max(0, cy - h), max(0, cx - h)
    sl = (slice(y0, y0 + 2 * h), slice(x0, x0 + 2 * h))
    nucc_img, sigc, nucc, ringc = nuc_img[sl], sig_img[sl], nuc[sl], ring[sl]

    def norm(a):
        a = a.astype(float)
        lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
        return np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1)

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))
    ax[0].imshow(norm(nucc_img), cmap="gray"); ax[0].set_title("C1 nuclei + segmentation")
    nb = find_boundaries(nucc, mode="outer")
    ov = np.zeros((*nucc.shape, 4)); ov[nb] = [1, 0, 0, 1]
    ax[0].imshow(ov)
    ax[1].imshow(norm(sigc), cmap="gray"); ax[1].set_title("C0 signal + ring (mask3)")
    rb = ringc > 0
    ov2 = np.zeros((*ringc.shape, 4)); ov2[rb] = [0, 1, 1, 0.6]
    ax[1].imshow(ov2)
    ax[2].imshow(norm(sigc), cmap="magma"); ax[2].set_title("C0 signal")
    # fixed-length white scale bar (SCALEBAR_UM microns) in the bottom-right corner
    from matplotlib.patches import Rectangle
    ph, pw = sigc.shape
    bar_len = SCALEBAR_UM / PIXEL_SIZE_UM          # length in pixels
    bar_h = max(3.0, 0.015 * ph)
    mx, my = 0.05 * pw, 0.07 * ph
    bx, by = pw - mx - bar_len, ph - my - bar_h
    ax[2].add_patch(Rectangle((bx, by), bar_len, bar_h,
                              color="white", ec="none", zorder=5))
    # no text label: the "50 µm" annotation is added manually
    for a in ax:
        a.axis("off")
    fig.suptitle(name)
    fig.tight_layout()
    fig.savefig(QC_DIR / f"{name}_qc.png", dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--one", type=str, help="single nuclei (C1) tiff path")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--qc", action="store_true", help="save QC overlays")
    ap.add_argument("--no-masks", action="store_true", help="do not save mask tiffs")
    ap.add_argument("--tiff-dir", help="override input TIFF directory")
    ap.add_argument("--results-dir", help="override output results directory")
    args = ap.parse_args()

    global TIFF_DIR, RESULTS_DIR, MASK_DIR, QC_DIR
    if args.tiff_dir:
        TIFF_DIR = Path(args.tiff_dir).resolve()
    if args.results_dir:
        RESULTS_DIR = Path(args.results_dir).resolve()
        MASK_DIR = RESULTS_DIR / "masks"
        QC_DIR = RESULTS_DIR / "qc"

    save_masks = not args.no_masks
    summaries = []
    if args.one:
        s = process_pair(Path(args.one), save_masks=save_masks, save_qc=args.qc)
        if s:
            summaries.append(s)
    elif args.all:
        nucs = sorted(TIFF_DIR.rglob(f"*{NUCLEI_SUFFIX}"))
        print(f"Processing {len(nucs)} image pairs", flush=True)
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        for i, p in enumerate(nucs, 1):
            print(f"[{i}/{len(nucs)}] {p.relative_to(TIFF_DIR)}", flush=True)
            try:
                s = process_pair(p, save_masks=save_masks, save_qc=args.qc)
                if s:
                    summaries.append(s)
                    # incremental save so partial progress survives a crash
                    pd.DataFrame(summaries).to_csv(
                        RESULTS_DIR / "per_image_summary.csv", index=False
                    )
            except Exception as e:
                print(f"  [ERROR] {e}", file=sys.stderr, flush=True)
    else:
        ap.error("use --one PATH or --all")

    if summaries:
        RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        out = RESULTS_DIR / "per_image_summary.csv"
        df = pd.DataFrame(summaries)
        if out.exists() and args.all is False:
            prev = pd.read_csv(out)
            df = pd.concat([prev[~prev.image.isin(df.image)], df], ignore_index=True)
        df.to_csv(out, index=False)
        print(f"\nWrote {out} ({len(df)} rows)")


if __name__ == "__main__":
    main()
