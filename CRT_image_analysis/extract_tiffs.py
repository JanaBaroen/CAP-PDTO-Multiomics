"""
TIF extraction script for organoid imaging analysis.

Reads all CZI files under source_file/, keeps only the "large" scenes
(area >= SIZE_THRESHOLD_RATIO * max scene area per file), and saves
each channel as a 16-bit TIFF under output_tiffs/.

Output naming:
  output_tiffs/{organoid_type}/{stem}_S{scene_idx:02d}_C{channel}.tif
"""

import aicspylibczi
import tifffile
import numpy as np
from pathlib import Path
import sys

# ── Configuration ────────────────────────────────────────────────────────────

# CHANGE THIS: set to your project root directory
ROOT = Path("/path/to/your/project")
SOURCE_DIR = ROOT / "source_file"
OUTPUT_DIR = ROOT / "output_tiffs"

# Keep scenes whose area >= this fraction of the largest scene in the file.
# 0.5 means "at least half the size of the biggest scene".
SIZE_THRESHOLD_RATIO = 0.5

N_CHANNELS = 2  # always 2 per file

# ── Helpers ──────────────────────────────────────────────────────────────────

def get_large_scenes(czi: aicspylibczi.CziFile, ratio: float) -> list:
    """Return [(scene_idx, (x, y, w, h)), ...] for scenes above the size threshold."""
    n_scenes = czi.get_dims_shape()[0]["S"][1] - czi.get_dims_shape()[0]["S"][0]
    bboxes = []
    for s in range(n_scenes):
        bb = czi.get_mosaic_scene_bounding_box(index=s)
        bboxes.append((s, bb, bb.w * bb.h))

    max_area = max(area for _, _, area in bboxes)
    threshold = ratio * max_area

    return [
        (s, (bb.x, bb.y, bb.w, bb.h))
        for s, bb, area in bboxes
        if area >= threshold
    ]


def extract_nonmosaic(czi, czi_path: Path, out_dir: Path, ratio: float) -> int:
    """Extract scenes from a non-mosaic (single-tile) multi-scene CZI."""
    dims = czi.get_dims_shape()[0]
    s0, s1 = dims.get("S", (0, 1))
    stem = czi_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)

    scenes = {}
    areas = {}
    for s in range(s0, s1):
        chans = []
        for ch in range(N_CHANNELS):
            data, _shp = czi.read_image(S=s, C=ch)
            arr = np.squeeze(data)
            chans.append(arr)
        scenes[s] = chans
        areas[s] = chans[0].shape[-1] * chans[0].shape[-2]

    threshold = ratio * max(areas.values())
    written = 0
    for s in range(s0, s1):
        if areas[s] < threshold:
            continue
        for ch in range(N_CHANNELS):
            img2d = scenes[s][ch]
            out_path = out_dir / f"{stem}_S{s:02d}_C{ch}.tif"
            tifffile.imwrite(str(out_path), img2d, photometric="minisblack")
            print(f"  wrote {out_path.name}  {img2d.shape}  {img2d.dtype}  [non-mosaic]")
            written += 1
    return written


def extract_file(czi_path: Path, out_dir: Path, ratio: float) -> int:
    """Extract large scenes from one CZI file. Returns number of TIFFs written."""
    czi = aicspylibczi.CziFile(str(czi_path))

    if not czi.is_mosaic():
        return extract_nonmosaic(czi, czi_path, out_dir, ratio)

    large_scenes = get_large_scenes(czi, ratio)
    if not large_scenes:
        print(f"  [SKIP] no scenes passed size filter: {czi_path.name}")
        return 0

    stem = czi_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    written = 0

    for scene_idx, region in large_scenes:
        for ch in range(N_CHANNELS):
            img = czi.read_mosaic(region=region, scale_factor=1.0, C=ch)
            # img shape: (1, H, W)  — drop the leading singleton
            img2d = img[0]

            out_name = f"{stem}_S{scene_idx:02d}_C{ch}.tif"
            out_path = out_dir / out_name
            tifffile.imwrite(str(out_path), img2d, photometric="minisblack")
            print(f"  wrote {out_path.relative_to(ROOT)}  {img2d.shape}  {img2d.dtype}")
            written += 1

    return written


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(SOURCE_DIR),
                    help="directory of .czi files (searched recursively)")
    ap.add_argument("--output", default=str(OUTPUT_DIR),
                    help="output directory for extracted TIFFs")
    args = ap.parse_args()
    source_dir = Path(args.source).resolve()
    output_dir = Path(args.output).resolve()

    if not source_dir.exists():
        print(f"ERROR: source directory not found: {source_dir}", file=sys.stderr)
        sys.exit(1)

    czi_files = sorted(source_dir.rglob("*.czi"))
    if not czi_files:
        print("No CZI files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(czi_files)} CZI files in {source_dir}. "
          f"Size threshold ratio: {SIZE_THRESHOLD_RATIO}\n")

    total_written = 0
    total_skipped = 0

    for czi_path in czi_files:
        # Mirror the source subfolder structure (e.g. PDTO4/, Controls/)
        rel_parent = czi_path.parent.relative_to(source_dir)
        out_dir = output_dir / rel_parent

        print(f"Processing: {czi_path.relative_to(ROOT)}")
        try:
            n = extract_file(czi_path, out_dir, SIZE_THRESHOLD_RATIO)
            if n == 0:
                total_skipped += 1
            else:
                total_written += n
        except Exception as e:
            print(f"  [ERROR] {e}", file=sys.stderr)
            total_skipped += 1

    print(f"\nDone. {total_written} TIFFs written, {total_skipped} files skipped/errored.")


if __name__ == "__main__":
    main()
