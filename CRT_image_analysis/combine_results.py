"""
Combine dataset-1 and dataset-2 ring-analysis results into results_final/combined/.

Reads the two per-dataset outputs already staged under results_final/ and writes:
  combined/per_image_summary.csv    all images, with a 'dataset' column
  combined/per_nucleus_all.csv      every nucleus from both datasets
  combined/condition_means.csv      mean +/- sd per (dataset, group, condition)
  combined/figures/ring_signal_by_condition.png
  combined/figures/ring_vs_nucleus_signal.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FINAL = ROOT / "results_final"
COMBINED = FINAL / "combined"
FIGDIR = COMBINED / "figures"

# CHANGE THIS: labels and subfolder names for the datasets you're combining
# (each subfolder should already contain per_image_summary.csv and a
# nuclei_csv/ directory, as written by nuclei_ring_analysis.py)
DATASETS = {
    "dataset1 (source_file)": FINAL / "dataset1_source_file",
    "dataset2 (source_file_2)": FINAL / "dataset2_source_file_2",
}
COND_ORDER = ["control", "untreated", "treated"]


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)

    # ---- combined per-image ----
    per_img = []
    per_nuc = []
    for label, d in DATASETS.items():
        s = pd.read_csv(d / "per_image_summary.csv")
        s.insert(0, "dataset", label)
        per_img.append(s)
        for f in sorted((d / "nuclei_csv").glob("*_nuclei.csv")):
            n = pd.read_csv(f)
            n.insert(0, "dataset", label)
            per_nuc.append(n)

    img = pd.concat(per_img, ignore_index=True)
    img.to_csv(COMBINED / "per_image_summary.csv", index=False)

    nuc = pd.concat(per_nuc, ignore_index=True)
    nuc.to_csv(COMBINED / "per_nucleus_all.csv", index=False)

    # ---- condition means ----
    means = (
        img.groupby(["dataset", "group", "condition"])
        .agg(
            n_images=("image", "count"),
            n_nuclei_total=("n_nuclei", "sum"),
            ring_signal_mean=("ring_signal_mean", "mean"),
            ring_signal_sd=("ring_signal_mean", "std"),
            nucleus_signal_mean=("nucleus_signal_mean", "mean"),
        )
        .reset_index()
    )
    means.to_csv(COMBINED / "condition_means.csv", index=False)
    print(means.to_string(index=False))

    # ---- Figure 1: per-image ring signal by condition, per dataset ----
    conds_all = [c for c in COND_ORDER if c in img["condition"].unique()]
    fig, axes = plt.subplots(1, len(DATASETS), figsize=(14, 6), sharey=False)
    cmap = plt.get_cmap("tab10")
    groups_all = sorted(img["group"].unique())
    gcolor = {g: cmap(i % 10) for i, g in enumerate(groups_all)}
    for ax, (label, _d) in zip(np.atleast_1d(axes), DATASETS.items()):
        sub_ds = img[img["dataset"] == label]
        conds = [c for c in conds_all if c in sub_ds["condition"].unique()]
        for ci, cond in enumerate(conds):
            sub = sub_ds[sub_ds["condition"] == cond]
            vals = sub["ring_signal_mean"].dropna()
            if len(vals):
                ax.boxplot(vals, positions=[ci], widths=0.6, showfliers=False,
                           patch_artist=True, boxprops=dict(facecolor="#eeeeee"))
            for grp, pts in sub.groupby("group"):
                y = pts["ring_signal_mean"].dropna()
                jit = (np.random.RandomState(hash(grp) % 999).rand(len(y)) - 0.5) * 0.3
                ax.scatter(np.full(len(y), ci) + jit, y, s=30,
                           color=gcolor[grp], label=grp, zorder=3)
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(conds)
        ax.set_title(label)
        ax.set_ylabel("Per-image mean C0 signal in ring (mask3)")
    # single legend
    handles = [plt.Line2D([0], [0], marker="o", ls="", color=gcolor[g], label=g)
               for g in groups_all]
    fig.legend(handles=handles, title="organoid", loc="upper right", fontsize=8)
    fig.suptitle("Perinuclear ring C0 signal by condition")
    fig.tight_layout(rect=[0, 0, 0.92, 0.96])
    fig.savefig(FIGDIR / "ring_signal_by_condition.png", dpi=130)
    plt.close(fig)

    # ---- Figure 2: ring vs nucleus C0 (per image), colored by dataset ----
    fig, ax = plt.subplots(figsize=(7, 7))
    for i, (label, _d) in enumerate(DATASETS.items()):
        sub = img[img["dataset"] == label]
        ax.scatter(sub["nucleus_signal_mean"], sub["ring_signal_mean"],
                   s=32, color=cmap(i), label=label)
    lim = [0, np.nanmax([img["nucleus_signal_mean"].max(),
                         img["ring_signal_mean"].max()]) * 1.05]
    ax.plot(lim, lim, "k--", lw=1, alpha=0.6)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("nucleus mean C0"); ax.set_ylabel("ring mean C0")
    ax.set_title("Ring vs nucleus C0 signal (per image)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGDIR / "ring_vs_nucleus_signal.png", dpi=130)
    plt.close(fig)

    print(f"\nImages: {len(img)}  Nuclei: {len(nuc)}")
    print(f"Wrote combined outputs to {COMBINED}")


if __name__ == "__main__":
    main()
