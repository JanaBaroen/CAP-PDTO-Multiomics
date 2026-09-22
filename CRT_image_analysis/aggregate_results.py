"""
Aggregate per-image ring-intensity results into condition comparisons + figures.

Reads results/per_image_summary.csv (written by nuclei_ring_analysis.py) and the
per-nucleus CSVs, then produces:
  results/master_summary.csv         per-image table (tidy)
  results/condition_means.csv        mean +/- sd per (group, condition)
  results/fig_ring_C1_by_condition.png
  results/fig_ring_vs_nucleus_C1.png
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"

COND_ORDER = ["control", "untreated", "treated"]


def main():
    global RESULTS
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", help="results directory to aggregate")
    args = ap.parse_args()
    if args.results_dir:
        RESULTS = Path(args.results_dir).resolve()

    summ = pd.read_csv(RESULTS / "per_image_summary.csv")
    summ.to_csv(RESULTS / "master_summary.csv", index=False)

    # condition means per group
    g = (
        summ.groupby(["group", "condition"])
        .agg(
            n_images=("image", "count"),
            n_nuclei_total=("n_nuclei", "sum"),
            ring_signal_mean=("ring_signal_mean", "mean"),
            ring_C1_sd=("ring_signal_mean", "std"),
            nucleus_signal_mean=("nucleus_signal_mean", "mean"),
        )
        .reset_index()
    )
    g.to_csv(RESULTS / "condition_means.csv", index=False)
    print(g.to_string(index=False))

    # ---- Figure 1: ring C1 per image, grouped by condition, colored by group ----
    summ["condition"] = pd.Categorical(
        summ["condition"], [c for c in COND_ORDER if c in summ["condition"].unique()]
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    conds = list(summ["condition"].cat.categories)
    groups = sorted(summ["group"].unique())
    cmap = plt.get_cmap("tab10")
    for ci, cond in enumerate(conds):
        sub = summ[summ["condition"] == cond]
        ax.boxplot(
            sub["ring_signal_mean"].dropna(), positions=[ci], widths=0.6,
            showfliers=False, patch_artist=True,
            boxprops=dict(facecolor="#dddddd"),
        )
        for gi, grp in enumerate(groups):
            pts = sub[sub["group"] == grp]["ring_signal_mean"].dropna()
            jitter = (np.random.RandomState(gi).rand(len(pts)) - 0.5) * 0.3
            ax.scatter(np.full(len(pts), ci) + jitter, pts, s=28,
                       color=cmap(gi % 10), label=grp if ci == 0 else None, zorder=3)
    ax.set_xticks(range(len(conds)))
    ax.set_xticklabels(conds)
    ax.set_ylabel("Per-image mean C0 signal in ring (mask3)")
    ax.set_title("Perinuclear ring C0 signal by condition")
    ax.legend(title="organoid", fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_ring_C1_by_condition.png", dpi=130)
    plt.close(fig)

    # ---- Figure 2: ring vs nucleus C1 (sanity / enrichment view) ----
    fig, ax = plt.subplots(figsize=(7, 7))
    for gi, grp in enumerate(groups):
        sub = summ[summ["group"] == grp]
        ax.scatter(sub["nucleus_signal_mean"], sub["ring_signal_mean"],
                   color=cmap(gi % 10), label=grp, s=30)
    lim = [0, np.nanmax([summ["nucleus_signal_mean"].max(), summ["ring_signal_mean"].max()]) * 1.05]
    ax.plot(lim, lim, "k--", lw=1, alpha=0.6)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("nucleus mean C0"); ax.set_ylabel("ring mean C0")
    ax.set_title("Ring vs nucleus C0 (per image)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_ring_vs_nucleus_C1.png", dpi=130)
    plt.close(fig)

    print(f"\nWrote master_summary.csv, condition_means.csv and 2 figures to {RESULTS}")


if __name__ == "__main__":
    main()
