#!/usr/bin/env python
"""Assemble the final manuscript figures from versioned source panels.

The analysis scripts write source panels to ``figures/source``. This script
only lays out those panels; it does not refit models or recalculate statistics.
Figure 5 is generated separately by ``make_virtual_knockout_figure.R``.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "figures" / "source"
FINAL = ROOT / "figures" / "final"
FINAL.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {"font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 150}
)


def source(name: str) -> Path:
    path = SOURCE / name
    if not path.exists():
        raise FileNotFoundError(
            f"Missing source panel {path}. Run workflow/run_core_analysis.ps1 first."
        )
    return path


def panel_label(ax, letter: str) -> None:
    ax.text(
        -0.018,
        1.025,
        letter,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


def add_panel(ax, filename: str, letter: str, title: str, crop_top: float = 1.0) -> None:
    image = mpimg.imread(source(filename))
    if crop_top < 1.0:
        image = image[: int(image.shape[0] * crop_top), :]
    ax.imshow(image)
    ax.axis("off")
    panel_label(ax, letter)
    ax.set_title(title, fontsize=9, loc="left", pad=5)


def save_both(fig, stem: str) -> None:
    fig.savefig(FINAL / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FINAL / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("saved", stem)


def grid_figure(items, nrows: int, ncols: int, stem: str, figsize) -> None:
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_1d(axes).ravel()
    for ax in axes:
        ax.axis("off")
    for ax, item in zip(axes, items):
        add_panel(ax, *item)
    fig.tight_layout(pad=1.1)
    save_both(fig, stem)


def ensure_expression_context_panel() -> None:
    """Build the Figure 4 expression-context source panel if needed."""
    output = SOURCE / "fig4a_sting_by_group.png"
    if output.exists():
        return
    modules = pd.read_csv(ROOT / "analysis" / "module_scores_by_library_v4.csv")
    stats = pd.read_csv(ROOT / "analysis" / "module_scores_liao_stats_v4.csv")
    stats = stats.set_index("score")
    order = ["Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH"]
    palette = dict(zip(order, sns.color_palette("Set2", len(order))))
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.2))
    for ax, column, title in [
        (axes[0], "EC_inflammation", "Endothelial inflammation"),
        (axes[1], "cGAS_STING", "cGAS-STING expression"),
    ]:
        sns.boxplot(data=modules, x="group", y=column, order=order, hue="group", palette=palette, legend=False, width=0.58, showfliers=False, ax=ax)
        sns.stripplot(data=modules, x="group", y=column, order=order, color="#222222", size=4.8, jitter=0.12, ax=ax)
        fdr = float(stats.loc[column, "Liao_KW_fdr"])
        ax.set_title(f"{title} module\nLiao exact KW FDR = {fdr:.2g}")
        ax.set_xlabel("")
        ax.set_ylabel("Mean module score")
        ax.tick_params(axis="x", rotation=35)
        ax.grid(axis="x", visible=False)
        sns.despine(ax=ax)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.86, bottom=0.24, wspace=0.30)
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


ensure_expression_context_panel()


grid_figure(
    [
        ("umap_clusters.png", "A", "Leiden clusters"),
        ("umap_celltype.png", "B", "Marker-panel cell-type labels"),
        ("umap_group.png", "C", "Clinical/source groups"),
        ("umap_dataset.png", "D", "Source datasets (descriptive visualization)"),
        ("fig1e_v4_ec_fraction.png", "E", "EC fraction by sampling unit"),
        ("fig1f_v4_retained_cells.png", "F", "Retained cells per library"),
    ],
    2,
    3,
    "Figure1",
    (16, 9.5),
)

grid_figure(
    [
        ("fig2a_ec_umap.png", "A", "EC clusters and panel-defined states"),
        ("fig2b_subtype_dotplot.png", "B", "Subtype marker panels"),
        ("fig2c_v4_subtype_composition.png", "C", "Subtype composition by sampling unit"),
        ("fig2d_v4_module_scores.png", "D", "Endothelial module scores by sampling unit"),
    ],
    2,
    2,
    "Figure2",
    (15, 11),
)

fig = plt.figure(figsize=(16, 10.5))
gs = GridSpec(2, 2, figure=fig, height_ratios=[0.95, 1.2], hspace=0.15, wspace=0.08)
add_panel(fig.add_subplot(gs[0, 0]), "fig3a_sonfh_effects_descriptive_v4.png", "A", "SONFH library-level endothelial effects")
add_panel(fig.add_subplot(gs[0, 1]), "fig3b_key_gene_effects_v4.png", "B", "Mitochondrial stress and selective-clearance genes")
add_panel(fig.add_subplot(gs[1, :]), "fig2d_v4_module_scores.png", "C", "Sampling-unit mitochondrial module scores", crop_top=0.49)
fig.subplots_adjust(left=0.025, right=0.985, top=0.97, bottom=0.03, hspace=0.14, wspace=0.08)
save_both(fig, "Figure3")

grid_figure(
    [
        ("fig3c_liao_gsea_v4.png", "A", "Within-cohort Hallmark enrichment"),
        ("fig4a_sting_by_group.png", "B", "cGAS-STING-related expression context"),
        ("fig5b_v4_ulm_liao.png", "C", "Signed TF activities in independent Liao donors"),
        ("fig4c_v4_liao_lr_points.png", "D", "EC-centred communication effects by donor"),
    ],
    2,
    2,
    "Figure4",
    (16, 11.5),
)

classifier = source("fig6_v4_nested_cv.png")
shutil.copy2(classifier, FINAL / "Figure6.png")
image = mpimg.imread(classifier)
fig, ax = plt.subplots(figsize=(image.shape[1] / 300, image.shape[0] / 300), dpi=300)
ax.imshow(image)
ax.axis("off")
fig.subplots_adjust(0, 0, 1, 1)
fig.savefig(FINAL / "Figure6.pdf", bbox_inches="tight", pad_inches=0)
plt.close(fig)
print("saved Figure6")

if not (FINAL / "Figure5.pdf").exists():
    print("Figure5 is pending: run plotting/make_virtual_knockout_figure.R")
