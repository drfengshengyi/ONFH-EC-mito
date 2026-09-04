#!/usr/bin/env python
"""Build reviewer-corrected Figures 2, 3, and 6 from versioned result tables.

The script intentionally reads plot labels and FDR values from machine-readable
analysis outputs.  It does not contain manually transcribed statistical values.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpecFromSubplotSpec


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "figure_inputs"
OUTPUT = ROOT / "figures" / "final"

GROUP_ORDER = ["Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH"]
GROUP_LABELS = ["Healthy", "HOA", "FNF", "ARCO 3A", "ARCO 4", "SONFH\n(libraries)"]
GROUP_COLORS = {
    "Healthy": "#2c7fb8",
    "HOA": "#56b4e9",
    "FNF": "#1b9e77",
    "ONFH_3A": "#e6a51a",
    "ONFH_4": "#d95f0e",
    "SONFH": "#cc79a7",
}
SUBTYPE_ORDER = [
    "Lymphatic",
    "Type H / EMCN-KDR",
    "Type R / bone-remodeling",
    "Venous / ACKR1",
]
SUBTYPE_COLUMNS = {
    "Lymphatic": "lymphatic",
    "Type H / EMCN-KDR": "typeH_EMCN_KDR",
    "Type R / bone-remodeling": "typeR_bone_remodel",
    "Venous / ACKR1": "venous_ACKR1",
}
SUBTYPE_COLORS = {
    "Lymphatic": "#cc79a7",
    "Type H / EMCN-KDR": "#0072b2",
    "Type R / bone-remodeling": "#009e73",
    "Venous / ACKR1": "#e69f00",
}
MODULE_LABELS = {
    "Mito_fission": "Mito fission",
    "Mito_fusion": "Mito fusion",
    "Mitophagy_core": "Mitophagy core",
    "cGAS_STING": "cGAS-STING",
    "EC_inflammation": "EC inflammation",
    "YAP_mTOR": "YAP-mTOR",
}


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(OUTPUT / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str, title: str, *, y: float = 1.06) -> None:
    ax.text(-0.10, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")
    ax.text(0.00, y, title, transform=ax.transAxes, fontsize=10.5, va="bottom")


def ordered_jitter(n: int, center: float) -> np.ndarray:
    if n <= 1:
        return np.array([center])
    return center + np.linspace(-0.11, 0.11, n)


def draw_group_points(
    ax: plt.Axes,
    data: pd.DataFrame,
    value_column: str,
    ylabel: str,
    title: str,
    *,
    proportion: bool = False,
) -> None:
    for index, group in enumerate(GROUP_ORDER):
        values = pd.to_numeric(data.loc[data["group"] == group, value_column], errors="coerce").dropna().to_numpy()
        if values.size == 0:
            continue
        ax.scatter(
            ordered_jitter(values.size, index),
            values,
            s=19,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.35,
            zorder=3,
        )
        mean = float(values.mean())
        sd = float(values.std(ddof=1)) if values.size > 1 else 0.0
        lower = max(0.0, mean - sd) if proportion else mean - sd
        upper = mean + sd
        ax.vlines(index, lower, upper, color="#444444", lw=0.9, zorder=2)
        ax.hlines([lower, upper], index - 0.06, index + 0.06, color="#444444", lw=0.9, zorder=2)
    ax.set_xticks(range(len(GROUP_ORDER)), GROUP_LABELS, rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title, pad=6)
    if proportion:
        ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#e7e7e7", linewidth=0.5, zorder=0)


def make_figure2() -> None:
    umap = pd.read_csv(INPUT / "figure2_umap.csv.gz", compression="gzip")
    dot = pd.read_csv(INPUT / "figure2_marker_dotplot.csv")
    composition = pd.read_csv(INPUT / "ec_subtype_composition_v4.csv")
    modules = pd.read_csv(INPUT / "module_scores_by_library_v4.csv")
    stats = pd.read_csv(INPUT / "module_scores_liao_stats_v4.csv").set_index("module")

    fig = plt.figure(figsize=(14.6, 10.2))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.86, 1.14],
        width_ratios=[0.97, 1.23],
        left=0.055,
        right=0.985,
        bottom=0.075,
        top=0.95,
        wspace=0.27,
        hspace=0.40,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    for subtype in SUBTYPE_ORDER:
        subset = umap.loc[umap["subtype"] == subtype]
        ax_a.scatter(subset["UMAP1"], subset["UMAP2"], s=2.3, alpha=0.72, color=SUBTYPE_COLORS[subtype], label=subtype)
    ax_a.set(xticks=[], yticks=[], xlabel="", ylabel="")
    ax_a.spines[["left", "bottom"]].set_visible(False)
    ax_a.legend(
        frameon=True,
        facecolor="white",
        framealpha=0.88,
        edgecolor="none",
        loc="upper left",
        fontsize=7.4,
        markerscale=3,
    )
    panel_label(ax_a, "A", "EC subtypes (marker-panel defined)")

    ax_b = fig.add_subplot(outer[0, 1])
    gene_order = list(dict.fromkeys(dot["gene"].astype(str)))
    subtype_y = {name: idx for idx, name in enumerate(SUBTYPE_ORDER[::-1])}
    gene_x = {name: idx for idx, name in enumerate(gene_order)}
    scatter = ax_b.scatter(
        dot["gene"].map(gene_x),
        dot["subtype"].map(subtype_y),
        s=18 + 80 * pd.to_numeric(dot["fraction"], errors="coerce").fillna(0),
        c=pd.to_numeric(dot["scaled_mean"], errors="coerce").fillna(0),
        cmap="viridis",
        vmin=0,
        vmax=1,
        edgecolor="#777777",
        linewidth=0.25,
    )
    ax_b.set_xticks(range(len(gene_order)), gene_order, rotation=90)
    ax_b.set_yticks(range(len(SUBTYPE_ORDER)), SUBTYPE_ORDER[::-1])
    ax_b.set_xlim(-0.7, len(gene_order) - 0.3)
    for boundary in (3.5, 7.5, 10.5):
        ax_b.axvline(boundary, color="#d9d9d9", lw=0.7)
    cbar = fig.colorbar(scatter, ax=ax_b, orientation="horizontal", fraction=0.09, pad=0.20, aspect=30)
    cbar.set_label("Scaled mean expression", fontsize=7.5)
    for fraction, x in zip((0.25, 0.50, 0.75, 1.0), (1.04, 1.10, 1.16, 1.22)):
        ax_b.scatter(x, 0.80, s=18 + 80 * fraction, transform=ax_b.transAxes, color="#777777", clip_on=False)
        ax_b.text(x, 0.69, f"{int(fraction * 100)}", transform=ax_b.transAxes, ha="center", fontsize=6.6)
    ax_b.text(1.13, 0.94, "Cells expressing\nmarker (%)", transform=ax_b.transAxes, ha="center", va="top", fontsize=7.2)
    panel_label(ax_b, "B", "Subtype marker panels")

    grid_c = GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[1, 0], hspace=0.52, wspace=0.30)
    axes_c: list[plt.Axes] = []
    for idx, subtype in enumerate(SUBTYPE_ORDER):
        ax = fig.add_subplot(grid_c[idx // 2, idx % 2])
        draw_group_points(ax, composition, SUBTYPE_COLUMNS[subtype], "% of ECs", subtype, proportion=True)
        axes_c.append(ax)
    panel_label(axes_c[0], "C", "Subtype composition by sampling unit (mean ± SD)", y=1.30)

    grid_d = GridSpecFromSubplotSpec(2, 3, subplot_spec=outer[1, 1], hspace=0.62, wspace=0.38)
    axes_d: list[plt.Axes] = []
    for idx, module in enumerate(MODULE_LABELS):
        ax = fig.add_subplot(grid_d[idx // 3, idx % 3])
        fdr = float(stats.loc[module, "Liao_KW_fdr"])
        draw_group_points(
            ax,
            modules,
            module,
            "Mean module score",
            f"{MODULE_LABELS[module]}\nLiao exact KW FDR = {fdr:.2f}",
        )
        axes_d.append(ax)
    panel_label(axes_d[0], "D", "Module scores by sampling unit (mean ± SD)", y=1.30)
    save(fig, "Figure2")


def make_figure3() -> None:
    de = pd.read_csv(INPUT / "de_ec_SONFH_vs_HOA_descriptive_v4.csv", index_col=0)
    effects = pd.read_csv(INPUT / "fig3_key_gene_effects_v4.csv", index_col=0)
    modules = pd.read_csv(INPUT / "module_scores_by_library_v4.csv")
    stats = pd.read_csv(INPUT / "module_scores_liao_stats_v4.csv").set_index("module")

    # Use an explicit layout here rather than constrained_layout. The latter
    # can assign a disproportionate inter-panel gutter to nested GridSpecs,
    # which previously left a large visual gap in panel C.
    fig = plt.figure(figsize=(12.4, 8.4))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.04, 0.96],
        width_ratios=[1.42, 0.78],
        left=0.065,
        right=0.985,
        bottom=0.105,
        top=0.945,
        wspace=0.20,
        hspace=0.42,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    x = np.log10(pd.to_numeric(de["baseMean"], errors="coerce").clip(lower=0) + 1)
    y = pd.to_numeric(de["log2FoldChange"], errors="coerce")
    # Panel A supplies transcriptome-wide context; panel B carries the complete
    # gene-by-gene readout.  Restricting direct labels here to the narrative
    # anchors prevents redundant callout ladders from obscuring the data cloud.
    ax_a.scatter(x, y, s=5.5, color="#aebdc6", alpha=0.24, linewidth=0, rasterized=True)
    ax_a.axhline(0, color="#33434d", lw=0.85, zorder=1)
    ax_a.grid(axis="y", color="#edf1f3", linewidth=0.55, zorder=0)

    targets = [gene for gene in effects.index if gene in de.index]
    focal = {"BAK1", "BAX", "SQSTM1", "CALCOCO2", "OPTN"}
    contextual = [gene for gene in targets if gene not in focal]
    ax_a.scatter(
        x.loc[contextual],
        y.loc[contextual],
        s=27,
        color="#4f8094",
        alpha=0.82,
        edgecolor="white",
        linewidth=0.55,
        zorder=4,
        label="Other prespecified genes",
    )

    marker_specs = {
        "BAK1": {"marker": "o", "face": "#d95f59", "edge": "white", "size": 38},
        "BAX": {"marker": "o", "face": "#d95f59", "edge": "white", "size": 38},
        "SQSTM1": {"marker": "o", "face": "#c94f4f", "edge": "#8f3333", "size": 48},
        "CALCOCO2": {"marker": "D", "face": "#dfa62f", "edge": "white", "size": 52},
        "OPTN": {"marker": "s", "face": "white", "edge": "#7562a3", "size": 48},
    }
    label_specs = {
        "BAK1": {"offset": (-8, 12), "ha": "right", "color": "#a94442"},
        "BAX": {"offset": (9, 12), "ha": "left", "color": "#a94442"},
        # Place OPTN above-left at approximately the same visual height as BAX.
        # This keeps its label and leader clear of the adjacent CALCOCO2 marker.
        "OPTN": {"offset": (-10, 25), "ha": "right", "color": "#69538f"},
        "CALCOCO2": {"offset": (-10, -17), "ha": "right", "color": "#a56f00"},
        "SQSTM1": {"offset": (10, -2), "ha": "left", "color": "#9f3d3d"},
    }
    for gene, spec in marker_specs.items():
        if gene not in de.index:
            continue
        ax_a.scatter(
            [x.loc[gene]],
            [y.loc[gene]],
            s=spec["size"],
            marker=spec["marker"],
            facecolor=spec["face"],
            edgecolor=spec["edge"],
            linewidth=1.15 if gene in {"SQSTM1", "OPTN"} else 0.65,
            zorder=6,
        )
        label = label_specs[gene]
        ax_a.annotate(
            gene,
            (float(x.loc[gene]), float(y.loc[gene])),
            xytext=label["offset"],
            textcoords="offset points",
            ha=label["ha"],
            va="center",
            fontsize=7.3,
            color=label["color"],
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.90},
            arrowprops={"arrowstyle": "-", "color": label["color"], "lw": 0.55, "shrinkA": 1.5, "shrinkB": 2},
            zorder=7,
        )
    ax_a.set_xlabel("Mean endothelial pseudobulk abundance, log10(baseMean + 1)")
    ax_a.set_ylabel("Descriptive log2 fold change")
    ax_a.legend(loc="lower right", frameon=False, fontsize=7.0, handletextpad=0.35, borderaxespad=0.6)
    ax_a.text(
        0.98,
        0.96,
        "DESCRIPTIVE ONLY\nSONFH participant mapping unavailable",
        transform=ax_a.transAxes,
        ha="right",
        va="top",
        fontsize=6.9,
        color="#65517f",
        linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.38", "facecolor": "#f3eff7", "edgecolor": "#d8cee4", "linewidth": 0.6},
    )
    ax_a.text(
        0.02,
        0.035,
        "Complete prespecified-gene effects are shown in B",
        transform=ax_a.transAxes,
        fontsize=6.9,
        color="#66747c",
    )
    panel_label(ax_a, "A", "Transcriptome-wide context of SONFH library effects")

    ax_b = fig.add_subplot(outer[0, 1])
    matrix = effects.to_numpy(dtype=float)
    limit = max(1.3, float(np.nanmax(np.abs(matrix))))
    heat = ax_b.imshow(matrix, cmap="coolwarm", vmin=-limit, vmax=limit, aspect="auto")
    ax_b.set_yticks(range(len(effects.index)), effects.index)
    for tick in ax_b.get_yticklabels():
        if tick.get_text() == "CALCOCO2":
            tick.set_color("#b97800")
            tick.set_fontweight("bold")
    ax_b.set_xticks(
        [0, 1],
        [
            "Liao ONFH 3A vs HOA\n(independent donors)",
            "SONFH vs HOA\n(descriptive cross-cohort;\nlibrary-level effects)",
        ],
        fontsize=7.0,
    )
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            ax_b.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=6.5, color="white" if abs(value) > 0.75 else "#333333")
    cbar = fig.colorbar(heat, ax=ax_b, fraction=0.045, pad=0.04)
    cbar.set_label("log2 fold change")
    ax_b.spines[["left", "bottom"]].set_visible(False)
    panel_label(ax_b, "B", "Heterogeneous stress and selective-clearance effects")

    grid_c = GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1, :], wspace=0.24)
    axes_c: list[plt.Axes] = []
    for idx, module in enumerate(("Mito_fission", "Mito_fusion", "Mitophagy_core")):
        ax = fig.add_subplot(grid_c[0, idx])
        fdr = float(stats.loc[module, "Liao_KW_fdr"])
        display_label = (
            "Selective-clearance transcripts"
            if module == "Mitophagy_core"
            else MODULE_LABELS[module]
        )
        draw_group_points(
            ax,
            modules,
            module,
            "Mean module score",
            f"{display_label}\nLiao exact KW FDR = {fdr:.2f}",
        )
        axes_c.append(ax)
    # One shared y-axis label is sufficient for these aligned small multiples
    # and keeps the three panels visually continuous.
    for ax in axes_c[1:]:
        ax.set_ylabel("")
    panel_label(axes_c[0], "C", "Sampling-unit mitochondrial module scores (mean ± SD)", y=1.22)
    save(fig, "Figure3")


def make_figure6() -> None:
    permutation = pd.read_csv(INPUT / "diag_permutation_v4.csv")
    stability = pd.read_csv(INPUT / "diag_feature_stability_v4.csv", index_col=0).iloc[:15]
    nested = pd.read_csv(INPUT / "diag_nested_cv_performance_v7.csv")
    comparator = pd.read_csv(INPUT / "diag_ma_comparator_repeat_performance_v7.csv")
    summary = json.loads((INPUT / "diag_summary_v4.json").read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(11.4, 8.2), constrained_layout=True)
    outer = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0])

    ax_a = fig.add_subplot(outer[0, 0])
    null = pd.to_numeric(permutation["perm_auc"], errors="coerce").dropna()
    observed = float(summary["aggregate_oof_auc"])
    p_value = float(summary["permutation_empirical_p"])
    ax_a.hist(null, bins=26, color="#56b4e9", edgecolor="white", linewidth=0.5)
    ax_a.axvline(observed, color="#d55e00", lw=2.2, label=f"Observed aggregate OOF AUC = {observed:.3f}")
    ax_a.set_xlabel("Aggregate repeated-OOF AUC under label permutation")
    ax_a.set_ylabel("Permutation count")
    ax_a.legend(frameon=False, fontsize=7.4, loc="upper left")
    panel_label(ax_a, "A", f"Label-permutation null ({len(null):,} runs; empirical p = {p_value:.4f})")

    ax_b = fig.add_subplot(outer[0, 1])
    ordered = stability.sort_values("selection_frequency", ascending=True)
    ax_b.barh(ordered.index, ordered["selection_frequency"], color="#0072b2")
    ax_b.set_xlim(0, 1.0)
    ax_b.set_xlabel("Selection frequency across 25 outer folds")
    panel_label(ax_b, "B", "Feature-selection stability")

    labels = [
        "EC mitochondrial\ncandidates",
        f"Ma et al. comparator\n({len(summary['ma2024_available_genes'])}/7 genes available)",
    ]

    def performance_panel(
        ax: plt.Axes,
        first: np.ndarray,
        second: np.ndarray,
        *,
        metric_label: str,
        baseline: float,
        panel: str,
        title: str,
        ylim: tuple[float, float],
    ) -> None:
        for idx, (values, color) in enumerate(((first, "#0072b2"), (second, "#cc79a7"))):
            ax.scatter(
                ordered_jitter(len(values), idx),
                values,
                s=44,
                color=color,
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
            mean = float(values.mean())
            sd = float(values.std(ddof=1))
            ax.errorbar(idx, mean, yerr=sd, color="#333333", capsize=6, lw=1.2, zorder=2)
        ax.axhline(baseline, color="#888888", lw=0.9, ls="--")
        ax.text(
            0.02,
            0.05,
            f"Dashed line: baseline = {baseline:.2f}",
            transform=ax.transAxes,
            fontsize=7.0,
            color="#666666",
        )
        ax.set_xlim(-0.42, 1.42)
        ax.set_ylim(*ylim)
        ax.set_xticks([0, 1], labels)
        ax.set_ylabel(metric_label)
        panel_label(ax, panel, title)

    ax_c = fig.add_subplot(outer[1, 0])
    nested_auc = pd.to_numeric(nested["AUC"], errors="coerce").dropna().to_numpy()
    comparator_auc = pd.to_numeric(comparator["AUC"], errors="coerce").dropna().to_numpy()
    performance_panel(
        ax_c,
        nested_auc,
        comparator_auc,
        metric_label="AUC per outer-CV repeat",
        baseline=0.50,
        panel="C",
        title="Repeat-level AUC (mean ± SD)",
        ylim=(0.45, 1.02),
    )

    ax_d = fig.add_subplot(outer[1, 1])
    nested_ap = pd.to_numeric(nested["average_precision"], errors="coerce").dropna().to_numpy()
    comparator_ap = pd.to_numeric(comparator["average_precision"], errors="coerce").dropna().to_numpy()
    performance_panel(
        ax_d,
        nested_ap,
        comparator_ap,
        metric_label="Average precision per outer-CV repeat",
        baseline=float(summary["average_precision_prevalence_baseline"]),
        panel="D",
        title="Repeat-level average precision (mean ± SD)",
        ylim=(0.70, 1.02),
    )

    save(fig, "Figure7")


def main() -> None:
    configure_style()
    make_figure2()
    make_figure3()
    make_figure6()
    print("Reviewer-corrected Figures 2, 3, and 6 written to figures/final/.")


if __name__ == "__main__":
    main()
