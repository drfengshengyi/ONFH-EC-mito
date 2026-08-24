#!/usr/bin/env python
"""Build Figure 4 entirely from versioned statistical tables.

Every displayed FDR is read from a frozen analysis output.  The figure no
longer depends on a manually assembled base PDF or post-hoc text overlays.
"""
from __future__ import annotations

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
LIAO_GROUPS = ["HOA", "ONFH_3A", "ONFH_4"]
LIAO_LABELS = ["HOA", "ARCO 3A", "ARCO 4"]
FOCUS_TFS = ["RELA", "NFKB1", "ATF4", "FOXO3", "STAT3", "HIF1A"]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.3,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.3,
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.75,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str, title: str, *, x: float = -0.12, y: float = 1.08) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")
    ax.text(0.0, y, title, transform=ax.transAxes, fontsize=11, va="bottom")


def fdr_text(value: float) -> str:
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10**exponent)
        return rf"FDR {coefficient:.2g} x $10^{{{exponent}}}$"
    if value < 0.01:
        return f"FDR {value:.4f}"
    return f"FDR {value:.3f}" if value < 0.1 else f"FDR {value:.2f}"


def ordered_jitter(n: int, center: float, width: float = 0.12) -> np.ndarray:
    if n <= 1:
        return np.array([center])
    return center + np.linspace(-width, width, n)


def draw_enrichment(ax: plt.Axes) -> None:
    gsea = pd.read_csv(INPUT / "gsea_ONFH3A_vs_HOA_H.csv")
    gsea = gsea.dropna(subset=["NES", "padj"]).copy()
    gsea["label"] = gsea["pathway"].str.replace("HALLMARK_", "", regex=False).str.replace("_", " ")
    significant = gsea.loc[gsea["padj"] < 0.05]
    positive = significant.loc[significant["NES"] > 0].nlargest(6, "NES")
    negative = significant.loc[significant["NES"] < 0].nsmallest(6, "NES")
    plot = pd.concat([negative, positive]).drop_duplicates("pathway").sort_values("NES")

    colors = np.where(plot["NES"] >= 0, "#df7621", "#2b8cbe")
    y = np.arange(len(plot))
    ax.barh(y, plot["NES"], color=colors, height=0.72)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks(y, plot["label"])
    ax.set_xlabel("Normalized enrichment score")
    ax.grid(axis="x", color="#e6e6e6", linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

    for yi, (_, row) in enumerate(plot.iterrows()):
        value = float(row["NES"])
        ax.text(
            value - 0.055 if value >= 0 else value + 0.055,
            yi,
            fdr_text(float(row["padj"])),
            va="center",
            ha="right" if value >= 0 else "left",
            color="white",
            fontsize=7.1,
            fontweight="bold",
        )
    panel_label(ax, "A", "Within-cohort Hallmark enrichment", x=-0.11)


def draw_module_axis(ax: plt.Axes, module: str, title: str, modules: pd.DataFrame, stats: pd.DataFrame) -> None:
    arrays: list[np.ndarray] = []
    for group in GROUP_ORDER:
        values = pd.to_numeric(modules.loc[modules["group"] == group, module], errors="coerce").dropna().to_numpy()
        arrays.append(values)
    box = ax.boxplot(
        arrays,
        positions=np.arange(len(GROUP_ORDER)),
        widths=0.58,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#333333", "linewidth": 1.0},
        boxprops={"linewidth": 0.8, "color": "#555555"},
        whiskerprops={"linewidth": 0.8, "color": "#555555"},
        capprops={"linewidth": 0.8, "color": "#555555"},
    )
    for patch, group in zip(box["boxes"], GROUP_ORDER):
        patch.set_facecolor(GROUP_COLORS[group])
        patch.set_alpha(0.82)
    for index, (group, values) in enumerate(zip(GROUP_ORDER, arrays)):
        ax.scatter(
            ordered_jitter(len(values), index, 0.09),
            values,
            s=13,
            color=GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.25,
            zorder=3,
        )
    fdr = float(stats.loc[module, "Liao_KW_fdr"])
    ax.set_title(f"{title}\nLiao exact KW FDR = {fdr:.2f}", pad=6, fontsize=8.5)
    ax.set_xticks(range(len(GROUP_ORDER)), GROUP_LABELS, rotation=25, ha="right")
    ax.set_ylabel("Mean module score")
    ax.grid(axis="y", color="#e7e7e7", linewidth=0.5, zorder=0)


def draw_modules(fig: plt.Figure, slot) -> None:
    modules = pd.read_csv(INPUT / "module_scores_by_library_v4.csv")
    stats = pd.read_csv(INPUT / "module_scores_liao_stats_v4.csv").set_index("module")
    sub = GridSpecFromSubplotSpec(1, 2, subplot_spec=slot, wspace=0.44)
    axes = [fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])]
    draw_module_axis(axes[0], "EC_inflammation", "EC inflammation", modules, stats)
    draw_module_axis(axes[1], "cGAS_STING", "cGAS-STING", modules, stats)
    panel_label(axes[0], "B", "Sampling-unit inflammatory scores", x=-0.24)


def draw_tf_activities(fig: plt.Figure, slot) -> None:
    activity = pd.read_csv(INPUT / "sample_level_tf_ulm_v4.csv")
    stats = pd.read_csv(INPUT / "tf_stats_ulm_v4.csv").set_index("TF")
    liao = activity.loc[activity["group"].isin(LIAO_GROUPS) & activity["independent_for_inference"].astype(bool)]
    sub = GridSpecFromSubplotSpec(2, 3, subplot_spec=slot, hspace=0.48, wspace=0.35)
    axes: list[plt.Axes] = []
    for index, tf in enumerate(FOCUS_TFS):
        ax = fig.add_subplot(sub[index // 3, index % 3])
        for group_index, group in enumerate(LIAO_GROUPS):
            values = pd.to_numeric(liao.loc[liao["group"] == group, tf], errors="coerce").dropna().to_numpy()
            ax.scatter(
                ordered_jitter(len(values), group_index, 0.08),
                values,
                s=19,
                color=GROUP_COLORS[group],
                edgecolor="white",
                linewidth=0.3,
                zorder=3,
            )
            if len(values):
                mean = float(values.mean())
                sd = float(values.std(ddof=1)) if len(values) > 1 else 0.0
                ax.errorbar(group_index, mean, yerr=sd, color="#444444", capsize=2.5, lw=0.9, zorder=2)
        ax.set_xticks(range(3), LIAO_LABELS, rotation=20, ha="right")
        ax.set_ylabel(tf)
        ax.set_title(f"{tf} (Liao FDR = {float(stats.loc[tf, 'Liao_KW_fdr']):.2f})", fontsize=7.7, pad=5)
        ax.grid(axis="y", color="#e7e7e7", linewidth=0.45, zorder=0)
        axes.append(ax)
    panel_label(axes[0], "C", "Signed TF activities in independent Liao donors", x=-0.29, y=1.27)


def draw_communication(ax: plt.Axes) -> None:
    top = pd.read_csv(INPUT / "comm_top_ONFH_4_vs_ONFH_3A_v4.csv").head(15)
    long = pd.read_csv(INPUT / "comm_scores_v4_long.csv.gz", compression="gzip")
    independent = long["independent_for_inference"].astype(str).str.lower().eq("true")
    long = long.loc[long["group"].isin(["ONFH_3A", "ONFH_4"]) & independent]
    keys = top[["sender", "receiver", "pair"]].drop_duplicates()
    plot = long.merge(keys, on=["sender", "receiver", "pair"], how="inner")
    plot["label"] = plot["sender"] + "->" + plot["receiver"] + ": " + plot["pair"]
    labels = (top["sender"] + "->" + top["receiver"] + ": " + top["pair"]).tolist()

    for index, label in enumerate(labels):
        center = len(labels) - 1 - index
        for group, offset in (("ONFH_3A", 0.13), ("ONFH_4", -0.13)):
            values = pd.to_numeric(
                plot.loc[(plot["label"] == label) & (plot["group"] == group), "score"], errors="coerce"
            ).dropna().to_numpy()
            y = ordered_jitter(len(values), center + offset, 0.045)
            ax.scatter(
                values,
                y,
                s=18,
                color=GROUP_COLORS[group],
                edgecolor="white",
                linewidth=0.3,
                label=("ARCO 3A" if group == "ONFH_3A" else "ARCO 4") if index == 0 else None,
                zorder=3,
            )
    ax.set_yticks(np.arange(len(labels)), labels[::-1])
    ax.set_xlabel("Ligand-receptor score in independent Liao donors")
    ax.grid(axis="x", color="#e7e7e7", linewidth=0.5, zorder=0)
    ax.legend(frameon=False, loc="lower right", fontsize=7.2)
    panel_label(ax, "D", "EC-centered communication effects by donor", x=-0.11)


def main() -> None:
    configure_style()
    fig = plt.figure(figsize=(16.0, 11.2))
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=[0.92, 1.18],
        width_ratios=[1.03, 0.97],
        left=0.105,
        right=0.985,
        bottom=0.075,
        top=0.94,
        hspace=0.37,
        wspace=0.30,
    )
    draw_enrichment(fig.add_subplot(outer[0, 0]))
    draw_modules(fig, outer[0, 1])
    draw_tf_activities(fig, outer[1, 0])
    draw_communication(fig.add_subplot(outer[1, 1]))

    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / "Figure4.pdf", bbox_inches="tight", pad_inches=0.10)
    fig.savefig(OUTPUT / "Figure4.png", dpi=300, bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    print("Figure 4 rebuilt from versioned GSEA, module, TF, and communication tables.")


if __name__ == "__main__":
    main()
