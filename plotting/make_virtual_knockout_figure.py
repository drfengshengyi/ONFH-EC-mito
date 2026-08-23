"""Publication Figure 5 with official scTenifoldKnk network and manifold outputs.

This figure borrows the visual grammar of established in-silico perturbation
papers while preserving the meaning of the present analysis:

* the network panel contains actual SQSTM1 outgoing edges from the official
  scTenifoldKnk wild-type tensor networks;
* the vector panels show WT-to-virtual-KO displacement of genes in the official
  manifold alignment, not cell-state trajectories;
* the pathway curves display cumulative recovery along perturbation ranks, with
  p values/FDR taken from the prespecified matched-null analysis.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = ROOT / "figures" / "final"
OUT.mkdir(parents=True, exist_ok=True)
FIGURE_DATA = RESULTS / "official_r_vko_figure_data"
PRIMARY = RESULTS / "official_r_vko_manuscript"
DEFAULT = RESULTS / "official_r_vko_official_default"
TARGET = "SQSTM1"

COLORS = {
    "EC inflammation": "#E45756",
    "OXPHOS": "#4C78A8",
    "Angiogenesis": "#59A14F",
    "Mitophagy": "#F2CF5B",
    "Mito stress": "#B279A2",
    "Other": "#D9D9D9",
}

sns.set_theme(style="white", context="paper")
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.titleweight": "normal",
        "figure.dpi": 150,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label.upper(),
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        va="top",
        ha="left",
    )


def rounded_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    text: str,
    facecolor: str,
    fontsize: float = 8.0,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        transform=ax.transAxes,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=facecolor,
        edgecolor="#4B5563",
        linewidth=1.0,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize,
        linespacing=1.22,
    )


def workflow_panel(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Official-package virtual-knockout design", loc="left", pad=5)

    rounded_box(
        ax,
        (0.02, 0.76),
        0.25,
        0.12,
        "HOA1: 26 ECs\nExcluded from network",
        "#E5E7EB",
        fontsize=7.7,
    )
    rounded_box(
        ax,
        (0.02, 0.51),
        0.25,
        0.17,
        "HOA2 control ECs\nn = 1,043\nSQSTM1+ = 70.8%",
        "#DDEBF7",
        fontsize=7.7,
    )
    rounded_box(
        ax,
        (0.02, 0.26),
        0.25,
        0.17,
        "HOA3 control ECs\nn = 760\nSQSTM1+ = 72.2%",
        "#E2F0D9",
        fontsize=7.7,
    )
    rounded_box(
        ax,
        (0.38, 0.40),
        0.26,
        0.27,
        "Official R packages\nscTenifoldKnk v1.1\nscTenifoldNet v1.4\n300 shared genes",
        "#FFF2CC",
    )
    rounded_box(
        ax,
        (0.73, 0.40),
        0.25,
        0.27,
        "Virtual SQSTM1 KO\nSQSTM1 network row = 0\nWT--KO manifold alignment\nDifferential regulation",
        "#FCE4D6",
    )

    arrow = dict(arrowstyle="-|>", lw=1.5, color="#4B5563", mutation_scale=12)
    ax.annotate("", (0.38, 0.57), (0.27, 0.595), xycoords=ax.transAxes, arrowprops=arrow)
    ax.annotate("", (0.38, 0.49), (0.27, 0.345), xycoords=ax.transAxes, arrowprops=arrow)
    ax.annotate("", (0.73, 0.535), (0.64, 0.535), xycoords=ax.transAxes, arrowprops=arrow)
    ax.text(0.295, 0.82, "not modeled", transform=ax.transAxes, fontsize=6.6, color="#6B7280", va="center")

    ax.text(
        0.01,
        0.08,
        "Parameter audit",
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        0.01,
        0.015,
        (
            "Primary: 20 networks, q=0.95, 30 manifold dimensions   |   "
            "Package-default sensitivity: 10 networks, q=0.90, 2 dimensions\n"
            "No downstream gene or pathway was FDR-significant in both donors under both profiles."
        ),
        transform=ax.transAxes,
        fontsize=7.5,
        va="bottom",
    )
    panel_label(ax, "A")


def read_primary_results() -> dict[str, pd.DataFrame]:
    return {
        donor: pd.read_csv(PRIMARY / f"vko_sqstm1_{donor}_official_r.csv").set_index("Gene")
        for donor in ("hoa2", "hoa3")
    }


def pathway_sets() -> dict[str, set[str]]:
    table = pd.read_csv(PRIMARY / "vko_sqstm1_pathway_enrichment_official_r.csv")
    return {
        row.pathway: set(str(row.genes).split(";"))
        for row in table.itertuples(index=False)
    }


def gene_category(gene: str, pathways: dict[str, set[str]]) -> str:
    if gene in pathways["EC_inflammation"]:
        return "EC inflammation"
    if gene in pathways["OXPHOS"]:
        return "OXPHOS"
    if gene in pathways["Angiogenesis"]:
        return "Angiogenesis"
    if gene in pathways["Mitophagy_core"]:
        return "Mitophagy"
    stress_sets = (
        pathways["ROS_defense"]
        | pathways["mtDNA_release"]
        | pathways["cGAS_STING"]
        | pathways["Mito_fission"]
        | pathways["Mito_proteostasis"]
    )
    if gene in stress_sets:
        return "Mito stress"
    return "Other"


def network_table() -> pd.DataFrame:
    edges = []
    for donor in ("hoa2", "hoa3"):
        frame = pd.read_csv(FIGURE_DATA / f"sqstm1_outgoing_edges_{donor}.csv")
        frame = frame[["gene", "wt_outgoing_weight", "absolute_wt_weight"]].rename(
            columns={
                "wt_outgoing_weight": f"weight_{donor}",
                "absolute_wt_weight": f"abs_{donor}",
            }
        )
        edges.append(frame)
    merged = edges[0].merge(edges[1], on="gene", how="inner")
    for donor in ("hoa2", "hoa3"):
        maximum = merged[f"abs_{donor}"].max()
        merged[f"scaled_{donor}"] = merged[f"abs_{donor}"] / maximum
    merged["consensus_edge_score"] = merged[["scaled_hoa2", "scaled_hoa3"]].mean(axis=1)
    merged["edge_sign"] = np.select(
        [
            (merged["weight_hoa2"] > 0) & (merged["weight_hoa3"] > 0),
            (merged["weight_hoa2"] < 0) & (merged["weight_hoa3"] < 0),
        ],
        ["positive", "negative"],
        default="discordant",
    )
    return merged.sort_values("consensus_edge_score", ascending=False)


def selected_network_nodes(table: pd.DataFrame, pathways: dict[str, set[str]]) -> list[str]:
    stable = pd.read_csv(RESULTS / "official_r_vko_cross_profile_gene_summary.csv").head(12)["Gene"].tolist()
    selected: list[str] = []

    def add(gene: str) -> None:
        if gene != TARGET and gene not in selected:
            selected.append(gene)

    for gene in stable:
        add(gene)
    for gene in table.head(30)["gene"]:
        add(gene)
    category_candidates: dict[str, list[str]] = defaultdict(list)
    for gene in table["gene"]:
        category_candidates[gene_category(gene, pathways)].append(gene)
    for category in COLORS:
        for gene in category_candidates[category][:4]:
            add(gene)
    return selected[:42]


def network_panel(ax: plt.Axes) -> None:
    pathways = pathway_sets()
    table = network_table().set_index("gene")
    nodes = selected_network_nodes(table.reset_index(), pathways)
    results = read_primary_results()
    stable = set(
        pd.read_csv(RESULTS / "official_r_vko_cross_profile_gene_summary.csv").head(12)["Gene"]
    )
    significant = {
        gene
        for donor in ("hoa2", "hoa3")
        for gene, row in results[donor].iterrows()
        if gene != TARGET and float(row["adjusted_p_value"]) < 0.05
    }

    graph = nx.DiGraph()
    graph.add_node(TARGET)
    for gene in nodes:
        graph.add_edge(TARGET, gene)

    ordered = sorted(
        nodes,
        key=lambda gene: (
            list(COLORS).index(gene_category(gene, pathways)),
            -float(table.loc[gene, "consensus_edge_score"]),
        ),
    )
    positions = {TARGET: np.array([0.0, 0.0])}
    for index, gene in enumerate(ordered):
        angle = 2 * np.pi * index / len(ordered) + np.pi / 2
        radius = 0.80 if index % 2 == 0 else 1.00
        positions[gene] = np.array([radius * np.cos(angle), radius * np.sin(angle)])

    edge_colors = {
        "positive": "#D95F59",
        "negative": "#4C78A8",
        "discordant": "#9CA3AF",
    }
    nx.draw_networkx_edges(
        graph,
        positions,
        ax=ax,
        edge_color=[edge_colors[table.loc[gene, "edge_sign"]] for gene in nodes],
        width=[0.55 + 2.1 * float(table.loc[gene, "consensus_edge_score"]) for gene in nodes],
        alpha=0.62,
        arrows=True,
        arrowstyle="-|>",
        arrowsize=7,
        node_size=1,
        min_target_margin=8,
    )

    for gene in nodes:
        score = float(table.loc[gene, "consensus_edge_score"])
        category = gene_category(gene, pathways)
        nx.draw_networkx_nodes(
            graph,
            positions,
            nodelist=[gene],
            node_color=[COLORS[category]],
            node_size=65 + 230 * score,
            edgecolors="#1F2937" if gene in significant else "#FFFFFF",
            linewidths=1.5 if gene in significant else 0.6,
            ax=ax,
        )
    nx.draw_networkx_nodes(
        graph,
        positions,
        nodelist=[TARGET],
        node_color=["#F6D6D3"],
        node_size=1500,
        edgecolors="#B23A33",
        linewidths=2.0,
        ax=ax,
    )
    ax.text(0, 0, "SQSTM1\nWT source", ha="center", va="center", fontsize=9, fontweight="bold")

    label_genes = {
        "NFKBIA", "ICAM1", "TNFAIP3", "SELE", "MT-ND1", "MT-ATP6",
        "MT-ND4", "MT-CO2", "MT-CO1", "VWF", "ENG", "SERPINE1",
        "CXCL2", "MCL1", "KDM6B", "TIPARP",
    } | significant
    label_offsets = {
        "MT-ATP6": (-0.055, 0.055),
        "MT-ND1": (0.045, -0.035),
        "MCL1": (-0.050, 0.055),
        "TNFAIP3": (-0.025, 0.030),
        "SELE": (0.025, -0.015),
    }
    for gene in nodes:
        if gene not in label_genes:
            continue
        x, y = positions[gene]
        label_radius = 1.10 if np.hypot(x, y) > 0.9 else 0.91
        angle = np.arctan2(y, x)
        tx, ty = label_radius * np.cos(angle), label_radius * np.sin(angle)
        tx += label_offsets.get(gene, (0.0, 0.0))[0]
        ty += label_offsets.get(gene, (0.0, 0.0))[1]
        ax.text(
            tx,
            ty,
            gene,
            fontsize=6.4,
            fontstyle="italic",
            fontweight="bold" if gene in stable else "normal",
            ha="left" if tx >= 0 else "right",
            va="center",
        )

    category_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color, markeredgecolor="none", markersize=6, label=label)
        for label, color in COLORS.items()
    ]
    sign_handles = [
        Line2D([0], [0], color=edge_colors[sign], lw=1.7, label=f"{sign.capitalize()} edge sign")
        for sign in ("positive", "negative", "discordant")
    ]
    ax.legend(
        handles=category_handles + sign_handles,
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.085),
        frameon=False,
        fontsize=6.4,
        handlelength=1.5,
        columnspacing=1.0,
    )
    ax.set_title("SQSTM1 WT neighborhood removed by virtual KO", loc="left", pad=5)
    ax.set_xlim(-1.30, 1.30)
    ax.set_ylim(-1.38, 1.24)
    ax.set_aspect("equal")
    ax.axis("off")
    panel_label(ax, "B", x=-0.05)


def manifold_subpanel(ax: plt.Axes, donor: str) -> None:
    manifold = pd.read_csv(FIGURE_DATA / f"wt_ko_manifold_{donor}.csv")
    wide = manifold.pivot(index="gene", columns="state", values=["PC1", "PC2"])
    results = pd.read_csv(PRIMARY / f"vko_sqstm1_{donor}_official_r.csv").set_index("Gene")
    genes = wide.index.intersection(results.index).difference([TARGET])
    wide = wide.loc[genes]
    dx = wide[("PC1", "KO")] - wide[("PC1", "WT")]
    dy = wide[("PC2", "KO")] - wide[("PC2", "WT")]

    # The official WT--KO coordinates differ at a much smaller numerical scale
    # than the manifold itself.  As in vector-field perturbation figures, a
    # single donor-specific display factor is applied to every arrow; direction
    # and relative magnitude are preserved.  Robust axis limits prevent a few
    # manifold outliers from flattening the central vector field.
    x0 = wide[("PC1", "WT")]
    y0 = wide[("PC2", "WT")]
    x_low, x_high = x0.quantile([0.01, 0.99])
    y_low, y_high = y0.quantile([0.01, 0.99])
    x_pad = max((x_high - x_low) * 0.10, 1e-4)
    y_pad = max((y_high - y_low) * 0.10, 1e-4)
    robust_span = max(x_high - x_low, y_high - y_low)
    displacement = np.hypot(dx, dy)
    reference = max(float(displacement.quantile(0.95)), 1e-12)
    arrow_factor = 0.045 * robust_span / reference
    display_dx = dx * arrow_factor
    display_dy = dy * arrow_factor

    visible = (
        x0.between(x_low - x_pad, x_high + x_pad)
        & y0.between(y_low - y_pad, y_high + y_pad)
    )
    ax.scatter(x0[visible], y0[visible], s=7, c="#D1D5DB", alpha=0.52, linewidths=0)
    ax.quiver(
        x0[visible],
        y0[visible],
        display_dx[visible],
        display_dy[visible],
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.0028,
        headwidth=3.2,
        headlength=4.2,
        color="#9CA3AF",
        alpha=0.35,
    )
    top = results.drop(index=TARGET, errors="ignore").sort_values("rank").head(10)
    top_genes = [gene for gene in top.index if gene in wide.index]
    fdr_genes = [
        gene
        for gene in top_genes
        if gene != TARGET and float(results.loc[gene, "adjusted_p_value"]) < 0.05
    ]
    ax.quiver(
        wide.loc[top_genes, ("PC1", "WT")],
        wide.loc[top_genes, ("PC2", "WT")],
        display_dx.loc[top_genes],
        display_dy.loc[top_genes],
        angles="xy",
        scale_units="xy",
        scale=1,
        width=0.007,
        headwidth=4.0,
        headlength=5.0,
        color="#4C78A8" if donor == "hoa2" else "#59A14F",
        alpha=0.9,
    )
    if fdr_genes:
        ax.scatter(
            x0.loc[fdr_genes] + display_dx.loc[fdr_genes],
            y0.loc[fdr_genes] + display_dy.loc[fdr_genes],
            s=38,
            facecolor="none",
            edgecolor="#D62728",
            linewidth=1.1,
            zorder=5,
        )
    ax.set_xlim(x_low - x_pad, x_high + x_pad)
    ax.set_ylim(y_low - y_pad, y_high + y_pad)

    # Label only FDR-significant downstream genes and place their names in a
    # dedicated column with leader lines.  This prevents the dense central
    # manifold from obscuring gene names.
    label_genes = [gene for gene in results.sort_values("rank").index if gene in fdr_genes]
    if donor == "hoa2":
        label_y = np.linspace(0.91, 0.70, max(len(label_genes), 1))
    else:
        label_y = np.linspace(0.54, 0.45, max(len(label_genes), 1))
    for gene, y_text in zip(label_genes, label_y):
        x = float(x0.loc[gene] + display_dx.loc[gene])
        y = float(y0.loc[gene] + display_dy.loc[gene])
        ax.annotate(
            gene,
            xy=(x, y),
            xycoords="data",
            xytext=(0.63, float(y_text)),
            textcoords=ax.transAxes,
            fontsize=6.3,
            fontstyle="italic",
            ha="left",
            va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=0.35),
            arrowprops=dict(arrowstyle="-", color="#6B7280", lw=0.55, shrinkA=1, shrinkB=1),
        )
    ax.text(
        0.02,
        0.02,
        f"Arrows enlarged {arrow_factor:,.0f}x for display",
        transform=ax.transAxes,
        fontsize=6.0,
        color="#4B5563",
        ha="left",
        va="bottom",
    )
    ax.set_title(donor.upper(), fontsize=9, pad=2)
    ax.set_xlabel("Gene-manifold PC1", fontsize=7)
    ax.set_ylabel("Gene-manifold PC2", fontsize=7)
    ax.tick_params(labelsize=6)
    sns.despine(ax=ax)


def manifold_panel(fig: plt.Figure, spec) -> None:
    subgrid = spec.subgridspec(1, 2, wspace=0.30)
    axes = [fig.add_subplot(subgrid[0, index]) for index in range(2)]
    manifold_subpanel(axes[0], "hoa2")
    manifold_subpanel(axes[1], "hoa3")
    axes[0].text(
        0.0,
        1.18,
        "WT-to-virtual-KO displacement in the aligned gene manifold",
        transform=axes[0].transAxes,
        fontsize=11,
        ha="left",
        va="bottom",
    )
    axes[0].text(
        0.0,
        -0.18,
        (
            "Gray arrows: downstream genes; colored arrows: top 10 ranks; red rings/labels: FDR<0.05.\n"
            "Uniform magnification preserves direction/relative size; these are gene-manifold displacements, not cell-fate vectors."
        ),
        transform=axes[0].transAxes,
        fontsize=6.8,
        ha="left",
        va="top",
    )
    panel_label(axes[0], "C", x=-0.20, y=1.20)


def recovery_curve(ax: plt.Axes, pathway: str, donor: str, display: str, color: str) -> None:
    results = pd.read_csv(PRIMARY / f"vko_sqstm1_{donor}_official_r.csv").sort_values("rank")
    pathway_table = pd.read_csv(PRIMARY / "vko_sqstm1_pathway_enrichment_official_r.csv").set_index("pathway")
    genes = set(str(pathway_table.loc[pathway, "genes"]).split(";"))
    genes.discard(TARGET)
    member_ranks = sorted(results.loc[results["Gene"].isin(genes), "rank"].astype(int).tolist())
    n_total = len(results)
    x = np.concatenate(([0.0], np.asarray(member_ranks) / n_total, [1.0]))
    y = np.concatenate(([0.0], np.arange(1, len(member_ranks) + 1) / len(member_ranks), [1.0]))
    ax.step(x, y, where="post", color=color, lw=1.8)
    ax.fill_between(x, y, x, step="post", where=y >= x, color=color, alpha=0.12)
    ax.plot([0, 1], [0, 1], color="#9CA3AF", ls="--", lw=0.9)
    ax.vlines(np.asarray(member_ranks) / n_total, -0.055, 0.0, color="#111827", lw=0.55)
    p_value = float(pathway_table.loc[pathway, f"empirical_p_{donor}"])
    fdr = float(pathway_table.loc[pathway, f"empirical_fdr_{donor}"])
    ax.set_title(f"{display}\n{donor.upper()} primary", fontsize=9, fontweight="bold", pad=4)
    ax.text(
        0.03,
        0.95,
        f"empirical p={p_value:.3f}\nBH FDR={fdr:.3f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.07, 1.02)
    ax.set_xlabel("Perturbation rank percentile", fontsize=7)
    ax.tick_params(labelsize=6)
    sns.despine(ax=ax)


def pathway_panel(fig: plt.Figure, spec) -> None:
    subgrid = spec.subgridspec(1, 3, wspace=0.30)
    axes = [fig.add_subplot(subgrid[0, index]) for index in range(3)]
    recovery_curve(axes[0], "EC_inflammation", "hoa3", "EC inflammation", COLORS["EC inflammation"])
    recovery_curve(axes[1], "OXPHOS", "hoa2", "Oxidative phosphorylation", COLORS["OXPHOS"])
    recovery_curve(axes[2], "Mitophagy_core", "hoa3", "Mitophagy core", COLORS["Mitophagy"])
    axes[0].set_ylabel("Cumulative fraction of pathway genes", fontsize=7)
    axes[1].set_ylabel("")
    axes[2].set_ylabel("")
    axes[0].text(
        0.0,
        1.18,
        "Matched-null pathway rank recovery",
        transform=axes[0].transAxes,
        fontsize=11,
        ha="left",
        va="bottom",
    )
    axes[0].text(
        0.0,
        -0.18,
        (
            "Curves show pathway-member recovery along official perturbation ranks; p/FDR use 20,000 matched null sets.\n"
            "These descriptive curves are not GSEA enrichment-score curves."
        ),
        transform=axes[0].transAxes,
        fontsize=6.8,
        ha="left",
        va="top",
    )
    panel_label(axes[0], "D", x=-0.24, y=1.20)


def main() -> None:
    fig = plt.figure(figsize=(16.0, 10.4), facecolor="white")
    grid = fig.add_gridspec(
        2,
        5,
        width_ratios=[1, 1, 1, 1, 1],
        height_ratios=[1.08, 0.92],
        left=0.055,
        right=0.98,
        top=0.97,
        bottom=0.12,
        hspace=0.30,
        wspace=0.36,
    )
    workflow_panel(fig.add_subplot(grid[0, :2]))
    network_panel(fig.add_subplot(grid[0, 2:]))
    manifold_panel(fig, grid[1, :2])
    pathway_panel(fig, grid[1, 2:])

    png_out = OUT / "Figure5.png"
    pdf_out = OUT / "Figure5.pdf"
    fig.savefig(png_out, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png_out)
    print(pdf_out)


if __name__ == "__main__":
    main()
