# -*- coding: utf-8 -*-
"""Figure 3 v4: separate valid Liao-donor inference from SONFH-library effects.

The four public SONFH sequencing libraries represent three reported patients,
but their participant mapping is not public.  Consequently, SONFH contrasts
are plotted as descriptive effect sizes without p values or FDR labels.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import seaborn as sns

from v4_common import ANALYSIS, FIGS, save_figure, setup_plot

setup_plot()
import matplotlib.pyplot as plt


KEY_GENES = [
    "BAX", "BAK1", "SQSTM1", "CALCOCO2", "BNIP3L", "MAP1LC3B",
    "GABARAP", "TBK1", "TMEM173", "CGAS", "EIF4EBP1", "DNM1L",
    "OPA1", "MFN1", "MFN2",
]


def read_de(name: str) -> pd.DataFrame:
    df = pd.read_csv(ANALYSIS / name, index_col=0)
    df.index = df.index.astype(str).str.upper()
    return df


# Panel a: descriptive SONFH-library effect-size plot.  No p-value geometry.
sonfh = read_de("de_ec_SONFH_vs_HOA_descriptive_v4.csv")
sonfh["abundance"] = np.log10(sonfh["baseMean"].clip(lower=0) + 1)
sonfh["highlight"] = sonfh.index.isin(KEY_GENES)
fig, ax = plt.subplots(figsize=(7.4, 5.4))
ax.scatter(
    sonfh.loc[~sonfh["highlight"], "abundance"],
    sonfh.loc[~sonfh["highlight"], "log2FoldChange"],
    s=9, alpha=0.22, color="#8c96c6", linewidth=0,
)
ax.scatter(
    sonfh.loc[sonfh["highlight"], "abundance"],
    sonfh.loc[sonfh["highlight"], "log2FoldChange"],
    s=30, alpha=0.9, color="#d7301f", linewidth=0.3, edgecolor="white",
)
# The selected genes occupy a narrow band around zero.  A deterministic
# two-column layout is more legible and reproducible than initializing every
# label at its point and asking an automatic solver to untangle the cluster.
label_layout = {
    # Left label column, ordered from top to bottom.
    "BAK1": (1.72, 2.40, "right"),
    "EIF4EBP1": (1.72, 1.65, "right"),
    "MFN1": (1.72, 0.90, "right"),
    "OPA1": (1.72, 0.15, "right"),
    "MFN2": (1.72, -0.60, "right"),
    "GABARAP": (1.72, -1.35, "right"),
    "CGAS": (1.72, -2.10, "right"),
    # Right label column, ordered from top to bottom.
    "TBK1": (3.62, 2.65, "left"),
    "BAX": (3.62, 1.90, "left"),
    "TMEM173": (3.62, 1.15, "left"),
    "DNM1L": (3.62, 0.40, "left"),
    "MAP1LC3B": (3.62, -0.35, "left"),
    "BNIP3L": (3.62, -1.10, "left"),
    "SQSTM1": (3.62, -1.85, "left"),
    "CALCOCO2": (3.62, -2.60, "left"),
}
for gene, row in sonfh.loc[sonfh["highlight"]].iterrows():
    label_x, label_y, alignment = label_layout[gene]
    ax.annotate(
        gene,
        xy=(row["abundance"], row["log2FoldChange"]),
        xytext=(label_x, label_y),
        textcoords="data",
        ha=alignment,
        va="center",
        fontsize=7.2,
        color="#333333",
        zorder=4,
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white", alpha=0.90,
                  edgecolor="none"),
        arrowprops=dict(arrowstyle="-", lw=0.45, color="#6B7280",
                        shrinkA=1.5, shrinkB=2.5),
    )
ax.axhline(0, color="#333333", lw=0.8)
ax.set_xlabel("Mean endothelial pseudobulk abundance, log10(baseMean + 1)")
ax.set_ylabel("Descriptive log2 fold change")
ax.set_title("SONFH libraries versus Liao HOA donors\n(effect sizes only; participant mapping unavailable)")
save_figure(fig, FIGS / "fig3a_sonfh_effects_descriptive_v4.png")
plt.close(fig)


# Panel b: focused mitochondrial-quality-control effects across one valid and
# one descriptive contrast. Significance marks apply only to the Liao contrast.
liao = read_de("de_ec_ONFH_3A_vs_HOA_v4.csv")
genes = [g for g in KEY_GENES if g in liao.index and g in sonfh.index]
heat = pd.DataFrame({
    "Liao ONFH 3A vs HOA\n(independent donors)": liao.loc[genes, "log2FoldChange"],
    "SONFH vs HOA\n(descriptive libraries)": sonfh.loc[genes, "log2FoldChange"],
}, index=genes)
heat.to_csv(ANALYSIS / "fig3_key_gene_effects_v4.csv")
fig, ax = plt.subplots(figsize=(6.4, 6.2))
lim = max(1.0, float(np.nanpercentile(np.abs(heat.values), 95)))
sns.heatmap(
    heat, cmap="vlag", center=0, vmin=-lim, vmax=lim, annot=True, fmt=".2f",
    annot_kws={"fontsize": 7}, linewidths=0.5, linecolor="white",
    cbar_kws={"label": "log2 fold change", "shrink": 0.75}, ax=ax,
)
for row_i, gene in enumerate(genes):
    if pd.notna(liao.loc[gene, "padj"]) and liao.loc[gene, "padj"] < 0.05:
        ax.text(0.92, row_i + 0.5, "*", ha="center", va="center", fontsize=12, color="black")
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("Selected endothelial gene effects\nNo selected gene reached FDR < 0.05 in the Liao contrast")
save_figure(fig, FIGS / "fig3b_key_gene_effects_v4.png")
plt.close(fig)


# Panel c: pathway inference is restricted to independent Liao participants.
gsea = pd.read_csv(ANALYSIS / "gsea_ONFH3A_vs_HOA_H.csv")
gsea = gsea.dropna(subset=["NES", "padj"]).copy()
gsea["label"] = gsea["pathway"].str.replace("HALLMARK_", "", regex=False).str.replace("_", " ")
sig = gsea[gsea["padj"] < 0.05]
if len(sig) < 10:
    sig = gsea
pos = sig.sort_values("NES", ascending=False).head(6)
neg = sig.sort_values("NES", ascending=True).head(6)
plot = pd.concat([neg, pos]).drop_duplicates("pathway").sort_values("NES")
plot.to_csv(ANALYSIS / "gsea_ONFH3A_vs_HOA_H_selected_v4.csv", index=False)
fig, ax = plt.subplots(figsize=(7.2, 5.8))
colors = np.where(plot["NES"] >= 0, "#cb181d", "#2171b5")
ax.barh(plot["label"], plot["NES"], color=colors, alpha=0.88)
ax.axvline(0, color="#333333", lw=0.8)
def fdr_label(value):
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10 ** exponent)
        return rf"FDR ${coefficient:.1f}\times10^{{{exponent}}}$"
    return f"FDR {value:.2g}"
for y, (_, row) in enumerate(plot.iterrows()):
    x = row["NES"]
    ax.text(x - 0.06 if x >= 0 else x + 0.06, y, fdr_label(row["padj"]),
            va="center", ha="right" if x >= 0 else "left", fontsize=8.2,
            color="white", fontweight="bold")
ax.set_xlabel("Normalized enrichment score")
ax.set_ylabel("")
ax.set_title("Preranked competitive Hallmark enrichment: ARCO 3A versus HOA")
ax.text(
    0.01, 0.01,
    "fgseaMultilevel gene-set null; participant labels were not permuted",
    transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color="#444444",
)
save_figure(fig, FIGS / "fig3c_liao_gsea_v4.png")
plt.close(fig)

print("Wrote Figure 3 v4 panels and source tables.")
