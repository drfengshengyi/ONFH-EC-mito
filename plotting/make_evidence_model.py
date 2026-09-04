# -*- coding: utf-8 -*-
"""Create an evidence-constrained summary schematic for manuscript v4."""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT_PNG = ROOT / "figures" / "final" / "SupplementaryFigureS1.png"
OUT_PDF = ROOT / "figures" / "final" / "SupplementaryFigureS1.pdf"
OUT_PNG.parent.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 14, "figure.dpi": 150,
})


def box(ax, xy, wh, text, face, edge="#2d3436", fontsize=9, lw=1.1):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
        facecolor=face, edgecolor=edge, linewidth=lw,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, wrap=True, linespacing=1.25)
    return patch


def arrow(ax, start, end, label=None, dashed=False, color="#636e72"):
    patch = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=12, linewidth=1.2,
        linestyle="--" if dashed else "-", color=color,
        connectionstyle="arc3,rad=0.0",
    )
    ax.add_patch(patch)
    if label:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.02,
                label, ha="center", va="bottom", fontsize=7.5, color=color,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5})


fig, ax = plt.subplots(figsize=(13.2, 7.2))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
ax.set_title("Evidence-constrained model of endothelial mitochondrial stress and selective-autophagy receptors in ONFH",
             pad=18, fontweight="bold")

ax.text(0.17, 0.91, "OBSERVED IN THE REANALYSED DATA", ha="center", va="center",
        fontsize=10, fontweight="bold", color="#1b4f72")
ax.text(0.50, 0.91, "NOT ESTABLISHED BY THESE DATA", ha="center", va="center",
        fontsize=10, fontweight="bold", color="#7d6608")
ax.text(0.83, 0.91, "TESTABLE NEXT STEPS", ha="center", va="center",
        fontsize=10, fontweight="bold", color="#196f3d")

box(ax, (0.04, 0.64), (0.26, 0.18),
    "SONFH libraries versus HOA\n(descriptive only)\n\nBAX/BAK1 and EIF4EBP1 higher;\nSQSTM1/CALCOCO2 lower; OPTN not\nsimilarly prioritized",
    "#d6eaf8", edge="#2874a6")
box(ax, (0.04, 0.38), (0.26, 0.18),
    "Independent Liao donors\n\nNo FDR-significant TF activity,\nYAP/TAZ signature, type H fraction,\nor communication contrast",
    "#eaf2f8", edge="#5dade2")
box(ax, (0.04, 0.12), (0.26, 0.18),
    "Peripheral serum transcriptomes\n\nRepeated nested-CV aggregate OOF\nAUC = 0.870; no external validation\nor fixed diagnostic signature",
    "#e8f8f5", edge="#17a589")

box(ax, (0.39, 0.64), (0.23, 0.18),
    "A receptor-specific causal mechanism\n\nSQSTM1: primary computational case\nCALCOCO2: secondary candidate\nOPTN: mechanistic-context control",
    "#fcf3cf", edge="#b7950b")
box(ax, (0.39, 0.38), (0.23, 0.18),
    "cGAS-STING activation\n\nCGAS and TMEM173 expression alone\ndoes not establish pathway activity",
    "#fdebd0", edge="#ca6f1e")
box(ax, (0.39, 0.12), (0.23, 0.18),
    "Stage progression or causal\nordering across cohorts\n\nSONFH libraries cannot be placed\non the Liao ARCO trajectory",
    "#fce4ec", edge="#c2185b")

box(ax, (0.70, 0.64), (0.26, 0.18),
    "Test the receptor system experimentally:\nparallel SQSTM1/CALCOCO2 loss-rescue,\nOPTN comparator, and single versus\ncombined perturbations",
    "#e9f7ef", edge="#239b56")
box(ax, (0.70, 0.38), (0.26, 0.18),
    "Measure mitochondrial respiration,\nmtDNA release and phospho-TBK1/IRF3;\nperturb cGAS/STING experimentally",
    "#e9f7ef", edge="#239b56")
box(ax, (0.70, 0.12), (0.26, 0.18),
    "Lock the classifier and threshold,\nthen test in a prospectively collected\nexternal steroid-exposed cohort",
    "#e9f7ef", edge="#239b56")

arrow(ax, (0.30, 0.73), (0.39, 0.73), "hypothesis", dashed=True, color="#a04000")
arrow(ax, (0.62, 0.73), (0.70, 0.73), "requires validation", dashed=True, color="#a04000")
arrow(ax, (0.30, 0.47), (0.39, 0.47), "not supported", dashed=True, color="#a04000")
arrow(ax, (0.62, 0.47), (0.70, 0.47), "direct test", dashed=True, color="#a04000")
arrow(ax, (0.30, 0.21), (0.39, 0.21), "cannot infer", dashed=True, color="#a04000")
arrow(ax, (0.62, 0.21), (0.70, 0.21), "validation design", dashed=True, color="#a04000")

ax.text(0.5, 0.035,
        "Solid statements summarize observed analysis outputs; dashed arrows denote hypotheses or missing evidence.",
        ha="center", va="center", fontsize=8, color="#555555")

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Wrote {OUT_PNG}")
print(f"Wrote {OUT_PDF}")
