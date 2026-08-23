# -*- coding: utf-8 -*-
# W3-4: EC refined annotation + final gene-set scoring + sample-level stats + Fig 2
import time, warnings, json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as st
import scanpy as sc

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "analysis/py_fig2.log"
FIGS = ROOT / "figures" / "source"
FIGS.mkdir(parents=True, exist_ok=True)
if LOG.exists():
    LOG.unlink()

def lg(*a):
    s = time.strftime("%H:%M:%S") + " | " + " ".join(str(x) for x in a)
    print(s, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(s + "\n")

from v4_common import setup_plot
setup_plot()
import matplotlib.pyplot as plt
import seaborn as sns

lg("start")
ec = sc.read_h5ad(str(ROOT / "analysis/ec_annotated.h5ad"))
lg("EC:", ec.shape)

# ---------- 1. DEG-driven cluster markers ----------
sc.tl.rank_genes_groups(ec, groupby="leiden", method="wilcoxon")
mk = sc.get.rank_genes_groups_df(ec, group=None)
top = mk.groupby("group").head(8)
top.to_csv(ROOT / "analysis/ec_cluster_markers_top8.csv", index=False)
lg("cluster markers saved")

# ---------- 2. rescore with final gene sets ----------
GS = json.load(open(ROOT / "analysis/genesets_final.json", encoding="utf-8"))
for name, genes in GS.items():
    sc.tl.score_genes(ec, genes, score_name=name)
lg("final scores done:", len(GS), "sets")

# derived indices
ec.obs["Dyn_imbalance"] = ec.obs["Mito_fission"] - ec.obs["Mito_fusion"]
ec.obs["MQC_failure"] = ec.obs["Dyn_imbalance"] - ec.obs["Mitophagy_core"]

# ---------- 3. sample-level aggregation ----------
score_cols = ["Mito_fission", "Mito_fusion", "Dyn_imbalance", "Mitophagy_core",
              "MQC_failure", "Mito_biogenesis", "Mito_proteostasis", "OXPHOS",
              "ROS_defense", "mtDNA_release", "cGAS_STING", "EC_inflammation",
              "Angiogenesis", "YAP_mTOR"]
samp = ec.obs.groupby(["sample", "dataset", "group"], observed=True)[score_cols].mean().reset_index()
samp.to_csv(ROOT / "analysis/sample_level_scores.csv", index=False)
lg("sample-level table:", samp.shape)

# ---------- 4. statistics ----------
order = ["Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH"]
pairs = [("ONFH_3A", "HOA"), ("ONFH_4", "HOA"), ("ONFH_4", "ONFH_3A"),
         ("SONFH", "HOA"), ("SONFH", "ONFH_4"), ("ONFH_3A", "Healthy")]
stat_rows = []
for sc_col in score_cols:
    groups = [samp.loc[samp.group == g, sc_col].values for g in order]
    groups = [g for g in groups if len(g) > 0]
    try:
        kw = st.kruskal(*groups).pvalue
    except Exception:
        kw = np.nan
    rec = {"score": sc_col, "KW_p": round(kw, 4)}
    for a, b in pairs:
        xa = samp.loc[samp.group == a, sc_col].values
        xb = samp.loc[samp.group == b, sc_col].values
        if len(xa) >= 2 and len(xb) >= 2:
            p = st.mannwhitneyu(xa, xb, alternative="two-sided").pvalue
            eff = np.median(xa) - np.median(xb)
        else:
            p, eff = np.nan, np.nan
        rec[f"{a}_vs_{b}_p"] = round(p, 4) if p == p else ""
        rec[f"{a}_vs_{b}_eff"] = round(eff, 3) if eff == eff else ""
    stat_rows.append(rec)
stats_df = pd.DataFrame(stat_rows)
stats_df.to_csv(ROOT / "analysis/sample_level_stats.csv", index=False)
lg("stats saved")
print(stats_df[["score", "KW_p", "ONFH_3A_vs_HOA_p", "ONFH_4_vs_ONFH_3A_p",
                "SONFH_vs_HOA_p", "SONFH_vs_ONFH_4_p"]].to_string())

# ---------- 5. Fig 2 panels ----------
# Fig2a: EC UMAP subtypes (already have, regenerate cleaner)
sc.pl.umap(ec, color="EC_subtype", show=False, frameon=False, size=14)
plt.savefig(FIGS / "fig2a_ec_umap.png", dpi=300, bbox_inches="tight"); plt.close()

# Fig2b: dotplot of subtype markers
marker_genes = {"typeH": ["EMCN", "KDR", "CDH5", "DACH1"],
                "typeR": ["SMAD1", "PPARG", "NOTCH4", "COL4A1"],
                "venous_ACKR1": ["ACKR1", "SELE", "VCAM1", "PLAT"],
                "lymphatic": ["PROX1", "LYVE1", "PDPN"]}
marker_genes = {k: [g for g in v if g in ec.var_names] for k, v in marker_genes.items()}
marker_genes = {k: v for k, v in marker_genes.items() if v}
sc.pl.dotplot(ec, marker_genes, groupby="EC_subtype", show=False,
              standard_scale="var")
plt.savefig(FIGS / "fig2b_subtype_dotplot.png", dpi=300, bbox_inches="tight"); plt.close()
lg("fig2a/b done")

# Fig2c: sample-level boxplots for key scores
key_scores = ["Mito_fission", "Mito_fusion", "Dyn_imbalance", "Mitophagy_core",
              "cGAS_STING", "EC_inflammation"]
palette = dict(zip(order, sns.color_palette("Set2", len(order))))
fig, axes = plt.subplots(2, 3, figsize=(13, 8))
for ax, sc_col in zip(axes.ravel(), key_scores):
    sns.boxplot(data=samp, x="group", y=sc_col, order=order, ax=ax,
                palette=palette, width=0.55, showfliers=False)
    sns.stripplot(data=samp, x="group", y=sc_col, order=order, ax=ax,
                  color="black", size=4, jitter=0.15)
    ax.set_title(sc_col); ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)
    # annotate key p values
    row = stats_df[stats_df.score == sc_col].iloc[0]
    txt = []
    for pcol, lab in [("ONFH_3A_vs_HOA_p", "3A vs HOA"),
                      ("ONFH_4_vs_ONFH_3A_p", "4 vs 3A"),
                      ("SONFH_vs_HOA_p", "SONFH vs HOA")]:
        v = row.get(pcol)
        if v != "" and v == v:
            txt.append(f"{lab}: p={v}")
    ax.text(0.02, 0.98, "\n".join(txt), transform=ax.transAxes,
            fontsize=6.5, va="top")
fig.suptitle("EC mitochondrial dynamics / mitophagy / STING module scores (sample-level)")
fig.tight_layout()
fig.savefig(FIGS / "fig2c_scores_sample_level.png", dpi=300, bbox_inches="tight")
plt.close(fig)
lg("fig2c done")

# Fig2d: subtype x group mean score heatmap (key sets)
piv = ec.obs.groupby(["EC_subtype", "group"], observed=True)[
    ["Dyn_imbalance", "Mitophagy_core", "cGAS_STING", "OXPHOS"]].mean()
heat = piv.unstack("group")
heat.columns = [f"{s}|{g}" for s, g in heat.columns]
heat.to_csv(ROOT / "analysis/ec_subtype_group_scores.csv")
h2 = ec.obs.groupby(["EC_subtype", "group"], observed=True)["Dyn_imbalance"].mean().unstack()
h2 = h2.reindex(columns=[g for g in order if g in h2.columns])
fig, ax = plt.subplots(figsize=(7.5, 4.2))
sns.heatmap(h2, cmap="RdBu_r", center=0, annot=True, fmt=".2f", ax=ax,
            cbar_kws={"label": "Dynamics imbalance (fission-fusion)"})
ax.set_title("EC subtype × group: dynamics imbalance")
fig.savefig(FIGS / "fig2d_subtype_group_heatmap.png", dpi=300, bbox_inches="tight")
plt.close(fig)
lg("fig2d done")

# Fig2e: per-sample correlation: imbalance vs mitophagy, imbalance vs STING
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, ycol, ylab in [(axes[0], "Mitophagy_core", "Mitophagy score"),
                       (axes[1], "cGAS_STING", "cGAS-STING score")]:
    for g in order:
        d = samp[samp.group == g]
        ax.scatter(d["Dyn_imbalance"], d[ycol], label=g, s=42,
                   color=palette[g], edgecolor="k", linewidth=0.4)
    r, p = st.spearmanr(samp["Dyn_imbalance"], samp[ycol])
    ax.set_xlabel("Dynamics imbalance (fission-fusion)")
    ax.set_ylabel(ylab)
    ax.set_title(f"spearman r={r:.2f}, p={p:.3g}")
axes[1].legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
fig.tight_layout()
fig.savefig(FIGS / "fig2e_correlations.png", dpi=300, bbox_inches="tight")
plt.close(fig)
lg("fig2e done")

# Fig2f: subtype composition by group with sample points
comp_s = pd.crosstab([ec.obs.group, ec.obs["sample"]], ec.obs.EC_subtype,
                     normalize="index") * 100
comp_s = comp_s.reset_index().melt(id_vars=["group", "sample"],
                                   var_name="subtype", value_name="pct")
comp_s.to_csv(ROOT / "analysis/ec_subtype_composition_samples.csv", index=False)
fig, ax = plt.subplots(figsize=(9, 5.4))
sns.boxplot(data=comp_s, x="subtype", y="pct", hue="group", ax=ax,
            hue_order=order, palette=palette, width=0.6, showfliers=False)
sns.stripplot(data=comp_s, x="subtype", y="pct", hue="group", ax=ax,
              hue_order=order, dodge=True, palette="dark:k", size=3,
              jitter=0.12, legend=False)
ax.set_ylabel("% of EC"); ax.set_xlabel("")
ax.tick_params(axis="x", rotation=15)
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
fig.savefig(FIGS / "fig2f_subtype_composition.png", dpi=300, bbox_inches="tight")
plt.close(fig)
lg("fig2f done")

ec.write(str(ROOT / "analysis/ec_final.h5ad"))
lg("saved ec_final.h5ad; DONE")
