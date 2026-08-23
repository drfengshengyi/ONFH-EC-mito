# -*- coding: utf-8 -*-
"""Figure 5 v4: correctly signed DoRothEA ULM activity.

Key safeguards
--------------
* Interaction signs come from ``consensus_stimulation`` and
  ``consensus_inhibition``; ``consensus_direction`` is never used as a sign.
* Activity is the t statistic from decoupler's actual univariate linear model.
* Formal tests are restricted to the nine independent Liao HOA/ARCO donors.
* GSE290411's four sequencing libraries are retained for visualization only,
  because the public record reports three participants without a library map.
"""
from __future__ import annotations

import time
from pathlib import Path

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from v4_common import (
    ANALYSIS,
    FIGS,
    GROUP_ORDER,
    LIAO_INFERENCE_GROUPS,
    ROOT,
    attach_v4_metadata,
    bh_fdr,
    exact_kruskal_p,
    exact_mwu_summary,
    save_figure,
    setup_plot,
)

setup_plot()
import matplotlib.pyplot as plt
import seaborn as sns

LOG = ANALYSIS / "py_fig5_v4.log"


def lg(*parts) -> None:
    line = time.strftime("%H:%M:%S") + " | " + " ".join(map(str, parts))
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


if LOG.exists():
    LOG.unlink()
lg("start")

ec = attach_v4_metadata(sc.read_h5ad(ANALYSIS / "ec_final.h5ad"))
lg("EC", ec.shape)

# Build an auditable, non-ambiguous signed network.
dor = pd.read_csv(ROOT / "data" / "dorothea_ABC.tsv", sep="\t")
stim = as_bool(dor["consensus_stimulation"])
inhib = as_bool(dor["consensus_inhibition"])
directed = as_bool(dor["consensus_direction"])
dor["weight"] = np.select([stim & ~inhib, inhib & ~stim], [1.0, -1.0], default=np.nan)
dor["sign_status"] = np.select(
    [stim & ~inhib, inhib & ~stim, stim & inhib],
    ["stimulation", "inhibition", "ambiguous_both"],
    default="unknown",
)

audit = (
    dor.groupby(["consensus_direction", "sign_status"], dropna=False)
    .size()
    .rename("n_edges")
    .reset_index()
)
audit.to_csv(ANALYSIS / "dorothea_sign_audit_v4.csv", index=False)

net0 = dor.loc[
    directed & dor["weight"].notna(),
    ["source_genesymbol", "target_genesymbol", "weight"],
].rename(columns={"source_genesymbol": "source", "target_genesymbol": "target"})
net0 = net0[net0["target"].isin(ec.var_names)].drop_duplicates()

# If the same TF-target edge has conflicting consensus signs, discard it.
edge_nsign = net0.groupby(["source", "target"])["weight"].nunique()
conflicting = edge_nsign[edge_nsign > 1].index
if len(conflicting):
    conflict_frame = pd.DataFrame(conflicting.tolist(), columns=["source", "target"])
    conflict_frame.to_csv(ANALYSIS / "dorothea_conflicting_edges_v4.csv", index=False)
    net0 = net0.merge(conflict_frame.assign(_conflict=True), how="left", on=["source", "target"])
    net0 = net0.loc[net0["_conflict"].isna(), ["source", "target", "weight"]]

target_counts = net0.groupby("source")["target"].nunique()
keep_sources = target_counts[target_counts >= 10].index
net = net0[net0["source"].isin(keep_sources)].copy()
net.sort_values(["source", "target"]).to_csv(ANALYSIS / "dorothea_signed_network_v4.csv", index=False)
lg("signed network", len(net), "edges", len(keep_sources), "regulons")

# decoupler ULM: score = t statistic of the slope in expression ~ edge weight.
dc.mt.ulm(ec, net, tmin=10, raw=False, empty=True, bsize=2000, verbose=True)
ulm = ec.obsm["score_ulm"].copy()
ulm.index = ec.obs_names
ulm.to_csv(ANALYSIS / "ec_tf_ulm_v4.csv.gz", compression="gzip")
lg("ULM score matrix", ulm.shape)

obs_cols = [
    "sample",
    "dataset",
    "group",
    "participant_id",
    "inferential_unit",
    "independent_for_inference",
]
cell = pd.concat([ec.obs[obs_cols].reset_index(drop=True), ulm.reset_index(drop=True)], axis=1)
sample_tf = (
    cell.groupby(obs_cols, observed=True, dropna=False)[ulm.columns]
    .mean()
    .reset_index()
)
sample_tf.to_csv(ANALYSIS / "sample_level_tf_ulm_v4.csv", index=False)

# Formal inference: only HOA, ARCO 3A and ARCO 4 donors from the Liao cohort.
liao = sample_tf[
    (sample_tf["dataset"] == "liao_alcohol")
    & sample_tf["group"].isin(LIAO_INFERENCE_GROUPS)
    & sample_tf["independent_for_inference"]
].copy()
if liao["sample"].nunique() != 9:
    raise RuntimeError(f"expected 9 independent Liao inference samples, got {liao['sample'].nunique()}")

pairs = [("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA")]
rows = []
for tf in ulm.columns:
    rec = {"TF": tf}
    arrays = [liao.loc[liao["group"] == g, tf].to_numpy() for g in LIAO_INFERENCE_GROUPS]
    rec["Liao_KW_p"] = exact_kruskal_p(arrays)
    rec["Liao_KW_method"] = "exhaustive permutation, tie-corrected"
    for a, b in pairs:
        xa = liao.loc[liao["group"] == a, tf].to_numpy()
        xb = liao.loc[liao["group"] == b, tf].to_numpy()
        key = f"{a}_vs_{b}"
        summary = exact_mwu_summary(xa, xb)
        rec[f"{key}_effect"] = float(np.median(xa) - np.median(xb))
        rec[f"{key}_hodges_lehmann"] = summary["hodges_lehmann"]
        rec[f"{key}_rank_biserial"] = summary["rank_biserial"]
        rec[f"{key}_p"] = summary["p"]
        rec[f"{key}_minimum_p_no_ties"] = summary["minimum_p_no_ties"]
    # Cross-cohort library-level effect only; no p value is calculated.
    xs = sample_tf.loc[sample_tf["group"] == "SONFH", tf].to_numpy()
    xh = sample_tf.loc[sample_tf["group"] == "HOA", tf].to_numpy()
    rec["SONFH_libraries_vs_HOA_effect_descriptive"] = float(np.median(xs) - np.median(xh))
    rows.append(rec)

stats = pd.DataFrame(rows)
stats["Liao_KW_fdr"] = bh_fdr(stats["Liao_KW_p"])
for a, b in pairs:
    key = f"{a}_vs_{b}"
    stats[f"{key}_fdr"] = bh_fdr(stats[f"{key}_p"])
stats.to_csv(ANALYSIS / "tf_stats_ulm_v4.csv", index=False)

# Pre-specified YAP/TAZ-TEAD signature: score and test it instead of reporting coverage alone.
yap_genes = [
    "CTGF", "CYR61", "ANKRD1", "AXL", "AMOTL2", "BIRC5", "BCL2", "FOSL1",
    "NT5E", "CRIM1", "ITGB2", "SERPINE1", "TGFB2", "DIAPH3", "MYL9", "GADD45A",
]
present = [g for g in yap_genes if g in ec.var_names]
idx = [ec.var_names.get_loc(g) for g in present]
x = ec.X[:, idx]
x = x.toarray() if sp.issparse(x) else np.asarray(x)
x = x.astype(np.float64)
sd = x.std(axis=0, ddof=1)
sd[sd == 0] = 1
yap_cell = ((x - x.mean(axis=0)) / sd).mean(axis=1)
yap = ec.obs[obs_cols].copy()
yap["YAPTAZ_TEAD_signature"] = yap_cell
yap_sample = yap.groupby(obs_cols, observed=True, dropna=False)["YAPTAZ_TEAD_signature"].mean().reset_index()
yap_sample.to_csv(ANALYSIS / "yaptaz_signature_v4.csv", index=False)
yap_liao = yap_sample[
    (yap_sample["dataset"] == "liao_alcohol")
    & yap_sample["group"].isin(LIAO_INFERENCE_GROUPS)
]
yap_arrays = [
    yap_liao.loc[yap_liao["group"] == g, "YAPTAZ_TEAD_signature"].to_numpy()
    for g in LIAO_INFERENCE_GROUPS
]
yap_result = {
    "n_genes_present": len(present),
    "n_genes_defined": len(yap_genes),
    "genes_present": ";".join(present),
    "Liao_KW_p": exact_kruskal_p(yap_arrays),
    "Liao_KW_method": "exhaustive permutation, tie-corrected",
}
for a, b in pairs:
    xa = yap_liao.loc[yap_liao["group"] == a, "YAPTAZ_TEAD_signature"].to_numpy()
    xb = yap_liao.loc[yap_liao["group"] == b, "YAPTAZ_TEAD_signature"].to_numpy()
    key = f"{a}_vs_{b}"
    summary = exact_mwu_summary(xa, xb)
    yap_result[f"{key}_effect"] = float(np.median(xa) - np.median(xb))
    yap_result[f"{key}_hodges_lehmann"] = summary["hodges_lehmann"]
    yap_result[f"{key}_rank_biserial"] = summary["rank_biserial"]
    yap_result[f"{key}_p"] = summary["p"]
    yap_result[f"{key}_minimum_p_no_ties"] = summary["minimum_p_no_ties"]
pd.DataFrame([yap_result]).to_csv(ANALYSIS / "yaptaz_signature_stats_v4.csv", index=False)
lg("YAP/TAZ signature", f"{len(present)}/{len(yap_genes)} genes", "Liao KW", yap_result["Liao_KW_p"])

# ------------------------------ figures ------------------------------
palette = dict(zip(GROUP_ORDER, sns.color_palette("Set2", len(GROUP_ORDER))))

# A: all libraries, clearly labelled as library-level visualization.
tf_cols = list(ulm.columns)
top = sample_tf[tf_cols].var().sort_values(ascending=False).head(24).index
heat = sample_tf.set_index("sample")[list(top)]
heat_group = sample_tf.set_index("sample")["group"]
heat = heat.loc[heat_group.sort_values(key=lambda s: s.map({g: i for i, g in enumerate(GROUP_ORDER)})).index]
heat_z = (heat - heat.mean(axis=0)) / heat.std(axis=0, ddof=1)
fig, ax = plt.subplots(figsize=(9.2, 7.2))
sns.heatmap(heat_z.T, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "row z-score"})
ax.set_xlabel("sequencing libraries grouped by cohort/group")
ax.set_ylabel("")
ax.set_title("Top variable TF activities (signed DoRothEA ULM; visualization only)")
save_figure(fig, FIGS / "fig5a_v4_ulm_heatmap.png")
plt.close(fig)

# B: formal Liao-only donor display.
focus = [
    t for t in ["TEAD1", "TEAD2", "TEAD4", "YAP1", "WWTR1", "RELA", "NFKB1",
                "IRF3", "ATF4", "FOXO3", "STAT3", "HIF1A"] if t in tf_cols
]
fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.2))
for ax, tf in zip(axes.ravel(), focus):
    sns.stripplot(
        data=liao, x="group", y=tf, order=LIAO_INFERENCE_GROUPS, hue="group",
        palette=palette, legend=False, size=6, jitter=0.12, ax=ax,
    )
    for j, group in enumerate(LIAO_INFERENCE_GROUPS):
        values = liao.loc[liao["group"] == group, tf]
        ax.errorbar(j, values.mean(), yerr=values.std(ddof=1), color="black", capsize=3, lw=1.2)
    row = stats.loc[stats["TF"] == tf].iloc[0]
    ax.set_title(f"{tf} (Liao KW FDR={row['Liao_KW_fdr']:.2f})")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=22)
for ax in axes.ravel()[len(focus):]:
    ax.axis("off")
fig.suptitle("TF activity in independent Liao donors (points; mean +/- SD)", y=1.01)
fig.tight_layout()
save_figure(fig, FIGS / "fig5b_v4_ulm_liao.png")
plt.close(fig)

# C: one valid within-cohort contrast; exact MWU p values have a discrete floor at n=3 vs 3.
contrast = "ONFH_4_vs_ONFH_3A"
plot_df = stats.dropna(subset=[f"{contrast}_p"]).copy()
plot_df["mlog10p"] = -np.log10(plot_df[f"{contrast}_p"])
fig, ax = plt.subplots(figsize=(7.6, 5.8))
ax.scatter(plot_df[f"{contrast}_effect"], plot_df["mlog10p"], s=20, c="#9ecae1", edgecolor="white")
for tf in ["TEAD1", "TEAD2", "TEAD4", "RELA", "IRF3", "ATF4", "STAT3", "HIF1A"]:
    row = plot_df[plot_df["TF"] == tf]
    if row.empty:
        continue
    r = row.iloc[0]
    ax.scatter(r[f"{contrast}_effect"], r["mlog10p"], s=48, c="#cb181d", edgecolor="black", lw=0.4)
    ax.annotate(tf, (r[f"{contrast}_effect"], r["mlog10p"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
ax.set_xlabel("median ULM activity difference (ARCO 4 - ARCO 3A)")
ax.set_ylabel("-log10 exact two-sided MWU p")
ax.set_title("Within-cohort TF activity effects (Liao; n=3 vs 3)\nno regulon significant after BH correction")
ax.text(0.99, 0.02, "Discrete exact-test p values; effect sizes are emphasized", transform=ax.transAxes,
        ha="right", color="#555555", fontsize=8, style="italic")
fig.tight_layout()
save_figure(fig, FIGS / "fig5c_v4_ulm_liao_effects.png")
plt.close(fig)

lg("done")
