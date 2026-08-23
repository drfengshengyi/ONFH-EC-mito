# -*- coding: utf-8 -*-
"""Composition summaries with correct sampling-unit labels and valid contrasts."""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
import scanpy as sc

from v4_common import (
    ANALYSIS,
    FIGS,
    GROUP_ORDER,
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

LOG = ANALYSIS / "py_comp_v4.log"
LIAO_GROUPS = ["FNF", "HOA", "ONFH_3A", "ONFH_4"]
PAIRWISE = [("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA")]


def lg(*parts) -> None:
    line = time.strftime("%H:%M:%S") + " | " + " ".join(map(str, parts))
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def add_stats(
    table: pd.DataFrame,
    feature_columns: list[str],
    feature_name: str,
    *,
    kw_groups: list[str] | None = None,
) -> pd.DataFrame:
    kw_groups = LIAO_GROUPS if kw_groups is None else kw_groups
    rows = []
    for feature in feature_columns:
        rec = {feature_name: feature}
        liao = table[(table["dataset"] == "liao_alcohol") & table["group"].isin(LIAO_GROUPS)]
        arrays = [liao.loc[liao["group"] == group, feature].dropna().to_numpy() for group in kw_groups]
        rec["Liao_KW_p"] = exact_kruskal_p(arrays)
        rec["Liao_KW_method"] = "exhaustive permutation, tie-corrected"
        for a, b in PAIRWISE:
            va = liao.loc[liao["group"] == a, feature].dropna().to_numpy()
            vb = liao.loc[liao["group"] == b, feature].dropna().to_numpy()
            key = f"{a}_vs_{b}"
            summary = exact_mwu_summary(va, vb)
            rec[f"{key}_median_difference"] = float(np.median(va) - np.median(vb))
            rec[f"{key}_effect"] = rec[f"{key}_median_difference"]
            rec[f"{key}_hodges_lehmann"] = summary["hodges_lehmann"]
            rec[f"{key}_rank_biserial"] = summary["rank_biserial"]
            rec[f"{key}_p"] = summary["p"]
            rec[f"{key}_minimum_p_no_ties"] = summary["minimum_p_no_ties"]
        sonfh = table.loc[table["group"] == "SONFH", feature].dropna().to_numpy()
        hoa = table.loc[table["group"] == "HOA", feature].dropna().to_numpy()
        rec["SONFH_libraries_vs_HOA_effect_descriptive"] = float(np.median(sonfh) - np.median(hoa))
        rows.append(rec)
    result = pd.DataFrame(rows)
    result["Liao_KW_fdr"] = bh_fdr(result["Liao_KW_p"])
    for a, b in PAIRWISE:
        key = f"{a}_vs_{b}"
        result[f"{key}_fdr"] = bh_fdr(result[f"{key}_p"])
    return result


def summarize(table: pd.DataFrame, feature_columns: list[str], feature_name: str) -> pd.DataFrame:
    rows = []
    for feature in feature_columns:
        for group in GROUP_ORDER:
            values = table.loc[table["group"] == group, feature].dropna()
            rows.append(
                {
                    feature_name: feature,
                    "group": group,
                    "unit": "sequencing library" if group == "SONFH" else "participant",
                    "n_units": len(values),
                    "mean": values.mean(),
                    "sd": values.std(ddof=1),
                    "median": values.median(),
                    "min": values.min(),
                    "max": values.max(),
                }
            )
    return pd.DataFrame(rows)


if LOG.exists():
    LOG.unlink()
atlas = attach_v4_metadata(sc.read_h5ad(ANALYSIS / "atlas_annotated.h5ad"))
ec = attach_v4_metadata(sc.read_h5ad(ANALYSIS / "ec_final.h5ad"))
lg("atlas", atlas.shape, "EC", ec.shape)

metadata = (
    atlas.obs[["sample", "dataset", "group", "inferential_unit", "independent_for_inference"]]
    .drop_duplicates("sample")
    .set_index("sample")
)

# All-cell composition.
composition = pd.crosstab(atlas.obs["sample"], atlas.obs["cell_type"], normalize="index") * 100
composition = metadata.join(composition, how="right").reset_index()
composition.to_csv(ANALYSIS / "composition_by_library_v4.csv", index=False)
cell_types = [c for c in composition.columns if c not in metadata.columns and c != "sample"]
summarize(composition, cell_types, "cell_type").to_csv(
    ANALYSIS / "composition_group_summary_v4.csv", index=False
)
add_stats(composition, cell_types, "cell_type").to_csv(
    ANALYSIS / "composition_stats_v4.csv", index=False
)

# EC-subtype composition.
ec_metadata = (
    ec.obs[["sample", "dataset", "group", "inferential_unit", "independent_for_inference"]]
    .drop_duplicates("sample")
    .set_index("sample")
)
subtype = pd.crosstab(ec.obs["sample"], ec.obs["EC_subtype"], normalize="index") * 100
subtype = ec_metadata.join(subtype, how="right").reset_index()
subtype.to_csv(ANALYSIS / "ec_subtype_composition_v4.csv", index=False)
subtypes = [c for c in subtype.columns if c not in ec_metadata.columns and c != "sample"]
summarize(subtype, subtypes, "subtype").to_csv(
    ANALYSIS / "ec_subtype_group_summary_v4.csv", index=False
)
subtype_stats = add_stats(subtype, subtypes, "subtype")
subtype_stats.to_csv(ANALYSIS / "ec_subtype_stats_v4.csv", index=False)

palette = dict(zip(GROUP_ORDER, sns.color_palette("Set2", len(GROUP_ORDER))))
display_order = GROUP_ORDER
display_labels = ["Healthy", "HOA", "FNF", "ARCO 3A", "ARCO 4", "SONFH\n(libraries)"]

# Figure 1e: points plus actual mean +/- SD.
fig, ax = plt.subplots(figsize=(6.1, 4.5))
sns.stripplot(
    data=composition, x="group", y="EC", order=display_order, hue="group",
    palette=palette, legend=False, size=6, jitter=0.12, ax=ax,
)
for i, group in enumerate(display_order):
    values = composition.loc[composition["group"] == group, "EC"].dropna()
    ax.errorbar(i, values.mean(), yerr=values.std(ddof=1), color="black", capsize=4, lw=1.3)
ax.set_xticks(range(len(display_order)), display_labels, rotation=20, ha="right")
ax.set_xlabel("")
ax.set_ylabel("EC fraction (% of all retained cells)")
ax.set_title("EC composition by sampling unit (points; mean +/- SD)")
ax.text(
    0.99, 0.98, "SONFH points are sequencing libraries; no participant-level test",
    transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#555555",
)
fig.tight_layout()
save_figure(fig, FIGS / "fig1e_v4_ec_fraction.png")
plt.close(fig)

# Figure 1f: retained-cell counts expose unequal sampling depth across libraries.
sample_counts = atlas.obs["sample"].astype(str).value_counts().rename("retained_cells").to_frame()
sample_counts = metadata[["dataset", "group", "inferential_unit"]].join(sample_counts, how="inner")
sample_counts.to_csv(ANALYSIS / "retained_cell_counts_v4.csv")
sample_counts = sample_counts.sort_values(
    "group", key=lambda s: s.map({g: i for i, g in enumerate(GROUP_ORDER)})
)
fig, ax = plt.subplots(figsize=(7.0, 4.5))
ax.bar(
    sample_counts.index,
    sample_counts["retained_cells"],
    color=[palette[g] for g in sample_counts["group"]],
    edgecolor="white",
    linewidth=0.4,
)
ax.set_ylabel("Retained cells")
ax.set_xlabel("Sequencing library")
ax.set_title("Unequal retained-cell yield across sequencing libraries")
ax.tick_params(axis="x", rotation=55, labelsize=7)
ax.text(
    0.99, 0.98, "Counts describe sequencing libraries, not biological replication",
    transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color="#555555",
)
fig.tight_layout()
save_figure(fig, FIGS / "fig1f_v4_retained_cells.png")
plt.close(fig)

# Figure 2c: one readable facet per retained subtype.
fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.8), sharex=True)
for ax, subtype_name in zip(axes.ravel(), subtypes):
    sns.stripplot(
        data=subtype, x="group", y=subtype_name, order=display_order, hue="group",
        palette=palette, legend=False, size=5.5, jitter=0.12, ax=ax,
    )
    for i, group in enumerate(display_order):
        values = subtype.loc[subtype["group"] == group, subtype_name].dropna()
        ax.errorbar(i, values.mean(), yerr=values.std(ddof=1), color="black", capsize=3, lw=1.1)
    ax.set_xticks(range(len(display_order)), display_labels, rotation=25, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("% of ECs")
    label = subtype_name.replace("typeH_EMCN_KDR", "type H / EMCN-KDR")
    label = label.replace("typeR_bone_remodel", "type R / bone-remodelling")
    label = label.replace("venous_ACKR1", "venous / ACKR1")
    ax.set_title(label)
for ax in axes.ravel()[len(subtypes):]:
    ax.axis("off")
fig.suptitle("EC subtype composition by sampling unit (points; mean +/- SD)", y=1.01)
fig.tight_layout()
save_figure(fig, FIGS / "fig2c_v4_subtype_composition.png")
plt.close(fig)

# Figure 2d: sampling-unit module means, with inference restricted to Liao HOA/ONFH.
module_columns = [
    "Mito_fission", "Mito_fusion", "Mitophagy_core",
    "cGAS_STING", "EC_inflammation", "YAP_mTOR",
]
module_scores = ec.obs.groupby("sample", observed=True)[module_columns].mean()
module_scores = ec_metadata.join(module_scores, how="right").reset_index()
module_scores.to_csv(ANALYSIS / "module_scores_by_library_v4.csv", index=False)
module_stats = add_stats(
    module_scores,
    module_columns,
    "module",
    kw_groups=["HOA", "ONFH_3A", "ONFH_4"],
)
module_stats.to_csv(ANALYSIS / "module_scores_liao_stats_v4.csv", index=False)

fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.8))
for ax, module in zip(axes.ravel(), module_columns):
    sns.stripplot(
        data=module_scores, x="group", y=module, order=display_order, hue="group",
        palette=palette, legend=False, size=5, jitter=0.12, ax=ax,
    )
    for i, group in enumerate(display_order):
        values = module_scores.loc[module_scores["group"] == group, module].dropna()
        ax.errorbar(i, values.mean(), yerr=values.std(ddof=1), color="black", capsize=3, lw=1.0)
    fdr = module_stats.loc[module_stats["module"] == module, "Liao_KW_fdr"].iloc[0]
    ax.set_title(f"{module} (Liao HOA/ONFH FDR={fdr:.2g})")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30, labelsize=7)
fig.suptitle("Endothelial module scores by sampling unit (points; mean +/- SD)")
fig.tight_layout()
save_figure(fig, FIGS / "fig2d_v4_module_scores.png")
plt.close(fig)

type_h = subtype_stats[subtype_stats["subtype"] == "typeH_EMCN_KDR"]
if len(type_h):
    lg("type H Liao KW", float(type_h.iloc[0]["Liao_KW_p"]), "FDR", float(type_h.iloc[0]["Liao_KW_fdr"]))
lg("done")
