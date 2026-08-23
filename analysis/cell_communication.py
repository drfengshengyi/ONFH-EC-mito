# -*- coding: utf-8 -*-
"""Figure 4 v4: conservative EC-centred ligand-receptor scoring.

Multi-subunit complexes require every annotated subunit to be present and to
pass the prevalence threshold. Complex expression is summarized by the
geometric mean of subunit means. Formal tests use only independent Liao donors;
GSE290411 libraries are reported descriptively without participant-level p
values.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
import scanpy as sc

from v4_common import (
    ANALYSIS,
    FIGS,
    LIAO_INFERENCE_GROUPS,
    ROOT,
    attach_v4_metadata,
    bh_fdr,
    exact_mwu_summary,
    save_figure,
    setup_plot,
)

setup_plot()
import matplotlib.pyplot as plt
import seaborn as sns

MIN_CELLS = 30
MIN_PREVALENCE = 0.10
LOG = ANALYSIS / "py_fig4_v4.log"


def lg(*parts) -> None:
    line = time.strftime("%H:%M:%S") + " | " + " ".join(map(str, parts))
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def geometric_mean(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return values
    positive = np.all(values > 0, axis=1)
    out = np.zeros(values.shape[0], dtype=float)
    out[positive] = np.exp(np.mean(np.log(values[positive]), axis=1))
    return out


if LOG.exists():
    LOG.unlink()
lg("start")

adata = attach_v4_metadata(sc.read_h5ad(ANALYSIS / "atlas_annotated.h5ad"))
adata.obs["sample"] = adata.obs["sample"].astype(str)
adata.obs["group"] = adata.obs["group"].astype(str)
adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
lg("atlas", adata.shape)

# ------------------------- interaction database -------------------------
interactions = pd.read_csv(ROOT / "data" / "cellchatdb_interactions.csv")
complex_table = pd.read_csv(ROOT / "data" / "cellchatdb_complex_named.csv")
complex_map: dict[str, tuple[str, ...]] = {}
for _, row in complex_table.iterrows():
    components = tuple(
        str(row[col]).strip()
        for col in complex_table.columns
        if col.startswith("subunit_") and pd.notna(row[col]) and str(row[col]).strip()
    )
    if components:
        complex_map[str(row["complex_name"]).strip()] = components

universe = set(adata.var_names)


def resolve(name: str) -> tuple[str, ...]:
    name = str(name).strip()
    return complex_map.get(name, (name,))


pair_map = {}
skipped_missing_subunit = 0
missing_subunit_rows = []
for _, row in interactions.iterrows():
    ligand = resolve(row["ligand"])
    receptor = resolve(row["receptor"])
    # Never silently shrink a complex when one essential subunit is missing.
    if not set(ligand).issubset(universe) or not set(receptor).issubset(universe):
        skipped_missing_subunit += 1
        missing_subunit_rows.append({
            "interaction_name": row.get("interaction_name", ""),
            "ligand": row["ligand"],
            "receptor": row["receptor"],
            "missing_ligand_subunits": ";".join(sorted(set(ligand).difference(universe))),
            "missing_receptor_subunits": ";".join(sorted(set(receptor).difference(universe))),
        })
        continue
    key = (tuple(sorted(ligand)), tuple(sorted(receptor)))
    if key in pair_map:
        continue
    label = row.get("interaction_name_2")
    if not isinstance(label, str) or not label.strip():
        label = row["interaction_name"]
    pair_map[key] = {
        "pair": label,
        "pathway": row["pathway_name"],
        "ligand": ligand,
        "receptor": receptor,
        "ligand_is_complex": len(ligand) > 1,
        "receptor_is_complex": len(receptor) > 1,
    }
pairs = list(pair_map.values())
pd.DataFrame(missing_subunit_rows).to_csv(
    ANALYSIS / "comm_missing_subunit_v4.csv", index=False
)
pd.DataFrame([
    {"category": "retained_complete_unique_interactions", "n": len(pairs)},
    {"category": "excluded_rows_missing_essential_subunit", "n": skipped_missing_subunit},
]).to_csv(ANALYSIS / "comm_complex_filter_audit_v4.csv", index=False)
lg("usable unique pairs", len(pairs), "skipped for missing essential subunit", skipped_missing_subunit)

needed_genes = sorted({g for pair in pairs for g in pair["ligand"] + pair["receptor"]})
gene_to_needed = {g: i for i, g in enumerate(needed_genes)}
ligand_idx = [np.array([gene_to_needed[g] for g in pair["ligand"]], dtype=int) for pair in pairs]
receptor_idx = [np.array([gene_to_needed[g] for g in pair["receptor"]], dtype=int) for pair in pairs]

# ---------------------- per-library cell-type summaries ----------------------
x = adata.X.tocsr()
columns = np.array([adata.var_names.get_loc(g) for g in needed_genes], dtype=int)
samples = adata.obs["sample"].to_numpy()
groups = adata.obs["group"].to_numpy()
cell_types = adata.obs["cell_type"].to_numpy()
types = sorted(np.unique(cell_types))
sample_list = sorted(np.unique(samples))
sample_group = {sample: groups[samples == sample][0] for sample in sample_list}
sample_inferential = {
    sample: bool(adata.obs.loc[adata.obs["sample"] == sample, "independent_for_inference"].iloc[0])
    for sample in sample_list
}


def sample_type_stats(sample: str):
    rows = np.flatnonzero(samples == sample)
    local_types = cell_types[rows]
    mean = np.full((len(types), len(needed_genes)), np.nan, dtype=np.float32)
    prevalence = np.full_like(mean, np.nan)
    n_cells = np.zeros(len(types), dtype=int)
    for type_idx, cell_type in enumerate(types):
        selected = rows[local_types == cell_type]
        n_cells[type_idx] = len(selected)
        if len(selected) == 0:
            continue
        sub = x[selected][:, columns]
        mean[type_idx] = np.asarray(sub.mean(axis=0)).ravel()
        prevalence[type_idx] = np.asarray((sub > 0).mean(axis=0)).ravel()
    return mean, prevalence, n_cells


means, prevalences, cell_counts = {}, {}, {}
for sample in sample_list:
    means[sample], prevalences[sample], cell_counts[sample] = sample_type_stats(sample)
lg("per-library summaries complete", len(sample_list))


def scores_for_sample(sample: str) -> np.ndarray:
    mean = means[sample]
    prevalence = prevalences[sample]
    n_cells = cell_counts[sample]
    out = np.full((len(pairs), len(types), len(types)), np.nan, dtype=np.float32)
    for pair_idx in range(len(pairs)):
        li = ligand_idx[pair_idx]
        ri = receptor_idx[pair_idx]
        ligand_mean = geometric_mean(mean[:, li])
        receptor_mean = geometric_mean(mean[:, ri])
        # Every essential subunit must clear the prevalence threshold.
        ligand_ok = np.all(prevalence[:, li] >= MIN_PREVALENCE, axis=1)
        receptor_ok = np.all(prevalence[:, ri] >= MIN_PREVALENCE, axis=1)
        ligand_ok &= (n_cells >= MIN_CELLS) & np.isfinite(ligand_mean)
        receptor_ok &= (n_cells >= MIN_CELLS) & np.isfinite(receptor_mean)
        valid = np.outer(ligand_ok, receptor_ok)
        score = np.outer(np.nan_to_num(ligand_mean), np.nan_to_num(receptor_mean))
        out[pair_idx] = np.where(valid, score, np.nan)
    return out


scores = {sample: scores_for_sample(sample) for sample in sample_list}
ec_index = types.index("EC")

# Store individual library/donor scores for transparent missingness auditing.
long_rows = []
for sample in sample_list:
    for pair_idx, pair in enumerate(pairs):
        for sender_idx, sender in enumerate(types):
            for receiver_idx, receiver in enumerate(types):
                if sender != "EC" and receiver != "EC":
                    continue
                value = scores[sample][pair_idx, sender_idx, receiver_idx]
                long_rows.append(
                    {
                        "sample": sample,
                        "group": sample_group[sample],
                        "independent_for_inference": sample_inferential[sample],
                        "sender": sender,
                        "receiver": receiver,
                        "pair": pair["pair"],
                        "pathway": pair["pathway"],
                        "ligand_is_complex": pair["ligand_is_complex"],
                        "receptor_is_complex": pair["receptor_is_complex"],
                        "score": value,
                    }
                )
long = pd.DataFrame(long_rows)
long.to_csv(ANALYSIS / "comm_scores_v4_long.csv.gz", index=False, compression="gzip")
lg("long EC-centred score table", long.shape)

# ------------------------ valid within-cohort inference ------------------------
formal_contrasts = [("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA")]
rows = []
keys = ["sender", "receiver", "pair", "pathway", "ligand_is_complex", "receptor_is_complex"]
for identity, data in long.groupby(keys, observed=True, dropna=False):
    rec = dict(zip(keys, identity))
    for a, b in formal_contrasts:
        va = data.loc[
            (data["group"] == a) & data["independent_for_inference"], "score"
        ].dropna().to_numpy()
        vb = data.loc[
            (data["group"] == b) & data["independent_for_inference"], "score"
        ].dropna().to_numpy()
        prefix = f"{a}_vs_{b}"
        rec[f"{prefix}_n_a"] = len(va)
        rec[f"{prefix}_n_b"] = len(vb)
        if len(va) >= 2 and len(vb) >= 2:
            rec[f"{prefix}_median_a"] = float(np.median(va))
            rec[f"{prefix}_median_b"] = float(np.median(vb))
            rec[f"{prefix}_delta"] = float(np.median(va) - np.median(vb))
            summary = exact_mwu_summary(va, vb)
            rec[f"{prefix}_hodges_lehmann"] = summary["hodges_lehmann"]
            rec[f"{prefix}_rank_biserial"] = summary["rank_biserial"]
            rec[f"{prefix}_p"] = summary["p"]
            rec[f"{prefix}_minimum_p_no_ties"] = summary["minimum_p_no_ties"]

    # SONFH libraries versus HOA: effect only, explicitly non-inferential.
    sonfh = data.loc[data["group"] == "SONFH", "score"].dropna().to_numpy()
    hoa = data.loc[
        (data["group"] == "HOA") & data["independent_for_inference"], "score"
    ].dropna().to_numpy()
    rec["SONFH_libraries_n"] = len(sonfh)
    rec["HOA_donors_n"] = len(hoa)
    if len(sonfh) and len(hoa):
        rec["SONFH_libraries_vs_HOA_delta_descriptive"] = float(np.median(sonfh) - np.median(hoa))
    rows.append(rec)

result = pd.DataFrame(rows)
for a, b in formal_contrasts:
    prefix = f"{a}_vs_{b}"
    result[f"{prefix}_fdr"] = bh_fdr(result[f"{prefix}_p"])
result.to_csv(ANALYSIS / "comm_donor_level_v4.csv", index=False)

# ------------------------------- figures -------------------------------
contrast = "ONFH_4_vs_ONFH_3A"
delta_col = f"{contrast}_delta"
top = result.dropna(subset=[delta_col]).copy()
top["abs_delta"] = top[delta_col].abs()
top = top.sort_values("abs_delta", ascending=False).head(15)
top.to_csv(ANALYSIS / "comm_top_ONFH_4_vs_ONFH_3A_v4.csv", index=False)
top_labels = (top["sender"] + "->" + top["receiver"] + ": " + top["pair"]).tolist()
top_key = top[["sender", "receiver", "pair"]].drop_duplicates()
plot_long = long[
    long["group"].isin(["ONFH_3A", "ONFH_4"])
    & long["independent_for_inference"]
].merge(top_key, on=["sender", "receiver", "pair"], how="inner")
plot_long["label"] = plot_long["sender"] + "->" + plot_long["receiver"] + ": " + plot_long["pair"]

fig, ax = plt.subplots(figsize=(9.5, 7.5))
sns.stripplot(
    data=plot_long, x="score", y="label", order=top_labels, hue="group",
    hue_order=["ONFH_3A", "ONFH_4"], palette=["#3690c0", "#cb181d"],
    dodge=True, size=5.5, jitter=0.13, ax=ax,
)
ax.set_xlabel("ligand-receptor score in individual Liao donors")
ax.set_ylabel("")
ax.set_title("Largest within-cohort EC-centred effects (ARCO 4 vs 3A)\nall donor points shown; FDR values in source table")
ax.legend(title="", loc="lower right")
fig.tight_layout()
save_figure(fig, FIGS / "fig4c_v4_liao_lr_points.png")
plt.close(fig)

sender_ec = result[result["sender"] == "EC"].dropna(subset=[delta_col]).copy()
pathway = sender_ec.groupby(["pathway", "receiver"], observed=True)[delta_col].median().unstack()
pathway = pathway.fillna(0)
pathway = pathway.loc[pathway.abs().sum(axis=1).sort_values(ascending=False).head(18).index]
fig, ax = plt.subplots(figsize=(8.4, 6.5))
sns.heatmap(pathway, cmap="RdBu_r", center=0, ax=ax, cbar_kws={"label": "median score difference"})
ax.set_xlabel("receiver cell type")
ax.set_ylabel("")
ax.set_title("EC-sender pathway effects: ARCO 4 minus ARCO 3A (Liao cohort)")
fig.tight_layout()
save_figure(fig, FIGS / "fig4d_v4_liao_pathway.png")
plt.close(fig)

for a, b in formal_contrasts:
    prefix = f"{a}_vs_{b}"
    n_sig = int((result[f"{prefix}_fdr"] < 0.05).sum())
    lg(prefix, "BH-FDR<0.05", n_sig)
lg("done")
