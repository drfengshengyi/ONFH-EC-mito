# -*- coding: utf-8 -*-
"""Audit endothelial subtype labels and quantify an ambiguity-aware sensitivity.

This does not overwrite the published panel-based labels.  It records the top
two marker-panel scores for each EC cluster and repeats the type-H composition
summary after clusters with a small top-two margin are labelled Ambiguous.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
from v4_common import (
    ANALYSIS,
    LIAO_INFERENCE_GROUPS,
    exact_kruskal_p,
    exact_mwu_summary,
    read_sample_metadata,
)

MARGIN_THRESHOLD = 0.15

scores = pd.read_csv(ANALYSIS / "ec_subtype_scores.csv", index_col=0)
scores.index = scores.index.astype(str)
ranked = np.argsort(-scores.to_numpy(), axis=1)
labels = scores.columns.to_numpy()
audit = pd.DataFrame({
    "cluster": scores.index,
    "top_panel": labels[ranked[:, 0]],
    "second_panel": labels[ranked[:, 1]],
    "top_score": scores.to_numpy()[np.arange(len(scores)), ranked[:, 0]],
    "second_score": scores.to_numpy()[np.arange(len(scores)), ranked[:, 1]],
})
audit["top_two_margin"] = audit["top_score"] - audit["second_score"]
audit["ambiguous_margin_lt_0.15"] = audit["top_two_margin"] < MARGIN_THRESHOLD

adata = sc.read_h5ad(ANALYSIS / "ec_final.h5ad", backed="r")
obs = adata.obs[["sample", "leiden", "EC_subtype"]].copy()
obs["leiden"] = obs["leiden"].astype(str)
assigned = (
    obs.groupby("leiden", observed=True)["EC_subtype"]
    .agg(lambda x: x.astype(str).value_counts().index[0])
    .rename("assigned_subtype")
)
audit = audit.merge(assigned, left_on="cluster", right_index=True, how="left")
audit["assigned_matches_top_panel"] = audit["assigned_subtype"] == audit["top_panel"]
audit["strict_subtype"] = np.where(
    audit["ambiguous_margin_lt_0.15"], "Ambiguous", audit["top_panel"]
)
audit.to_csv(ANALYSIS / "ec_annotation_confidence_v4.csv", index=False)

strict_map = audit.set_index("cluster")["strict_subtype"]
obs["strict_subtype"] = obs["leiden"].map(strict_map).fillna("Ambiguous")
counts = obs.groupby(["sample", "strict_subtype"], observed=True).size().unstack(fill_value=0)
composition = counts.div(counts.sum(axis=1), axis=0) * 100
metadata = read_sample_metadata().set_index("sample")
composition = metadata[["dataset", "group", "independent_for_inference"]].join(
    composition, how="inner"
).reset_index()
composition.to_csv(ANALYSIS / "ec_subtype_composition_strict_v4.csv", index=False)

type_h = "typeH_EMCN_KDR"
liao = composition[composition["group"].isin(LIAO_INFERENCE_GROUPS)].copy()
groups = [liao.loc[liao["group"] == g, type_h].dropna() for g in LIAO_INFERENCE_GROUPS]
kw_p = exact_kruskal_p(groups)
rows = [{"test": "Liao_KW", "effect_percentage_points": np.nan, "p": kw_p}]
for a, b in [("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA")]:
    xa = liao.loc[liao["group"] == a, type_h].dropna()
    xb = liao.loc[liao["group"] == b, type_h].dropna()
    summary = exact_mwu_summary(xa, xb)
    rows.append({
        "test": f"{a}_vs_{b}",
        "effect_percentage_points": float(xa.median() - xb.median()),
        "hodges_lehmann": summary["hodges_lehmann"],
        "rank_biserial": summary["rank_biserial"],
        "minimum_p_no_ties": summary["minimum_p_no_ties"],
        "p": summary["p"],
    })
pd.DataFrame(rows).to_csv(ANALYSIS / "ec_subtype_strict_stats_v4.csv", index=False)

print(
    f"Audited {len(audit)} EC clusters; "
    f"{audit['ambiguous_margin_lt_0.15'].sum()} had top-two margin < {MARGIN_THRESHOLD}. "
    f"Strict type-H Liao KW p={kw_p:.4g}."
)
