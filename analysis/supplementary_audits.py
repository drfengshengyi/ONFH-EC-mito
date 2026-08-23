# -*- coding: utf-8 -*-
"""Reviewer-requested small-sample and cross-cohort audit tables.

This script does not create additional biological replicates. It exposes the
attainable resolution of the prespecified tests, adds participant-level effect
sizes with explicitly exploratory bootstrap intervals, and summarizes EC QC
metrics that may differ across source cohorts.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

from v4_common import ANALYSIS, ROOT, bh_fdr, exact_mwu_summary


OUT = ROOT / "results" / "supplementary_tables"
OUT.mkdir(parents=True, exist_ok=True)
CONTRASTS = [("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA")]


def exhaustive_hl_bootstrap_ci(a, b, confidence: float = 0.95) -> tuple[float, float]:
    """Percentile interval from every within-group participant resample.

    At n=3 per group this is a finite 27 x 27 enumeration. The interval is
    reported as descriptive because conventional asymptotics are unreliable
    and an exact two-sided 95% shift interval is unbounded at this resolution.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    estimates = []
    for ia in product(range(len(a)), repeat=len(a)):
        aa = a[list(ia)]
        for ib in product(range(len(b)), repeat=len(b)):
            bb = b[list(ib)]
            estimates.append(float(np.median(aa[:, None] - bb[None, :])))
    alpha = (1 - confidence) / 2
    return tuple(np.quantile(estimates, [alpha, 1 - alpha]).astype(float))


def effect_rows(table: pd.DataFrame, family: str, endpoints: list[str], unit: str) -> list[dict]:
    rows = []
    for endpoint in endpoints:
        for group_a, group_b in CONTRASTS:
            a = table.loc[table["group"] == group_a, endpoint].dropna().to_numpy(dtype=float)
            b = table.loc[table["group"] == group_b, endpoint].dropna().to_numpy(dtype=float)
            summary = exact_mwu_summary(a, b)
            ci_low, ci_high = exhaustive_hl_bootstrap_ci(a, b)
            rows.append(
                {
                    "family": family,
                    "endpoint": endpoint,
                    "unit": unit,
                    "contrast": f"{group_a}_vs_{group_b}",
                    "n_a": len(a),
                    "n_b": len(b),
                    "median_a": np.median(a) if len(a) else np.nan,
                    "median_b": np.median(b) if len(b) else np.nan,
                    "median_difference_a_minus_b": (
                        float(np.median(a) - np.median(b)) if len(a) and len(b) else np.nan
                    ),
                    "hodges_lehmann_shift_a_minus_b": summary["hodges_lehmann"],
                    "rank_biserial_a_gt_b": summary["rank_biserial"],
                    "hl_participant_bootstrap_ci95_low": ci_low,
                    "hl_participant_bootstrap_ci95_high": ci_high,
                    "exact_permutation_mwu_p": summary["p"],
                    "minimum_two_sided_p_no_ties": summary["minimum_p_no_ties"],
                    "interval_note": (
                        "Exhaustive within-group participant-resampling percentile interval; "
                        "descriptive at n=3/group. The exact two-sided 95% shift interval is "
                        "unbounded because the smallest no-tie p value is 0.10."
                    ),
                }
            )
    return rows


def build_small_sample_table() -> pd.DataFrame:
    rows: list[dict] = []

    composition = pd.read_csv(ANALYSIS / "composition_by_library_v4.csv")
    composition = composition.loc[composition["dataset"] == "liao_alcohol"]
    rows.extend(effect_rows(composition, "major-cell composition", ["EC"], "percentage points"))

    subtype = pd.read_csv(ANALYSIS / "ec_subtype_composition_v4.csv")
    subtype = subtype.loc[subtype["dataset"] == "liao_alcohol"]
    subtype_endpoints = [
        "lymphatic", "typeH_EMCN_KDR", "typeR_bone_remodel", "venous_ACKR1"
    ]
    rows.extend(effect_rows(subtype, "EC-state composition", subtype_endpoints, "percentage points"))

    modules = pd.read_csv(ANALYSIS / "module_scores_by_library_v4.csv")
    modules = modules.loc[modules["dataset"] == "liao_alcohol"]
    module_endpoints = [
        "Mito_fission", "Mito_fusion", "Mitophagy_core",
        "cGAS_STING", "EC_inflammation", "YAP_mTOR",
    ]
    rows.extend(effect_rows(modules, "EC module score", module_endpoints, "mean module-score units"))

    tf = pd.read_csv(ANALYSIS / "sample_level_tf_ulm_v4.csv")
    tf = tf.loc[
        (tf["dataset"] == "liao_alcohol")
        & tf["group"].isin(["HOA", "ONFH_3A", "ONFH_4"])
        & tf["independent_for_inference"].astype(bool)
    ]
    tf_endpoints = [name for name in ["RELA", "NFKB1", "ATF4", "FOXO3", "STAT3", "HIF1A"] if name in tf]
    rows.extend(effect_rows(tf, "signed TF activity", tf_endpoints, "ULM activity units"))

    top = pd.read_csv(ANALYSIS / "comm_top_ONFH_4_vs_ONFH_3A_v4.csv")
    long = pd.read_csv(ANALYSIS / "comm_scores_v4_long.csv.gz")
    keys = top[["sender", "receiver", "pair"]].drop_duplicates()
    lr = long.loc[
        long["group"].isin(["HOA", "ONFH_3A", "ONFH_4"])
        & long["independent_for_inference"].astype(bool)
    ].merge(keys, on=["sender", "receiver", "pair"], how="inner")
    lr["endpoint"] = lr["sender"] + "->" + lr["receiver"] + ": " + lr["pair"]
    for endpoint, data in lr.groupby("endpoint", observed=True):
        wide = data[["group", "score"]].rename(columns={"score": endpoint})
        rows.extend(effect_rows(wide, "selected EC-centered ligand-receptor score", [endpoint], "score units"))

    result = pd.DataFrame(rows)
    result["bh_fdr_within_family_and_contrast"] = np.nan
    for _, indices in result.groupby(["family", "contrast"], observed=True).groups.items():
        result.loc[indices, "bh_fdr_within_family_and_contrast"] = bh_fdr(
            result.loc[indices, "exact_permutation_mwu_p"]
        )
    result.to_csv(ANALYSIS / "small_sample_effects_v7.csv", index=False)
    result.to_csv(OUT / "TableS4_small_sample_effects_v7.csv", index=False)
    return result


def build_cross_cohort_qc() -> tuple[pd.DataFrame, pd.DataFrame]:
    atlas = sc.read_h5ad(ANALYSIS / "atlas_annotated.h5ad", backed="r")
    obs = atlas.obs.loc[atlas.obs["cell_type"].astype(str) == "EC"].copy()
    qc = (
        obs.groupby(["sample", "dataset", "group"], observed=True)
        .agg(
            retained_ecs=("cell_type", "size"),
            median_genes_per_ec=("nFeature", "median"),
            median_umis_per_ec=("nCount", "median"),
            median_mitochondrial_percent_per_ec=("percent.mt", "median"),
        )
        .reset_index()
    )
    qc["sampling_unit"] = np.where(qc["group"] == "SONFH", "unresolved sequencing library", "participant")
    qc.to_csv(ANALYSIS / "ec_qc_by_sampling_unit_v7.csv", index=False)
    qc.to_csv(OUT / "TableS3_ec_qc_by_sampling_unit_v7.csv", index=False)

    rows = []
    for dataset, data in qc.groupby("dataset", observed=True):
        rec = {
            "dataset": dataset,
            "n_public_sampling_units": len(data),
            "sampling_unit": "; ".join(sorted(data["sampling_unit"].unique())),
        }
        for metric in [
            "retained_ecs", "median_genes_per_ec", "median_umis_per_ec",
            "median_mitochondrial_percent_per_ec",
        ]:
            rec[f"{metric}_median_across_units"] = data[metric].median()
            rec[f"{metric}_min_across_units"] = data[metric].min()
            rec[f"{metric}_max_across_units"] = data[metric].max()
        rows.append(rec)
    cohort = pd.DataFrame(rows)
    cohort.to_csv(ANALYSIS / "ec_qc_by_cohort_v7.csv", index=False)
    cohort.to_csv(OUT / "TableS3_ec_qc_by_cohort_v7.csv", index=False)
    return qc, cohort


def build_direction_and_mapping_tables() -> None:
    effects = pd.read_csv(ANALYSIS / "fig3_key_gene_effects_v4.csv", index_col=0)
    effects.columns = ["liao_onfh3a_vs_hoa_log2fc", "sonfh_libraries_vs_hoa_log2fc_descriptive"]
    effects["direction_concordant"] = (
        np.sign(effects["liao_onfh3a_vs_hoa_log2fc"])
        == np.sign(effects["sonfh_libraries_vs_hoa_log2fc_descriptive"])
    )
    effects["interpretation"] = np.where(
        effects["direction_concordant"],
        "directionally concordant; magnitude remains vulnerable to cross-cohort technical confounding",
        "directionally discordant; no cross-cohort replication",
    )
    effects.to_csv(ANALYSIS / "cross_cohort_direction_concordance_v7.csv")
    effects.to_csv(OUT / "TableS6_cross_cohort_direction_concordance_v7.csv")

    mapping = pd.read_csv(ANALYSIS / "ec_annotation_confidence_v4.csv")
    mapping.to_csv(OUT / "TableS5_cluster_to_state_mapping_v7.csv", index=False)


def main() -> None:
    effects = build_small_sample_table()
    _, cohort = build_cross_cohort_qc()
    build_direction_and_mapping_tables()
    print(f"Wrote {len(effects)} small-sample effect rows")
    print(cohort.to_string(index=False))
    print(f"Supplementary source tables: {OUT}")


if __name__ == "__main__":
    main()
