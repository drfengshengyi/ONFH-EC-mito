"""Postprocess donor-separated official-R scTenifoldKnk results.

The script applies the same prespecified, expression/prevalence-matched pathway
null used in the earlier Python-port analysis. It processes both the
manuscript-matched and official-default R profiles, writes profile-specific
pathway tables, and creates cross-profile gene/pathway sensitivity summaries.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import combine_pvalues, spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
VKO_INPUTS = ROOT / "virtual_knockout"
ANALYSIS = ROOT / "analysis"
TARGET = "SQSTM1"
DONORS = ("hoa2", "hoa3")
SEED = 20260820
N_PERMUTATIONS = 20_000
PROFILES = {
    "manuscript": RESULTS / "official_r_vko_manuscript",
    "official_default": RESULTS / "official_r_vko_official_default",
}


def read_profile(directory: Path) -> dict[str, pd.DataFrame]:
    results: dict[str, pd.DataFrame] = {}
    for donor in DONORS:
        path = directory / f"vko_sqstm1_{donor}_official_r.csv"
        frame = pd.read_csv(path)
        required = {"Gene", "rank", "p_value", "adjusted_p_value"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        results[donor] = frame
    return results


def matched_pathway_enrichment(
    results: dict[str, pd.DataFrame],
    features: pd.DataFrame,
    gene_sets: dict[str, list[str]],
) -> pd.DataFrame:
    universe = sorted(
        (set(results["hoa2"]["Gene"]) & set(results["hoa3"]["Gene"]))
        - {TARGET}
    )
    feature_meta = features.set_index("gene").reindex(universe).copy()
    matching_columns = [
        "prevalence_hoa2",
        "prevalence_hoa3",
        "mean_log1p_cpm_hoa2",
        "mean_log1p_cpm_hoa3",
    ]
    if feature_meta[matching_columns].isna().any().any():
        raise ValueError("Feature metadata are incomplete for the R-result universe")

    matching_score = pd.concat(
        [
            feature_meta[column].rank(method="average", pct=True)
            for column in matching_columns
        ],
        axis=1,
    ).mean(axis=1)
    feature_meta["matching_stratum"] = pd.qcut(
        matching_score.rank(method="first"), q=5, labels=False
    )

    rows: list[dict[str, object]] = []
    for pathway_idx, (pathway, members) in enumerate(gene_sets.items()):
        genes = sorted(set(members) & set(universe))
        if len(genes) < 5:
            continue
        row: dict[str, object] = {
            "pathway": pathway,
            "n_genes": len(genes),
            "genes": ";".join(genes),
            "null_model": (
                "20,000 random sets matched by quintiles of control detection "
                "prevalence and mean log1p CPM"
            ),
        }
        stratum_counts = (
            feature_meta.loc[genes, "matching_stratum"].value_counts().to_dict()
        )
        donor_p: list[float] = []
        for donor_idx, donor in enumerate(DONORS):
            ranks = results[donor].set_index("Gene")["rank"].reindex(universe)
            observed = float(ranks.reindex(genes).mean())
            rank_values = ranks.to_numpy()
            strata = feature_meta["matching_stratum"].to_numpy()
            rng = np.random.default_rng(SEED + pathway_idx * 10 + donor_idx)
            permuted_sums = np.zeros(N_PERMUTATIONS, dtype=float)
            for stratum, count in stratum_counts.items():
                pool = rank_values[strata == stratum]
                random_matrix = rng.random((N_PERMUTATIONS, len(pool)))
                sampled = np.argpartition(random_matrix, count - 1, axis=1)[
                    :, :count
                ]
                permuted_sums += pool[sampled].sum(axis=1)
            permuted = permuted_sums / len(genes)
            empirical_p = float(
                (1 + np.count_nonzero(permuted <= observed))
                / (N_PERMUTATIONS + 1)
            )
            row[f"mean_rank_{donor}"] = observed
            row[f"empirical_p_{donor}"] = empirical_p
            donor_p.append(empirical_p)
        row["replicated_nominal_0_05"] = bool(max(donor_p) < 0.05)
        row["fisher_p"] = float(combine_pvalues(donor_p, method="fisher")[1])
        rows.append(row)

    output = pd.DataFrame(rows)
    for donor in DONORS:
        output[f"empirical_fdr_{donor}"] = multipletests(
            output[f"empirical_p_{donor}"], method="fdr_bh"
        )[1]
    output["replicated_fdr_0_05"] = (
        output[[f"empirical_fdr_{donor}" for donor in DONORS]] < 0.05
    ).all(axis=1)
    output["fisher_fdr"] = multipletests(
        output["fisher_p"], method="fdr_bh"
    )[1]
    return output.sort_values(
        ["fisher_fdr", "fisher_p", "pathway"]
    ).reset_index(drop=True)


def gene_sensitivity(
    profile_results: dict[str, dict[str, pd.DataFrame]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged: pd.DataFrame | None = None
    correlations: list[dict[str, object]] = []
    for profile, results in profile_results.items():
        for donor, frame in results.items():
            keep = frame[
                ["Gene", "rank", "p_value", "adjusted_p_value"]
            ].rename(
                columns={
                    "rank": f"rank_{profile}_{donor}",
                    "p_value": f"p_value_{profile}_{donor}",
                    "adjusted_p_value": f"adjusted_p_value_{profile}_{donor}",
                }
            )
            merged = keep if merged is None else merged.merge(
                keep, on="Gene", how="inner"
            )
    assert merged is not None

    for donor in DONORS:
        rho, p_value = spearmanr(
            merged[f"rank_manuscript_{donor}"],
            merged[f"rank_official_default_{donor}"],
        )
        correlations.append(
            {
                "comparison": "official-R manuscript vs official-default",
                "donor": donor,
                "spearman_rho": rho,
                "p_value": p_value,
            }
        )
    rank_columns = [column for column in merged if column.startswith("rank_")]
    merged["mean_rank_across_profiles_and_donors"] = merged[rank_columns].mean(
        axis=1
    )
    merged["max_rank_across_profiles_and_donors"] = merged[rank_columns].max(
        axis=1
    )
    fdr_columns = [
        column for column in merged if column.startswith("adjusted_p_value_")
    ]
    merged["fdr_significant_all_profiles_and_donors"] = (
        merged[fdr_columns] < 0.05
    ).all(axis=1)
    merged["top10_all_profiles_and_donors"] = (merged[rank_columns] <= 10).all(
        axis=1
    )
    merged = merged.sort_values(
        ["mean_rank_across_profiles_and_donors", "Gene"]
    ).reset_index(drop=True)
    return merged, pd.DataFrame(correlations)


def pathway_sensitivity(
    pathway_results: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    merged: pd.DataFrame | None = None
    for profile, frame in pathway_results.items():
        keep_columns = [
            "pathway",
            "n_genes",
            "genes",
            *[f"mean_rank_{donor}" for donor in DONORS],
            *[f"empirical_p_{donor}" for donor in DONORS],
            *[f"empirical_fdr_{donor}" for donor in DONORS],
            "replicated_nominal_0_05",
            "replicated_fdr_0_05",
            "fisher_p",
            "fisher_fdr",
        ]
        keep = frame[keep_columns].copy()
        keep = keep.rename(
            columns={
                column: f"{column}_{profile}"
                for column in keep.columns
                if column not in {"pathway", "n_genes", "genes"}
            }
        )
        merged = keep if merged is None else merged.merge(
            keep, on=["pathway", "n_genes", "genes"], how="inner"
        )
    assert merged is not None
    fdr_columns = [
        f"empirical_fdr_{donor}_{profile}"
        for profile in PROFILES
        for donor in DONORS
    ]
    merged["fdr_significant_all_profiles_and_donors"] = (
        merged[fdr_columns] < 0.05
    ).all(axis=1)
    merged["minimum_empirical_p"] = merged[
        [
            f"empirical_p_{donor}_{profile}"
            for profile in PROFILES
            for donor in DONORS
        ]
    ].min(axis=1)
    return merged.sort_values(
        ["minimum_empirical_p", "pathway"]
    ).reset_index(drop=True)


def main() -> None:
    features = pd.read_csv(VKO_INPUTS / "vko_selected_features_v5.csv")
    with (ANALYSIS / "genesets_final.json").open(encoding="utf-8") as handle:
        gene_sets = json.load(handle)

    profile_results = {
        profile: read_profile(directory)
        for profile, directory in PROFILES.items()
    }
    pathway_results: dict[str, pd.DataFrame] = {}
    for profile, directory in PROFILES.items():
        table = matched_pathway_enrichment(
            profile_results[profile], features, gene_sets
        )
        table.insert(0, "profile", profile)
        output_path = (
            directory / "vko_sqstm1_pathway_enrichment_official_r.csv"
        )
        table.to_csv(output_path, index=False)
        pathway_results[profile] = table
        print(output_path)

    gene_table, correlations = gene_sensitivity(profile_results)
    gene_table.to_csv(
        RESULTS / "official_r_vko_cross_profile_gene_summary.csv", index=False
    )
    correlations.to_csv(
        RESULTS / "official_r_vko_cross_profile_correlations.csv", index=False
    )
    pathway_sensitivity(pathway_results).to_csv(
        RESULTS / "official_r_vko_pathway_sensitivity.csv", index=False
    )
    print(RESULTS / "official_r_vko_cross_profile_gene_summary.csv")
    print(RESULTS / "official_r_vko_pathway_sensitivity.csv")


if __name__ == "__main__":
    main()
