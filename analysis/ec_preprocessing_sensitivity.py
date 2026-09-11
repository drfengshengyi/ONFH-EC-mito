# -*- coding: utf-8 -*-
"""Sensitivity audits for EC integration, annotation, and low-depth pseudobulk.

This script complements the frozen primary analysis without overwriting it. It
performs three prespecified checks requested during Scientific Reports revision:

1. Rebuild the EC PCA/neighbour graph/Leiden solution without ComBat.
2. Scan top-two annotation-margin thresholds of 0.10, 0.15, and 0.20.
3. Audit participant-level EC counts and raw UMI depth, including the low-EC
   HOA1 pseudobulk, and extract the existing omit-HOA1 pathway refit.

The input EC object contains normalized log1p expression plus per-cell ``nCount``
values. No serum analysis is run here.
"""
from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score

from v4_common import LIAO_INFERENCE_GROUPS, exact_kruskal_p, exact_mwu_summary


SUBTYPE_PANELS = {
    "typeH_EMCN_KDR": ["EMCN", "KDR", "CDH5", "DACH1", "VWF"],
    "arterial": ["SOX17", "HEY1", "GJA5", "DLL4", "EFNB2", "SOX13"],
    "venous_ACKR1": ["ACKR1", "SELE", "VCAM1", "PLAT"],
    "tip_angiogenic": ["ESM1", "ANGPT2", "APLN", "DLL4"],
    "lymphatic": ["PROX1", "LYVE1", "PDPN"],
    "typeR_bone_remodel": ["SMAD1", "PPARG", "NOTCH4", "COL4A1"],
}
MAJOR_CELL_PANELS = {
    "EC": ["PECAM1", "VWF", "CDH5", "EMCN", "KDR"],
    "Pericyte_SMC": ["RGS5", "ACTA2", "MYH11", "PDGFRB"],
    "Osteoblast": ["RUNX2", "SP7", "ALPL", "IBSP", "BGLAP", "COL1A1", "SPP1"],
    "MSC_stromal": ["CXCL12", "NT5E", "THY1", "ENG", "APOD", "CFD"],
    "Adipo_lineage": ["ADIPOQ", "PLIN1", "FABP4", "LPL"],
    "Chondrocyte": ["COL2A1", "ACAN", "SOX9"],
    "Osteoclast": ["ACP5", "CTSK", "MMP9", "TNFRSF11A", "CALCR"],
    "Myeloid": ["CD68", "CD14", "LYZ", "CSF1R", "AIF1"],
    "T_cell": ["CD3D", "CD3E", "CD8A", "CD4"],
    "NK": ["GNLY", "NKG7", "NCAM1"],
    "B_cell": ["CD19", "MS4A1", "CD79A"],
    "Plasma": ["JCHAIN", "MZB1", "XBP1"],
    "Granulocyte": ["S100A8", "S100A9", "FCGR3B", "MPO"],
    "Mast": ["TPSAB1", "TPSB2", "KIT"],
    "Platelet": ["PPBP", "PF4", "GP9"],
    "RBC": ["HBB", "HBA1", "HBA2"],
    "pDC": ["LILRA4", "IRF7", "CLEC4C"],
}
THRESHOLDS = (0.10, 0.15, 0.20)
PAIRWISE = (("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-analysis",
        required=True,
        type=Path,
        help="Directory containing ec_final.h5ad and ec_subtype_scores.csv.",
    )
    parser.add_argument(
        "--project-root",
        default=Path(__file__).resolve().parents[1],
        type=Path,
        help="Repository root containing results/participant_fgsea_stability.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        type=Path,
        help="Output directory (default: <project-root>/results/sr_preprocessing_sensitivity).",
    )
    return parser.parse_args()


def panel_scores(
    adata: sc.AnnData, cluster_key: str, panel_definition: dict[str, list[str]]
) -> pd.DataFrame:
    panels = {
        name: [gene for gene in genes if gene in adata.var_names]
        for name, genes in panel_definition.items()
    }
    empty = [name for name, genes in panels.items() if not genes]
    if empty:
        raise ValueError(f"No measured genes for panels: {empty}")

    cluster = adata.obs[cluster_key].astype(str).to_numpy()
    cluster_ids = sorted(np.unique(cluster), key=lambda value: int(value))
    gene_index = {gene: i for i, gene in enumerate(adata.var_names)}
    matrix = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    scores = pd.DataFrame(index=cluster_ids, columns=panels, dtype=float)
    for cluster_id in cluster_ids:
        rows = np.flatnonzero(cluster == cluster_id)
        for panel, genes in panels.items():
            columns = [gene_index[gene] for gene in genes]
            scores.loc[cluster_id, panel] = float(
                np.asarray(matrix[rows][:, columns].mean(axis=0)).ravel().mean()
            )
    return scores


def score_audit(scores: pd.DataFrame) -> pd.DataFrame:
    values = scores.to_numpy()
    ranked = np.argsort(-values, axis=1)
    labels = scores.columns.to_numpy()
    audit = pd.DataFrame(
        {
            "cluster": scores.index.astype(str),
            "top_panel": labels[ranked[:, 0]],
            "second_panel": labels[ranked[:, 1]],
            "top_score": values[np.arange(len(scores)), ranked[:, 0]],
            "second_score": values[np.arange(len(scores)), ranked[:, 1]],
        }
    )
    audit["top_two_margin"] = audit["top_score"] - audit["second_score"]
    return audit


def composition(obs: pd.DataFrame, state_key: str) -> pd.DataFrame:
    counts = obs.groupby(["sample", state_key], observed=True).size().unstack(fill_value=0)
    fractions = counts.div(counts.sum(axis=1), axis=0) * 100
    metadata = (
        obs[["sample", "dataset", "group"]]
        .drop_duplicates("sample")
        .set_index("sample")
    )
    return metadata.join(fractions, how="inner").reset_index()


def type_h_statistics(comp: pd.DataFrame, analysis: str) -> pd.DataFrame:
    type_h = "typeH_EMCN_KDR"
    if type_h not in comp.columns:
        comp = comp.copy()
        comp[type_h] = 0.0
    liao = comp[comp["group"].isin(LIAO_INFERENCE_GROUPS)].copy()
    grouped = [liao.loc[liao["group"] == group, type_h].to_numpy() for group in LIAO_INFERENCE_GROUPS]
    rows: list[dict[str, object]] = [
        {
            "analysis": analysis,
            "test": "Liao_exact_Kruskal_Wallis",
            "group_a": "HOA_ONFH_3A_ONFH_4",
            "group_b": "",
            "n_a": len(liao),
            "n_b": np.nan,
            "median_a": np.nan,
            "median_b": np.nan,
            "median_difference_percentage_points": np.nan,
            "hodges_lehmann": np.nan,
            "rank_biserial": np.nan,
            "exact_p": exact_kruskal_p(grouped),
        }
    ]
    for group_a, group_b in PAIRWISE:
        a = liao.loc[liao["group"] == group_a, type_h].to_numpy()
        b = liao.loc[liao["group"] == group_b, type_h].to_numpy()
        summary = exact_mwu_summary(a, b)
        rows.append(
            {
                "analysis": analysis,
                "test": "exact_Mann_Whitney",
                "group_a": group_a,
                "group_b": group_b,
                "n_a": len(a),
                "n_b": len(b),
                "median_a": float(np.median(a)),
                "median_b": float(np.median(b)),
                "median_difference_percentage_points": float(np.median(a) - np.median(b)),
                "hodges_lehmann": summary["hodges_lehmann"],
                "rank_biserial": summary["rank_biserial"],
                "exact_p": summary["p"],
            }
        )
    return pd.DataFrame(rows)


def run_no_combat(adata: sc.AnnData, output_dir: Path) -> dict[str, object]:
    hvg = adata.var.get("highly_variable")
    if hvg is None or int(np.asarray(hvg).sum()) != 1500:
        sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor="seurat", subset=False)
        hvg = adata.var["highly_variable"]

    reduced = adata[:, np.asarray(hvg, dtype=bool)].copy()
    reduced.X = reduced.X.toarray() if sp.issparse(reduced.X) else np.asarray(reduced.X).copy()
    sc.pp.scale(reduced, max_value=10)
    sc.pp.pca(reduced, n_comps=20, svd_solver="arpack")
    adata.obsm["X_pca_no_combat"] = reduced.obsm["X_pca"]
    sc.pp.neighbors(adata, use_rep="X_pca_no_combat", n_neighbors=15, n_pcs=20)
    sc.tl.leiden(
        adata,
        key_added="leiden_no_combat",
        resolution=0.4,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=0,
    )

    scores = panel_scores(adata, "leiden_no_combat", SUBTYPE_PANELS)
    audit = score_audit(scores)
    subtype_map = audit.set_index("cluster")["top_panel"]
    adata.obs["EC_subtype_no_combat"] = (
        adata.obs["leiden_no_combat"].astype(str).map(subtype_map).astype(str)
    )

    original_cluster = adata.obs["leiden"].astype(str)
    no_combat_cluster = adata.obs["leiden_no_combat"].astype(str)
    original_subtype = adata.obs["EC_subtype"].astype(str)
    no_combat_subtype = adata.obs["EC_subtype_no_combat"].astype(str)
    comp = composition(adata.obs, "EC_subtype_no_combat")
    stats = type_h_statistics(comp, "no_ComBat")

    audit.to_csv(output_dir / "no_combat_cluster_annotation_scores.csv", index=False)
    comp.to_csv(output_dir / "no_combat_ec_state_composition.csv", index=False)
    stats.to_csv(output_dir / "no_combat_typeH_exact_statistics.csv", index=False)

    summary = {
        "n_cells": int(adata.n_obs),
        "n_hvg": int(np.asarray(hvg).sum()),
        "n_original_clusters": int(original_cluster.nunique()),
        "n_no_combat_clusters": int(no_combat_cluster.nunique()),
        "cluster_adjusted_rand_index": float(adjusted_rand_score(original_cluster, no_combat_cluster)),
        "subtype_adjusted_rand_index": float(adjusted_rand_score(original_subtype, no_combat_subtype)),
        "cell_level_subtype_concordance": float((original_subtype == no_combat_subtype).mean()),
        "no_combat_typeH_exact_omnibus_p": float(stats.iloc[0]["exact_p"]),
    }
    pd.DataFrame([summary]).to_csv(output_dir / "no_combat_summary.csv", index=False)
    return summary


def run_atlas_no_combat(
    atlas_path: Path, output_dir: Path
) -> dict[str, object]:
    """Repeat atlas neighbours/clustering from the archived uncorrected PCA."""
    atlas = sc.read_h5ad(atlas_path)
    if atlas.n_obs != 115572:
        raise ValueError(
            f"Expected the frozen 115,572-cell atlas; found {atlas.n_obs}"
        )
    if "X_pca" not in atlas.obsm:
        raise ValueError("atlas_annotated.h5ad lacks the archived uncorrected X_pca")

    sc.pp.neighbors(
        atlas,
        use_rep="X_pca",
        n_neighbors=15,
        n_pcs=30,
        random_state=0,
    )
    sc.tl.leiden(
        atlas,
        key_added="leiden_no_combat",
        resolution=0.5,
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=0,
    )
    scores = panel_scores(atlas, "leiden_no_combat", MAJOR_CELL_PANELS)
    audit = score_audit(scores)
    type_map = audit.set_index("cluster")["top_panel"]
    atlas.obs["cell_type_no_combat"] = (
        atlas.obs["leiden_no_combat"].astype(str).map(type_map).astype(str)
    )

    original_cluster = atlas.obs["leiden"].astype(str)
    no_combat_cluster = atlas.obs["leiden_no_combat"].astype(str)
    original_type = atlas.obs["cell_type"].astype(str)
    no_combat_type = atlas.obs["cell_type_no_combat"].astype(str)
    original_ec = original_type.eq("EC")
    no_combat_ec = no_combat_type.eq("EC")
    intersection = int((original_ec & no_combat_ec).sum())
    union = int((original_ec | no_combat_ec).sum())

    per_sample = (
        atlas.obs.assign(
            original_ec=original_ec.to_numpy(dtype=int),
            no_combat_ec=no_combat_ec.to_numpy(dtype=int),
        )
        .groupby("sample", observed=True)
        .agg(
            dataset=("dataset", "first"),
            group=("group", "first"),
            retained_cells=("sample", "size"),
            original_combat_ec_cells=("original_ec", "sum"),
            no_combat_ec_cells=("no_combat_ec", "sum"),
        )
        .reset_index()
    )
    per_sample["original_combat_ec_percent"] = (
        100 * per_sample["original_combat_ec_cells"] / per_sample["retained_cells"]
    )
    per_sample["no_combat_ec_percent"] = (
        100 * per_sample["no_combat_ec_cells"] / per_sample["retained_cells"]
    )
    per_sample["difference_percentage_points"] = (
        per_sample["no_combat_ec_percent"] - per_sample["original_combat_ec_percent"]
    )

    audit.to_csv(output_dir / "atlas_no_combat_cluster_annotation_scores.csv", index=False)
    per_sample.to_csv(output_dir / "atlas_no_combat_ec_composition.csv", index=False)
    summary = {
        "n_cells": int(atlas.n_obs),
        "n_original_clusters": int(original_cluster.nunique()),
        "n_no_combat_clusters": int(no_combat_cluster.nunique()),
        "cluster_adjusted_rand_index": float(adjusted_rand_score(original_cluster, no_combat_cluster)),
        "cell_type_adjusted_rand_index": float(adjusted_rand_score(original_type, no_combat_type)),
        "cell_level_type_concordance": float((original_type == no_combat_type).mean()),
        "original_combat_ec_cells": int(original_ec.sum()),
        "no_combat_ec_cells": int(no_combat_ec.sum()),
        "ec_overlap_cells": intersection,
        "ec_jaccard": float(intersection / union),
    }
    pd.DataFrame([summary]).to_csv(output_dir / "atlas_no_combat_summary.csv", index=False)
    return summary


def run_threshold_scan(
    adata: sc.AnnData, source_scores: Path, output_dir: Path
) -> pd.DataFrame:
    scores = pd.read_csv(source_scores, index_col=0)
    scores.index = scores.index.astype(str)
    audit = score_audit(scores)
    cluster = adata.obs["leiden"].astype(str)
    outputs: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []

    for threshold in THRESHOLDS:
        label = f"margin_lt_{threshold:.2f}"
        audit[label] = audit["top_two_margin"] < threshold
        strict_map = audit.set_index("cluster").apply(
            lambda row: "Ambiguous" if row[label] else row["top_panel"], axis=1
        )
        state_key = f"strict_subtype_{threshold:.2f}"
        adata.obs[state_key] = cluster.map(strict_map).fillna("Ambiguous")
        comp = composition(adata.obs, state_key)
        stats = type_h_statistics(comp.rename(columns={state_key: "strict_subtype"}), f"margin_{threshold:.2f}")
        stats.insert(1, "margin_threshold", threshold)
        outputs.append(stats)
        ambiguous_cells = int((adata.obs[state_key] == "Ambiguous").sum())
        summaries.append(
            {
                "margin_threshold": threshold,
                "ambiguous_clusters": int(audit[label].sum()),
                "total_clusters": len(audit),
                "ambiguous_cells": ambiguous_cells,
                "ambiguous_cells_percent": 100 * ambiguous_cells / adata.n_obs,
                "typeH_exact_omnibus_p": float(stats.iloc[0]["exact_p"]),
            }
        )
        comp.to_csv(
            output_dir / f"annotation_margin_{threshold:.2f}_composition.csv", index=False
        )

    audit.to_csv(output_dir / "annotation_margin_threshold_cluster_audit.csv", index=False)
    pd.concat(outputs, ignore_index=True).to_csv(
        output_dir / "annotation_margin_threshold_exact_statistics.csv", index=False
    )
    summary = pd.DataFrame(summaries)
    summary.to_csv(output_dir / "annotation_margin_threshold_summary.csv", index=False)
    return summary


def run_low_depth_audit(
    adata: sc.AnnData, project_root: Path, output_dir: Path
) -> dict[str, object]:
    obs = adata.obs.copy()
    if "nCount" not in obs:
        raise ValueError("ec_final.h5ad lacks per-cell nCount required for the depth audit")
    sample = (
        obs.groupby("sample", observed=True)
        .agg(
            dataset=("dataset", "first"),
            group=("group", "first"),
            retained_ec_cells=("sample", "size"),
            stored_ec_nCount_sum=("nCount", "sum"),
            median_stored_ec_nCount=("nCount", "median"),
        )
        .reset_index()
    )
    liao_groups = set(LIAO_INFERENCE_GROUPS)
    sample["included_in_liao_pseudobulk"] = sample["group"].isin(liao_groups)
    sample["pseudobulk_inclusion_rule"] = np.where(
        sample["included_in_liao_pseudobulk"],
        "all independent Liao participants with at least one retained EC; no post hoc minimum EC threshold",
        "not part of the formal Liao HOA/ONFH_3A/ONFH_4 pseudobulk contrasts",
    )
    sample["included_in_SQSTM1_virtual_KO"] = sample["sample"].isin(["hoa2", "hoa3"])
    sample["virtual_KO_note"] = np.where(
        sample["sample"].eq("hoa1"),
        "excluded from donor-network perturbation because only 26 ECs remained",
        "separate from the low-depth pseudobulk inclusion decision",
    )
    sample.to_csv(output_dir / "pseudobulk_sampling_unit_depth_audit.csv", index=False)

    fgsea_path = (
        project_root
        / "results"
        / "participant_fgsea_stability"
        / "fgsea_selected_pathways_by_omission.csv"
    )
    fgsea = pd.read_csv(fgsea_path)
    hoa1_refit = fgsea.loc[fgsea["fit_id"].eq("omit_hoa1")].copy()
    hoa1_refit.to_csv(output_dir / "omit_hoa1_hallmark_refit.csv", index=False)
    key_pathways = {
        "HALLMARK_ALLOGRAFT_REJECTION",
        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
        "HALLMARK_INFLAMMATORY_RESPONSE",
    }
    key = hoa1_refit.loc[hoa1_refit["pathway"].isin(key_pathways)]
    hoa1 = sample.loc[sample["sample"].eq("hoa1")].iloc[0]
    summary = {
        "hoa1_retained_ec_cells": int(hoa1["retained_ec_cells"]),
        "hoa1_stored_ec_nCount_sum": int(round(float(hoa1["stored_ec_nCount_sum"]))),
        "hoa1_key_pathways_positive": int((key["NES"] > 0).sum()),
        "hoa1_key_pathways_fdr_lt_0_05": int((key["padj"] < 0.05).sum()),
        "hoa1_key_pathways_total": int(len(key)),
        "hoa1_key_pathway_NES_min": float(key["NES"].min()),
        "hoa1_key_pathway_NES_max": float(key["NES"].max()),
    }
    pd.DataFrame([summary]).to_csv(output_dir / "hoa1_low_depth_summary.csv", index=False)
    return summary


def main() -> None:
    args = parse_args()
    source_analysis = args.source_analysis.resolve()
    project_root = args.project_root.resolve()
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else project_root / "results" / "sr_preprocessing_sensitivity"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = source_analysis / "ec_final.h5ad"
    atlas_path = source_analysis / "atlas_annotated.h5ad"
    scores_path = source_analysis / "ec_subtype_scores.csv"
    if not h5ad_path.exists() or not atlas_path.exists() or not scores_path.exists():
        raise FileNotFoundError(
            "source-analysis must contain atlas_annotated.h5ad, ec_final.h5ad, "
            "and ec_subtype_scores.csv"
        )
    adata = sc.read_h5ad(h5ad_path)
    if adata.n_obs != 13426:
        raise ValueError(
            f"Expected the frozen 13,426-cell EC object used by the manuscript; found {adata.n_obs}"
        )

    threshold_summary = run_threshold_scan(adata, scores_path, output_dir)
    low_depth_summary = run_low_depth_audit(adata, project_root, output_dir)
    no_combat_summary = run_no_combat(adata, output_dir)
    atlas_no_combat_summary = run_atlas_no_combat(atlas_path, output_dir)

    provenance = {
        "source_h5ad": str(h5ad_path),
        "source_h5ad_bytes": h5ad_path.stat().st_size,
        "source_h5ad_mtime": h5ad_path.stat().st_mtime,
        "source_scores": str(scores_path),
        "source_atlas_h5ad": str(atlas_path),
        "output_dir": str(output_dir),
        "thresholds": list(THRESHOLDS),
        "scanpy": version("scanpy"),
        "anndata": version("anndata"),
        "scikit_learn": version("scikit-learn"),
        "random_state": 0,
        "serum_analysis_rerun": False,
        "no_combat_summary": no_combat_summary,
        "atlas_no_combat_summary": atlas_no_combat_summary,
        "threshold_summary": threshold_summary.to_dict(orient="records"),
        "low_depth_summary": low_depth_summary,
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(provenance, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
