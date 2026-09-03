# -*- coding: utf-8 -*-
"""Build participant-level pseudobulk ranks for Hallmark stability analysis.

The primary ONFH ARCO 3A-versus-HOA comparison contains three independent
participants per group. This script reconstructs the six EC pseudobulks from
the archived count matrix, fits the full contrast, and repeats the model after
omitting each participant once. The resulting rank files are consumed by
``participant_fgsea_stability.R``. Leave-one-participant-out analyses quantify
sensitivity to a single participant; they are not an independent validation.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
import scanpy as sc
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


BASE_DIR = Path(__file__).resolve().parents[1]
TARGET_GROUPS = ("ONFH_3A", "HOA")


def log(message: str) -> None:
    print(f"{time.strftime('%H:%M:%S')} | {message}", flush=True)


def fit_contrast(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    metadata = pd.DataFrame({"group": groups.astype(str)}, index=counts.index)
    matrix = counts.loc[:, counts.sum(axis=0) > 0]
    dds = DeseqDataSet(
        counts=matrix.astype(int),
        metadata=metadata,
        design="~group",
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["group", "ONFH_3A", "HOA"], quiet=True)
    stats.summary()
    result = stats.results_df.copy()
    result.index.name = "gene"
    return result.reset_index().sort_values("stat", ascending=False)


def build_six_pseudobulks(source: Path, metadata: pd.DataFrame) -> pd.DataFrame:
    h5ad = source / "ec_final.h5ad"
    rds_candidates = [source / "counts_qc.rds", source / "counts_qc_v4.rds"]
    rds = next((path for path in rds_candidates if path.exists()), None)
    if not h5ad.exists() or rds is None:
        raise FileNotFoundError(
            "The source analysis directory must contain ec_final.h5ad and counts_qc.rds"
        )

    target_samples = metadata.loc[metadata["group"].isin(TARGET_GROUPS), "sample"].tolist()
    if len(target_samples) != 6:
        raise RuntimeError(f"Expected six inferential participants, found {len(target_samples)}")

    log("reading endothelial cell identities")
    ec = sc.read_h5ad(h5ad, backed="r")
    cell_to_sample = ec.obs["sample"].astype(str)
    all_samples = cell_to_sample.drop_duplicates().tolist()

    log("reading archived sparse count matrix")
    r_object = rdata.read_rds(str(rds))
    genes = np.asarray(r_object.Dimnames[0]).astype(str)
    cells = np.asarray(r_object.Dimnames[1]).astype(str)
    matrix = sp.csc_matrix(
        (
            np.asarray(r_object.x, dtype=np.float64),
            np.asarray(r_object.i),
            np.asarray(r_object.p),
        ),
        shape=tuple(r_object.Dim),
    )
    del r_object

    sample_of = cell_to_sample.reindex(cells)
    keep = sample_of.isin(all_samples).to_numpy()
    if int(keep.sum()) == 0:
        raise RuntimeError("No target endothelial cells aligned to the archived count matrix")
    matrix = matrix[:, keep]
    sample_of = sample_of.iloc[np.flatnonzero(keep)].astype(str).to_numpy()

    pseudobulk = np.zeros((len(all_samples), len(genes)), dtype=np.int64)
    for row, sample in enumerate(all_samples):
        columns = np.flatnonzero(sample_of == sample)
        if len(columns) == 0:
            raise RuntimeError(f"No endothelial cells found for {sample}")
        pseudobulk[row] = np.asarray(matrix[:, columns].sum(axis=1)).ravel().astype(np.int64)
        log(f"{sample}: {len(columns):,} endothelial cells")

    counts = pd.DataFrame(pseudobulk, index=all_samples, columns=genes)
    # Match the primary pipeline exactly: filter genes across every EC sampling
    # unit before selecting the six independent participants used for inference.
    counts = counts.loc[:, counts.sum(axis=0) >= 10]
    return counts.loc[target_samples]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-analysis",
        type=Path,
        default=Path(os.environ.get("ONFH_SOURCE_ANALYSIS", BASE_DIR / "analysis")),
        help="Directory containing the archived EC h5ad and raw-count RDS.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=BASE_DIR / "results" / "participant_fgsea_stability",
    )
    args = parser.parse_args()
    source = args.source_analysis.resolve()
    output = args.output_dir.resolve()
    ranks_dir = output / "de_ranks"
    ranks_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(BASE_DIR / "analysis" / "sample_metadata_v4.csv")
    inference_flag = metadata["independent_for_inference"].astype(str).str.lower().eq("true")
    mask = (
        metadata["dataset"].eq("liao_alcohol")
        & inference_flag
        & metadata["group"].isin(TARGET_GROUPS)
    )
    metadata = metadata.loc[mask, ["sample", "group", "participant_id"]].copy()
    counts = build_six_pseudobulks(source, metadata)
    counts.T.rename_axis("gene").reset_index().to_csv(
        output / "participant_ec_pseudobulk_counts.csv.gz", index=False, compression="gzip"
    )

    groups = metadata.set_index("sample").loc[counts.index, "group"]
    fits: list[dict[str, object]] = []
    designs = [("full", None)] + [(f"omit_{sample}", sample) for sample in counts.index]
    for fit_id, omitted in designs:
        keep_samples = counts.index if omitted is None else counts.index[counts.index != omitted]
        group_counts = groups.loc[keep_samples].value_counts().to_dict()
        if min(group_counts.values()) < 2:
            raise RuntimeError(f"{fit_id} leaves fewer than two participants in one group")
        log(f"fitting {fit_id}: {group_counts}")
        result = fit_contrast(counts.loc[keep_samples], groups.loc[keep_samples])
        result.to_csv(ranks_dir / f"{fit_id}.csv.gz", index=False, compression="gzip")
        fits.append(
            {
                "fit_id": fit_id,
                "omitted_sample": "" if omitted is None else omitted,
                "omitted_group": "" if omitted is None else groups.loc[omitted],
                "n_onfh_3a": int(group_counts.get("ONFH_3A", 0)),
                "n_hoa": int(group_counts.get("HOA", 0)),
                "n_genes": int(len(result)),
            }
        )

    pd.DataFrame(fits).to_csv(output / "participant_lopo_design.csv", index=False)
    provenance = {
        "analysis": "leave-one-participant-out DESeq2 rank stability",
        "contrast": "Liao ONFH ARCO 3A versus HOA",
        "biological_unit": "independent participant pseudobulk",
        "n_primary": {"ONFH_3A": 3, "HOA": 3},
        "n_leave_one_out_fits": 6,
        "interpretive_scope": (
            "sensitivity to removal of one participant; not an independent cohort validation"
        ),
        "source_analysis_directory": "external preprocessed analysis directory supplied at runtime",
        "source_input_files": ["ec_final.h5ad", "counts_qc.rds or counts_qc_v4.rds"],
    }
    (output / "participant_lopo_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    log(f"wrote {output}")


if __name__ == "__main__":
    main()
