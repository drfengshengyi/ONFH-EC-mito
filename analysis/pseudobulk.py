# -*- coding: utf-8 -*-
"""Version-4 EC pseudobulk analysis with explicit inferential units.

Formal p values/FDR are written only for independent Liao participants.  The
SONFH-versus-HOA model is used only to estimate descriptive library-level
log2 fold changes; its p values are deliberately discarded because four SONFH
libraries cannot be mapped to the three reported patients.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd
import rdata
import scanpy as sc
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from v4_common import ANALYSIS, attach_v4_metadata, read_sample_metadata

LOG = ANALYSIS / "py_pseudobulk_v4.log"


def lg(*parts) -> None:
    line = time.strftime("%H:%M:%S") + " | " + " ".join(map(str, parts))
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def fit_contrast(counts: pd.DataFrame, metadata: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    use = metadata[metadata["group"].isin([a, b])].copy()
    matrix = counts.loc[use.index]
    matrix = matrix.loc[:, matrix.sum(axis=0) > 0]
    dds = DeseqDataSet(
        counts=matrix.astype(int), metadata=use[["group"]], design="~group",
        refit_cooks=True, quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["group", a, b], quiet=True)
    stats.summary()
    return stats.results_df.sort_values("padj")


if LOG.exists():
    LOG.unlink()

ec = attach_v4_metadata(sc.read_h5ad(ANALYSIS / "ec_final.h5ad"))
ec_cells = set(ec.obs_names.astype(str))
ec_samples = ec.obs[["sample"]].copy()

r_object = rdata.read_rds(str(ANALYSIS / "counts_qc_v4.rds"))
genes = np.asarray(r_object.Dimnames[0]).astype(str)
cells = np.asarray(r_object.Dimnames[1]).astype(str)
matrix = sp.csc_matrix(
    (np.asarray(r_object.x, dtype=np.float64), np.asarray(r_object.i), np.asarray(r_object.p)),
    shape=tuple(r_object.Dim),
)
del r_object

keep = np.fromiter((cell in ec_cells for cell in cells), dtype=bool, count=len(cells))
cells = cells[keep]
matrix = matrix[:, keep]
sample_of = ec_samples["sample"].astype(str).reindex(cells)
if sample_of.isna().any():
    raise RuntimeError("EC cell names do not align between counts_qc.rds and ec_final.h5ad")

samples = sample_of.unique().tolist()
pseudobulk = np.zeros((len(samples), len(genes)), dtype=np.int64)
for row, sample in enumerate(samples):
    columns = np.flatnonzero(sample_of.to_numpy() == sample)
    pseudobulk[row] = np.asarray(matrix[:, columns].sum(axis=1)).ravel().astype(np.int64)
counts = pd.DataFrame(pseudobulk, index=samples, columns=genes)
counts = counts.loc[:, counts.sum(axis=0) >= 10]
metadata = read_sample_metadata().set_index("sample").loc[samples]
lg("pseudobulk", counts.shape)

summary = []
for a, b in [("ONFH_3A", "HOA"), ("ONFH_4", "ONFH_3A"), ("ONFH_4", "HOA")]:
    use = metadata[
        (metadata["dataset"] == "liao_alcohol")
        & metadata["independent_for_inference"]
        & metadata["group"].isin([a, b])
    ]
    result = fit_contrast(counts, use, a, b)
    result["inferential_status"] = "formal_independent_liao_participants"
    tag = f"{a}_vs_{b}"
    result.to_csv(ANALYSIS / f"de_ec_{tag}_v4.csv")
    summary.append({
        "contrast": tag,
        "n_a": int((use["group"] == a).sum()),
        "n_b": int((use["group"] == b).sum()),
        "status": "formal_independent_liao_participants",
    })
    lg(tag, "formal", len(result), "genes")

# Estimate but do not test the cross-cohort SONFH library effect.
a, b = "SONFH", "HOA"
use = metadata[metadata["group"].isin([a, b])]
descriptive = fit_contrast(counts, use, a, b)[["baseMean", "log2FoldChange"]].copy()
descriptive["inferential_status"] = (
    "descriptive_only_four_SONFH_libraries_from_three_reported_patients_mapping_unavailable"
)
descriptive.to_csv(ANALYSIS / "de_ec_SONFH_vs_HOA_descriptive_v4.csv")
summary.append({
    "contrast": "SONFH_vs_HOA",
    "n_a": int((use["group"] == a).sum()),
    "n_b": int((use["group"] == b).sum()),
    "status": "descriptive_library_effect_only_no_p_or_fdr",
})
pd.DataFrame(summary).to_csv(ANALYSIS / "pseudobulk_contrast_audit_v4.csv", index=False)
lg("SONFH_vs_HOA descriptive", len(descriptive), "genes; p/FDR omitted")
