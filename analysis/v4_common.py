# -*- coding: utf-8 -*-
"""Shared helpers for the v4 revision pipeline.

The v4 scripts are intentionally self-contained: project paths are resolved
relative to this file, plotting does not depend on ``daimon_runtime``, and
sample/library metadata are validated before any inferential analysis.
"""
from __future__ import annotations

import os
from functools import lru_cache
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


ROOT = Path(os.environ.get("ONFH_ROOT", Path(__file__).resolve().parents[1])).resolve()
ANALYSIS = ROOT / "analysis"
FIGS = ROOT / "figures" / "source"
MANUSCRIPT_FIGS = ROOT / "figures" / "final"
FIGS.mkdir(parents=True, exist_ok=True)
MANUSCRIPT_FIGS.mkdir(parents=True, exist_ok=True)

GROUP_ORDER = ["Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH"]
LIAO_INFERENCE_GROUPS = ["HOA", "ONFH_3A", "ONFH_4"]


def setup_plot() -> None:
    """Configure a deterministic, headless manuscript plotting style."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "savefig.facecolor": "white",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def bh_fdr(pvalues) -> np.ndarray:
    """Benjamini-Hochberg adjustment preserving missing values."""
    p = np.asarray(pvalues, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan, dtype=float)
    if not ok.any():
        return out
    values = p[ok]
    order = np.argsort(values)
    ranked = values[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out[np.where(ok)[0][order]] = q
    return out


def minimum_two_sided_mwu_p(n_a: int, n_b: int) -> float:
    """No-tie floor for an exhaustive two-sided Mann-Whitney test.

    With three observations per group, only ``C(6, 3)=20`` allocations exist,
    and the two most extreme allocations give a minimum two-sided p value of
    0.10. Ties can make the attainable p value still larger.
    """
    if n_a < 1 or n_b < 1:
        return np.nan
    return min(1.0, 2.0 / comb(n_a + n_b, n_a))


@lru_cache(maxsize=None)
def _group_a_allocations(n_total: int, n_a: int) -> tuple[tuple[int, ...], ...]:
    return tuple(combinations(range(n_total), n_a))


def exact_mwu_summary(a, b) -> dict[str, float]:
    """Tie-aware exhaustive two-sided Mann-Whitney summary.

    The p value is calculated by enumerating every allocation of the pooled
    observed values to groups of the original sizes. The rank-biserial effect
    is positive when values in ``a`` tend to exceed values in ``b``. The
    Hodges-Lehmann shift is the median of all pairwise ``a - b`` differences.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    n_a, n_b = len(a), len(b)
    if n_a < 2 or n_b < 2:
        return {
            "n_a": float(n_a), "n_b": float(n_b), "u": np.nan,
            "rank_biserial": np.nan, "hodges_lehmann": np.nan,
            "p": np.nan, "minimum_p_no_ties": minimum_two_sided_mwu_p(n_a, n_b),
        }

    pooled = np.concatenate([a, b])
    ranks = rankdata(pooled, method="average")
    u_observed = float(ranks[:n_a].sum() - n_a * (n_a + 1) / 2)
    center = n_a * n_b / 2
    observed_distance = abs(u_observed - center)

    extreme = 0
    allocations = _group_a_allocations(n_a + n_b, n_a)
    for indices in allocations:
        u = float(ranks[list(indices)].sum() - n_a * (n_a + 1) / 2)
        if abs(u - center) >= observed_distance - 1e-12:
            extreme += 1

    return {
        "n_a": float(n_a),
        "n_b": float(n_b),
        "u": u_observed,
        "rank_biserial": float(2 * u_observed / (n_a * n_b) - 1),
        "hodges_lehmann": float(np.median(a[:, None] - b[None, :])),
        "p": float(extreme / len(allocations)),
        "minimum_p_no_ties": minimum_two_sided_mwu_p(n_a, n_b),
    }


@lru_cache(maxsize=None)
def _label_partitions(group_sizes: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], ...], ...]:
    """Enumerate labelled partitions for small exact Kruskal-Wallis tests."""
    remaining = tuple(range(sum(group_sizes)))
    partitions: list[tuple[tuple[int, ...], ...]] = []

    def allocate(available: tuple[int, ...], group_index: int, chosen: list[tuple[int, ...]]) -> None:
        if group_index == len(group_sizes) - 1:
            partitions.append(tuple(chosen + [available]))
            return
        size = group_sizes[group_index]
        available_set = set(available)
        for group in combinations(available, size):
            rest = tuple(sorted(available_set.difference(group)))
            allocate(rest, group_index + 1, chosen + [tuple(group)])

    allocate(remaining, 0, [])
    return tuple(partitions)


def exact_kruskal_p(groups) -> float:
    """Exhaustive, tie-corrected Kruskal-Wallis p value for small groups."""
    arrays = [np.asarray(group, dtype=float) for group in groups]
    arrays = [array[np.isfinite(array)] for array in arrays]
    if len(arrays) < 2 or any(len(array) == 0 for array in arrays):
        return np.nan
    pooled = np.concatenate(arrays)
    ranks = rankdata(pooled, method="average")
    n_total = len(pooled)
    _, tie_counts = np.unique(pooled, return_counts=True)
    correction = 1.0 - np.sum(tie_counts**3 - tie_counts) / (n_total**3 - n_total)
    if correction <= 0:
        return 1.0

    def statistic(partition) -> float:
        h = 12.0 / (n_total * (n_total + 1)) * sum(
            ranks[list(indices)].sum() ** 2 / len(indices) for indices in partition
        ) - 3.0 * (n_total + 1)
        return float(h / correction)

    sizes = tuple(len(array) for array in arrays)
    observed_partition = []
    start = 0
    for size in sizes:
        observed_partition.append(tuple(range(start, start + size)))
        start += size
    observed = statistic(tuple(observed_partition))
    partitions = _label_partitions(sizes)
    extreme = sum(statistic(partition) >= observed - 1e-12 for partition in partitions)
    return float(extreme / len(partitions))


def read_sample_metadata() -> pd.DataFrame:
    path = ANALYSIS / "sample_metadata_v4.csv"
    meta = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {
        "sample",
        "dataset",
        "group",
        "participant_id",
        "inferential_unit",
        "independent_for_inference",
        "source_note",
    }
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"sample_metadata_v4.csv missing columns: {sorted(missing)}")
    if meta["sample"].duplicated().any():
        dup = meta.loc[meta["sample"].duplicated(), "sample"].tolist()
        raise ValueError(f"duplicate sample rows: {dup}")
    meta["independent_for_inference"] = (
        meta["independent_for_inference"].str.lower().map({"true": True, "false": False})
    )
    if meta["independent_for_inference"].isna().any():
        raise ValueError("independent_for_inference must contain only true/false")
    return meta


def attach_v4_metadata(adata, *, copy: bool = False):
    """Attach validated participant/inference metadata to an AnnData object."""
    if copy:
        adata = adata.copy()
    meta = read_sample_metadata().set_index("sample")
    observed = set(adata.obs["sample"].astype(str).unique())
    absent = observed.difference(meta.index)
    if absent:
        raise ValueError(f"samples missing from sample_metadata_v4.csv: {sorted(absent)}")
    sample = adata.obs["sample"].astype(str)
    for col in [
        "participant_id",
        "inferential_unit",
        "independent_for_inference",
        "source_note",
    ]:
        adata.obs[col] = sample.map(meta[col]).to_numpy()
    return adata


def liao_inference_samples(meta: pd.DataFrame | None = None) -> list[str]:
    """Return confirmed independent Liao samples used for formal inference."""
    if meta is None:
        meta = read_sample_metadata()
    mask = (
        (meta["dataset"] == "liao_alcohol")
        & meta["group"].isin(LIAO_INFERENCE_GROUPS)
        & meta["independent_for_inference"]
    )
    return meta.loc[mask, "sample"].tolist()


def save_figure(fig, path: Path, **kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white", **kwargs)
