#!/usr/bin/env python3
"""Descriptive spatial contextualization of prespecified ONFH EC programs.

This analysis uses the single-section osteoarthritic femoral-head Visium
CytAssist reference in GEO series GSE284089 (sample GSM8677818). It is an
anatomical contextualization only: spots are not biological replicates and the
output must not be interpreted as disease-matched validation.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp
from PIL import Image
from scipy.io import mmread
from scipy.stats import spearmanr


ACCESSION = "GSE284089"
SAMPLE = "GSM8677818"
SOURCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE284nnn/"
    "GSE284089/suppl/GSE284089_RAW.tar"
)

# The EC panel expands the five markers used for atlas annotation with
# canonical vascular-endothelial markers to improve robustness in sparse FFPE
# Visium data. State programs are frozen from the manuscript gene-set file.
EC_IDENTITY = [
    "PECAM1", "VWF", "CDH5", "CLDN5", "EMCN", "KDR", "FLT1", "ENG",
    "RAMP2", "RAMP3", "PLVAP", "KLF2", "KLF4", "ESAM", "EGFL7",
    "ECSCR", "EPAS1", "ESM1", "TEK", "ERG", "CLEC14A", "CA4", "RGCC",
]


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=root / "data" / "gse284089")
    parser.add_argument(
        "--gene-sets", type=Path, default=root / "analysis" / "genesets_final.json"
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=root / "results" / "spatial_contextualization",
    )
    parser.add_argument(
        "--figure-dir", type=Path, default=root / "figures" / "source"
    )
    parser.add_argument(
        "--final-figure-dir", type=Path, default=root / "figures" / "final"
    )
    parser.add_argument("--minimum-umi", type=int, default=100)
    parser.add_argument("--control-genes-per-target", type=int, default=20)
    parser.add_argument("--seed", type=int, default=284089)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gzip_json(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    sd = values.std(ddof=1)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values)
    return (values - values.mean()) / sd


def build_expression_matched_controls(
    mean_expression: np.ndarray,
    feature_names: np.ndarray,
    targets: list[str],
    all_signature_genes: set[str],
    controls_per_target: int,
    rng: np.random.Generator,
) -> list[str]:
    """Choose deterministic expression-bin-matched background genes."""
    frame = pd.DataFrame(
        {"gene": feature_names, "mean_expression": mean_expression}
    )
    frame = frame.loc[frame["mean_expression"] > 0].copy()
    # Rank-first qcut remains well defined despite many tied low-expression genes.
    frame["bin"] = pd.qcut(
        frame["mean_expression"].rank(method="first"),
        q=min(24, len(frame)),
        labels=False,
        duplicates="drop",
    )
    gene_to_bin = dict(zip(frame["gene"], frame["bin"]))
    chosen: list[str] = []
    for gene in targets:
        if gene not in gene_to_bin:
            continue
        pool = frame.loc[
            (frame["bin"] == gene_to_bin[gene])
            & (~frame["gene"].isin(all_signature_genes)),
            "gene",
        ].to_numpy()
        if len(pool) == 0:
            continue
        take = min(controls_per_target, len(pool))
        chosen.extend(rng.choice(pool, size=take, replace=False).tolist())
    return sorted(set(chosen))


def mean_rows(matrix: sp.csr_matrix, rows: list[int]) -> np.ndarray:
    if not rows:
        return np.full(matrix.shape[1], np.nan)
    return np.asarray(matrix[rows, :].mean(axis=0)).ravel()


def add_panel_label(ax: mpl.axes.Axes, label: str) -> None:
    ax.text(
        -0.055,
        1.035,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="right",
        va="bottom",
    )


def draw_spatial_panel(
    ax: mpl.axes.Axes,
    image: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray | None,
    title: str,
    crop: tuple[float, float, float, float],
    cmap: str = "RdBu_r",
    norm: mpl.colors.Normalize | None = None,
    point_size: float = 9.5,
    tissue_alpha: float = 0.58,
    colorbar_label: str | None = None,
) -> mpl.collections.PathCollection | None:
    height, width = image.shape[:2]
    ax.imshow(image, extent=(0, width, height, 0), interpolation="none", alpha=tissue_alpha)
    artist = None
    if values is None:
        ax.scatter(
            x,
            y,
            s=point_size,
            facecolors="none",
            edgecolors="#177E7B",
            linewidths=0.28,
            alpha=0.65,
            rasterized=True,
        )
    else:
        artist = ax.scatter(
            x,
            y,
            c=values,
            s=point_size,
            cmap=cmap,
            norm=norm,
            linewidths=0,
            alpha=0.92,
            rasterized=True,
        )
    xmin, xmax, ymin, ymax = crop
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(title, fontsize=9.2, fontweight="bold", pad=5)
    if artist is not None and colorbar_label:
        bar = ax.figure.colorbar(artist, ax=ax, fraction=0.035, pad=0.018)
        bar.ax.tick_params(labelsize=6.7, length=2)
        bar.set_label(colorbar_label, fontsize=7.2)
        bar.outline.set_linewidth(0.4)
    return artist


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.final_figure_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "matrix": args.data_dir / f"{SAMPLE}_matrix.mtx.gz",
        "features": args.data_dir / f"{SAMPLE}_features.tsv.gz",
        "barcodes": args.data_dir / f"{SAMPLE}_barcodes.tsv.gz",
        "positions": args.data_dir / f"{SAMPLE}_tissue_positions.csv.gz",
        "scales": args.data_dir / f"{SAMPLE}_scalefactors_json.json.gz",
        "image": args.data_dir / f"{SAMPLE}_tissue_hires_image.png.gz",
        "archive": args.data_dir / f"{ACCESSION}_RAW.tar",
    }
    missing = [str(path) for path in files.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required spatial inputs:\n" + "\n".join(missing))

    features = pd.read_csv(files["features"], sep="\t", header=None)
    barcodes = pd.read_csv(files["barcodes"], sep="\t", header=None)[0].astype(str)
    positions = pd.read_csv(files["positions"])
    matrix = mmread(files["matrix"], spmatrix=True).tocsr().astype(np.float64)
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"Matrix shape {matrix.shape} does not match features/barcodes "
            f"({len(features)}, {len(barcodes)})"
        )
    gene_names = features[1].astype(str).to_numpy()
    duplicate_symbol_rows = int(pd.Index(gene_names).duplicated(keep=False).sum())
    if duplicate_symbol_rows:
        # GEO supplies Ensembl IDs and display symbols. Repeated display symbols
        # are summed deterministically before any normalization or scoring.
        codes, unique_symbols = pd.factorize(gene_names, sort=False)
        aggregator = sp.csr_matrix(
            (np.ones(len(codes)), (codes, np.arange(len(codes)))),
            shape=(len(unique_symbols), len(codes)),
        )
        matrix = (aggregator @ matrix).tocsr()
        gene_names = np.asarray(unique_symbols, dtype=str)

    total_umi = np.asarray(matrix.sum(axis=0)).ravel()
    detected_genes = np.asarray((matrix > 0).sum(axis=0)).ravel()
    keep = total_umi >= args.minimum_umi
    matrix = matrix[:, keep]
    kept_barcodes = barcodes.to_numpy()[keep]
    total_umi = total_umi[keep]
    detected_genes = detected_genes[keep]

    position_index = positions.set_index("barcode")
    if not pd.Index(kept_barcodes).isin(position_index.index).all():
        raise ValueError("Not all retained barcodes have spatial coordinates")
    positions = position_index.loc[kept_barcodes].reset_index()

    scale_factors = read_gzip_json(files["scales"])
    scale = float(scale_factors["tissue_hires_scalef"])
    x = positions["pxl_col_in_fullres"].to_numpy(float) * scale
    y = positions["pxl_row_in_fullres"].to_numpy(float) * scale
    with gzip.open(files["image"], "rb") as handle:
        image = np.asarray(Image.open(handle).convert("RGB"))

    # Library-size normalization to 10,000 followed by log1p. Raw counts remain
    # untouched on disk and are not redistributed by the repository.
    inv_depth = np.divide(1.0e4, total_umi, out=np.zeros_like(total_umi), where=total_umi > 0)
    lognorm = matrix @ sp.diags(inv_depth)
    lognorm.data = np.log1p(lognorm.data)
    lognorm = lognorm.tocsr()

    with args.gene_sets.open("r", encoding="utf-8") as handle:
        frozen_sets = json.load(handle)
    programs = {
        "EC marker enrichment": EC_IDENTITY,
        "Stress/permeabilization": sorted(
            set(frozen_sets["mtDNA_release"] + ["EIF4EBP1"])
        ),
        "Selective clearance": frozen_sets["Mitophagy_core"],
        "EC inflammation": frozen_sets["EC_inflammation"],
    }
    all_signature_genes = set().union(*[set(v) for v in programs.values()])
    gene_index = {gene: idx for idx, gene in enumerate(gene_names)}
    mean_expression = np.asarray(lognorm.mean(axis=1)).ravel()
    rng = np.random.default_rng(args.seed)

    scores: dict[str, np.ndarray] = {}
    controls: dict[str, list[str]] = {}
    coverage_rows: list[dict] = []
    for program, requested in programs.items():
        present = [gene for gene in requested if gene in gene_index]
        detected = [
            gene
            for gene in present
            if matrix[gene_index[gene], :].getnnz() > 0
        ]
        control_genes = build_expression_matched_controls(
            mean_expression,
            gene_names,
            detected,
            all_signature_genes,
            args.control_genes_per_target,
            rng,
        )
        if not detected or not control_genes:
            raise RuntimeError(f"Insufficient genes to score {program}")
        target_rows = [gene_index[g] for g in detected]
        control_rows = [gene_index[g] for g in control_genes]
        raw_score = mean_rows(lognorm, target_rows) - mean_rows(lognorm, control_rows)
        scores[program] = zscore(raw_score)
        controls[program] = control_genes
        for gene in requested:
            idx = gene_index.get(gene)
            n_spots = int(matrix[idx, :].getnnz()) if idx is not None else 0
            coverage_rows.append(
                {
                    "program": program,
                    "gene": gene,
                    "present_in_feature_table": idx is not None,
                    "detected_in_qc_spots": n_spots,
                    "percent_qc_spots_detected": 100.0 * n_spots / matrix.shape[1],
                    "included_in_score": gene in detected,
                }
            )

    # Expression-matched controls are versioned so every score is replayable.
    control_rows = [
        {"program": program, "control_gene": gene}
        for program, genes in controls.items()
        for gene in genes
    ]
    pd.DataFrame(control_rows).to_csv(
        args.results_dir / "spatial_expression_matched_controls.csv", index=False
    )
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(args.results_dir / "spatial_module_gene_coverage.csv", index=False)

    spot_table = positions.copy()
    spot_table["total_umi"] = total_umi.astype(int)
    spot_table["detected_genes"] = detected_genes.astype(int)
    for name, values in scores.items():
        spot_table[name] = values
    ec_threshold = float(np.quantile(scores["EC marker enrichment"], 0.80))
    spot_table["ec_enriched_top_quintile"] = (
        spot_table["EC marker enrichment"] >= ec_threshold
    )
    sqstm1_idx = gene_index.get("SQSTM1")
    sqstm1 = (
        np.asarray(lognorm[sqstm1_idx, :].toarray()).ravel()
        if sqstm1_idx is not None
        else np.zeros(matrix.shape[1])
    )
    spot_table["SQSTM1_log1p_CPM10k"] = sqstm1
    spot_table.to_csv(
        args.results_dir / "GSE284089_spatial_spot_scores.csv.gz",
        index=False,
        compression="gzip",
    )

    correlation_rows = []
    ec_values = scores["EC marker enrichment"]
    ec_mask = spot_table["ec_enriched_top_quintile"].to_numpy(bool)
    for program in ["Stress/permeabilization", "Selective clearance", "EC inflammation"]:
        for subset, mask in [("all QC spots", np.ones(matrix.shape[1], bool)),
                             ("top-quintile EC-marker spots", ec_mask)]:
            rho = float(spearmanr(ec_values[mask], scores[program][mask]).statistic)
            correlation_rows.append(
                {
                    "subset": subset,
                    "program": program,
                    "n_spots": int(mask.sum()),
                    "spearman_rho_descriptive": rho,
                    "p_value": "not_tested_single_section",
                }
            )
    correlations = pd.DataFrame(correlation_rows)
    correlations.to_csv(
        args.results_dir / "spatial_within_section_correlations.csv", index=False
    )

    summary = pd.DataFrame(
        [
            {
                "accession": ACCESSION,
                "sample": SAMPLE,
                "disease_context": "osteoarthritis",
                "biological_samples": 1,
                "sections": 1,
                "spots_before_qc": int(len(keep)),
                "minimum_umi_inclusive": args.minimum_umi,
                "spots_after_qc": int(matrix.shape[1]),
                "median_umi_after_qc": float(np.median(total_umi)),
                "median_detected_genes_after_qc": float(np.median(detected_genes)),
                "ec_enriched_spots": int(ec_mask.sum()),
                "sqstm1_detected_spots": int((sqstm1 > 0).sum()),
                "sqstm1_detected_percent": float(100 * (sqstm1 > 0).mean()),
                "claim_boundary": (
                    "descriptive anatomical contextualization; not ONFH validation; "
                    "spots are not biological replicates"
                ),
            }
        ]
    )
    summary.to_csv(args.results_dir / "spatial_run_summary.csv", index=False)

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 9.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    pad = 55
    crop = (
        max(0.0, x.min() - pad),
        min(float(image.shape[1]), x.max() + pad),
        max(0.0, y.min() - pad),
        min(float(image.shape[0]), y.max() + pad),
    )
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(7.2, 5.25),
        layout="constrained",
        facecolor="white",
    )
    draw_spatial_panel(
        axes[0, 0], image, x, y, None,
        "H&E and retained spots",
        crop,
        tissue_alpha=1.0,
    )
    axes[0, 0].text(
        0.03,
        0.04,
        "GSE284089 · one section\n2,947 spots after UMI ≥100",
        transform=axes[0, 0].transAxes,
        fontsize=7.1,
        color="#1E2A35",
        ha="left",
        va="bottom",
        bbox={"boxstyle": "round,pad=0.28", "fc": "white", "ec": "#8FA3B5", "lw": 0.5, "alpha": 0.9},
    )
    draw_spatial_panel(
        axes[0, 1], image, x, y, scores["EC marker enrichment"],
        "EC marker enrichment",
        crop,
        cmap="viridis",
        norm=mpl.colors.Normalize(vmin=-2.0, vmax=2.5),
        colorbar_label="Matched-control score (z)",
    )
    draw_spatial_panel(
        axes[0, 2], image, x, y, scores["Stress/permeabilization"],
        "Stress program",
        crop,
        cmap="RdBu_r",
        norm=mpl.colors.TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5),
        colorbar_label="Matched-control score (z)",
    )
    draw_spatial_panel(
        axes[1, 0], image, x, y, scores["Selective clearance"],
        "Selective-clearance program",
        crop,
        cmap="RdBu_r",
        norm=mpl.colors.TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5),
        colorbar_label="Matched-control score (z)",
    )
    draw_spatial_panel(
        axes[1, 1], image, x, y, scores["EC inflammation"],
        "EC-inflammation program",
        crop,
        cmap="RdBu_r",
        norm=mpl.colors.TwoSlopeNorm(vmin=-2.5, vcenter=0.0, vmax=2.5),
        colorbar_label="Matched-control score (z)",
    )

    ax = axes[1, 2]
    height, width = image.shape[:2]
    ax.imshow(image, extent=(0, width, height, 0), interpolation="none", alpha=0.58)
    ax.scatter(x, y, s=9.5, c="#D0D4D8", linewidths=0, alpha=0.55, rasterized=True)
    detected_mask = sqstm1 > 0
    detected_artist = ax.scatter(
        x[detected_mask],
        y[detected_mask],
        c=sqstm1[detected_mask],
        s=12,
        cmap="Blues",
        norm=mpl.colors.Normalize(vmin=0, vmax=max(1.0, float(np.quantile(sqstm1[detected_mask], 0.98)))),
        linewidths=0.18,
        edgecolors="#1C3C5C",
        alpha=0.94,
        rasterized=True,
    )
    ax.set_xlim(crop[0], crop[1])
    ax.set_ylim(crop[3], crop[2])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        f"SQSTM1 transcript detection\n{detected_mask.sum():,}/{matrix.shape[1]:,} QC spots",
        fontsize=8.9,
        fontweight="bold",
        pad=5,
    )
    bar = fig.colorbar(detected_artist, ax=ax, fraction=0.035, pad=0.018)
    bar.ax.tick_params(labelsize=6.7, length=2)
    bar.set_label("log1p CPM10k", fontsize=7.2)
    bar.outline.set_linewidth(0.4)

    for label, ax in zip("ABCDEF", axes.ravel()):
        add_panel_label(ax, label)
    fig.suptitle(
        "Spatial contextualization in one osteoarthritic femoral-head section",
        fontsize=10.7,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.002,
        "Single-section descriptive maps; spots are not biological replicates and do not provide ONFH validation.",
        ha="center",
        va="bottom",
        fontsize=7.2,
        color="#4D5965",
    )
    figure_base = args.figure_dir / "fig5_spatial_contextualization"
    fig.savefig(figure_base.with_suffix(".pdf"), dpi=600, facecolor="white")
    fig.savefig(
        figure_base.with_suffix(".png"),
        dpi=400,
        facecolor="white",
        pil_kwargs={"compress_level": 6},
    )
    fig.savefig(figure_base.with_suffix(".svg"), facecolor="white")
    final_base = args.final_figure_dir / "Figure5"
    fig.savefig(final_base.with_suffix(".pdf"), dpi=600, facecolor="white")
    fig.savefig(
        final_base.with_suffix(".png"),
        dpi=400,
        facecolor="white",
        pil_kwargs={"compress_level": 6},
    )
    plt.close(fig)

    provenance = {
        "analysis": "GSE284089 external femoral-head spatial contextualization",
        "accession": ACCESSION,
        "sample": SAMPLE,
        "source_url": SOURCE_URL,
        "source_archive_sha256": sha256(files["archive"]),
        "input_sha256": {key: sha256(path) for key, path in files.items()},
        "minimum_umi_inclusive": args.minimum_umi,
        "normalization": "per-spot CPM to 10,000 followed by log1p",
        "duplicate_feature_handling": (
            f"summed {duplicate_symbol_rows} feature-table rows carrying duplicated "
            "display symbols before normalization"
        ),
        "score": (
            "mean target log1p CPM10k minus mean expression-bin-matched control "
            "log1p CPM10k; standardized across retained spots"
        ),
        "control_genes_per_target": args.control_genes_per_target,
        "random_seed": args.seed,
        "ec_enriched_definition": "top quintile of EC marker enrichment score",
        "statistical_inference": "none; within-section correlations are descriptive",
        "image_processing": "unaltered deposited high-resolution H&E image; display crop only",
        "claim_boundary": (
            "one osteoarthritic femoral-head section used as an anatomical spatial "
            "scaffold; not disease-matched validation or independent replication"
        ),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": mpl.__version__,
        },
    }
    with (args.results_dir / "spatial_provenance.json").open("w", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2, ensure_ascii=False)
    with (args.results_dir / "RUN_SUMMARY.txt").open("w", encoding="utf-8") as handle:
        handle.write(summary.to_string(index=False) + "\n\n")
        handle.write("Descriptive correlations (no inferential p values):\n")
        handle.write(correlations.to_string(index=False) + "\n")

    print(summary.to_string(index=False))
    print("\nDescriptive correlations:\n", correlations.to_string(index=False))
    print(f"\nFigure written to {figure_base.with_suffix('.pdf')}")
    print(f"Final figure written to {final_base.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
