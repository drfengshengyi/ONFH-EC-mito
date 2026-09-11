# -*- coding: utf-8 -*-
"""Hard-stop quality checks for the isolated enhanced GSE123568 analysis."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
ROOT = DEFAULT_ROOT
OUT = ROOT / "results" / "serum_enhanced_20260910"
PLAN = OUT / "00_frozen_plan"


def configure_paths(project_root: Path, output_dir: Path | None) -> None:
    global ROOT, OUT, PLAN
    ROOT = project_root.resolve()
    OUT = (output_dir if output_dir is not None else ROOT / "results" / "serum_enhanced_20260910").resolve()
    PLAN = OUT / "00_frozen_plan"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError("CRITICAL QC FAILURE: " + message)


def check_preflight() -> list[str]:
    notes = []
    required = ["analysis_plan_frozen.json", "input_manifest.csv", "sample_manifest.csv", "candidate_gene_manifest.csv", "ma_comparator_gene_manifest.csv", "probe_annotation_audit.csv", "file_checksums_sha256.txt", "probe_selection_leakage_unit_test.json", "protected_artifact_manifest_before.csv"]
    for name in required:
        require((PLAN / name).is_file(), f"missing preflight file {name}")
    samples = pd.read_csv(PLAN / "sample_manifest.csv")
    require(len(samples) == 40 and samples.sample_id.nunique() == 40, "sample manifest must contain 40 unique samples")
    require(int(samples.label_numeric.sum()) == 30 and int((samples.label_numeric == 0).sum()) == 10, "sample class counts must be 30/10")
    require(samples.sample_id.tolist() == [f"GSM{i}" for i in range(3507251, 3507291)], "sample order mismatch")
    candidate = pd.read_csv(PLAN / "candidate_gene_manifest.csv")
    ma = pd.read_csv(PLAN / "ma_comparator_gene_manifest.csv")
    require(len(candidate) == 77 and candidate.final_available.astype(bool).all(), "candidate manifest must have 77 available genes")
    require(ma.final_available.astype(bool).sum() == 4, "Ma comparator must have exactly four available genes")
    test = json.loads((PLAN / "probe_selection_leakage_unit_test.json").read_text(encoding="utf-8"))
    require(test["status"] == "PASS", "probe-selection leakage unit test did not pass")
    manifest = pd.read_csv(PLAN / "input_manifest.csv")
    for row in manifest.itertuples(index=False):
        path = Path(row.path)
        require(path.is_file(), f"input disappeared: {path}")
        require(sha256(path) == row.sha256, f"input hash changed: {path}")
    notes.append("preflight manifests, labels, genes, hashes and leakage unit test PASS")
    return notes


def check_result(directory: Path, expected: int, *, smoke: bool) -> list[str]:
    notes = []
    required = ["outer_fold_assignments.csv", "oof_predictions_long.csv", "oof_predictions_aggregated.csv", "repeat_performance.csv", "aggregate_performance.csv", "selected_hyperparameters.csv", "feature_coefficients_long.csv", "feature_stability_summary.csv", "features_per_outer_fit.csv", "label_permutation_null.csv", "label_permutation_summary.csv", "refit_oob_bootstrap_all.csv", "refit_oob_bootstrap_invalid.csv", "refit_oob_bootstrap_summary.csv", "matched_random_space_auc_null.csv", "matched_random_space_gene_map.csv.gz", "matched_random_space_qc.csv", "matched_random_space_summary.csv"]
    for name in required:
        require((directory / name).is_file(), f"missing result {directory.name}/{name}")
    assignments = pd.read_csv(directory / "outer_fold_assignments.csv")
    require(len(assignments) == 200, "outer assignments must have 200 rows")
    require((assignments.groupby(["repeat", "sample_id"]).size() == 1).all(), "each sample must be tested once per repeat")
    oof = pd.read_csv(directory / "oof_predictions_long.csv")
    require(len(oof) == 200 and oof[["candidate_probability", "ma_probability"]].notna().all().all(), "OOF long coverage/NA failure")
    require(((oof[["candidate_probability", "ma_probability"]] >= 0) & (oof[["candidate_probability", "ma_probability"]] <= 1)).all().all(), "probability outside [0,1]")
    aggregate = pd.read_csv(directory / "oof_predictions_aggregated.csv")
    require(len(aggregate) == 40 and aggregate.sample_id.nunique() == 40, "aggregate OOF must contain 40 unique samples")
    coefficients = pd.read_csv(directory / "feature_coefficients_long.csv")
    require(len(coefficients) == 25 * 77, "candidate coefficient table must have 25x77 rows")
    require(pd.read_csv(directory / "label_permutation_null.csv").task_id.nunique() == expected, "permutation task count mismatch")
    boot = pd.read_csv(directory / "refit_oob_bootstrap_all.csv")
    require(len(boot) == expected and boot.valid.astype(bool).all(), "valid refit bootstrap count mismatch")
    require(((boot.n_oob_case >= 2) & (boot.n_oob_control >= 2)).all(), "invalid OOB class count retained")
    random = pd.read_csv(directory / "matched_random_space_auc_null.csv")
    require(random.task_id.nunique() == expected, "matched random task count mismatch")
    maps = pd.read_csv(directory / "matched_random_space_gene_map.csv.gz")
    require(len(maps) == expected * 25 * 77, "matched gene-map row count mismatch")
    uniqueness = maps.groupby(["task_id", "repeat", "outer_fold"]).matched_gene.nunique()
    require(len(uniqueness) == expected * 25 and (uniqueness == 77).all(), "every matched task/fold must have 77 unique genes")
    require(set(maps.matched_gene).isdisjoint(set(maps.candidate_gene)), "matched background contains a candidate gene")
    if smoke:
        replay = directory / "determinism_replay.json.gz"
        require(replay.is_file(), "smoke replay record missing")
        with gzip.open(replay, "rt", encoding="utf-8") as handle:
            require(bool(json.load(handle)["pass"]), "same-seed deterministic replay failed")
    else:
        extra = ["fixed_oof_paired_bootstrap_20000.csv.gz", "fixed_oof_paired_bootstrap_summary.csv", "calibration_summary.csv", "calibration_bins.csv", "brier_comparison.csv", "serum_enhanced_summary.json"]
        for name in extra:
            require((directory / name).is_file(), f"missing formal result {name}")
        fixed = pd.read_csv(directory / "fixed_oof_paired_bootstrap_20000.csv.gz")
        require(len(fixed) == 20000, "fixed-prediction paired bootstrap must have 20,000 rows")
    notes.append(f"{directory.name}: OOF, {expected} permutations, {expected} valid bootstraps and {expected} matched spaces PASS")
    return notes


def check_protected() -> str:
    before = pd.read_csv(PLAN / "protected_artifact_manifest_before.csv")
    allowed_release_replacements = {
        (ROOT / "figures" / "final" / "Figure7.pdf").resolve(): OUT / "03_figures" / "Figure7_enhanced_review.pdf",
        (ROOT / "figures" / "final" / "Figure7.png").resolve(): OUT / "03_figures" / "Figure7_enhanced_review.png",
    }
    for row in before.itertuples(index=False):
        path = Path(row.path)
        replacement = allowed_release_replacements.get(path.resolve())
        if replacement is not None and replacement.is_file() and path.is_file() and sha256(path) == sha256(replacement):
            continue
        require(path.is_file(), f"protected artifact missing: {path}")
        require(path.stat().st_size == int(row.size_bytes) and sha256(path) == row.sha256, f"protected artifact changed: {path}")
    return "protected artifacts match frozen baselines or the explicitly allowed enhanced Figure 7 release"


def check_figures() -> list[str]:
    figdir = OUT / "03_figures"
    png = figdir / "Figure7_enhanced_review.png"
    pdf = figdir / "Figure7_enhanced_review.pdf"
    require(png.is_file() and pdf.is_file(), "review Figure 7 PNG/PDF missing")
    with Image.open(png) as image:
        require(image.width >= 2200, f"PNG width is only {image.width}px")
        dpi = image.info.get("dpi", (0, 0))
        require(float(dpi[0]) >= 299, f"PNG resolution below 300 dpi: {dpi}")
    require(len(PdfReader(str(pdf)).pages) == 1, "review PDF must have exactly one page")
    require((OUT / "02_formal" / "delong_comparison.csv").is_file(), "paired DeLong output missing")
    return ["review Figure 7 PNG/PDF dimensions and paired DeLong output PASS"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["preflight", "smoke", "formal", "final"], required=True)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()
    configure_paths(args.project_root, args.output_dir)
    notes = check_preflight()
    if args.stage in {"smoke", "formal", "final"}:
        notes += check_result(OUT / "01_smoke", 20, smoke=True)
    if args.stage in {"formal", "final"}:
        notes += check_result(OUT / "02_formal", 1000, smoke=False)
    if args.stage == "final":
        notes += check_figures()
        notes.append(check_protected())
        require((OUT / "04_report" / "GSE123568_Figure7_enhanced_analysis_report_zh.md").is_file(), "Chinese analysis report missing")
        require((OUT / "04_report" / "output_manifest.csv").is_file(), "output manifest missing")
    print("QC PASS")
    for note in notes:
        print("-", note)


if __name__ == "__main__":
    main()
