#!/usr/bin/env python
"""Static and release-integrity checks for the public repository."""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import math
import re
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CODE_SUFFIXES = {".py", ".R", ".r", ".ps1", ".md"}
FORBIDDEN = (
    re.compile(r"[A-Za-z]:[\\/](?:Users|onfh-pilot|ONFH2026dry)[\\/]", re.I),
    re.compile(r"/home/[^/]+/"),
)

REQUIRED = [
    "README.md",
    "analysis/sample_metadata_v4.csv",
    "analysis/v4_common.py",
    "analysis/prepare_matrices.R",
    "virtual_knockout/run_official_vko.R",
    "virtual_knockout/run_matched_control_vko.R",
    "virtual_knockout/postprocess_mt_exclusion.R",
    "virtual_knockout/vko_selected_features_v5.csv",
    "plotting/assemble_manuscript_figures.py",
    "plotting/make_reviewed_figures.py",
    "plotting/make_figure4.py",
    "plotting/make_virtual_knockout_figure.R",
    "results/figure_inputs/module_scores_liao_stats_v4.csv",
    "results/figure_inputs/gsea_ONFH3A_vs_HOA_H.csv",
    "results/figure_inputs/sample_level_tf_ulm_v4.csv",
    "results/figure_inputs/tf_stats_ulm_v4.csv",
    "results/figure_inputs/comm_scores_v4_long.csv.gz",
    "results/figure_inputs/comm_top_ONFH_4_vs_ONFH_3A_v4.csv",
    "results/figure_inputs/diag_summary_v4.json",
    "results/official_r_vko_manuscript_no_mt_encoded/vko_provenance_official_r.json",
    "results/official_r_vko_official_default_no_mt_encoded/vko_provenance_official_r.json",
    "results/official_r_vko_no_mt_gene_comparison.csv",
    "results/official_r_vko_no_mt_summary.csv",
    "results/official_r_vko_no_mt_cross_donor_audit.csv",
    "results/official_r_vko_no_mt_hoa3_fdr_audit.csv",
    "results/official_r_vko_no_mt_pathway_by_donor.csv",
    "results/official_r_vko_no_mt_pathway_summary.csv",
    "results/vko_mt_encoded_exclusion_manifest.csv",
    "results/official_r_vko_matched_controls/matched_control_selection.csv",
    "results/official_r_vko_matched_controls/matched_control_vko_summary.csv",
    "results/official_r_vko_matched_controls/sqstm1_matched_control_calibration.csv",
    "results/official_r_vko_matched_controls/frozen_wt_reuse_validation.csv",
    "results/official_r_vko_matched_controls/sqstm1_selection_rationale.csv",
    "workflow/run_core_analysis.ps1",
    "workflow/run_virtual_knockout.ps1",
    "workflow/run_figures.ps1",
    "results/Supplementary_Tables_S1-S10.xlsx",
]

DATA_REQUIRED = [
    *[f"data/liao2022/{sample}.rds" for sample in (
        "onfh1", "onfh2", "onfh3", "onfh4", "onfh5", "onfh6",
        "hoa1", "hoa2", "hoa3", "fnf1", "fnf2"
    )],
    "data/GSE123568_series_matrix.txt.gz",
    "data/GSE123568_family.soft.gz",
    "data/cellchatdb_interactions.csv",
    "data/cellchatdb_complex_named.csv",
    "data/dorothea_ABC.tsv",
]


def fail(errors: list[str], message: str) -> None:
    errors.append(message)
    print(f"FAIL  {message}")


def pass_(message: str) -> None:
    print(f"PASS  {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-data", action="store_true", help="also require downloaded primary data")
    parser.add_argument(
        "--rscript",
        type=Path,
        help="optional Rscript executable used to parse every versioned R file",
    )
    parser.add_argument(
        "--submission-dir",
        type=Path,
        help="optional submission figure directory whose PDFs must byte-match figures/final",
    )
    args = parser.parse_args()
    errors: list[str] = []

    missing = [name for name in REQUIRED if not (ROOT / name).exists()]
    if missing:
        fail(errors, "missing release files: " + ", ".join(missing))
    else:
        pass_(f"{len(REQUIRED)} required release files")

    python_files = sorted(
        path
        for directory in ("analysis", "virtual_knockout", "plotting", "qa")
        for path in (ROOT / directory).rglob("*.py")
    )
    syntax_errors = []
    for path in python_files:
        try:
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (SyntaxError, UnicodeDecodeError) as exc:
            syntax_errors.append(f"{path.relative_to(ROOT)}: {exc}")
    if syntax_errors:
        fail(errors, "Python syntax: " + " | ".join(syntax_errors))
    else:
        pass_(f"Python syntax ({len(python_files)} files)")

    if args.rscript:
        r_files = sorted(
            path
            for directory in ("analysis", "virtual_knockout", "plotting", "environment")
            for path in (ROOT / directory).rglob("*.R")
        )
        r_syntax_errors = []
        for path in r_files:
            r_path = str(path).replace("\\", "/")
            completed = subprocess.run(
                [str(args.rscript), "-e", f"invisible(parse(file={r_path!r}))"],
                cwd=ROOT,
                capture_output=True,
                text=True,
                check=False,
            )
            if completed.returncode != 0:
                detail = (completed.stderr or completed.stdout).strip().replace("\n", " | ")
                r_syntax_errors.append(f"{path.relative_to(ROOT)}: {detail}")
        if r_syntax_errors:
            fail(errors, "R syntax: " + " | ".join(r_syntax_errors))
        else:
            pass_(f"R syntax ({len(r_files)} files)")

    local_paths = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix not in CODE_SUFFIXES or ".git" in path.parts:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for line_number, line in enumerate(text.splitlines(), 1):
            if any(pattern.search(line) for pattern in FORBIDDEN):
                local_paths.append(f"{path.relative_to(ROOT)}:{line_number}")
    if local_paths:
        fail(errors, "machine-specific absolute paths: " + ", ".join(local_paths))
    else:
        pass_("no machine-specific absolute paths")

    with (ROOT / "analysis/sample_metadata_v4.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    participants = {row["participant_id"] for row in rows if row["participant_id"]}
    sonfh = [row for row in rows if row["dataset"] == "sonfh_cystic"]
    metadata_ok = (
        len(rows) == 19
        and len(participants) == 15
        and len(sonfh) == 4
        and all(row["independent_for_inference"].strip().lower() in {"false", "0"} for row in sonfh)
    )
    if metadata_ok:
        pass_("sampling-unit guardrail: 19 libraries, 15 mapped participants, 4 descriptive SONFH libraries")
    else:
        fail(errors, "sampling-unit metadata no longer matches the frozen inferential design")

    with (ROOT / "results/figure_inputs/gsea_ONFH3A_vs_HOA_H.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        gsea = {row["pathway"]: float(row["padj"]) for row in csv.DictReader(handle)}
    expected_gsea = {
        "HALLMARK_ALLOGRAFT_REJECTION": 1.12899319297421e-17,
        "HALLMARK_INTERFERON_GAMMA_RESPONSE": 2.11625060949194e-14,
        "HALLMARK_INTERFERON_ALPHA_RESPONSE": 8.32810112872156e-10,
        "HALLMARK_INFLAMMATORY_RESPONSE": 1.15308484547667e-9,
    }
    if all(
        name in gsea and math.isclose(gsea[name], value, rel_tol=1e-12)
        for name, value in expected_gsea.items()
    ):
        pass_("Figure 4A GSEA FDR values match the frozen manuscript results")
    else:
        fail(errors, "Figure 4A GSEA table no longer matches the frozen manuscript results")

    final_expected = [f"figures/final/Figure{i}.{ext}" for i in range(1, 7) for ext in ("pdf", "png")]
    final_expected += [f"figures/final/SupplementaryFigureS1.{ext}" for ext in ("pdf", "png")]
    final_missing = [name for name in final_expected if not (ROOT / name).exists()]
    if final_missing:
        fail(errors, "missing final figures: " + ", ".join(final_missing))
    else:
        pass_("all submission figures are present as PDF and PNG")

    if args.submission_dir:
        names = [f"Figure{i}.pdf" for i in range(1, 7)]
        submission_names = {name: name for name in names}
        supplement_candidates = ("SupplementaryFigureS1.pdf", "Supplementary_Figure_S1.pdf")
        supplement_name = next(
            (name for name in supplement_candidates if (args.submission_dir / name).exists()),
            supplement_candidates[0],
        )
        submission_names["SupplementaryFigureS1.pdf"] = supplement_name
        missing_submission = [
            submitted
            for submitted in submission_names.values()
            if not (args.submission_dir / submitted).exists()
        ]
        if missing_submission:
            fail(errors, "missing PDFs in --submission-dir: " + ", ".join(missing_submission))
        else:
            mismatched = [
                canonical
                for canonical, submitted in submission_names.items()
                if sha256(ROOT / "figures" / "final" / canonical)
                != sha256(args.submission_dir / submitted)
            ]
            if mismatched:
                fail(errors, "repository/submission figure hash mismatch: " + ", ".join(mismatched))
            else:
                pass_("repository and submission figure PDFs are byte-identical")

    if args.check_data:
        data_missing = [name for name in DATA_REQUIRED if not (ROOT / name).exists()]
        gse169 = list((ROOT / "data/gse169396").glob("*")) if (ROOT / "data/gse169396").exists() else []
        gse290 = list((ROOT / "data/gse290411").glob("*")) if (ROOT / "data/gse290411").exists() else []
        if len(gse169) < 12:
            data_missing.append("data/gse169396/ (12 expected 10x files)")
        if len(gse290) < 12:
            data_missing.append("data/gse290411/ (12 expected 10x files)")
        if data_missing:
            fail(errors, "missing primary-data inputs: " + ", ".join(data_missing))
        else:
            pass_("primary-data inputs")

    if errors:
        print(f"\nRepository check failed ({len(errors)} issue groups).", file=sys.stderr)
        return 1
    print("\nRepository check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
