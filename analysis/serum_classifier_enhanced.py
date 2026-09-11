# -*- coding: utf-8 -*-
"""Leakage-controlled enhanced GSE123568 serum analysis for Figure 7 review.

This is an exploratory internal evaluation.  It deliberately writes only to
results/serum_enhanced_20260910 and never updates manuscript-facing artifacts.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import traceback
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import joblib
import numpy as np
import pandas as pd
import scipy
import sklearn
from joblib import Parallel, delayed
from scipy.special import expit, logit
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

DEFAULT_ROOT = Path(__file__).resolve().parents[1]
ROOT = DEFAULT_ROOT
DATA = ROOT / "data"
OUT = ROOT / "results" / "serum_enhanced_20260910"
PLAN_DIR = OUT / "00_frozen_plan"
SMOKE_DIR = OUT / "01_smoke"
FORMAL_DIR = OUT / "02_formal"
FIG_DIR = OUT / "03_figures"
REPORT_DIR = OUT / "04_report"
CHECKPOINT_DIR = OUT / "checkpoints"
LOG = REPORT_DIR / "run_log.txt"

BASE_SEED = 20260820
OUTER_REPEATS = 5
OUTER_FOLDS = 5
INNER_FOLDS = 4
CS = np.logspace(-3, 2, 16)
CHECKPOINT_EVERY = 20

CANDIDATE_MODULES = [
    "Mito_fission", "Mito_fusion", "Mitophagy_core", "mtDNA_release",
    "cGAS_STING", "YAP_mTOR",
]
CANDIDATES = """AKT1 AMBRA1 ATG5 ATG7 BAK1 BAX BBC3 BCL2 BCL2L1 BCL2L13 BID BNIP3 BNIP3L CALCOCO2 CCL5 CHUK CXCL10 DNM1 DNM1L DNM2 EIF4EBP1 FIS1 FKBP8 FUNDC1 GABARAP GABARAPL2 GDAP1 GIGYF1 GIGYF2 HIF1A IKBKB IL6 INF2 IRF3 IRF7 LATS1 LATS2 MAP1LC3B MCL1 MFF MFN1 MFN2 MIEF1 MIEF2 MST1 MTCH2 MTFR1 MTFR1L MTOR NBR1 NFKB1 OPA1 OPTN PHB2 PINK1 PLD6 PMAIP1 RELA RHEB RPS6KB1 RPTOR SQSTM1 TAX1BP1 TBK1 TEAD1 TEAD4 TMEM173 TNF TSC1 TSC2 ULK1 VCP VDAC1 VDAC2 VDAC3 WWTR1 YAP1""".split()
MA_ALL = "BID FTH1 LACTB PDK3 RAB5IF SOD2 SQOR".split()
MA_AVAILABLE_EXPECTED = "BID FTH1 LACTB PDK3".split()
EXPECTED_SAMPLES = [f"GSM{i}" for i in range(3507251, 3507291)]
EXPECTED_CONTROLS = set(EXPECTED_SAMPLES[:10])
EXPECTED_CASES = set(EXPECTED_SAMPLES[10:])
PROTECTED = [ROOT / "results" / "figure_inputs", ROOT / "figures" / "final"]


def configure_paths(project_root: Path, data_dir: Path | None, output_dir: Path | None) -> None:
    """Configure repository-relative paths before any analysis phase runs."""
    global ROOT, DATA, OUT, PLAN_DIR, SMOKE_DIR, FORMAL_DIR, FIG_DIR
    global REPORT_DIR, CHECKPOINT_DIR, LOG, PROTECTED

    ROOT = project_root.resolve()
    DATA = (data_dir if data_dir is not None else ROOT / "data").resolve()
    OUT = (output_dir if output_dir is not None else ROOT / "results" / "serum_enhanced_20260910").resolve()
    PLAN_DIR = OUT / "00_frozen_plan"
    SMOKE_DIR = OUT / "01_smoke"
    FORMAL_DIR = OUT / "02_formal"
    FIG_DIR = OUT / "03_figures"
    REPORT_DIR = OUT / "04_report"
    CHECKPOINT_DIR = OUT / "checkpoints"
    LOG = REPORT_DIR / "run_log.txt"
    PROTECTED = [ROOT / "results" / "figure_inputs", ROOT / "figures" / "final"]


def log_message(*parts: object) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    line = time.strftime("%Y-%m-%d %H:%M:%S") + " | " + " ".join(map(str, parts))
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json_gz(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    with gzip.open(temporary, "wt", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, separators=(",", ":"))
    os.replace(temporary, path)


def atomic_csv(path: Path, table: pd.DataFrame, *, compression=None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    table.to_csv(temporary, index=False, compression=compression)
    os.replace(temporary, path)


def quantile_bin(values: np.ndarray, n: int = 10) -> np.ndarray:
    """Deterministic approximate quantile bins in 0..n-1, tolerating ties."""
    ranks = pd.Series(values).rank(method="average", pct=True).to_numpy()
    return np.minimum(n - 1, np.floor(np.maximum(0.0, ranks - 1e-12) * n)).astype(int)


@dataclass
class DataBundle:
    expression: np.ndarray  # probes x samples
    probe_ids: np.ndarray
    raw_symbols: np.ndarray
    symbols: np.ndarray
    sample_ids: list[str]
    source_metadata: list[str]
    y: np.ndarray
    gene_to_rows: dict[str, np.ndarray]


def parse_symbol(raw: str) -> tuple[str, str]:
    raw = str(raw).strip()
    if not raw or raw in {"---", "NA", "nan"}:
        return "", "missing_symbol"
    first = raw.split(" /// ", 1)[0].split("//", 1)[0].strip()
    if not first or first == "---":
        return "", "missing_symbol"
    status = "single_symbol"
    if " /// " in raw or "//" in raw:
        status = "multi_symbol_first_valid_used"
    return first, status


def load_probe_data() -> DataBundle:
    matrix_path = DATA / "GSE123568_series_matrix.txt.gz"
    soft_path = DATA / "GSE123568_family.soft.gz"
    if not matrix_path.is_file() or not soft_path.is_file():
        raise FileNotFoundError("Required archived GSE123568 inputs are missing")

    disease_line = None
    header = None
    rows: list[list[float]] = []
    probes: list[str] = []
    in_table = False
    with gzip.open(matrix_path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_characteristics_ch1") and "disease:" in line:
                if disease_line is not None:
                    raise RuntimeError("Expected exactly one disease metadata line")
                disease_line = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                header = None
            elif in_table and header is None:
                header = [x.strip('"') for x in line.rstrip().split("\t")]
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                parts = line.rstrip().split("\t")
                probes.append(parts[0].strip('"'))
                rows.append([float(x) if x not in {"", "NA", "null"} else np.nan for x in parts[1:]])
    if disease_line is None or header is None:
        raise RuntimeError("Series matrix metadata/table could not be parsed")
    samples = header[1:]
    if samples != EXPECTED_SAMPLES:
        raise RuntimeError(f"Sample IDs/order differ from frozen expectation: {samples}")
    if len(disease_line) != 40:
        raise RuntimeError(f"Expected 40 disease metadata values, found {len(disease_line)}")
    expression = np.asarray(rows, dtype=np.float64)
    if expression.shape != (len(probes), 40):
        raise RuntimeError(f"Unexpected probe matrix shape {expression.shape}")

    annotation: dict[str, str] = {}
    in_platform = False
    platform_header = None
    symbol_column = None
    with gzip.open(soft_path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_platform = True
                platform_header = None
                continue
            if line.startswith("!platform_table_end") and in_platform:
                break
            if not in_platform:
                continue
            cells = line.rstrip("\n\r").split("\t")
            if platform_header is None:
                platform_header = cells
                symbol_column = next((i for i, x in enumerate(cells) if "Gene Symbol" in x), None)
                if symbol_column is None:
                    raise RuntimeError("Gene Symbol annotation column not found")
                continue
            if len(cells) > symbol_column:
                annotation[cells[0].strip()] = cells[symbol_column].strip()

    raw_symbols = np.asarray([annotation.get(p, "") for p in probes], dtype=object)
    parsed = [parse_symbol(x) for x in raw_symbols]
    symbols = np.asarray([x[0] for x in parsed], dtype=object)
    gene_rows: dict[str, list[int]] = {}
    for i, gene in enumerate(symbols):
        if gene:
            gene_rows.setdefault(str(gene), []).append(i)
    gene_to_rows = {g: np.asarray(v, dtype=int) for g, v in gene_rows.items()}

    y = np.asarray([1 if s in EXPECTED_CASES else 0 for s in samples], dtype=int)
    parsed_from_metadata = np.asarray([1 if "disease: SONFH" in x else 0 for x in disease_line], dtype=int)
    if not np.array_equal(y, parsed_from_metadata):
        raise RuntimeError("GEO disease metadata does not match the frozen GSM labels")
    if int(y.sum()) != 30 or int((1 - y).sum()) != 10:
        raise RuntimeError("Expected 30 SONFH and 10 steroid-exposed non-SONFH samples")
    return DataBundle(expression, np.asarray(probes), raw_symbols, symbols, samples, disease_line, y, gene_to_rows)


def select_probes(bundle: DataBundle, genes: Sequence[str], train_indices: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    selected = []
    records = []
    for gene in genes:
        rows = bundle.gene_to_rows.get(gene, np.asarray([], dtype=int))
        if len(rows) == 0:
            raise RuntimeError(f"No platform probe for frozen gene {gene}")
        means = np.nanmean(bundle.expression[np.ix_(rows, train_indices)], axis=1)
        if not np.isfinite(means).all():
            raise RuntimeError(f"Nonfinite training mean during probe selection for {gene}")
        best = float(np.max(means))
        ties = rows[np.isclose(means, best, rtol=0.0, atol=0.0)]
        chosen = int(sorted(ties, key=lambda z: str(bundle.probe_ids[z]))[0])
        selected.append(chosen)
        records.append({
            "gene": gene, "selected_probe": str(bundle.probe_ids[chosen]),
            "n_candidate_probes": int(len(rows)), "training_mean_expression": best,
            "tie_count": int(len(ties)),
            "selection_rule": "highest outer-training mean; exact ties by lexicographic probe ID",
        })
    return np.asarray(selected, dtype=int), records


def model_cv_splits(y: np.ndarray, seed: int, groups: np.ndarray | None = None) -> tuple[list[tuple[np.ndarray, np.ndarray]], int]:
    if groups is None:
        n_splits = INNER_FOLDS
        splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        splits = list(splitter.split(np.zeros(len(y)), y))
    else:
        unique_counts = [len(np.unique(groups[y == cls])) for cls in (0, 1)]
        n_splits = min(INNER_FOLDS, *unique_counts)
        if n_splits < 2:
            raise ValueError("fewer than two unique original samples in an inner-CV class")
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        splits = list(splitter.split(np.zeros(len(y)), y, groups))
        for tr, va in splits:
            if set(groups[tr]).intersection(groups[va]):
                raise RuntimeError("duplicate GSM crossed an inner-CV split")
            if len(np.unique(y[tr])) < 2 or len(np.unique(y[va])) < 2:
                raise ValueError("inner group-CV fold lacks one class")
    return splits, int(n_splits)


def fit_predict(
    x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, *,
    penalty: str, seed: int, groups: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, dict]:
    splits, n_splits = model_cv_splits(y_train, seed, groups)
    pipe = Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            penalty=penalty, solver="liblinear", class_weight="balanced",
            max_iter=5000, random_state=seed,
        )),
    ])
    search = GridSearchCV(
        pipe, {"clf__C": CS}, scoring="roc_auc", cv=splits,
        refit=True, n_jobs=1, error_score="raise", return_train_score=False,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        search.fit(x_train, y_train)
    convergence = sum(issubclass(w.category, ConvergenceWarning) for w in caught)
    fitted = search.best_estimator_.named_steps["clf"]
    coef = fitted.coef_[0].astype(float)
    prediction = search.predict_proba(x_test)[:, 1]
    if not np.isfinite(prediction).all() or np.any((prediction < 0) | (prediction > 1)):
        raise RuntimeError("Invalid probability generated")
    details = {
        "best_C": float(search.best_params_["clf__C"]),
        "best_inner_auc": float(search.best_score_), "inner_splits": n_splits,
        "convergence_warnings": int(convergence), "n_iter": int(fitted.n_iter_[0]),
        "n_nonzero": int(np.sum(np.abs(coef) > 1e-12)),
    }
    return prediction, coef, details


def outer_splits(y: np.ndarray) -> Iterable[tuple[int, int, np.ndarray, np.ndarray]]:
    for repeat in range(OUTER_REPEATS):
        split = StratifiedKFold(OUTER_FOLDS, shuffle=True, random_state=repeat)
        for fold, (train, test) in enumerate(split.split(np.zeros(len(y)), y)):
            yield repeat, fold, train, test


def metric_row(model: str, y: np.ndarray, probability: np.ndarray) -> dict:
    return {
        "model": model,
        "auc": float(roc_auc_score(y, probability)),
        "average_precision": float(average_precision_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def observed_analysis(bundle: DataBundle, destination: Path) -> dict:
    destination.mkdir(parents=True, exist_ok=True)
    n = len(bundle.y)
    probabilities = {"candidate": np.full((OUTER_REPEATS, n), np.nan), "ma": np.full((OUTER_REPEATS, n), np.nan)}
    fold_rows, assignment_rows, probe_rows, hyper_rows, coefficient_rows = [], [], [], [], []
    for repeat, fold, train, test in outer_splits(bundle.y):
        if set(train).intersection(test):
            raise RuntimeError("Outer train/test overlap")
        for idx in test:
            assignment_rows.append({"repeat": repeat, "outer_fold": fold, "sample_id": bundle.sample_ids[idx], "label": int(bundle.y[idx])})
        seed = BASE_SEED + repeat * 100 + fold
        for model, genes, penalty in (("candidate", CANDIDATES, "l1"), ("ma", MA_AVAILABLE_EXPECTED, "l2")):
            selected_rows, selection = select_probes(bundle, genes, train)
            for row in selection:
                probe_rows.append({"repeat": repeat, "outer_fold": fold, "model": model, **row})
            x_train = bundle.expression[np.ix_(selected_rows, train)].T
            x_test = bundle.expression[np.ix_(selected_rows, test)].T
            pred, coef, details = fit_predict(x_train, bundle.y[train], x_test, penalty=penalty, seed=seed)
            probabilities[model][repeat, test] = pred
            hyper_rows.append({"repeat": repeat, "outer_fold": fold, "model": model, **details})
            if model == "candidate":
                for gene, value in zip(genes, coef):
                    coefficient_rows.append({"repeat": repeat, "outer_fold": fold, "gene": gene, "coefficient": float(value), "nonzero": bool(abs(value) > 1e-12)})
        fold_rows.append({"repeat": repeat, "outer_fold": fold, "n_train": len(train), "n_test": len(test), "n_train_case": int(bundle.y[train].sum()), "n_test_case": int(bundle.y[test].sum())})

    for model in probabilities:
        if not np.isfinite(probabilities[model]).all():
            raise RuntimeError(f"Missing OOF probability for {model}")
    aggregate = {m: p.mean(axis=0) for m, p in probabilities.items()}
    long_rows, repeat_rows = [], []
    for repeat in range(OUTER_REPEATS):
        for i, sample in enumerate(bundle.sample_ids):
            long_rows.append({
                "repeat": repeat, "sample_id": sample, "label": int(bundle.y[i]),
                "candidate_probability": float(probabilities["candidate"][repeat, i]),
                "ma_probability": float(probabilities["ma"][repeat, i]),
            })
        for model in ("candidate", "ma"):
            repeat_rows.append({"repeat": repeat, **metric_row(model, bundle.y, probabilities[model][repeat])})

    aggregate_rows = [metric_row(m, bundle.y, aggregate[m]) for m in ("candidate", "ma")]
    aggregate_df = pd.DataFrame(aggregate_rows)
    candidate_metrics = aggregate_df.set_index("model").loc["candidate"]
    ma_metrics = aggregate_df.set_index("model").loc["ma"]
    comparison = pd.DataFrame([{
        "comparison": "candidate_minus_ma",
        "delta_auc": float(candidate_metrics.auc - ma_metrics.auc),
        "delta_average_precision": float(candidate_metrics.average_precision - ma_metrics.average_precision),
        "delta_brier": float(candidate_metrics.brier - ma_metrics.brier),
        "brier_direction": "negative favors candidate",
    }])

    coefficients = pd.DataFrame(coefficient_rows)
    stability_rows = []
    for gene in CANDIDATES:
        values = coefficients.loc[coefficients.gene == gene, "coefficient"].to_numpy(float)
        nz = values[np.abs(values) > 1e-12]
        positive, negative = int(np.sum(nz > 0)), int(np.sum(nz < 0))
        dominant = "not_selected" if len(nz) == 0 else ("positive" if positive >= negative else "negative")
        consistency = np.nan if len(nz) == 0 else max(positive, negative) / len(nz)
        record = {
            "gene": gene, "n_outer_fits": 25, "n_nonzero": int(len(nz)),
            "nonzero_frequency": float(len(nz) / 25), "n_positive": positive, "n_negative": negative,
            "positive_frequency_all_fits": float(positive / 25), "negative_frequency_all_fits": float(negative / 25),
            "dominant_direction": dominant, "direction_consistency_given_selected": consistency,
            "coefficient_median_all": float(np.median(values)), "coefficient_q1_all": float(np.quantile(values, .25)),
            "coefficient_q3_all": float(np.quantile(values, .75)),
            "coefficient_median_nonzero": np.nan if len(nz) == 0 else float(np.median(nz)),
            "coefficient_q1_nonzero": np.nan if len(nz) == 0 else float(np.quantile(nz, .25)),
            "coefficient_q3_nonzero": np.nan if len(nz) == 0 else float(np.quantile(nz, .75)),
            "stable_frequency_ge_0_60": bool(len(nz) / 25 >= .60),
            "stable_frequency_and_direction": bool(len(nz) / 25 >= .60 and np.isfinite(consistency) and consistency >= .80),
        }
        for repeat in range(OUTER_REPEATS):
            rv = coefficients.loc[(coefficients.gene == gene) & (coefficients.repeat == repeat), "nonzero"]
            record[f"repeat_{repeat}_nonzero_count"] = int(rv.sum())
        stability_rows.append(record)
    features_per_fit = coefficients.groupby(["repeat", "outer_fold"], as_index=False)["nonzero"].sum().rename(columns={"nonzero": "n_selected_genes"})

    atomic_csv(destination / "outer_fold_assignments.csv", pd.DataFrame(assignment_rows))
    atomic_csv(destination / "outer_fold_qc.csv", pd.DataFrame(fold_rows))
    atomic_csv(destination / "probe_selection_by_outer_fold.csv", pd.DataFrame(probe_rows))
    atomic_csv(destination / "oof_predictions_long.csv", pd.DataFrame(long_rows))
    atomic_csv(destination / "oof_predictions_aggregated.csv", pd.DataFrame({
        "sample_id": bundle.sample_ids, "label": bundle.y,
        "candidate_probability": aggregate["candidate"], "ma_probability": aggregate["ma"],
    }))
    atomic_csv(destination / "repeat_performance.csv", pd.DataFrame(repeat_rows))
    atomic_csv(destination / "aggregate_performance.csv", aggregate_df)
    atomic_csv(destination / "selected_hyperparameters.csv", pd.DataFrame(hyper_rows))
    atomic_csv(destination / "model_comparison_fixed_oof.csv", comparison)
    atomic_csv(destination / "feature_coefficients_long.csv", coefficients)
    atomic_csv(destination / "feature_stability_summary.csv", pd.DataFrame(stability_rows))
    atomic_csv(destination / "features_per_outer_fit.csv", features_per_fit)
    return {"probabilities": probabilities, "aggregate": aggregate, "aggregate_performance": aggregate_df, "comparison": comparison}


def full_candidate_oof(bundle: DataBundle, y: np.ndarray) -> float:
    probability = np.full((OUTER_REPEATS, len(y)), np.nan)
    for repeat, fold, train, test in outer_splits(y):
        selected, _ = select_probes(bundle, CANDIDATES, train)
        pred, _, _ = fit_predict(
            bundle.expression[np.ix_(selected, train)].T, y[train],
            bundle.expression[np.ix_(selected, test)].T,
            penalty="l1", seed=BASE_SEED + repeat * 100 + fold,
        )
        probability[repeat, test] = pred
    mean_probability = probability.mean(axis=0)
    return float(roc_auc_score(y, mean_probability))


def permutation_task(task_id: int, bundle: DataBundle) -> dict:
    rng = np.random.default_rng(BASE_SEED + 10000 + task_id)
    y_perm = rng.permutation(bundle.y)
    return {"task_id": task_id, "auc": full_candidate_oof(bundle, y_perm)}


def bootstrap_task(attempt_id: int, bundle: DataBundle) -> dict:
    rng = np.random.default_rng(BASE_SEED + 20000 + attempt_id)
    case = np.flatnonzero(bundle.y == 1)
    control = np.flatnonzero(bundle.y == 0)
    train = np.concatenate([rng.choice(case, 30, replace=True), rng.choice(control, 10, replace=True)])
    used = set(train.tolist())
    oob = np.asarray([i for i in range(len(bundle.y)) if i not in used], dtype=int)
    n_oob_case = int(bundle.y[oob].sum()) if len(oob) else 0
    n_oob_control = int(len(oob) - n_oob_case)
    base = {"attempt_id": attempt_id, "n_oob_case": n_oob_case, "n_oob_control": n_oob_control, "n_unique_train": len(used)}
    if n_oob_case < 2 or n_oob_control < 2:
        return {**base, "valid": False, "invalid_reason": "OOB requires at least two samples per class"}
    groups = np.asarray([bundle.sample_ids[i] for i in train])
    unique_by_class = [len(np.unique(groups[bundle.y[train] == cls])) for cls in (0, 1)]
    if min(unique_by_class) < 2:
        return {**base, "valid": False, "invalid_reason": "inner grouped CV has fewer than two unique samples in a class"}
    metrics = {}
    for model, genes, penalty in (("candidate", CANDIDATES, "l1"), ("ma", MA_AVAILABLE_EXPECTED, "l2")):
        selected, _ = select_probes(bundle, genes, train)
        pred, coef, details = fit_predict(
            bundle.expression[np.ix_(selected, train)].T, bundle.y[train],
            bundle.expression[np.ix_(selected, oob)].T,
            penalty=penalty, seed=BASE_SEED + 21000 + attempt_id, groups=groups,
        )
        metrics[model] = {**metric_row(model, bundle.y[oob], pred), "best_C": details["best_C"], "n_selected": int(np.sum(np.abs(coef) > 1e-12))}
    return {
        **base, "valid": True, "invalid_reason": "",
        "candidate_auc": metrics["candidate"]["auc"], "candidate_ap": metrics["candidate"]["average_precision"], "candidate_brier": metrics["candidate"]["brier"],
        "ma_auc": metrics["ma"]["auc"], "ma_ap": metrics["ma"]["average_precision"], "ma_brier": metrics["ma"]["brier"],
        "delta_auc": metrics["candidate"]["auc"] - metrics["ma"]["auc"],
        "delta_ap": metrics["candidate"]["average_precision"] - metrics["ma"]["average_precision"],
        "delta_brier": metrics["candidate"]["brier"] - metrics["ma"]["brier"],
        "candidate_best_C": metrics["candidate"]["best_C"], "ma_best_C": metrics["ma"]["best_C"],
        "candidate_n_selected": metrics["candidate"]["n_selected"],
    }


def all_gene_fold(bundle: DataBundle, train: np.ndarray, test: np.ndarray) -> dict:
    genes = sorted(bundle.gene_to_rows)
    selected, records = select_probes(bundle, genes, train)
    values = bundle.expression[np.ix_(selected, train)]
    means = np.mean(values, axis=1)
    variances = np.var(values, axis=1, ddof=0)
    finite = np.isfinite(means) & np.isfinite(variances) & (variances > 0)
    genes_arr = np.asarray(genes, dtype=object)
    mean_bins = np.full(len(genes), -1, dtype=int)
    var_bins = np.full(len(genes), -1, dtype=int)
    mean_bins[finite] = quantile_bin(means[finite])
    var_bins[finite] = quantile_bin(np.log1p(variances[finite]))
    return {
        "genes": genes_arr, "rows": selected, "means": means, "logvars": np.log1p(variances),
        "mean_bins": mean_bins, "var_bins": var_bins, "finite": finite,
        "train": train, "test": test,
    }


def make_gene_fold_cache(bundle: DataBundle) -> list[dict]:
    cache = []
    for repeat, fold, train, test in outer_splits(bundle.y):
        item = all_gene_fold(bundle, train, test)
        item.update({"repeat": repeat, "fold": fold})
        cache.append(item)
    return cache


def matched_mapping(task_id: int, cache: dict) -> tuple[list[int], list[dict]]:
    genes = cache["genes"]
    lookup = {g: i for i, g in enumerate(genes)}
    excluded = set(CANDIDATES).union(MA_ALL)
    eligible = cache["finite"] & np.asarray([g not in excluded for g in genes], dtype=bool)
    unused = eligible.copy()
    rng = np.random.default_rng(BASE_SEED + 30000 + task_id * 100 + cache["repeat"] * 10 + cache["fold"])
    chosen_rows, records = [], []
    for candidate in CANDIDATES:
        ci = lookup[candidate]
        cb = (int(cache["mean_bins"][ci]), int(cache["var_bins"][ci]))
        if min(cb) < 0:
            raise RuntimeError(f"Nonfinite candidate matching statistic for {candidate}")
        if not bool(unused.any()):
            raise RuntimeError("Background pool exhausted before 77 unique matches")
        # The same frozen Manhattan-expansion rule, vectorized across the
        # platform background.  genes is lexicographically sorted, so the
        # eligible tie pool has a deterministic order before seeded drawing.
        distances = np.abs(cache["mean_bins"] - cb[0]) + np.abs(cache["var_bins"] - cb[1])
        distances = np.where(unused, distances, 999)
        minimum = int(distances.min())
        pool = np.flatnonzero(distances == minimum)
        pick = int(pool[int(rng.integers(0, len(pool)))])
        unused[pick] = False
        chosen_rows.append(int(cache["rows"][pick]))
        records.append({
            "task_id": task_id, "repeat": cache["repeat"], "outer_fold": cache["fold"],
            "candidate_gene": candidate, "matched_gene": str(genes[pick]),
            "candidate_mean_bin": cb[0], "candidate_var_bin": cb[1],
            "matched_mean_bin": int(cache["mean_bins"][pick]), "matched_var_bin": int(cache["var_bins"][pick]),
            "expanded": bool(minimum > 0), "manhattan_distance": minimum,
            "training_mean_difference": float(cache["means"][pick] - cache["means"][ci]),
            "training_log1p_variance_difference": float(cache["logvars"][pick] - cache["logvars"][ci]),
        })
    if len(chosen_rows) != 77 or len(set(chosen_rows)) != 77:
        raise RuntimeError("Matched space is not 77 unique genes")
    return chosen_rows, records


def random_space_task(task_id: int, bundle: DataBundle, fold_cache: list[dict]) -> dict:
    probability = np.full((OUTER_REPEATS, len(bundle.y)), np.nan)
    maps: list[dict] = []
    for cache in fold_cache:
        rows, records = matched_mapping(task_id, cache)
        train, test = cache["train"], cache["test"]
        pred, _, _ = fit_predict(
            bundle.expression[np.ix_(rows, train)].T, bundle.y[train],
            bundle.expression[np.ix_(rows, test)].T,
            penalty="l1", seed=BASE_SEED + cache["repeat"] * 100 + cache["fold"],
        )
        probability[cache["repeat"], test] = pred
        maps.extend(records)
    aggregate = probability.mean(axis=0)
    if not np.isfinite(aggregate).all():
        raise RuntimeError("Matched random space has missing OOF predictions")
    return {
        "task_id": task_id, "auc": float(roc_auc_score(bundle.y, aggregate)),
        "average_precision": float(average_precision_score(bundle.y, aggregate)),
        "brier": float(brier_score_loss(bundle.y, aggregate)), "map": maps,
    }


def run_indexed_tasks(family: str, run_name: str, target: int, jobs: int, worker, *, resume: bool) -> list[dict]:
    directory = CHECKPOINT_DIR / run_name / family
    directory.mkdir(parents=True, exist_ok=True)
    existing = sorted(directory.glob("task_*.json.gz"))
    if existing and not resume:
        raise RuntimeError(f"Existing {family} checkpoints require --resume or a new output directory")
    completed = {int(x.stem.split("_")[-1].split(".")[0]) for x in existing}
    pending = [i for i in range(target) if i not in completed]
    log_message(run_name, family, len(completed), "complete;", len(pending), "pending")
    for start in range(0, len(pending), CHECKPOINT_EVERY):
        batch = pending[start:start + CHECKPOINT_EVERY]
        try:
            rows = Parallel(n_jobs=jobs, prefer="threads")(delayed(worker)(i) for i in batch)
        except Exception:
            log_message("TASK FAILURE", run_name, family, "batch", batch, traceback.format_exc())
            raise
        for task_id, row in zip(batch, rows):
            atomic_json_gz(directory / f"task_{task_id:05d}.json.gz", row)
        log_message(run_name, family, "checkpoint", min(start + len(batch) + len(completed), target), "/", target)
    results = []
    for path in sorted(directory.glob("task_*.json.gz")):
        task_id = int(path.stem.split("_")[-1].split(".")[0])
        if task_id < target:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                results.append(json.load(handle))
    results.sort(key=lambda x: int(x["task_id"]))
    if len(results) != target or len({int(x["task_id"]) for x in results}) != target:
        raise RuntimeError(f"{family} checkpoint coverage/uniqueness failure")
    return results


def run_bootstrap_tasks(run_name: str, target: int, jobs: int, bundle: DataBundle, *, resume: bool) -> tuple[list[dict], list[dict]]:
    directory = CHECKPOINT_DIR / run_name / "bootstrap"
    directory.mkdir(parents=True, exist_ok=True)
    existing = sorted(directory.glob("attempt_*.json.gz"))
    if existing and not resume:
        raise RuntimeError("Existing bootstrap checkpoints require --resume or a new output directory")
    attempts: dict[int, dict] = {}
    for path in existing:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            item = json.load(handle)
        attempts[int(item["attempt_id"])] = item
    next_id = 0 if not attempts else max(attempts) + 1
    while sum(bool(x["valid"]) for x in attempts.values()) < target:
        if next_id >= 5000:
            raise RuntimeError(f"Failed to obtain {target} valid bootstrap runs within 5000 attempts")
        ids = list(range(next_id, min(5000, next_id + CHECKPOINT_EVERY)))
        try:
            rows = Parallel(n_jobs=jobs, prefer="threads")(delayed(bootstrap_task)(i, bundle) for i in ids)
        except Exception:
            log_message("TASK FAILURE", run_name, "bootstrap", "attempts", ids, traceback.format_exc())
            raise
        for item in rows:
            aid = int(item["attempt_id"])
            attempts[aid] = item
            atomic_json_gz(directory / f"attempt_{aid:05d}.json.gz", item)
        next_id += len(ids)
        log_message(run_name, "bootstrap checkpoint", sum(bool(x["valid"]) for x in attempts.values()), "valid /", target, ";", len(attempts), "attempts")
    valid = sorted((x for x in attempts.values() if x["valid"]), key=lambda x: int(x["attempt_id"]))[:target]
    cutoff = int(valid[-1]["attempt_id"])
    invalid = sorted((x for x in attempts.values() if not x["valid"] and int(x["attempt_id"]) <= cutoff), key=lambda x: int(x["attempt_id"]))
    return valid, invalid


def fixed_prediction_bootstrap(y: np.ndarray, candidate: np.ndarray, ma: np.ndarray, n: int = 20000) -> pd.DataFrame:
    rng = np.random.default_rng(BASE_SEED + 71)
    cases, controls = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    rows = []
    for task_id in range(n):
        idx = np.concatenate([rng.choice(cases, len(cases), replace=True), rng.choice(controls, len(controls), replace=True)])
        rows.append({
            "task_id": task_id,
            "delta_auc": roc_auc_score(y[idx], candidate[idx]) - roc_auc_score(y[idx], ma[idx]),
            "delta_ap": average_precision_score(y[idx], candidate[idx]) - average_precision_score(y[idx], ma[idx]),
            "delta_brier": brier_score_loss(y[idx], candidate[idx]) - brier_score_loss(y[idx], ma[idx]),
        })
    return pd.DataFrame(rows)


def calibration_outputs(y: np.ndarray, predictions: dict[str, np.ndarray], destination: Path) -> None:
    summaries, bins = [], []
    baseline = float(brier_score_loss(y, np.full(len(y), y.mean())))
    for model, probability in predictions.items():
        brier = float(brier_score_loss(y, probability))
        clipped = np.clip(probability, 1e-6, 1 - 1e-6)
        z = logit(clipped).reshape(-1, 1)
        status = "ok"
        intercept = slope = np.nan
        try:
            calibration = LogisticRegression(penalty=None, solver="lbfgs", max_iter=5000).fit(z, y)
            intercept, slope = float(calibration.intercept_[0]), float(calibration.coef_[0, 0])
        except Exception as exc:
            status = f"failed: {type(exc).__name__}: {exc}"
        summaries.append({
            "model": model, "brier": brier, "prevalence_baseline_brier": baseline,
            "brier_skill_score": float(1 - brier / baseline),
            "calibration_intercept_exploratory": intercept, "calibration_slope_exploratory": slope,
            "calibration_status": status,
        })
        order = np.argsort(probability, kind="stable")
        for bin_id, idx in enumerate(np.array_split(order, 4), 1):
            bins.append({"model": model, "bin": bin_id, "n": len(idx), "mean_predicted": float(np.mean(probability[idx])), "observed_fraction": float(np.mean(y[idx]))})
    atomic_csv(destination / "calibration_summary.csv", pd.DataFrame(summaries))
    atomic_csv(destination / "calibration_bins.csv", pd.DataFrame(bins))
    atomic_csv(destination / "brier_comparison.csv", pd.DataFrame(summaries)[["model", "brier", "prevalence_baseline_brier", "brier_skill_score"]])


def write_preflight(bundle: DataBundle) -> None:
    for directory in (PLAN_DIR, SMOKE_DIR, FORMAL_DIR, FIG_DIR, REPORT_DIR, CHECKPOINT_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    inputs = [DATA / "GSE123568_series_matrix.txt.gz", DATA / "GSE123568_family.soft.gz", ROOT / "analysis" / "genesets_final.json"]
    manifest = []
    for path in inputs:
        stat = path.stat()
        manifest.append({"path": str(path), "size_bytes": stat.st_size, "modified_time": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(stat.st_mtime)), "sha256": sha256(path), "parse_status": "parsed_ok"})
    atomic_csv(PLAN_DIR / "input_manifest.csv", pd.DataFrame(manifest))

    values = bundle.expression
    quantiles = np.nanquantile(values, [0, .01, .25, .5, .75, .99, 1])
    scale = pd.DataFrame({"quantile": ["min", "q01", "q25", "median", "q75", "q99", "max"], "value": quantiles})
    scale["missing_values"] = int(np.isnan(values).sum())
    scale["interpretation"] = "GEO archived processed log2-scale expression; no additional log transformation applied"
    atomic_csv(PLAN_DIR / "expression_scale_qc.csv", scale)
    if np.isnan(values).any() or quantiles[-1] > 30 or quantiles[0] < -5:
        raise RuntimeError(f"Expression-scale critical QC failed: {dict(zip(scale.quantile, quantiles))}")

    sample_manifest = pd.DataFrame({
        "sample_id": bundle.sample_ids, "label_numeric": bundle.y,
        "label_text": np.where(bundle.y == 1, "SONFH", "steroid-exposed non-SONFH"),
        "source_metadata_text": bundle.source_metadata, "included": True, "exclusion_reason": "",
    })
    atomic_csv(PLAN_DIR / "sample_manifest.csv", sample_manifest)

    with (ROOT / "analysis" / "genesets_final.json").open(encoding="utf-8") as handle:
        genesets = json.load(handle)
    frozen_union = sorted({g for module in CANDIDATE_MODULES for g in genesets[module]})
    if any(g not in frozen_union for g in CANDIDATES):
        raise RuntimeError("Frozen 77-gene list contains gene outside prespecified six modules")
    module_map = {g: ";".join(m for m in CANDIDATE_MODULES if g in genesets[m]) for g in CANDIDATES}
    candidate_manifest = pd.DataFrame([{
        "gene": g, "source_module": module_map[g], "platform_has_probe": g in bundle.gene_to_rows,
        "n_probes": len(bundle.gene_to_rows.get(g, [])), "final_available": g in bundle.gene_to_rows,
        "unavailable_reason": "" if g in bundle.gene_to_rows else "no annotated platform probe",
    } for g in CANDIDATES])
    ma_manifest = pd.DataFrame([{
        "gene": g, "source": "Ma et al. original seven-gene comparator", "platform_has_probe": g in bundle.gene_to_rows,
        "n_probes": len(bundle.gene_to_rows.get(g, [])), "final_available": g in MA_AVAILABLE_EXPECTED and g in bundle.gene_to_rows,
        "unavailable_reason": "" if g in MA_AVAILABLE_EXPECTED and g in bundle.gene_to_rows else "not measurable under frozen original annotation rule",
    } for g in MA_ALL])
    atomic_csv(PLAN_DIR / "candidate_gene_manifest.csv", candidate_manifest)
    atomic_csv(PLAN_DIR / "ma_comparator_gene_manifest.csv", ma_manifest)
    if int(candidate_manifest.final_available.sum()) != 77:
        raise RuntimeError("Candidate-gene availability is not exactly 77")
    measured_ma = ma_manifest.loc[ma_manifest.platform_has_probe, "gene"].tolist()
    if measured_ma != MA_AVAILABLE_EXPECTED:
        raise RuntimeError(f"Ma measurable genes differ from frozen four: {measured_ma}")

    counts = pd.Series(bundle.symbols[bundle.symbols != ""]).value_counts()
    audit = []
    for probe, raw, gene in zip(bundle.probe_ids, bundle.raw_symbols, bundle.symbols):
        parsed_gene, status = parse_symbol(str(raw))
        audit.append({"probe_id": probe, "raw_gene_symbol": raw, "parsed_gene_symbol": parsed_gene, "mapping_status": status, "n_probes_for_gene": int(counts.get(gene, 0)) if gene else 0})
    atomic_csv(PLAN_DIR / "probe_annotation_audit.csv", pd.DataFrame(audit))

    plan = {
        "plan_version": "GSE123568_Figure7_enhanced_20260910_v1", "created_before_formal_models": True,
        "samples": [{"sample_id": s, "label": int(y)} for s, y in zip(bundle.sample_ids, bundle.y)],
        "candidate_genes": CANDIDATES, "candidate_modules": CANDIDATE_MODULES,
        "ma_original_genes": MA_ALL, "ma_measurable_genes": MA_AVAILABLE_EXPECTED,
        "base_seed": BASE_SEED, "outer": "StratifiedKFold 5 folds x 5 repeats; random_state 0..4",
        "inner": "StratifiedKFold 4 folds; bootstrap uses StratifiedGroupKFold min(4, unique class groups)",
        "C_grid": CS.tolist(), "candidate_model": "StandardScaler pipeline + L1 liblinear logistic regression",
        "ma_model": "StandardScaler pipeline + L2 liblinear logistic regression", "class_weight": "balanced", "max_iter": 5000,
        "formal_counts": {"label_permutations": 1000, "valid_refit_oob_bootstraps": 1000, "matched_random_spaces": 1000, "fixed_prediction_bootstraps": 20000},
        "smoke_counts": {"label_permutations": 20, "valid_refit_oob_bootstraps": 20, "matched_random_spaces": 20},
        "jobs": 4, "checkpoint_interval": 20,
        "primary_metric": "ROC-AUC of per-sample probabilities averaged across five repeated nested outer CV runs",
        "secondary_metrics": ["average precision", "Brier score", "repeat-specific metrics", "paired candidate-minus-Ma differences"],
        "empirical_p_formula": "(1 + number(null >= observed))/(1 + number(valid null))",
        "bootstrap_validity": "OOB contains at least 2 SONFH and 2 controls; inner group CV has at least 2 unique GSM per class; maximum 5000 attempts",
        "random_space_matching": "outer-training-only mean and log1p variance deciles; no replacement; Manhattan expansion",
        "figure_panels": ["A label permutation", "B matched random spaces", "C feature stability", "D paired model metrics", "E ROC and exploratory calibration"],
        "interpretation_boundary": "exploratory internal cross-fitted evaluation; no external validation, locked signature, clinical threshold, diagnostic or causal claim",
        "stop_conditions": ["input/sample/gene mismatch", "uncertain expression scale", "probe leakage unit test failure", "smoke QC failure", "fewer than requested valid bootstrap runs", "incomplete matched 77-gene spaces"],
    }
    plan_path = PLAN_DIR / "analysis_plan_frozen.json"
    serialized = json.dumps(plan, indent=2, ensure_ascii=False) + "\n"
    if plan_path.exists() and plan_path.read_text(encoding="utf-8") != serialized:
        raise RuntimeError("Frozen plan exists but differs; use a new versioned output directory")
    if not plan_path.exists():
        plan_path.write_text(serialized, encoding="utf-8")

    # Leakage unit test: altering only held-out samples must not change selected probes.
    _, _, train, test = next(iter(outer_splits(bundle.y)))
    before, _ = select_probes(bundle, CANDIDATES, train)
    modified = DataBundle(bundle.expression.copy(), bundle.probe_ids, bundle.raw_symbols, bundle.symbols, bundle.sample_ids, bundle.source_metadata, bundle.y, bundle.gene_to_rows)
    modified.expression[:, test] += np.arange(len(test), dtype=float)[None, :] * 1000 + 777
    after, _ = select_probes(modified, CANDIDATES, train)
    if not np.array_equal(before, after):
        raise RuntimeError("Critical leakage unit test failed")
    (PLAN_DIR / "probe_selection_leakage_unit_test.json").write_text(json.dumps({"status": "PASS", "test": "altered held-out expression does not change training-fold probe IDs", "n_genes": 77}, indent=2), encoding="utf-8")

    checksums = []
    for path in inputs + [plan_path, ROOT / "analysis" / "serum_classifier_enhanced.py"]:
        checksums.append(f"{sha256(path)}  {path}")
    (PLAN_DIR / "file_checksums_sha256.txt").write_text("\n".join(checksums) + "\n", encoding="utf-8")
    protected = []
    for directory in PROTECTED:
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                protected.append({"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_csv(PLAN_DIR / "protected_artifact_manifest_before.csv", pd.DataFrame(protected))
    log_message("PREFLIGHT PASS", "matrix", bundle.expression.shape, "range", float(np.min(values)), float(np.max(values)), "candidate", 77, "Ma", 4)


def run_analysis(bundle: DataBundle, run_name: str, count: int, jobs: int, *, resume: bool) -> None:
    destination = SMOKE_DIR if run_name == "smoke" else FORMAL_DIR
    destination.mkdir(parents=True, exist_ok=True)
    log_message(run_name, "observed nested CV start")
    observed = observed_analysis(bundle, destination)
    observed_auc = float(observed["aggregate_performance"].set_index("model").loc["candidate", "auc"])

    permutations = run_indexed_tasks("permutation", run_name, count, jobs, lambda i: permutation_task(i, bundle), resume=resume)
    permutation_df = pd.DataFrame(permutations)
    atomic_csv(destination / "label_permutation_null.csv", permutation_df)
    p_perm = (1 + int(np.sum(permutation_df.auc >= observed_auc))) / (1 + len(permutation_df))
    atomic_csv(destination / "label_permutation_summary.csv", pd.DataFrame([{"observed_auc": observed_auc, "n_valid": len(permutation_df), "n_ge_observed": int(np.sum(permutation_df.auc >= observed_auc)), "empirical_p": p_perm}]))

    valid, invalid = run_bootstrap_tasks(run_name, count, jobs, bundle, resume=resume)
    valid_df, invalid_df = pd.DataFrame(valid), pd.DataFrame(invalid)
    atomic_csv(destination / "refit_oob_bootstrap_all.csv", valid_df)
    invalid_columns = ["attempt_id", "n_oob_case", "n_oob_control", "n_unique_train", "valid", "invalid_reason"]
    atomic_csv(destination / "refit_oob_bootstrap_invalid.csv", invalid_df if len(invalid_df) else pd.DataFrame(columns=invalid_columns))
    summary_rows = []
    for metric in ["candidate_auc", "candidate_ap", "candidate_brier", "ma_auc", "ma_ap", "ma_brier", "delta_auc", "delta_ap", "delta_brier"]:
        summary_rows.append({"metric": metric, "n_valid": len(valid_df), "median": float(valid_df[metric].median()), "q025": float(valid_df[metric].quantile(.025)), "q975": float(valid_df[metric].quantile(.975)), "interval_name": "refit-aware OOB bootstrap sensitivity interval"})
    atomic_csv(destination / "refit_oob_bootstrap_summary.csv", pd.DataFrame(summary_rows))

    fold_cache = make_gene_fold_cache(bundle)
    random_results = run_indexed_tasks("matched_random", run_name, count, jobs, lambda i: random_space_task(i, bundle, fold_cache), resume=resume)
    random_metrics = pd.DataFrame([{k: v for k, v in item.items() if k != "map"} for item in random_results])
    maps = pd.DataFrame([row for item in random_results for row in item["map"]])
    atomic_csv(destination / "matched_random_space_auc_null.csv", random_metrics)
    atomic_csv(destination / "matched_random_space_gene_map.csv.gz", maps, compression="gzip")
    qc = maps.groupby("task_id").agg(n_rows=("matched_gene", "size"), n_unique_task_fold=("matched_gene", lambda x: int(x.groupby([maps.loc[x.index, "repeat"], maps.loc[x.index, "outer_fold"]]).nunique().min())), expanded_fraction=("expanded", "mean"), max_distance=("manhattan_distance", "max"), mean_abs_expression_difference=("training_mean_difference", lambda x: float(np.mean(np.abs(x)))), mean_abs_logvariance_difference=("training_log1p_variance_difference", lambda x: float(np.mean(np.abs(x))))).reset_index()
    qc["complete_25_folds_x_77"] = qc.n_rows == 25 * 77
    qc["all_folds_have_77_unique"] = qc.n_unique_task_fold == 77
    atomic_csv(destination / "matched_random_space_qc.csv", qc)
    if not bool((qc.complete_25_folds_x_77 & qc.all_folds_have_77_unique).all()):
        raise RuntimeError("Matched-space critical QC failed")
    p_specificity = (1 + int(np.sum(random_metrics.auc >= observed_auc))) / (1 + len(random_metrics))
    atomic_csv(destination / "matched_random_space_summary.csv", pd.DataFrame([{"observed_candidate_auc": observed_auc, "n_matched_spaces": len(random_metrics), "n_ge_observed": int(np.sum(random_metrics.auc >= observed_auc)), "empirical_p_specificity": p_specificity, "matching_scope": "outer-training-fold probe selection and mean/log1p-variance deciles"}]))

    if run_name == "smoke":
        replay_perm = permutation_task(0, bundle)
        replay_boot = bootstrap_task(int(valid[0]["attempt_id"]), bundle)
        replay_random = random_space_task(0, bundle, fold_cache)
        deterministic = replay_perm == permutations[0] and replay_boot == valid[0] and replay_random == random_results[0]
        atomic_json_gz(destination / "determinism_replay.json.gz", {"pass": deterministic, "permutation": replay_perm, "bootstrap": replay_boot, "random": replay_random})
        if not deterministic:
            raise RuntimeError("Smoke determinism replay failed")
    else:
        aggregated = pd.read_csv(destination / "oof_predictions_aggregated.csv")
        candidate_probability = aggregated.candidate_probability.to_numpy(float)
        ma_probability = aggregated.ma_probability.to_numpy(float)
        fixed = fixed_prediction_bootstrap(bundle.y, candidate_probability, ma_probability, 20000)
        atomic_csv(destination / "fixed_oof_paired_bootstrap_20000.csv.gz", fixed, compression="gzip")
        fixed_summary = []
        for metric in ["delta_auc", "delta_ap", "delta_brier"]:
            fixed_summary.append({"metric": metric, "median": float(fixed[metric].median()), "q025": float(fixed[metric].quantile(.025)), "q975": float(fixed[metric].quantile(.975)), "two_sided_zero_p": float(min(1, 2 * min((1 + (fixed[metric] <= 0).sum()) / (len(fixed) + 1), (1 + (fixed[metric] >= 0).sum()) / (len(fixed) + 1)))), "scope": "paired stratified bootstrap of fixed aggregate cross-fitted predictions"})
        atomic_csv(destination / "fixed_oof_paired_bootstrap_summary.csv", pd.DataFrame(fixed_summary))
        calibration_outputs(bundle.y, {"candidate": candidate_probability, "ma": ma_probability}, destination)
        summary = {
            "scope": "exploratory internal cross-fitted evaluation; not external validation",
            "candidate": observed["aggregate_performance"].set_index("model").loc["candidate"].to_dict(),
            "ma": observed["aggregate_performance"].set_index("model").loc["ma"].to_dict(),
            "comparison": observed["comparison"].iloc[0].to_dict(),
            "label_permutation": {"n": len(permutation_df), "empirical_p": p_perm},
            "matched_random_spaces": {"n": len(random_metrics), "empirical_p_specificity": p_specificity},
            "refit_oob_bootstrap": {"n_valid": len(valid_df), "n_invalid_before_target": len(invalid_df)},
        }
        (destination / "serum_enhanced_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float) + "\n", encoding="utf-8")
    log_message(run_name.upper(), "ANALYSIS COMPLETE", "observed AUC", observed_auc, "permutation p", p_perm, "specificity p", p_specificity)


def write_environment() -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    details = [
        f"timestamp={time.strftime('%Y-%m-%dT%H:%M:%S')}", f"command={' '.join(sys.argv)}",
        f"executable={sys.executable}", f"python={sys.version}", f"platform={platform.platform()}",
        f"cpu_count={os.cpu_count()}", "jobs_requested=4", f"numpy={np.__version__}",
        f"pandas={pd.__version__}", f"scipy={scipy.__version__}", f"scikit-learn={sklearn.__version__}", f"joblib={joblib.__version__}",
    ]
    try:
        frozen = subprocess.run([sys.executable, "-m", "pip", "freeze"], check=True, capture_output=True, text=True).stdout
        details.extend(["", "pip_freeze:", frozen.rstrip()])
    except Exception as exc:
        details.extend(["", f"pip_freeze_failed={type(exc).__name__}: {exc}"])
    (REPORT_DIR / "python_environment.txt").write_text("\n".join(details) + "\n", encoding="utf-8")


def finalize_outputs() -> None:
    summary_path = FORMAL_DIR / "serum_enhanced_summary.json"
    delong_path = FORMAL_DIR / "delong_comparison.csv"
    if not summary_path.is_file() or not delong_path.is_file():
        raise RuntimeError("Formal summary and R paired DeLong output are required before finalization")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    delong = pd.read_csv(delong_path).iloc[0]
    summary["paired_delong_p"] = float(delong.paired_delong_p)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    legacy = {
        "candidate_auc": 0.870, "candidate_average_precision": 0.927,
        "ma_auc": 0.813, "delta_auc": 0.057, "paired_delong_p": 0.192,
    }
    candidate = summary["candidate"]
    ma = summary["ma"]
    comparison = summary["comparison"]
    enhanced = {
        "candidate_auc": candidate["auc"], "candidate_average_precision": candidate["average_precision"],
        "ma_auc": ma["auc"], "delta_auc": comparison["delta_auc"], "paired_delong_p": float(delong.paired_delong_p),
    }
    reasons = {
        "candidate_auc": "outer-training-fold probe collapse and inner-fold scaler fitting",
        "candidate_average_precision": "outer-training-fold probe collapse and inner-fold scaler fitting",
        "ma_auc": "outer-training-fold probe collapse and inner-fold scaler fitting",
        "delta_auc": "both candidate and comparator were rebuilt under the stricter paired pipeline",
        "paired_delong_p": "paired DeLong applied to the newly aggregated cross-fitted predictions",
    }
    rows = []
    for metric in legacy:
        old, new = float(legacy[metric]), float(enhanced[metric])
        rows.append({
            "metric": metric, "legacy_value": old, "enhanced_value": new,
            "absolute_difference": new - old,
            "relative_difference": np.nan if old == 0 else (new - old) / abs(old),
            "plausible_pipeline_source": reasons[metric],
            "changes_cautious_conclusion": False,
        })
    atomic_csv(FORMAL_DIR / "legacy_vs_enhanced_comparison.csv", pd.DataFrame(rows))

    perm = pd.read_csv(FORMAL_DIR / "label_permutation_summary.csv").iloc[0]
    matched = pd.read_csv(FORMAL_DIR / "matched_random_space_summary.csv").iloc[0]
    refit = pd.read_csv(FORMAL_DIR / "refit_oob_bootstrap_summary.csv").set_index("metric")
    fixed = pd.read_csv(FORMAL_DIR / "fixed_oof_paired_bootstrap_summary.csv").set_index("metric")
    stability = pd.read_csv(FORMAL_DIR / "feature_stability_summary.csv")
    calibration = pd.read_csv(FORMAL_DIR / "calibration_summary.csv").set_index("model")
    invalid = pd.read_csv(FORMAL_DIR / "refit_oob_bootstrap_invalid.csv")
    input_manifest = pd.read_csv(PLAN_DIR / "input_manifest.csv")
    top = stability.sort_values(["nonzero_frequency", "direction_consistency_given_selected"], ascending=[False, False]).head(10)

    if float(perm.empirical_p) < .05:
        perm_text = "候选流程在本数据集内部优于完整标签随机化零分布；这不是外部验证。"
    else:
        perm_text = "完整流程标签置换未显示候选流程超出随机标签零分布。"
    if float(matched.empirical_p_specificity) < .05:
        match_text = "候选空间在本内部分析中优于表达量/方差匹配的同规模随机空间；这不等于外部特异性验证。"
    else:
        match_text = "候选空间未显示优于表达量/方差匹配的同规模随机空间，因此不能主张候选空间特异性。"
    fixed_auc = fixed.loc["delta_auc"]
    refit_auc = refit.loc["delta_auc"]
    superior = float(delong.paired_delong_p) < .05 and float(fixed_auc.q025) > 0
    compare_text = "候选模型在固定 OOF 配对比较中有统计优势。" if superior else "候选模型相对 Ma 比较器没有得到稳健的统计优势证据。"

    inputs_md = "\n".join(f"- `{row.path}`，SHA-256 `{row.sha256}`" for row in input_manifest.itertuples(index=False))
    top_md = "\n".join(
        f"- {row.gene}: 非零频率 {row.nonzero_frequency:.2f}，主方向 {row.dominant_direction}，入模后方向一致率 {row.direction_consistency_given_selected if np.isfinite(row.direction_consistency_given_selected) else 'NA'}"
        for row in top.itertuples(index=False)
    )
    report = f"""# GSE123568 Figure 7 增强版血清分析报告

## 分析目的与证据边界

本分析针对 40 例 GSE123568 外周血清微阵列样本，评估预先冻结的 77 基因线粒体候选空间。它是探索性的内部交叉拟合评估，不是独立外部验证、临床诊断模型、锁定签名或因果证据；未选择临床阈值。

## Material Passport

- 数据集：GSE123568（GPL15207）
- 样本：30 例 SONFH，10 例 steroid-exposed non-SONFH
- 主模型：训练折内探针折叠 + StandardScaler/L1 logistic regression Pipeline
- 比较器：Ma 原 7 基因中按冻结注释规则可测的 BID、FTH1、LACTB、PDK3；L2 logistic regression
- 外层：5 折 × 5 次重复；内层：4 折；C 网格：`numpy.logspace(-3, 2, 16)`
- 冻结种子：{BASE_SEED}
- 正式计算：1,000 次完整标签置换、1,000 个有效 refit-aware OOB bootstrap、1,000 个匹配随机 77 基因空间

## 输入与 SHA-256

{inputs_md}

输入矩阵是 GEO 存档的处理后表达矩阵，实际范围见 `00_frozen_plan/expression_scale_qc.csv`；没有再次 log2 转换或整矩阵重新标准化。

## 样本、基因和防泄漏核对

样本 ID 与冻结顺序 GSM3507251–GSM3507290 完全一致，标签为 10 个对照和 30 个 SONFH。77 个候选基因与 Ma 4 个可测比较基因均通过平台探针核对。每一外层训练折独立选择训练均值最高的代表探针，完全相同均值时按 probe ID 字典序选择；同一 probe 再应用于对应测试折。人为改变测试折表达值后，训练折所选 probe 完全不变，泄漏单元测试通过。`StandardScaler` 位于 `Pipeline` 内，内层每个训练子折独立拟合。

## 内部交叉拟合性能

| 模型 | Aggregate AUC | AP | Brier |
|---|---:|---:|---:|
| 77 基因候选空间 | {candidate['auc']:.4f} | {candidate['average_precision']:.4f} | {candidate['brier']:.4f} |
| Ma 四基因比较器 | {ma['auc']:.4f} | {ma['average_precision']:.4f} | {ma['brier']:.4f} |

配对差值（candidate − Ma）：AUC {comparison['delta_auc']:.4f}，AP {comparison['delta_average_precision']:.4f}，Brier {comparison['delta_brier']:.4f}（Brier 越低越好）。配对 DeLong p={float(delong.paired_delong_p):.4g}。固定 cross-fitted 预测的 20,000 次配对 bootstrap AUC 差区间为 {fixed_auc.q025:.4f} 至 {fixed_auc.q975:.4f}。{compare_text}

## 模型重建不确定性

1,000 个有效 refit-aware OOB bootstrap 的 candidate−Ma AUC 差重拟合敏感性区间为 {refit_auc.q025:.4f} 至 {refit_auc.q975:.4f}；在达到 1,000 个有效重采样前记录 {len(invalid)} 个无效尝试。它包含重新抽样、训练折探针选择、标准化、调参和变量选择，但 OOB 样本较少且仍来自同一数据集，不能称为外部泛化置信区间。

## 标签置换和候选空间特异性

- 完整流程标签置换：n={int(perm.n_valid)}，经验 p={float(perm.empirical_p):.6f}。{perm_text}
- 表达/方差匹配随机空间：n={int(matched.n_matched_spaces)}，经验 p={float(matched.empirical_p_specificity):.6f}。{match_text}

随机空间的代表探针、表达均值、log1p 方差以及二维十分位匹配均只使用对应外层训练折；每个任务的每个外层折均包含 77 个不重复非候选基因。完整匹配映射保存在压缩 CSV 中。

## 基因选择和方向稳定性

25 个外层最终拟合中非零频率最高的基因如下。非零频率和方向一致性是稳定性描述，不是显著性、因果效应或独立生物标志物证据。

{top_md}

## Brier 与探索性校准

样本内患病率常数 0.75 的 Brier 基线为 {calibration.loc['candidate','prevalence_baseline_brier']:.4f}。候选模型 Brier={calibration.loc['candidate','brier']:.4f}，Brier skill={calibration.loc['candidate','brier_skill_score']:.4f}；Ma 比较器 Brier={calibration.loc['ma','brier']:.4f}，Brier skill={calibration.loc['ma','brier_skill_score']:.4f}。校准截距和斜率及 4 个等频分箱仅作探索性展示；n=40 且类别 30:10，不能解释为临床校准验证。未做 Hosmer–Lemeshow 检验，未选择诊断阈值。

## 与旧 Figure 7 的比较

旧值只是只读复核参照。新值允许因训练折内探针折叠和内层 scaler 防泄漏而变化，详见 `02_formal/legacy_vs_enhanced_comparison.csv`。无论数值方向如何，原来的审慎结论不变：这是内部、探索性、候选优先级生成分析。

## 失败、警告与限制

- 无效 OOB bootstrap 尝试数：{len(invalid)}；原因逐条保留在 `refit_oob_bootstrap_invalid.csv`。
- 模型收敛警告保存在 `selected_hyperparameters.csv`。
- 没有独立外部队列、预注册锁定签名或临床阈值。
- 5 次 CV repeat 不是 5 个独立患者，未据此进行显著性检验。
- 固定 OOF bootstrap 与 refit-aware OOB bootstrap 分开报告，二者估计对象不同。

## 结果判定

{perm_text} {match_text} {compare_text} 最终只能支持“在 GSE123568 内部探索性评估中量化候选空间的判别与稳定性”，不能支持临床诊断、外部泛化、候选空间因果性或优于既有比较器的无条件主张。

## 输出

完整文件级清单和 SHA-256 见 `04_report/output_manifest.csv`。Figure 7 审阅版位于 `03_figures`；本任务没有替换 `figures/final/Figure7.*`、没有改写 `results/figure_inputs`、没有修改主稿或投稿包。
"""
    report_path = REPORT_DIR / "GSE123568_Figure7_enhanced_analysis_report_zh.md"
    report_path.write_text(report, encoding="utf-8")
    write_environment()

    before = pd.read_csv(PLAN_DIR / "protected_artifact_manifest_before.csv")
    for row in before.itertuples(index=False):
        path = Path(row.path)
        if not path.is_file() or path.stat().st_size != int(row.size_bytes) or sha256(path) != row.sha256:
            raise RuntimeError(f"Protected artifact changed: {path}")

    manifest_rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file() and path.name != "output_manifest.csv":
            manifest_rows.append({"relative_path": str(path.relative_to(OUT)), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    for path in [ROOT / "analysis" / "serum_classifier_enhanced.py", ROOT / "analysis" / "requirements_serum_enhanced.txt", ROOT / "plotting" / "make_figure7_enhanced.R", ROOT / "workflow" / "run_serum_enhanced.ps1", ROOT / "qa" / "check_serum_enhanced.py"]:
        manifest_rows.append({"relative_path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    atomic_csv(REPORT_DIR / "output_manifest.csv", pd.DataFrame(manifest_rows))
    log_message("FINALIZATION COMPLETE", "report", report_path)


def preflight_mode() -> None:
    if not OUT.exists():
        OUT.mkdir(parents=True)
    bundle = load_probe_data()
    write_preflight(bundle)
    write_environment()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["preflight", "smoke", "formal", "finalize"], required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--data-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    configure_paths(args.project_root, args.data_dir, args.output_dir)
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if args.phase == "preflight":
        preflight_mode()
        return
    if args.phase == "finalize":
        finalize_outputs()
        return
    plan = PLAN_DIR / "analysis_plan_frozen.json"
    if not plan.is_file():
        raise RuntimeError("Run preflight and freeze the plan before any model fitting")
    bundle = load_probe_data()
    run_analysis(bundle, args.phase, 20 if args.phase == "smoke" else 1000, args.jobs, resume=args.resume)
    write_environment()


if __name__ == "__main__":
    main()
