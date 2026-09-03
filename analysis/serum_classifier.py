# -*- coding: utf-8 -*-
"""Figure 7 v4: reproducible leakage-free blood classifier evaluation.

The primary statistic is the AUC of per-sample probabilities averaged across
five repeated outer five-fold CV runs. Every permutation repeats that exact
pipeline. The script checkpoints permutation results and can be resumed.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

warnings.simplefilter("ignore", FutureWarning)

from v4_common import ANALYSIS, FIGS, ROOT, save_figure, setup_plot

setup_plot()
import matplotlib.pyplot as plt
import seaborn as sns

BASE_SEED = 20260820
N_REPEATS = 5
N_OUTER_FOLDS = 5
N_INNER_FOLDS = 4
CS = np.logspace(-3, 2, 16)
DATA_ROOT = Path(os.environ.get("ONFH_DATA_ROOT", ROOT / "data")).resolve()
OUTPUT = Path(os.environ.get("ONFH_OUTPUT_DIR", ANALYSIS)).resolve()
LOG = OUTPUT / "py_fig7_v4.log"


def lg(*parts) -> None:
    line = time.strftime("%H:%M:%S") + " | " + " ".join(map(str, parts))
    print(line, flush=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def load_gse123568() -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    matrix_path = DATA_ROOT / "GSE123568_series_matrix.txt.gz"
    with gzip.open(matrix_path, "rt", errors="replace") as handle:
        lines = handle.readlines()

    disease_lines = [
        line for line in lines
        if line.startswith("!Sample_characteristics_ch1") and "disease:" in line
    ]
    if len(disease_lines) != 1:
        raise RuntimeError(f"expected one disease metadata line, got {len(disease_lines)}")
    disease = [token.strip().strip('"') for token in disease_lines[0].split("\t")[1:]]

    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    header = [x.strip('"') for x in lines[start + 1].strip().split("\t")]
    rows, probe_ids = [], []
    for line in lines[start + 2 :]:
        if line.startswith("!series_matrix_table_end"):
            break
        parts = line.strip().split("\t")
        probe_ids.append(parts[0].strip('"'))
        rows.append([float(x) for x in parts[1:]])
    expr = pd.DataFrame(rows, index=probe_ids, columns=header[1:])
    if expr.shape[1] != len(disease):
        raise RuntimeError("expression columns and disease labels are misaligned")

    with gzip.open(DATA_ROOT / "GSE123568_family.soft.gz", "rt", errors="replace") as handle:
        soft = handle.read()
    platform = soft.split("!platform_table_begin", 1)[1].split("!platform_table_end", 1)[0]
    platform_lines = platform.strip().split("\n")
    platform_header = platform_lines[0].split("\t")
    symbol_col = next(i for i, name in enumerate(platform_header) if "Gene Symbol" in name)
    probe_to_symbol: dict[str, str] = {}
    for row in platform_lines[1:]:
        cells = row.split("\t")
        if len(cells) <= symbol_col:
            continue
        symbol = cells[symbol_col].strip().split(" /// ")[0].split("//")[0].strip()
        if symbol:
            probe_to_symbol[cells[0]] = symbol

    expr["symbol"] = expr.index.map(probe_to_symbol)
    expr = expr.dropna(subset=["symbol"])
    expr = expr[expr["symbol"] != ""]
    sample_columns = [c for c in expr.columns if c != "symbol"]
    expr["rowmean"] = expr[sample_columns].mean(axis=1)
    expr = expr.sort_values("rowmean", ascending=False).drop_duplicates("symbol")
    expr = expr.set_index("symbol").drop(columns="rowmean")

    labels = np.array([1 if "disease: SONFH" in item else 0 for item in disease], dtype=int)
    if (int(labels.sum()), int((1 - labels).sum())) != (30, 10):
        raise RuntimeError("expected 30 SONFH and 10 steroid-exposed non-SONFH samples")
    return expr, labels, sample_columns


def candidate_matrix(expr: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    with (ANALYSIS / "genesets_final.json").open(encoding="utf-8") as handle:
        gene_sets = json.load(handle)
    modules = ["Mito_fission", "Mito_fusion", "Mitophagy_core", "mtDNA_release", "cGAS_STING", "YAP_mTOR"]
    candidates = sorted({gene for module in modules for gene in gene_sets[module] if gene in expr.index})
    return expr.loc[candidates].T.to_numpy(dtype=float), np.asarray(candidates)


def fit_l1_predict(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_test: np.ndarray,
    *,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    scaler = StandardScaler().fit(x_train)
    xtr = scaler.transform(x_train)
    xte = scaler.transform(x_test)
    inner = StratifiedKFold(N_INNER_FOLDS, shuffle=True, random_state=seed)
    model = LogisticRegressionCV(
        Cs=CS,
        cv=inner,
        penalty="l1",
        solver="liblinear",
        scoring="roc_auc",
        class_weight="balanced",
        max_iter=5000,
        random_state=seed,
        n_jobs=1,
        refit=True,
    )
    model.fit(xtr, y_train)
    selected = np.flatnonzero(np.abs(model.coef_[0]) > 1e-12)
    return model.predict_proba(xte)[:, 1], selected


def fit_l2_predict(x_train, y_train, x_test, *, seed: int) -> np.ndarray:
    scaler = StandardScaler().fit(x_train)
    xtr = scaler.transform(x_train)
    xte = scaler.transform(x_test)
    inner = StratifiedKFold(N_INNER_FOLDS, shuffle=True, random_state=seed)
    model = LogisticRegressionCV(
        Cs=CS,
        cv=inner,
        penalty="l2",
        solver="liblinear",
        scoring="roc_auc",
        class_weight="balanced",
        max_iter=5000,
        random_state=seed,
        n_jobs=1,
        refit=True,
    )
    model.fit(xtr, y_train)
    return model.predict_proba(xte)[:, 1]


def repeated_nested_oof(
    x: np.ndarray,
    y: np.ndarray,
    *,
    record_features: bool,
    l2: bool = False,
) -> dict:
    probabilities = np.full((N_REPEATS, len(y)), np.nan, dtype=float)
    selected: list[np.ndarray] = []
    repeat_aucs = []
    repeat_average_precision = []
    for repeat in range(N_REPEATS):
        outer = StratifiedKFold(N_OUTER_FOLDS, shuffle=True, random_state=repeat)
        for fold, (train, test) in enumerate(outer.split(x, y)):
            seed = BASE_SEED + repeat * 100 + fold
            if l2:
                pred = fit_l2_predict(x[train], y[train], x[test], seed=seed)
            else:
                pred, chosen = fit_l1_predict(x[train], y[train], x[test], seed=seed)
                if record_features:
                    selected.append(chosen)
            probabilities[repeat, test] = pred
        repeat_aucs.append(float(roc_auc_score(y, probabilities[repeat])))
        repeat_average_precision.append(float(average_precision_score(y, probabilities[repeat])))
    mean_probability = probabilities.mean(axis=0)
    return {
        "probabilities": probabilities,
        "mean_probability": mean_probability,
        "repeat_aucs": np.asarray(repeat_aucs),
        "repeat_average_precision": np.asarray(repeat_average_precision),
        "aggregate_auc": float(roc_auc_score(y, mean_probability)),
        "aggregate_average_precision": float(average_precision_score(y, mean_probability)),
        "selected": selected,
    }


def stratified_prediction_bootstrap(y, probability, n_boot: int = 5000) -> tuple[float, float]:
    """CI for the AUC of fixed, aggregated cross-fitted predictions."""
    rng = np.random.default_rng(BASE_SEED)
    case = np.flatnonzero(y == 1)
    control = np.flatnonzero(y == 0)
    aucs = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = np.concatenate(
            [rng.choice(case, len(case), replace=True), rng.choice(control, len(control), replace=True)]
        )
        aucs[i] = roc_auc_score(y[idx], probability[idx])
    return tuple(np.percentile(aucs, [2.5, 97.5]).astype(float))


def stratified_paired_auc_bootstrap(
    y: np.ndarray,
    probability_a: np.ndarray,
    probability_b: np.ndarray,
    n_boot: int = 20000,
) -> dict[str, float]:
    """Paired uncertainty for the fixed cross-fitted AUC difference.

    Both models are evaluated on identical stratified resamples of the same
    participants. This estimates uncertainty of the archived cross-fitted
    predictions and, deliberately, does not claim to include model-refitting
    variability.
    """
    rng = np.random.default_rng(BASE_SEED + 71)
    case = np.flatnonzero(y == 1)
    control = np.flatnonzero(y == 0)
    delta = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = np.concatenate(
            [rng.choice(case, len(case), replace=True), rng.choice(control, len(control), replace=True)]
        )
        delta[i] = roc_auc_score(y[idx], probability_a[idx]) - roc_auc_score(
            y[idx], probability_b[idx]
        )
    observed = roc_auc_score(y, probability_a) - roc_auc_score(y, probability_b)
    lower, upper = np.percentile(delta, [2.5, 97.5]).astype(float)
    lower_tail = (np.sum(delta <= 0) + 1) / (n_boot + 1)
    upper_tail = (np.sum(delta >= 0) + 1) / (n_boot + 1)
    return {
        "delta_auc": float(observed),
        "ci_low": float(lower),
        "ci_high": float(upper),
        "two_sided_p": float(min(1.0, 2 * min(lower_tail, upper_tail))),
        "n_bootstrap": int(n_boot),
    }


def one_permutation(perm_id: int, x: np.ndarray, y: np.ndarray) -> dict:
    rng = np.random.default_rng(BASE_SEED + 10_000 + perm_id)
    yp = rng.permutation(y)
    result = repeated_nested_oof(x, yp, record_features=False)
    return {
        "perm_id": perm_id,
        "perm_auc": result["aggregate_auc"],
        "perm_repeat_auc_mean": float(result["repeat_aucs"].mean()),
    }


def run_permutations(x, y, target: int, jobs: int, batch_size: int = 20) -> pd.DataFrame:
    out = OUTPUT / "diag_permutation_v4.csv"
    if out.exists():
        done = pd.read_csv(out)
        expected = {"perm_id", "perm_auc", "perm_repeat_auc_mean"}
        if not expected.issubset(done.columns):
            raise RuntimeError("existing diag_permutation_v4.csv has an incompatible schema")
        done = done.drop_duplicates("perm_id").sort_values("perm_id")
    else:
        done = pd.DataFrame(columns=["perm_id", "perm_auc", "perm_repeat_auc_mean"])

    completed = set(done["perm_id"].astype(int)) if len(done) else set()
    pending = [i for i in range(target) if i not in completed]
    lg("permutations", len(completed), "complete", len(pending), "pending", "jobs", jobs)
    for start in range(0, len(pending), batch_size):
        batch = pending[start : start + batch_size]
        rows = Parallel(n_jobs=jobs, prefer="processes")(
            delayed(one_permutation)(i, x, y) for i in batch
        )
        done = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
        done = done.drop_duplicates("perm_id").sort_values("perm_id")
        done.to_csv(out, index=False)
        lg("permutation checkpoint", len(done), "/", target)
    return done[done["perm_id"].astype(int) < target].sort_values("perm_id")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Skip the legacy matplotlib preview; manuscript figures are generated in R.",
    )
    args = parser.parse_args()
    if args.permutations < 0:
        raise ValueError("--permutations must be non-negative")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    if LOG.exists():
        LOG.unlink()
    expr, y, sample_names = load_gse123568()
    x, genes = candidate_matrix(expr)
    lg("GSE123568 serum", x.shape, "SONFH", int(y.sum()), "steroid-exposed non-SONFH", int((1-y).sum()))

    observed = repeated_nested_oof(x, y, record_features=True)
    repeat_aucs = observed["repeat_aucs"]
    aggregate_auc = observed["aggregate_auc"]
    aggregate_average_precision = observed["aggregate_average_precision"]
    ci_low, ci_high = stratified_prediction_bootstrap(y, observed["mean_probability"])
    lg("observed aggregate OOF AUC", f"{aggregate_auc:.3f}", "repeat mean", f"{repeat_aucs.mean():.3f}")

    feature_counts = pd.Series(0, index=genes, dtype=float)
    for chosen in observed["selected"]:
        feature_counts.iloc[chosen] += 1
    feature_frequency = (feature_counts / (N_REPEATS * N_OUTER_FOLDS)).sort_values(ascending=False)
    feature_frequency.rename("selection_frequency").to_csv(OUTPUT / "diag_feature_stability_v4.csv")
    pd.DataFrame({"repeat": range(N_REPEATS), "AUC": repeat_aucs}).to_csv(
        OUTPUT / "diag_nested_cv_auc_v4.csv", index=False
    )
    pd.DataFrame(
        {
            "repeat": range(N_REPEATS),
            "AUC": repeat_aucs,
            "average_precision": observed["repeat_average_precision"],
        }
    ).to_csv(OUTPUT / "diag_nested_cv_performance_v7.csv", index=False)

    ma_all = ["BID", "FTH1", "LACTB", "PDK3", "RAB5IF", "SOD2", "SQOR"]
    ma_available = [g for g in ma_all if g in expr.index]
    ma = repeated_nested_oof(expr.loc[ma_available].T.to_numpy(dtype=float), y, record_features=False, l2=True)
    paired = stratified_paired_auc_bootstrap(
        y, observed["mean_probability"], ma["mean_probability"]
    )

    oof_rows = []
    for repeat in range(N_REPEATS):
        for i, sample in enumerate(sample_names):
            oof_rows.append(
                {
                    "repeat": repeat,
                    "sample": sample,
                    "label": int(y[i]),
                    "oof_probability": float(observed["probabilities"][repeat, i]),
                    "candidate_oof_probability": float(observed["probabilities"][repeat, i]),
                    "ma_oof_probability": float(ma["probabilities"][repeat, i]),
                }
            )
    pd.DataFrame(oof_rows).to_csv(OUTPUT / "diag_oof_predictions_v4.csv", index=False)
    pd.DataFrame(
        {
            "sample": sample_names,
            "label": y.astype(int),
            "candidate_mean_oof_probability": observed["mean_probability"],
            "ma_mean_oof_probability": ma["mean_probability"],
        }
    ).to_csv(OUTPUT / "diag_oof_predictions_aggregated_v8.csv", index=False)
    pd.DataFrame(
        {
            "repeat": range(N_REPEATS),
            "AUC": ma["repeat_aucs"],
            "average_precision": ma["repeat_average_precision"],
        }
    ).to_csv(OUTPUT / "diag_ma_comparator_repeat_performance_v7.csv", index=False)
    pd.DataFrame(
        [
            {
                "comparison": "mitochondrial_candidate_space_minus_Ma_four_gene_model",
                **paired,
                "scope": "paired stratified bootstrap of fixed aggregated cross-fitted predictions; model-refitting variability excluded",
            }
        ]
    ).to_csv(OUTPUT / "diag_paired_model_comparison_v8.csv", index=False)

    if args.permutations:
        permutation = run_permutations(x, y, args.permutations, max(1, args.jobs))
        null = permutation["perm_auc"].to_numpy(dtype=float)
        permutation_p = float((np.sum(null >= aggregate_auc) + 1) / (len(null) + 1))
    else:
        permutation = pd.DataFrame(columns=["perm_id", "perm_auc", "perm_repeat_auc_mean"])
        null = np.array([], dtype=float)
        permutation_p = np.nan

    summary = {
        "primary_statistic": "AUC of per-sample OOF probabilities averaged over 5 repeated outer 5-fold CV runs",
        "model": "standardization plus inner-CV-tuned L1 logistic regression, fitted inside each outer training fold",
        "n_samples": int(len(y)),
        "n_sonfh": int(y.sum()),
        "n_steroid_exposed_non_sonfh": int((1 - y).sum()),
        "sample_source": "peripheral serum",
        "candidate_genes_on_array": int(len(genes)),
        "aggregate_oof_auc": round(aggregate_auc, 4),
        "aggregate_oof_average_precision": round(aggregate_average_precision, 4),
        "average_precision_prevalence_baseline": round(float(y.mean()), 4),
        "repeat_auc_mean": round(float(repeat_aucs.mean()), 4),
        "repeat_auc_sd": round(float(repeat_aucs.std(ddof=1)), 4),
        "repeat_average_precision_mean": round(float(observed["repeat_average_precision"].mean()), 4),
        "repeat_average_precision_sd": round(float(observed["repeat_average_precision"].std(ddof=1)), 4),
        "cross_fitted_prediction_bootstrap_ci95": [round(ci_low, 4), round(ci_high, 4)],
        "ci_scope": "stratified sample bootstrap of fixed aggregated cross-fitted predictions",
        "n_permutations": int(len(null)),
        "permutation_empirical_p": None if not np.isfinite(permutation_p) else round(permutation_p, 6),
        "permutation_statistic_matches_observed": True,
        "ma2024_available_genes": ma_available,
        "ma2024_aggregate_oof_auc": round(float(ma["aggregate_auc"]), 4),
        "ma2024_aggregate_oof_average_precision": round(float(ma["aggregate_average_precision"]), 4),
        "ma2024_repeat_auc_mean": round(float(ma["repeat_aucs"].mean()), 4),
        "ma2024_repeat_auc_sd": round(float(ma["repeat_aucs"].std(ddof=1)), 4),
        "ma2024_repeat_average_precision_mean": round(float(ma["repeat_average_precision"].mean()), 4),
        "ma2024_repeat_average_precision_sd": round(float(ma["repeat_average_precision"].std(ddof=1)), 4),
        "candidate_minus_ma_auc": round(float(paired["delta_auc"]), 4),
        "paired_fixed_prediction_bootstrap_ci95": [
            round(float(paired["ci_low"]), 4),
            round(float(paired["ci_high"]), 4),
        ],
        "paired_fixed_prediction_bootstrap_p": round(float(paired["two_sided_p"]), 6),
        "paired_comparison_scope": "paired stratified bootstrap of fixed aggregated cross-fitted predictions; model-refitting variability excluded",
        "feature_stability_top": {k: round(float(v), 3) for k, v in feature_frequency.head(15).items()},
    }
    with (OUTPUT / "diag_summary_v4.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)

    if args.no_plot:
        lg("legacy matplotlib preview skipped; use the R manuscript-figure workflow")
        lg("done", json.dumps(summary, ensure_ascii=False))
        return

    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.3))
    ax = axes[0]
    if len(null):
        ax.hist(null, bins=min(30, max(10, len(null) // 20)), color="#9ecae1", edgecolor="white")
    ax.axvline(aggregate_auc, color="#cb181d", lw=2, label=f"observed aggregate OOF AUC={aggregate_auc:.3f}")
    title_p = "not run" if not np.isfinite(permutation_p) else f"p={permutation_p:.4f}"
    ax.set_title(f"Matched full-pipeline permutation ({len(null)}), {title_p}")
    ax.set_xlabel("aggregate repeated-OOF AUC under label permutation")
    ax.legend(fontsize=7.5)

    ax = axes[1]
    top = feature_frequency[feature_frequency > 0].head(15)
    ax.barh(top.index[::-1], top.values[::-1], color="#4292c6")
    ax.set_xlabel("selection frequency across 25 outer folds")
    ax.set_title("Feature-selection stability")

    ax = axes[2]
    perf = pd.DataFrame(
        {
            "model": ["EC mitochondrial\ncandidates"] * N_REPEATS + [f"Ma et al.\n({len(ma_available)}/7 genes)"] * N_REPEATS,
            "repeat_auc": list(repeat_aucs) + list(ma["repeat_aucs"]),
        }
    )
    sns.stripplot(data=perf, x="model", y="repeat_auc", hue="model", legend=False, palette="Set2", size=7, jitter=0.08, ax=ax)
    for j, values in enumerate([repeat_aucs, ma["repeat_aucs"]]):
        ax.errorbar(j, values.mean(), yerr=values.std(ddof=1), color="black", capsize=4, lw=1.3)
    ax.set_ylim(0.35, 1.03)
    ax.set_xlabel("")
    ax.set_ylabel("AUC per outer-CV repeat")
    ax.set_title("Five repeats of nested outer five-fold CV")
    fig.suptitle("GSE123568 peripheral-serum classifier: exploratory internal validation", y=1.02)
    fig.tight_layout()
    save_figure(fig, FIGS / "fig6_v4_nested_cv.png")
    plt.close(fig)

    lg("done", json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
