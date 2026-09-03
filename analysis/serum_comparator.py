"""Export the deterministic Ma et al. comparator used in Figure 7.

This script performs the comparator fit once and writes the repeat-level
performance table.  Submission-figure assembly reads that saved table and does
not refit any model.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import serum_classifier as diagnostic


ANALYSIS = Path(__file__).resolve().parent
OUTPUT = ANALYSIS / "diag_ma_comparator_repeat_performance_v7.csv"
MA_GENES = ["BID", "FTH1", "LACTB", "PDK3", "RAB5IF", "SOD2", "SQOR"]


def main() -> None:
    expression, labels, _ = diagnostic.load_gse123568()
    available = [gene for gene in MA_GENES if gene in expression.index]
    result = diagnostic.repeated_nested_oof(
        expression.loc[available].T.to_numpy(dtype=float),
        labels,
        record_features=False,
        l2=True,
    )
    table = pd.DataFrame(
        {
            "repeat": range(len(result["repeat_aucs"])),
            "AUC": result["repeat_aucs"],
            "average_precision": result["repeat_average_precision"],
            "n_genes_available": len(available),
            "genes_available": ";".join(available),
        }
    )
    table.to_csv(OUTPUT, index=False)
    print(OUTPUT)


if __name__ == "__main__":
    main()
