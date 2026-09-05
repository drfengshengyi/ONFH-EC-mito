# SQSTM1 virtual knockout

This module uses the official `scTenifoldKnk` and `scTenifoldNet` R packages. HOA2 and HOA3 endothelial cells are modeled separately.

## Install

```powershell
Rscript environment\install_official_sctenifoldknk.R
```

The primary profile uses 20 networks, `q = 0.95`, and 30 manifold dimensions. The sensitivity profile uses 10 networks, `q = 0.90`, and 2 dimensions. Features, barcodes, input sizes, and SHA256 values are frozen in this directory.

## Run

```powershell
.\workflow\run_virtual_knockout.ps1 -Rscript Rscript -Python python -Cores 1
```

The workflow runs:

1. donor-separated primary and sensitivity profiles;
2. complete refits after excluding five mtDNA-encoded genes;
3. 20 expression- and topology-matched gene perturbations;
4. cross-donor and matched-null audits;
5. plot-data export and Figure 6 generation.

Outputs are stored under `results/official_r_vko_*`. Frozen result checks are summarized in `OFFICIAL_R_VKO_RESULTS.md`.

## Interpretation

The analysis sets the inferred `SQSTM1` outgoing network row to zero and ranks gene-manifold changes. It does not model cell-state trajectories or validate an experimental knockout. Matched genes are calibration comparators, not biologically inert controls.
