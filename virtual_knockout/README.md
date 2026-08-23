# Official-R SQSTM1 virtual knockout

This module uses the official `scTenifoldKnk`/`scTenifoldNet` R implementation. HOA2 and HOA3 control endothelial cells are analyzed separately; their networks are never pooled as if they were a single participant.

## Install

```powershell
Rscript environment\install_official_sctenifoldknk.R
```

The primary manuscript profile uses 20 networks, `q = 0.95`, and 30 manifold dimensions. A package-default sensitivity profile uses 10 networks, `q = 0.90`, and 2 dimensions. The frozen 300-gene feature list and cell-barcode manifest are versioned in this directory.

## Run

```powershell
.\workflow\run_virtual_knockout.ps1 -Rscript Rscript -Python python -Cores 1
```

The workflow performs:

1. donor-separated official-R perturbation for the manuscript profile;
2. donor-separated package-default sensitivity analysis;
3. matched-null pathway post-processing and cross-profile audit;
4. export of plot-ready R objects; and
5. generation of Figure 5.

Outputs are written under `results/official_r_vko_*`. Read `OFFICIAL_R_VKO_RESULTS.md` for the frozen result audit.

## Interpretation

Zeroing the inferred `SQSTM1` regulatory row simulates loss of its outgoing network influence. The output is a gene-manifold perturbation rank. It is not a cell-state velocity field, a differentiation trajectory, or experimental knockout validation.
