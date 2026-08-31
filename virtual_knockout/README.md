# Official-R SQSTM1 virtual knockout

This module uses the official `scTenifoldKnk`/`scTenifoldNet` R implementation. HOA2 and HOA3 control endothelial cells are analyzed separately; their networks are never pooled as if they were a single participant.

## Install

```powershell
Rscript environment\install_official_sctenifoldknk.R
```

The primary manuscript profile uses 20 networks, `q = 0.95`, and 30 manifold dimensions. A package-default sensitivity profile uses 10 networks, `q = 0.90`, and 2 dimensions. The frozen 300-gene feature list and cell-barcode manifest are versioned in this directory.
The byte sizes and SHA256 checksums of the exact HOA2/HOA3 RDS inputs used for the release are recorded in `vko_input_checksums.tsv`; primary matrices themselves remain excluded from Git.

## Run

```powershell
.\workflow\run_virtual_knockout.ps1 -Rscript Rscript -Python python -Cores 1
```

The workflow performs:

1. donor-separated official-R perturbation for the manuscript profile;
2. donor-separated package-default sensitivity analysis;
3. complete donor-specific refits of both profiles after removing the five mtDNA-encoded frozen features;
4. twenty expression- and WT-topology-matched comparator-gene perturbations under both profiles;
5. common-nuclear-universe rank and matched-null pathway audits;
6. export of plot-ready R objects; and
7. R-only generation of Figure 5.

Outputs are written under `results/official_r_vko_*`. Read `OFFICIAL_R_VKO_RESULTS.md` for the frozen result audit.
The exact cross-donor common-nuclear correlations and the primary-refit HOA3 FDR-positive identities used in Table S10g are frozen in `results/official_r_vko_no_mt_cross_donor_audit.csv` and `results/official_r_vko_no_mt_hoa3_fdr_audit.csv`, respectively.

## Interpretation

Zeroing the inferred `SQSTM1` regulatory row simulates loss of its outgoing network influence. The output is a gene-manifold perturbation rank. It is not a cell-state velocity field, a differentiation trajectory, or experimental knockout validation.

The matched genes are computational comparators selected before their perturbation outputs were inspected. They are not assumed to be biologically inert. The mtDNA-feature-exclusion analysis removes `MT-ATP6`, `MT-CO1`, `MT-CO2`, `MT-ND1`, and `MT-ND4` from the input expression matrices and completely rebuilds each network on the remaining 295 genes.
