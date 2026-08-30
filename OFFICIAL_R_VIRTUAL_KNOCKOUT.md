# Official-R scTenifoldKnk audit

This workflow reruns the donor-separated `SQSTM1` analysis with the authors' official `cailab-tamu/scTenifoldKnk` R package. It freezes the existing 300-gene feature list and audited HOA endothelial-cell barcodes.

## Install

From the repository root:

```powershell
Rscript environment/install_official_sctenifoldknk.R
```

If `Rscript` is not on `PATH`, pass its full executable path to the PowerShell workflow instead of editing any analysis script.

## Run

Run the manuscript profile and package-default sensitivity analysis together:

```powershell
.\workflow\run_virtual_knockout.ps1 -Rscript Rscript -Python python -Cores 1
```

The core and sensitivity calls include:

```powershell
Rscript virtual_knockout/run_official_vko.R --profile=manuscript --cores=1
Rscript virtual_knockout/run_official_vko.R --profile=official-default --cores=1
Rscript virtual_knockout/run_official_vko.R --profile=manuscript --cores=1 --exclude-mt-encoded
Rscript virtual_knockout/run_official_vko.R --profile=official-default --cores=1 --exclude-mt-encoded
Rscript virtual_knockout/run_matched_control_vko.R --cores=1
Rscript virtual_knockout/postprocess_mt_exclusion.R
```

Use one core for the cleanest reproducibility audit. Donor-specific CSV/RDS files, the consensus table, provenance JSON, session information, and run summary are written under `results/`.

## Important implementation details

The official R wrapper applies CPM normalization, samples cells with replacement through `scTenifoldNet::makeNetworks`, and internally calls `set.seed(1)`. Agreement between parameter profiles is treated as sensitivity evidence; disagreement is retained and reported.

Virtual-knockout rankings remain computational predictions. They are not wet-lab validation, a CellOracle trajectory perturbation, or causal proof of `SQSTM1` function.

The mtDNA-feature sensitivity completely refits the donor networks on 295 genes. Matched-gene perturbations reuse the frozen donor-specific WT networks because the official wrapper constructs WT before the target row is zeroed; an exact SQSTM1 rerun is retained as a validation check. Comparator genes are not biological negative controls.

See `OFFICIAL_R_VKO_RESULTS.md` for the completed result audit.
