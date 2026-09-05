# ONFH endothelial transcriptomic reanalysis

Code and processed outputs for:

> **Participant-aware transcriptomic benchmarking of endothelial states and mitochondrial hypotheses in osteonecrosis of the femoral head**

The repository separates analysis, virtual knockout, spatial contextualization, and figure generation. Raw public matrices are not redistributed. Compact results and final figures are included for audit and reuse.

## Repository map

| Directory | Contents |
|---|---|
| `analysis/` | Core Python and R analyses |
| `virtual_knockout/` | Donor-separated `scTenifoldKnk` analyses |
| `plotting/` | Figures 1–7 and Supplementary Figure S1 |
| `workflow/` | Reproducible PowerShell entry points |
| `results/` | Processed tables and supplementary results |
| `figures/` | Source panels and final figures |
| `data/` | Download instructions and dataset manifest |
| `environment/` | Frozen Python and R environments |
| `qa/` | Syntax, integrity, and portability checks |

## Public datasets

| Accession | Role |
|---|---|
| SRP361778 | Femoral-head scRNA-seq; participant-level analyses and HOA2/HOA3 virtual knockout |
| GSE169396 | Healthy femoral-head scRNA-seq context |
| GSE290411 | SONFH libraries; descriptive contrasts because participant mapping is unavailable |
| GSE284089 | One OA femoral-head spatial section; anatomical context only |
| GSE123568 | Serum expression; repeated nested cross-validation |

Expected files are listed in [`data/README.md`](data/README.md) and [`data/datasets.tsv`](data/datasets.tsv). Virtual-knockout input checksums are in [`virtual_knockout/vko_input_checksums.tsv`](virtual_knockout/vko_input_checksums.tsv).

## Quick start

Run from the repository root in PowerShell. Override `-Python` or `-Rscript` when needed.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r environment\requirements-python.lock.txt
Rscript environment\check_r_packages.R

.\workflow\run_core_analysis.ps1 -Python python -Rscript Rscript -Jobs 4
.\workflow\run_virtual_knockout.ps1 -Python python -Rscript Rscript -Cores 1
.\workflow\run_spatial_contextualization.ps1 -Python python
.\workflow\run_figures.ps1 -Python python -Rscript Rscript
python qa\check_repository.py --rscript Rscript
```

`workflow/run_all.ps1` runs the complete sequence. Install missing virtual-knockout R packages with `environment/install_official_sctenifoldknk.R`.

## Analysis limits

- Formal inference uses participants or sampling units when their identities are available.
- GSE290411 is descriptive because its participant-to-library map is unresolved.
- GSE284089 contains one OA spatial section. It provides anatomical context, not ONFH or causal validation.
- The `SQSTM1` analysis is an expression-derived network perturbation, not an experimental knockout or cell-fate simulation.
- MtDNA-feature-excluded refits and 20 matched-gene perturbations calibrate robustness and specificity.
- Virtual-knockout ranks remain exploratory. The interpretation emphasizes a heterogeneous selective-autophagy receptor system, with `CALCOCO2`/NDP52 as a secondary candidate and `OPTN` as a context control.

## Reproducibility

- `analysis/sample_metadata_v4.csv` is the authoritative sampling-unit map.
- Randomized analyses use fixed seeds recorded in code and provenance files.
- Python 3.12 is fixed in `.python-version`; R versions are in `environment/r-package-versions.tsv`.
- `results/README.md` links each figure and table to its producer and inputs.
- Final submission figures are in `figures/final/` as PDF and PNG.

## Data availability

All analyzed datasets are public. This repository provides code, frozen manifests, processed tables, and figure scripts. No new primary sequencing data were generated.

## Contacts

- Qiugen Wang: wangqiugen@126.com
- Jianguang Xu: xjg@shutcm.edu.cn

## Third-party data

Public datasets and databases remain subject to their original terms of use.
