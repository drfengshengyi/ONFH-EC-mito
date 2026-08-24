# ONFH endothelial and mitochondrial reanalysis

Reproducible code and processed outputs for the manuscript:

> **Participant-aware multi-cohort analysis characterizes heterogeneous endothelial and mitochondrial profiles in osteonecrosis of the femoral head**

This release separates data preparation, biological analyses, the official-R `SQSTM1` virtual-knockout audit, and manuscript plotting. Primary public matrices are not redistributed. Processed result tables and final figures are included so the reported outputs can be inspected without rebuilding the complete atlas.

## What is in this repository

| Directory | Contents |
|---|---|
| `analysis/` | Ordered Python/R analysis modules and frozen analysis manifests |
| `virtual_knockout/` | Donor-separated official `scTenifoldKnk` analysis and post-processing |
| `plotting/` | Figure 1–6 and Supplementary Figure S1 generation/assembly code |
| `workflow/` | PowerShell entry points for the three reproducible stages |
| `results/` | Processed tables, virtual-knockout results, and Supplementary Tables S1–S10 |
| `figures/source/` | Versioned source panels used by the layout script |
| `figures/final/` | Submission-ready PDF and PNG figures |
| `data/` | Download instructions and a machine-readable dataset manifest; raw data are ignored |
| `environment/` | Frozen Python dependencies and R installation helper |
| `qa/` | Repository integrity and portability checks |

## Public datasets

| Accession | Analysis role |
|---|---|
| SRP361778 | Liao femoral-head single-cell cohort; participant-level inference and HOA control-EC virtual knockout |
| GSE169396 | Independent healthy femoral-head single-cell atlas component |
| GSE290411 | Steroid-induced ONFH libraries; descriptive cross-cohort contrasts only because participant mapping was unavailable |
| GSE123568 | Peripheral-serum expression dataset; repeated nested cross-validation |

The exact local filenames expected by the scripts are listed in [`data/README.md`](data/README.md) and [`data/datasets.tsv`](data/datasets.tsv).

## Quick start

Run these commands from the repository root in PowerShell. Override `-Python` or `-Rscript` if the executables are not on `PATH`.

```powershell
# 1. Create an environment
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r environment\requirements-python.lock.txt

# 2. Rebuild the main analysis after placing public data under data/
.\workflow\run_core_analysis.ps1 -Python python -Rscript Rscript -Jobs 4

# 3. Run both official-R virtual-knockout profiles
.\workflow\run_virtual_knockout.ps1 -Python python -Rscript Rscript -Cores 1

# 4. Rebuild final figures from the versioned source panels/results
.\workflow\run_figures.ps1 -Python python

# 5. Audit repository structure, syntax, metadata, and portability
python qa\check_repository.py

# Optional release audit: also require byte-identical manuscript figure PDFs
python qa\check_repository.py --submission-dir <path-to-submission-figures>
```

For a complete run, use `workflow/run_all.ps1`. The virtual-knockout R packages are installed separately with `environment/install_official_sctenifoldknk.R`; see [`virtual_knockout/README.md`](virtual_knockout/README.md).

## Analysis boundary

- Statistical inference is performed at the participant/sampling-unit level where participant identity is available.
- The four GSE290411 SONFH libraries are retained for descriptive visualization and effect-size comparison, not treated as four confirmed independent participants.
- The `SQSTM1` virtual knockout is an expression-derived network perturbation performed separately in HOA2 and HOA3. It is **not** a wet-lab knockout and **not** a CellOracle cell-fate or differentiation-trajectory simulation.
- Virtual-knockout gene and pathway ranks are exploratory predictions. They do not establish causal `SQSTM1` dependence, mitochondrial function, or clinical utility.

## Reproducibility notes

- `analysis/sample_metadata_v4.csv` is the authoritative library-to-participant map.
- Randomized analyses use fixed seeds recorded in code and provenance files.
- `results/official_r_vko_manuscript/` and `results/official_r_vko_official_default/` retain parameter, package-version, and session information.
- Final figure PDFs in `figures/final/` are the submission versions; PNG counterparts are included for rapid review.

## Data availability

All datasets analyzed in this study are publicly available: SRP361778 in the NCBI Sequence Read Archive and GSE169396, GSE290411, and GSE123568 in NCBI GEO. Analysis code, frozen manifests, processed result tables, and figure-generation code are provided in this repository. No newly generated primary sequencing data are reported.

## Contacts

- Qiugen Wang: wangqiugen@126.com
- Jianguang Xu: xjg@shutcm.edu.cn

## License and third-party data

No license is asserted for the public datasets or third-party databases. Users should follow the terms of the original repositories and packages.
