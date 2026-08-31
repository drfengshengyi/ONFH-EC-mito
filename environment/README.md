# Environments

- `requirements-python.lock.txt`: complete frozen Python environment used for the analysis release.
- `requirements-vko-python.txt`: smaller Python dependency set for virtual-knockout post-processing.
- `install_official_sctenifoldknk.R`: installs/checks the official R packages required for the virtual knockout.
- `r-package-versions.tsv`: exact R runtime and directly used package versions archived for this release.
- `check_r_packages.R`: read-only verification of the installed R environment against that manifest.

R package versions used by each completed virtual-knockout run are recorded in the corresponding `results/official_r_vko_*/sessionInfo.txt` and provenance JSON, which is more precise than a generic package list.

Verify the manuscript R environment before a complete rebuild:

```powershell
Rscript environment\check_r_packages.R
```

The installation helper installs missing packages from the declared CRAN and CaiLab repositories; the version checker is the release guardrail and will stop on a missing or mismatched package.
