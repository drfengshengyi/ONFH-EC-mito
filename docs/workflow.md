# Reproducible workflow

```text
public matrices + frozen sampling-unit map
                |
                v
QC -> atlas integration -> cell-type annotation -> endothelial subset
                |
                +-> participant-aware composition / pseudobulk / pathways
                +-> ligand-receptor and signed-regulon summaries
                +-> enhanced serum nested-CV and matched-space calibration
                +-> no-ComBat / annotation-margin / HOA1 sensitivity audits
                |
                v
          source figure panels -> final figure assembly

HOA2 and HOA3 EC count matrices + frozen genes/barcodes
                |
                v
official scTenifoldKnk (donor separated; two parameter profiles)
                |
                v
matched-null audit + plot-data export -> Figure 6
```

## Entry points

| Command | Recalculates models? | Requires primary data? |
|---|---:|---:|
| `workflow/run_core_analysis.ps1` | Yes | Yes |
| `workflow/run_virtual_knockout.ps1` | Yes | HOA2/HOA3 RDS only |
| `workflow/run_figures.ps1` | No, except deterministic plot summaries | No, when versioned panels/results are present |
| `workflow/run_sr_preprocessing_sensitivity.ps1` | Yes | Frozen atlas and EC H5AD objects |
| `workflow/run_serum_enhanced.ps1` | Yes; formal phase is compute-intensive | GSE123568 matrix and SOFT files |
| `workflow/run_all.ps1` | Yes | Yes |

Every entry point resolves paths relative to its own repository location, so the repository can be cloned anywhere. `ONFH_ROOT` is set internally for R and Python modules.

The two Scientific Reports sensitivity workflows are kept outside `run_all.ps1` so a routine figure rebuild cannot trigger long resampling or overwrite the frozen primary analysis.
