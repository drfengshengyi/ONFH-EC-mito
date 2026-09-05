# Plotting

| Output | Canonical producer |
|---|---|
| Figure 1 | `assemble_manuscript_figures.py` |
| Figures 2, 3, and 7 | `make_reviewed_figures.py` |
| Figure 4 | `make_figure4.py` |
| Figure 5 | `analysis/spatial_contextualization.py` |
| Figure 6 | `make_genes_virtual_knockout_figure.R` |
| Supplementary Figure S1 | `make_evidence_model.py` |

Figure 6 uses donor-separated official-R outputs, nuclear-gene ranks, matched-gene controls, and matched-null pathway results. `CALCOCO2` is shown within the comparator distribution as a secondary candidate.

Legacy Figure 6 scripts are retained for provenance. They are not called by the current workflow and write only to legacy locations when applicable.

Run all canonical plotting steps with:

```powershell
.\workflow\run_figures.ps1 -Python python -Rscript Rscript
```

Final PDF and 300-dpi PNG files are written to `figures/final/`.
