# Plotting and figure assembly

- `assemble_manuscript_figures.py` assembles Figures 1–4 and 6 from versioned source panels.
- `make_virtual_knockout_figure.py` generates Figure 5 from the official-R outputs and matched-null post-processing tables.
- `make_evidence_model.py` generates Supplementary Figure S1.

Run all plotting steps with:

```powershell
.\workflow\run_figures.ps1 -Python python
```

The scripts write both 300-dpi PNG and vector/container PDF versions to `figures/final/`. Panel letters are added during assembly, outside the source-panel content, to prevent clipping. Final submitted PDFs are versioned alongside the plotting code.
