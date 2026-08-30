# Plotting and figure assembly

- `assemble_manuscript_figures.py` assembles the versioned source panels.
- `make_virtual_knockout_figure.R` generates Figure 5 from the official-R outputs and matched-null post-processing tables. R is the sole drawing and export backend for this figure.
- `make_evidence_model.py` generates Supplementary Figure S1.
- `make_reviewed_figures.py` rebuilds Figures 2, 3, and 6 from the versioned plotting inputs. It clips percentage error intervals at zero, uses publication-facing subtype labels, marks the SONFH contrast as descriptive, and reads displayed FDR values from the result table.
- `make_figure4.py` rebuilds all four Figure 4 panels from versioned GSEA, module-score, TF-activity, and communication tables. Every displayed FDR is table-derived; no post-hoc text overlay is used.
- `correct_figure4_fdr.py` is retained only as a legacy record and is not called by the current workflow.

Run all plotting steps with:

```powershell
.\workflow\run_figures.ps1 -Python python -Rscript Rscript
```

The scripts write both 300-dpi PNG and PDF versions to `figures/final/`. Panel letters are positioned outside the plotted data region to prevent clipping. Final submitted PDFs are versioned alongside the plotting code.
