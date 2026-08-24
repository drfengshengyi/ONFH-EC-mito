# Plotting and figure assembly

- `assemble_manuscript_figures.py` assembles the versioned source panels.
- `make_virtual_knockout_figure.py` generates Figure 5 from the official-R outputs and matched-null post-processing tables.
- `make_evidence_model.py` generates Supplementary Figure S1.
- `make_reviewed_figures.py` rebuilds Figures 2, 3, and 6 from the versioned plotting inputs. It clips percentage error intervals at zero, uses publication-facing subtype labels, marks the SONFH contrast as descriptive, and reads displayed FDR values from the result table.
- `correct_figure4_fdr.py` synchronizes the Figure 4 cGAS-STING title with the same machine-readable FDR table.

Run all plotting steps with:

```powershell
.\workflow\run_figures.ps1 -Python python
```

The scripts write both 300-dpi PNG and PDF versions to `figures/final/`. Panel letters are positioned outside the plotted data region to prevent clipping. Final submitted PDFs are versioned alongside the plotting code.
