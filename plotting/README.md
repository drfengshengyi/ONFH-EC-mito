# Plotting and figure assembly

- `assemble_manuscript_figures.py` assembles Figure 1 from versioned source panels.
- `analysis/spatial_contextualization.py` generates Figure 5 from the GSE284089 spatial matrix, coordinates, and deposited H&E image.
- `make_genes_virtual_knockout_figure.R` generates the Genes-ready Figure 6 from the official-R outputs, donor-separated nuclear-gene ranks, matched-gene controls, and matched-null pathway tables. It highlights `CALCOCO2` inside the comparator distribution to make its post-analysis secondary-candidate status auditable without implying specificity. R is the sole drawing and export backend for this figure.
- `make_virtual_knockout_figure.R` is the earlier Bone Reports/Life layout retained for provenance; it is not called by the canonical plotting workflow.
- `make_virtual_knockout_figure.py` is a legacy diagnostic retained for provenance. It writes only to `figures/legacy/` and cannot overwrite the canonical Figure 6.
- `make_evidence_model.py` generates Supplementary Figure S1.
- `make_reviewed_figures.py` rebuilds Figures 2, 3, and 7 from the versioned plotting inputs. Figure 3 distinguishes the post-analysis `CALCOCO2`/NDP52 secondary-candidate tier from `OPTN` as a mechanistic-context control; it otherwise preserves the frozen paired 15-gene display. The script clips percentage error intervals at zero, uses publication-facing subtype labels, marks the SONFH contrast as descriptive, and reads displayed FDR values from the result table.
- `make_figure4.py` rebuilds all four Figure 4 panels from versioned GSEA, module-score, TF-activity, and communication tables. Every displayed FDR is table-derived; no post-hoc text overlay is used.
- `correct_figure4_fdr.py` is retained only as a legacy record and is not called by the current workflow.

Each final figure has one canonical producer: Figure 1, `assemble_manuscript_figures.py`; Figures 2, 3, and 7, `make_reviewed_figures.py`; Figure 4, `make_figure4.py`; Figure 5, `analysis/spatial_contextualization.py`; Figure 6, `make_genes_virtual_knockout_figure.R`; and Supplementary Figure S1, `make_evidence_model.py`.

Run all plotting steps with:

```powershell
.\workflow\run_figures.ps1 -Python python -Rscript Rscript
```

The scripts write both 300-dpi PNG and PDF versions to `figures/final/`. Panel letters are positioned outside the plotted data region to prevent clipping. Final submitted PDFs are versioned alongside the plotting code.
