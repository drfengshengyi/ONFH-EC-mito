# Result-to-code map

The repository versions compact analysis results and the submission figures so
that every displayed value can be traced to a machine-readable source. Primary
public matrices are not redistributed.

| Deliverable | Canonical producer | Versioned direct inputs |
|---|---|---|
| Figure 1 | `plotting/assemble_manuscript_figures.py` | `figures/source/umap_*`, `fig1e_v4_ec_fraction.png`, `fig1f_v4_retained_cells.png` |
| Figure 2 | `plotting/make_reviewed_figures.py` | `results/figure_inputs/figure2_*`, `ec_subtype_composition_v4.csv`, `module_scores_*` |
| Figure 3 | `plotting/make_reviewed_figures.py` | `de_ec_SONFH_vs_HOA_descriptive_v4.csv`, `fig3_key_gene_effects_v4.csv`, `module_scores_*` |
| Figure 4 | `plotting/make_figure4.py` | GSEA, module-score, signed-TF and communication tables in `results/figure_inputs/` |
| Figure 5 | `analysis/spatial_contextualization.py` | GSE284089 spatial matrix, coordinates, deposited H&E image, and frozen gene sets |
| Figure 6 | `plotting/make_genes_virtual_knockout_figure.R` | official-R donor results, matched-gene controls, nuclear refits, and matched-null pathway tables |
| Figure 7 | `plotting/make_figure7_enhanced.R` | frozen plan and formal outputs in `results/serum_enhanced_20260910/` |
| Supplementary Figure S1 | `plotting/make_evidence_model.py` | prespecified evidence model encoded in the script |
| Supplementary Figure S2 | `plotting/make_figure7_enhanced.R` | refit-aware OOB, selected-feature count, matched-space QC, and calibration-bin diagnostics |
| Supplementary Tables S1-S13 | `tools/append_sr_sensitivity_tables.py` and `qa/build_supplementary_table_s13.py` | S10 contains virtual-knockout audits, S11 spatial context, S12 preprocessing/annotation/HOA1 sensitivity, and S13 the enhanced serum audit; the assembled submission workbook is kept local |

`workflow/run_figures.ps1` is the canonical figure entry point. Historical
scripts that are not part of that workflow must not write to `figures/final/`.
The optional `--submission-dir` repository audit requires byte-identical PDF
figures and a byte-identical supplementary workbook.

Large serum checkpoints, the 52.5-MB per-match gene map, rendered figure files, and manuscript/submission artifacts are intentionally not versioned in this release. The frozen design, compact formal tables, plot data, and code are retained; omitted artifacts are reproducible from `workflow/run_serum_enhanced.ps1`.
