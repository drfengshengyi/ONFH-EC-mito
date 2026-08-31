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
| Figure 5 | `plotting/make_virtual_knockout_figure.R` | official-R donor results, `results/official_r_vko_figure_data/`, matched-null pathway tables |
| Figure 6 | `plotting/make_reviewed_figures.py` | `results/figure_inputs/diag_*` |
| Supplementary Figure S1 | `plotting/make_evidence_model.py` | prespecified evidence model encoded in the script |
| Supplementary Tables S1-S10 | versioned workbook | machine-readable analysis tables; S10g-5 is reconciled against the two `official_r_vko_no_mt_*_audit.csv` files by `qa/check_repository.py` |

`workflow/run_figures.ps1` is the canonical figure entry point. Historical
scripts that are not part of that workflow must not write to `figures/final/`.
The optional `--submission-dir` repository audit requires byte-identical PDF
figures and a byte-identical supplementary workbook.
