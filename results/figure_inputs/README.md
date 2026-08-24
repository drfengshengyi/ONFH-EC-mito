# Versioned plotting inputs

These compact tables are the machine-readable inputs for `plotting/make_reviewed_figures.py` and `plotting/correct_figure4_fdr.py`.

- `figure2_umap.csv.gz` and `figure2_marker_dotplot.csv` contain the 13,426-cell EC embedding and marker-panel summaries used in Figure 2A-B.
- `ec_subtype_composition_v4.csv`, `module_scores_by_library_v4.csv`, and `module_scores_liao_stats_v4.csv` contain sampling-unit values and exact-permutation FDR values used in Figures 2-4.
- `de_ec_SONFH_vs_HOA_descriptive_v4.csv` and `fig3_key_gene_effects_v4.csv` contain the explicitly descriptive cross-cohort effects used in Figure 3.
- `diag_*` files contain the repeated nested-CV, label-permutation, feature-stability, and comparator results used in Figure 6.

The plotting scripts read statistical labels directly from these files; displayed FDR values are not manually transcribed.
