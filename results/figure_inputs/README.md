# Versioned plotting inputs

These compact tables are the machine-readable inputs for `plotting/make_reviewed_figures.py` and `plotting/make_figure4.py`.

- `figure2_umap.csv.gz` and `figure2_marker_dotplot.csv` contain the 13,426-cell EC embedding and marker-panel summaries used in Figure 2A-B.
- `ec_subtype_composition_v4.csv`, `module_scores_by_library_v4.csv`, and `module_scores_liao_stats_v4.csv` contain sampling-unit values and exact-permutation FDR values used in Figures 2-4.
- `de_ec_SONFH_vs_HOA_descriptive_v4.csv` and `fig3_key_gene_effects_v4.csv` contain the explicitly descriptive cross-cohort effects used in Figure 3.
- `gsea_ONFH3A_vs_HOA_H.csv` is the frozen Hallmark enrichment result used in Figure 4A and in the manuscript text.
- `sample_level_tf_ulm_v4.csv` and `tf_stats_ulm_v4.csv` contain donor-level signed TF activities and exact-permutation FDR values used in Figure 4C.
- `comm_scores_v4_long.csv.gz` and `comm_top_ONFH_4_vs_ONFH_3A_v4.csv` contain the donor-level communication scores and presorted largest effects used in Figure 4D.
- `diag_*` files contain the repeated nested-CV, label-permutation, feature-stability, and comparator results used in Figure 7.

The plotting scripts read statistical labels directly from these files; displayed FDR values are not manually transcribed.
