# Core analysis

The scripts are ordered by dependency rather than by historical filename. Run the complete sequence with `workflow/run_core_analysis.ps1`.

| Order | Script | Main purpose |
|---:|---|---|
| 1 | `prepare_matrices.R` | Load the four public datasets, harmonize gene symbols, apply cell QC, and attach sampling-unit metadata |
| 2 | `atlas_normalize.py` | Normalize the combined atlas and select variable genes |
| 3 | `atlas_batch_correct.py` | Construct the batch-aware latent representation |
| 4 | `atlas_cluster.py` | Neighborhood graph, clustering, and UMAP |
| 5 | `atlas_annotate.py` | Marker-panel cell-type annotation and atlas panels |
| 6 | `endothelial_subset.py` | Endothelial subset, state scoring, and EC UMAP |
| 7 | `build_genesets.py` | Freeze analysis gene sets against the observed feature space |
| 8 | `endothelial_panels.py` | Endothelial marker/state source panels |
| 9–15 | audit/inference modules | Annotation audit, composition, pseudobulk, enrichment, communication, regulons, and serum classifier |
| 16 | `supplementary_audits.py` | Supporting audit tables |
| 17 | `write_provenance.py` | Checksums and software/provenance record |

`v4_common.py` contains shared metadata, exact-test, multiple-testing, and plotting utilities. Analysis tables and large intermediate objects are written to `analysis/`; final manuscript deliverables are copied or generated under `results/` and `figures/`.

The SONFH entries with missing participant identifiers have `independent_for_inference = FALSE` in `sample_metadata_v4.csv`. Do not change that flag unless a verifiable participant-to-library map becomes available.
