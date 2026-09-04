# Official-R virtual-knockout result audit

Base-profile run date: 2026-08-22. Calibration run date: 2026-08-29.

Software: R 4.6.1, scTenifoldKnk 1.1, scTenifoldNet 1.4.

## Why SQSTM1 was perturbed

SQSTM1 was fixed before any perturbation output was inspected. It belonged to
the prespecified selective-clearance panel, was detected in 70.8% of HOA2 ECs
and 72.2% of HOA3 ECs, and had concordant negative effect estimates in the
descriptive SONFH-versus-HOA contrast (-1.265 log2 fold change) and the Liao
ARCO 3A-versus-HOA contrast (-0.450; not FDR significant). These criteria
support an exploratory loss hypothesis, not a causal target claim. The complete
decision trail is in `results/official_r_vko_matched_controls/sqstm1_selection_rationale.csv`
and Supplementary Table S10f.

## Original 300-gene profiles

### Manuscript profile

Parameters: 20 networks, 500 cells/network, 3 PCs, q=0.95, tensor rank 3,
maximum 500 tensor iterations, and 30 manifold dimensions.

- SQSTM1 ranked first in HOA2 and HOA3, confirming perturbation encoding.
- HOA2 downstream FDR < 0.05: MT-ND1, MT-ATP6, MT-ND4, and MT-CO2.
- HOA3 downstream FDR < 0.05: VWF and MT-ND1.
- MT-ND1 was the only downstream gene significant in both donors.
- Cross-donor rank Spearman rho was 0.4716 (p=5.09e-18).

### Package-default profile

Parameters: 10 networks, 500 cells/network, 3 PCs, q=0.90, tensor rank 3,
maximum 1,000 tensor iterations, and 2 manifold dimensions.

- SQSTM1 again ranked first in both donors.
- HOA2 downstream FDR < 0.05: MT-ATP6 and MT-ND4.
- HOA3 downstream FDR < 0.05: MT-ND1.
- No downstream gene was FDR significant in both donors.
- Cross-donor rank Spearman rho was 0.5365 (p=9.23e-24).

Within-donor rank agreement between the two original profiles was 0.7587 for
HOA2 and 0.8984 for HOA3.

## Complete mtDNA-feature-exclusion refits

The five mtDNA-encoded frozen features (MT-ATP6, MT-CO1, MT-CO2, MT-ND1, and
MT-ND4) were removed before network construction. Each donor/profile network
was rebuilt on 295 genes without replacement features. Original and refitted
models were compared on the same 294 non-target nuclear-encoded genes.

| Profile | Donor | Original-vs-refit rho | Top-20 overlap | Top-50 overlap | Original nuclear FDR count | Refit nuclear FDR count | FDR overlap |
|---|---|---:|---:|---:|---:|---:|---:|
| Manuscript | HOA2 | 0.609 | 10/20 | 27/50 | 0 | 0 | 0 |
| Manuscript | HOA3 | 0.751 | 16/20 | 39/50 | 1 | 3 | 1 |
| Package default | HOA2 | 0.738 | 10/20 | 32/50 | 0 | 0 | 0 |
| Package default | HOA3 | 0.782 | 16/20 | 37/50 | 0 | 0 | 0 |

In the manuscript-profile refit, HOA2 had no downstream FDR hit; HOA3 retained
VWF and added NFKBIA and C7, so none replicated across donors. Neither donor
had a downstream FDR hit under package defaults. Cross-donor common-nuclear
rank correlations were 0.555 and 0.632 for the manuscript and package-default
profiles, respectively.

The exact cross-donor coefficients, P values, overlap counts, FDR counts, and
gene identities are frozen in `results/official_r_vko_no_mt_cross_donor_audit.csv`.
The VWF/NFKBIA/C7 row-level statistics are frozen in
`results/official_r_vko_no_mt_hoa3_fdr_audit.csv`. That file keeps both the
official-output BH family (295 rows, including SQSTM1) and the common-nuclear
BH family (294 shared non-target nuclear genes) explicit.

EC-inflammation enrichment remained confined to HOA3 after refitting
(manuscript FDR=0.0320; package-default FDR=0.0115). HOA2 FDR values were 0.448
and 0.590, and no pathway reached FDR < 0.05 in both donors. The common-universe
rank, pathway, and provenance tables are under `results/official_r_vko_no_mt_*`
and `results/vko_mt_encoded_exclusion_manifest.csv`.

## Twenty matched-comparator perturbations

Twenty comparator genes were fixed using only control-EC expression,
detection prevalence, WT out-degree, and absolute WT out-strength. Disease
effects and perturbation outputs did not enter matching. Each comparator was
perturbed under both profiles in both donors. They are computational
comparators, not validated biological negative controls.

- Manuscript profile: SQSTM1 non-target cross-donor rho=0.466 versus comparator
  median 0.393 (range 0.263-0.569; empirical p=0.333).
- Manuscript profile: SQSTM1 top-20 overlap fraction=0.55 versus comparator
  median 0.525 (empirical p=0.524).
- Package default: SQSTM1 rho=0.532 versus comparator median 0.434
  (range 0.317-0.633; empirical p=0.238).
- Package default: SQSTM1 top-20 overlap fraction=0.65 versus comparator median
  0.425 (empirical p=0.0952).
- Downstream FDR-count endpoints did not show SQSTM1 specificity (all empirical
  p>=0.476 in the manuscript profile and p>=0.952 under package defaults).
- Re-running SQSTM1 from each frozen WT network reproduced all 300 ranks exactly
  with a maximum absolute p-value difference of zero.

`CALCOCO2`/NDP52 is one of the 20 matched comparator perturbations. Its
cross-donor correlations were 0.346 in the manuscript profile and 0.317 under
package defaults, with no cross-donor replicated downstream FDR gene. It is
highlighted in Figure 6 as a transcript-supported, post-analysis secondary candidate inside
the calibration distribution, not as a validated second target. `OPTN` is
retained only as a mechanistic-context control: its descriptive SONFH-versus-HOA
effect was positive (+0.469 log2 fold change), and its ranks after `SQSTM1`
perturbation were 227 in HOA2 and 165 in HOA3 in the mtDNA-feature-excluded
manuscript profile. The gene-level role and source audit is
`results/selective_autophagy_receptor_audit.csv` and Supplementary Table S1b.

## Interpretation

The official implementation yields donor-correlated computational rankings,
but the calibration analyses do not identify a donor-replicated nuclear-encoded
downstream hit or establish SQSTM1 as an exceptional perturbation target. The
combined interpretation is therefore a heterogeneous selective-autophagy
receptor context, not a single-gene mechanism. The leading mtDNA-encoded ranks
are sensitive to the modeled feature space, and the HOA3 EC-inflammation signal
does not replicate in HOA2. These outputs are a transparent
candidate-prioritization audit, not experimental knockout evidence, a cell-fate
trajectory simulation, or causal validation.
