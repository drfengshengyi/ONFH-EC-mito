# Official-R virtual-knockout result audit

Run date: 2026-08-22

Software: R 4.6.1, scTenifoldKnk 1.1, scTenifoldNet 1.4.

## Manuscript-matched profile

Parameters: 20 networks, 500 cells/network, 3 PCs, q=0.95, tensor rank 3,
maximum 500 tensor iterations, and 30 manifold dimensions.

- SQSTM1 ranked first in HOA2 and HOA3, confirming perturbation encoding.
- HOA2 downstream FDR < 0.05: MT-ND1, MT-ATP6, MT-ND4, and MT-CO2.
- HOA3 downstream FDR < 0.05: VWF and MT-ND1.
- MT-ND1 was the only downstream gene significant in both donors.
- Cross-donor rank Spearman rho was 0.4716 (p=5.09e-18).
- R-versus-Python rank Spearman rho was 0.5950 in both donors.

## Official-default sensitivity profile

Parameters: 10 networks, 500 cells/network, 3 PCs, q=0.90, tensor rank 3,
maximum 1,000 tensor iterations, and 2 manifold dimensions.

- SQSTM1 again ranked first in both donors.
- HOA2 downstream FDR < 0.05: MT-ATP6 and MT-ND4.
- HOA3 downstream FDR < 0.05: MT-ND1.
- No downstream gene was FDR significant in both donors.
- Cross-donor rank Spearman rho was 0.5365 (p=9.23e-24).
- R-versus-Python rank Spearman rho was 0.5860 in HOA2 and 0.6390 in HOA3.

## Cross-parameter robustness

Within-donor rank agreement between the two official-R profiles was 0.7587 for
HOA2 and 0.8984 for HOA3. SQSTM1, MT-ND1, MT-ATP6, MT-ND4, MT-CO2, VWF,
NFKBIA, ICAM1, ENG, and MCL1 remained among the leading cross-donor candidates.

## Interpretation

The official implementation supports a stable mitochondrial/endothelial
candidate-ranking pattern after SQSTM1 network perturbation. It does not support
a parameter-robust claim that a particular downstream gene is replicated at
FDR < 0.05 across both donors. MT-ND1 satisfies that criterion only in the
manuscript-matched profile and loses HOA2 FDR significance under the official
default profile.

These results should be reported as exploratory, parameter-sensitive network
predictions. They do not constitute experimental or causal validation of an
SQSTM1 mechanism.
