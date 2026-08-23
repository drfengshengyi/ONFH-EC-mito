# -*- coding: utf-8 -*-
# Stage 5: EC subset re-clustering + subtype annotation + module scores preview
import time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "analysis/py_atlas.log"
FIGS = ROOT / "figures" / "source"
FIGS.mkdir(parents=True, exist_ok=True)

def lg(*a):
    s = time.strftime("%H:%M:%S") + " | " + " ".join(str(x) for x in a)
    print(s, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(s + "\n")

from v4_common import setup_plot
setup_plot()
import matplotlib.pyplot as plt

lg("stage 5 start")
adata = sc.read_h5ad(str(ROOT / "analysis/atlas_annotated.h5ad"))
ec = adata[adata.obs.cell_type == "EC"].copy()
lg("EC cells:", ec.shape[0])
print(pd.crosstab(ec.obs.group, ec.obs["sample"]))

# re-embed EC
sc.pp.highly_variable_genes(ec, n_top_genes=1500, flavor="seurat", subset=False)
hvg = ec.var.highly_variable.values
sub = ec[:, hvg].copy()
sub.X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
sc.pp.combat(sub, key="sample")
if not sub.X.flags.writeable:
    sub.X = np.array(sub.X, dtype=np.float32)
sc.pp.scale(sub, max_value=10)
sc.pp.pca(sub, n_comps=20, svd_solver="arpack")
lg("EC combat+PCA done")

ec.obsm["X_pca_combat"] = sub.obsm["X_pca"]
sc.pp.neighbors(ec, use_rep="X_pca_combat", n_neighbors=15, n_pcs=20)
sc.tl.leiden(ec, resolution=0.4, flavor="igraph", n_iterations=2, directed=False)
lg("EC leiden:", ec.obs.leiden.nunique(), "clusters")
sc.tl.umap(ec)
lg("EC UMAP done")

# ---- subtype markers ----
SUBTYPE = {
    "typeH_EMCN_KDR": ["EMCN", "KDR", "CDH5", "DACH1", "VWF"],
    "arterial": ["SOX17", "HEY1", "GJA5", "DLL4", "EFNB2", "SOX13"],
    "venous_ACKR1": ["ACKR1", "SELE", "VCAM1", "PLAT"],
    "tip_angiogenic": ["ESM1", "ANGPT2", "APLN", "DLL4"],
    "lymphatic": ["PROX1", "LYVE1", "PDPN"],
    "typeR_bone_remodel": ["SMAD1", "PPARG", "NOTCH4", "COL4A1"],
}
SUBTYPE = {k: [g for g in v if g in ec.var_names] for k, v in SUBTYPE.items()}

# module scores (paper preview)
GENESETS = {
    "Mito_fission": ["DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "MTFR1L", "GIGYF1"],
    "Mito_fusion": ["OPA1", "MFN1", "MFN2"],
    "Mitophagy": ["PINK1", "PRKN", "BNIP3", "BNIP3L", "FUNDC1", "MAP1LC3B", "SQSTM1", "VCP", "ATG5", "ATG7"],
    "cGAS_STING": ["MB21D1", "TMEM173", "TBK1", "IRF3", "NFKB1", "CXCL10", "CCL5", "IFNB1"],
    "OXPHOS": ["NDUFA1", "NDUFB8", "UQCRB", "COX5A", "COX7A2", "ATP5F1A", "ATP5F1B", "MT-ND1", "MT-CO1", "MT-ATP6"],
    "FAO": ["CPT1A", "CPT2", "ACADM", "ACADVL", "HADHA", "HADHB", "PDK4", "SLC27A4"],
    "Angiogenesis": ["KDR", "FLT1", "ENG", "VWF", "EMCN", "ANGPT1", "TEK", "ESM1"],
}
for name, gs in {**GENESETS}.items():
    g2 = [g for g in gs if g in ec.var_names]
    sc.tl.score_genes(ec, g2, score_name=name)
lg("module scores done")

# subtype assignment by panel mean per cluster
cl_ids = sorted(ec.obs.leiden.unique(), key=int)
gix = {g: i for i, g in enumerate(ec.var_names)}
Xc = ec.X.tocsr() if sp.issparse(ec.X) else sp.csr_matrix(ec.X)
clusters = ec.obs.leiden.values
scores = pd.DataFrame(index=cl_ids, columns=list(SUBTYPE.keys()), dtype=float)
for cl in cl_ids:
    rsel = np.where(clusters == cl)[0]
    for pname, genes in SUBTYPE.items():
        if not genes:
            continue
        cols = [gix[g] for g in genes]
        scores.loc[cl, pname] = np.asarray(Xc[rsel][:, cols].mean(axis=0)).ravel().mean()
scores.to_csv(ROOT / "analysis/ec_subtype_scores.csv")
assign = scores.idxmax(axis=1)
ec.obs["EC_subtype"] = [assign[c] for c in clusters]
print(scores.round(3).to_string())
print(pd.crosstab(ec.obs.EC_subtype, ec.obs.group))

# ---- figures ----
sc.pl.umap(ec, color="EC_subtype", show=False, frameon=False)
plt.savefig(FIGS / "ec_umap_subtype.png", dpi=220, bbox_inches="tight"); plt.close()
sc.pl.umap(ec, color="group", show=False, frameon=False)
plt.savefig(FIGS / "ec_umap_group.png", dpi=220, bbox_inches="tight"); plt.close()
sc.pl.umap(ec, color="leiden", legend_loc="on data", show=False, frameon=False)
plt.savefig(FIGS / "ec_umap_leiden.png", dpi=220, bbox_inches="tight"); plt.close()
sc.pl.umap(ec, color=["EMCN", "KDR", "ACKR1", "SOX17"], show=False, frameon=False)
plt.savefig(FIGS / "ec_umap_markers.png", dpi=220, bbox_inches="tight"); plt.close()

# score boxplots by group
import seaborn as sns
order = [g for g in ["Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH"] if g in set(ec.obs.group)]
scoredf = ec.obs[["group", "sample", "EC_subtype"] + list(GENESETS.keys())].copy()
fig, axes = plt.subplots(2, 4, figsize=(16, 7.5))
for ax, gs in zip(axes.ravel(), GENESETS.keys()):
    sns.boxplot(data=scoredf, x="group", y=gs, order=order, ax=ax,
                showfliers=False, palette="Set2", width=0.6)
    ax.set_title(gs); ax.set_xlabel(""); ax.tick_params(axis="x", rotation=30)
for ax in axes.ravel()[len(GENESETS):]:
    ax.axis("off")
fig.suptitle("EC module scores by group (preview)")
fig.tight_layout()
fig.savefig(FIGS / "ec_module_scores_by_group.png", dpi=220, bbox_inches="tight")
plt.close(fig)
lg("EC figs done")

# EC subtype composition by group
comp = pd.crosstab(ec.obs.group, ec.obs.EC_subtype, normalize="index") * 100
comp.to_csv(ROOT / "analysis/ec_subtype_composition.csv")
fig, ax = plt.subplots(figsize=(8, 5))
comp.reindex(order).plot(kind="bar", stacked=True, ax=ax, colormap="Set3", width=0.75)
ax.set_ylabel("EC subtype composition (%)"); ax.set_xlabel("")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
plt.xticks(rotation=30, ha="right")
fig.savefig(FIGS / "ec_subtype_composition.png", dpi=220, bbox_inches="tight")
plt.close(fig)

ec.write(str(ROOT / "analysis/ec_annotated.h5ad"))
lg("saved ec_annotated.h5ad; DONE stage 5")
