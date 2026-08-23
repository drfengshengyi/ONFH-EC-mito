# -*- coding: utf-8 -*-
# Atlas stage C: marker-based cluster annotation + global figures
import time, warnings
from pathlib import Path
import numpy as np
import pandas as pd
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

lg("stage C start")
adata = sc.read_h5ad(str(ROOT / "analysis/atlas_umap.h5ad"))
lg("loaded:", adata.shape)

PANELS = {
    "EC": ["PECAM1", "VWF", "CDH5", "EMCN", "KDR"],
    "Pericyte_SMC": ["RGS5", "ACTA2", "MYH11", "PDGFRB"],
    "Osteoblast": ["RUNX2", "SP7", "ALPL", "IBSP", "BGLAP", "COL1A1", "SPP1"],
    "MSC_stromal": ["CXCL12", "NT5E", "THY1", "ENG", "APOD", "CFD"],
    "Adipo_lineage": ["ADIPOQ", "PLIN1", "FABP4", "LPL"],
    "Chondrocyte": ["COL2A1", "ACAN", "SOX9"],
    "Osteoclast": ["ACP5", "CTSK", "MMP9", "TNFRSF11A", "CALCR"],
    "Myeloid": ["CD68", "CD14", "LYZ", "CSF1R", "AIF1"],
    "T_cell": ["CD3D", "CD3E", "CD8A", "CD4"],
    "NK": ["GNLY", "NKG7", "NCAM1"],
    "B_cell": ["CD19", "MS4A1", "CD79A"],
    "Plasma": ["JCHAIN", "MZB1", "XBP1"],
    "Granulocyte": ["S100A8", "S100A9", "FCGR3B", "MPO"],
    "Mast": ["TPSAB1", "TPSB2", "KIT"],
    "Platelet": ["PPBP", "PF4", "GP9"],
    "RBC": ["HBB", "HBA1", "HBA2"],
    "pDC": ["LILRA4", "IRF7", "CLEC4C"],
}
panels = {k: [g for g in v if g in adata.var_names] for k, v in PANELS.items()}
panels = {k: v for k, v in panels.items() if v}

import scipy.sparse as sp
X = adata.X
Xc = X.tocsr()
clusters = adata.obs["leiden"].astype(str).values
cl_ids = sorted(np.unique(clusters), key=int)
gix = {g: i for i, g in enumerate(adata.var_names)}

# mean log-expression per cluster per panel
scores = pd.DataFrame(index=cl_ids, columns=list(panels.keys()), dtype=float)
for cl in cl_ids:
    mask = clusters == cl
    n = mask.sum()
    rowsel = np.where(mask)[0]
    for pname, genes in panels.items():
        cols = [gix[g] for g in genes]
        vals = np.asarray(Xc[rowsel][:, cols].mean(axis=0)).ravel()
        scores.loc[cl, pname] = vals.mean()
lg("panel scores computed")

scores.to_csv(ROOT / "analysis/cluster_marker_scores.csv")
anno = scores.idxmax(axis=1).rename("cell_type")
top2 = scores.apply(lambda r: r.nlargest(2).index.tolist(), axis=1)
anno_df = pd.DataFrame({"cluster": cl_ids,
                        "cell_type": anno.values,
                        "top2_panels": ["/".join(x) for x in top2.values]})
# tie-break: EC vs Pericyte confusion handled by raw EC score
anno_df.to_csv(ROOT / "analysis/cluster_annotation.csv", index=False)
print(anno_df.to_string())

# map annotation to cells
ct_map = dict(zip(anno_df.cluster, anno_df.cell_type))
adata.obs["cell_type"] = [ct_map[c] for c in clusters]

# composition per group
comp = pd.crosstab(adata.obs["group"], adata.obs["cell_type"], normalize="index") * 100
comp.to_csv(ROOT / "analysis/composition_by_group.csv")
comp_s = pd.crosstab(adata.obs["sample"], adata.obs["cell_type"], normalize="index") * 100
comp_s.to_csv(ROOT / "analysis/composition_by_sample.csv")
lg("composition tables saved")

# ---- figures ----
sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False,
           title="Leiden clusters", frameon=False)
plt.savefig(FIGS / "umap_clusters.png", dpi=220, bbox_inches="tight"); plt.close()
lg("fig umap_clusters")

sc.pl.umap(adata, color="group", show=False, frameon=False)
plt.savefig(FIGS / "umap_group.png", dpi=220, bbox_inches="tight"); plt.close()

sc.pl.umap(adata, color="dataset", show=False, frameon=False)
plt.savefig(FIGS / "umap_dataset.png", dpi=220, bbox_inches="tight"); plt.close()

sc.pl.umap(adata, color="cell_type", show=False, frameon=False)
plt.savefig(FIGS / "umap_celltype.png", dpi=220, bbox_inches="tight"); plt.close()

sc.pl.umap(adata, color=["PECAM1", "VWF", "EMCN", "KDR"], show=False, frameon=False)
plt.savefig(FIGS / "umap_ec_markers.png", dpi=220, bbox_inches="tight"); plt.close()
lg("umap figs done")

# composition stacked bar
import seaborn as sns
order = ["Healthy", "HOA", "FNF", "ONFH_3A", "ONFH_4", "SONFH"]
comp_o = comp.reindex([g for g in order if g in comp.index])
fig, ax = plt.subplots(figsize=(9, 5.2))
comp_o.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", width=0.75)
ax.set_ylabel("Composition (%)"); ax.set_xlabel("")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
plt.xticks(rotation=30, ha="right")
fig.savefig(FIGS / "composition_by_group.png", dpi=220, bbox_inches="tight")
plt.close(fig)
lg("composition fig done")

adata.write(str(ROOT / "analysis/atlas_annotated.h5ad"))
lg("saved atlas_annotated.h5ad; DONE stage C")
