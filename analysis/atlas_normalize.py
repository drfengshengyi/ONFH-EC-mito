# -*- coding: utf-8 -*-
# Atlas stage A: load rds -> anndata -> normalize -> HVG -> PCA -> harmony -> neighbors -> save
import time, warnings
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import rdata

warnings.filterwarnings("ignore")
ROOT = str(Path(__file__).resolve().parents[1])
LOG = ROOT + "/analysis/py_atlas.log"

def lg(*a):
    s = time.strftime("%H:%M:%S") + " | " + " ".join(str(x) for x in a)
    print(s, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(s + "\n")

lg("stage A start")
o = rdata.read_rds(ROOT + "/analysis/counts_qc_v4.rds")
lg("rds parsed")
genes = np.asarray(o.Dimnames[0]).astype(str)
cells = np.asarray(o.Dimnames[1]).astype(str)
mat = sp.csc_matrix((np.asarray(o.x, dtype=np.float32),
                     np.asarray(o.i), np.asarray(o.p)), shape=tuple(o.Dim))
lg("matrix:", mat.shape)
del o

meta = rdata.read_rds(ROOT + "/analysis/cell_meta_v4.rds")
lg("meta:", meta.shape, list(meta.columns))

import anndata
import scanpy as sc

X = mat.T.tocsr()  # cells x genes
del mat
adata = anndata.AnnData(X=X)
adata.var_names = genes
adata.obs_names = cells
for c in [
    "sample", "dataset", "group", "participant_id", "inferential_unit",
    "independent_for_inference", "source_note", "nFeature", "nCount",
    "percent.mt", "percent.hb",
]:
    adata.obs[c] = meta[c].values
lg("anndata:", adata.shape)
del X

sc.pp.normalize_total(adata, target_sum=1e4)
lg("normalized")
sc.pp.log1p(adata)
lg("log1p done")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=False)
lg("HVG:", int(adata.var.highly_variable.sum()))
sc.pp.pca(adata, n_comps=30, use_highly_variable=True, svd_solver="arpack")
lg("PCA done")
adata.write(ROOT + "/analysis/atlas_pca.h5ad")
lg("checkpoint saved (atlas_pca.h5ad); DONE stage A (batch correction handled in A2)")
