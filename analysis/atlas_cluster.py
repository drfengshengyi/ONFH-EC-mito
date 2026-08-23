# -*- coding: utf-8 -*-
# Atlas stage B: leiden clustering + UMAP
import time, warnings
from pathlib import Path
import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")
ROOT = str(Path(__file__).resolve().parents[1])
LOG = ROOT + "/analysis/py_atlas.log"

def lg(*a):
    s = time.strftime("%H:%M:%S") + " | " + " ".join(str(x) for x in a)
    print(s, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(s + "\n")

lg("stage B start")
adata = sc.read_h5ad(ROOT + "/analysis/atlas_neighbors.h5ad")
lg("loaded:", adata.shape)

sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)
lg("leiden done:", adata.obs["leiden"].nunique(), "clusters")

adata.write(ROOT + "/analysis/atlas_leiden.h5ad")
lg("leiden checkpoint saved")

sc.tl.umap(adata)
lg("UMAP done")

adata.write(ROOT + "/analysis/atlas_umap.h5ad")
lg("saved atlas_umap.h5ad; DONE stage B")
