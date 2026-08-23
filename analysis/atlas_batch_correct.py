# -*- coding: utf-8 -*-
# Atlas stage A2: from checkpoint -> combat on HVGs -> PCA -> neighbors -> save
import time, warnings
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scanpy as sc

warnings.filterwarnings("ignore")
ROOT = str(Path(__file__).resolve().parents[1])
LOG = ROOT + "/analysis/py_atlas.log"

def lg(*a):
    s = time.strftime("%H:%M:%S") + " | " + " ".join(str(x) for x in a)
    print(s, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(s + "\n")

lg("stage A2 start (combat route)")
adata = sc.read_h5ad(ROOT + "/analysis/atlas_pca.h5ad")
lg("loaded:", adata.shape)

hvg = adata.var.highly_variable.values
sub = adata[:, hvg].copy()
lg("HVG submatrix:", sub.shape)

# combat needs dense input
sub.X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
lg("densified:", round(sub.X.nbytes / 1e9, 2), "GB")

sc.pp.combat(sub, key="sample")
lg("combat done")
if not sub.X.flags.writeable:
    sub.X = np.array(sub.X, dtype=np.float32)
    lg("made X writable")
sub.write(ROOT + "/analysis/atlas_combat_hvg.h5ad")
lg("combat checkpoint saved")

import traceback
try:
    sc.pp.scale(sub, max_value=10)
    lg("scaled")
    sc.pp.pca(sub, n_comps=30, svd_solver="arpack")
    lg("PCA done")
except Exception:
    lg("SCALE/PCA FAILED:\n" + traceback.format_exc())
    raise

adata.obsm["X_pca_combat"] = sub.obsm["X_pca"]
sc.pp.neighbors(adata, use_rep="X_pca_combat", n_neighbors=15, n_pcs=30)
lg("neighbors done")

adata.write(ROOT + "/analysis/atlas_neighbors.h5ad")
lg("saved atlas_neighbors.h5ad; DONE stage A2")
