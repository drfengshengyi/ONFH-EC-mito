# -*- coding: utf-8 -*-
"""Write a compact machine-readable provenance manifest for this release."""
from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

from v4_common import ANALYSIS, ROOT


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


scripts = sorted(
    path
    for directory in ("analysis", "virtual_knockout", "plotting")
    for pattern in ("*.py", "*.R")
    for path in (ROOT / directory).glob(pattern)
)
key_inputs = [
    ANALYSIS / "atlas_annotated.h5ad",
    ANALYSIS / "ec_final.h5ad",
    ANALYSIS / "sample_metadata_v4.csv",
    ROOT / "data" / "cellchatdb_interactions.csv",
    ROOT / "data" / "cellchatdb_complex_named.csv",
    ROOT / "data" / "dorothea_ABC.tsv",
]
manifest = {
    "created_utc": datetime.now(timezone.utc).isoformat(),
    "python": sys.version,
    "platform": platform.platform(),
    "project_root": ".",
    "scripts": {str(p.relative_to(ROOT)): sha256(p) for p in scripts},
    "key_inputs": {
        str(p.relative_to(ROOT)): {"bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in key_inputs if p.exists()
    },
}
(ANALYSIS / "provenance_v4.json").write_text(
    json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
)
print("Wrote analysis/provenance_v4.json")
