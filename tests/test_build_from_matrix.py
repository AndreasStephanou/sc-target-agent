"""scripts/build_from_matrix.py: the common GEO 'one dense matrix + metadata' deposit -> .h5ad."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def test_build_from_matrix_genes_x_cells(tmp_path):
    rng = np.random.default_rng(0)
    cells = [f"S{i // 20 + 1}T_{j}" for i in range(60) for j in [i]]
    genes = [f"G{i}" for i in range(50)] + ["GPC3"]
    counts = pd.DataFrame(rng.poisson(1.0, size=(len(genes), len(cells))), index=genes, columns=cells)
    meta = pd.DataFrame({"Cell": cells, "sample": [c.split("_")[0] for c in cells],
                         "site": ["Tumor" if c.startswith(("S1", "S2")) else "PVTT" for c in cells], "celltype": "Hep"})
    counts.to_csv(tmp_path / "counts.txt.gz", sep="\t")
    meta.to_csv(tmp_path / "meta.txt.gz", sep="\t", index=False)
    out = tmp_path / "x.h5ad"
    r = subprocess.run([sys.executable, str(ROOT / "scripts/build_from_matrix.py"), "--counts", str(tmp_path / "counts.txt.gz"), "--meta", str(tmp_path / "meta.txt.gz"),
                        "--out", str(out), "--cell-col", "Cell", "--sample-col", "sample", "--site-col", "site", "--celltype-col", "celltype", "--keep-site", "Tumor"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "orientation: genes x cells" in r.stdout
    a = ad.read_h5ad(out)
    assert a.n_obs == 40 and a.n_vars == 51
    assert set(a.obs["condition"]) == {"Tumor"} and set(a.obs["sample"]) == {"S1T", "S2T"}
    assert "author_celltype" in a.obs
    assert type(a.X).__name__ == "csr_matrix"
