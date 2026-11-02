"""scripts/build_from_10x_samples.py: per-sample 10x triplets + sample sheet -> one .h5ad."""
from __future__ import annotations

import gzip
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]


def _write_10x(d: Path, n_cells: int, genes: list[str], prefix: str, seed: int):
    d.mkdir(parents=True)
    rng = np.random.default_rng(seed)
    X = sp.random(len(genes), n_cells, density=0.2, random_state=seed, data_rvs=lambda k: rng.poisson(2, k) + 1).tocoo()
    with gzip.open(d / f"{prefix}_matrix.mtx.gz", "wb") as fh:
        scipy.io.mmwrite(fh, X)
    pd.DataFrame({"id": [f"ENSG{i:05d}" for i in range(len(genes))], "sym": genes, "t": "Gene Expression"}).to_csv(d / f"{prefix}_features.tsv.gz", sep="\t", header=False, index=False)
    pd.Series([f"BC{i}-1" for i in range(n_cells)]).to_csv(d / f"{prefix}_barcodes.tsv.gz", sep="\t", header=False, index=False)


def test_build_from_sheet(tmp_path):
    genes = [f"G{i}" for i in range(40)] + ["FOLR1", "EPCAM"]
    rows = []
    for i, (smp, cond) in enumerate([("GSM1", "Tumor"), ("GSM2", "Omentum"), ("GSM3", "Omentum")]):
        d = tmp_path / smp
        _write_10x(d, 250 + 10 * i, genes if i < 2 else genes[:41], smp, i)  # third library lacks EPCAM
        rows.append({"sample": smp, "condition": cond, "patient": f"P{i}", "title": f"P{i}_{cond.lower()}", "dir": str(d)})
    sheet = tmp_path / "samples.csv"
    pd.DataFrame(rows).to_csv(sheet, index=False)
    out = tmp_path / "x.h5ad"
    r = subprocess.run([sys.executable, str(ROOT / "scripts/build_from_10x_samples.py"), "--sheet", str(sheet), "--out", str(out), "--min-cells", "100"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    a = ad.read_h5ad(out)
    assert a.n_obs == 250 + 260 + 270
    assert a.n_vars == 41 and "EPCAM" not in a.var_names and "FOLR1" in a.var_names  # intersection of gene sets
    assert set(a.obs["condition"]) == {"Tumor", "Omentum"} and a.obs["sample"].nunique() == 3
    assert a.obs_names[0].startswith("GSM1_BC")
    assert type(a.X).__name__ == "csr_matrix"


def test_empty_sheet_is_refused(tmp_path):
    sheet = tmp_path / "empty.csv"
    pd.DataFrame(columns=["sample", "condition", "dir"]).to_csv(sheet, index=False)
    r = subprocess.run([sys.executable, str(ROOT / "scripts/build_from_10x_samples.py"), "--sheet", str(sheet), "--out", str(tmp_path / "y.h5ad")], capture_output=True, text=True)
    assert r.returncode == 2 and "empty" in r.stderr
