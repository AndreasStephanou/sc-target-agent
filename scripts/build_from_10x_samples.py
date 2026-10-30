"""Build one .h5ad from per-sample 10x matrices listed in a sample sheet.

    python scripts/build_from_10x_samples.py --sheet data/gastric/raw/samples.csv --out data/gastric/gastric.h5ad

The sheet (from fetch_geo_10x.py, or hand-written) has columns: sample, condition, patient (optional),
dir. Each dir holds matrix.mtx(.gz) + features/genes.tsv(.gz) + barcodes.tsv(.gz), with any prefix
(GEO adds the GSM id). Cells are tagged <sample>_<barcode>; obs gets sample, condition, patient and
every other sheet column. Genes are matched by symbol across samples (intersection); X holds raw
integer counts. Prints cells per sample x condition at the end so the labels can be checked first.
"""
from __future__ import annotations

import argparse
import gzip
import io
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def _open(p: Path):
    return gzip.open(p, "rb") if p.suffix == ".gz" else open(p, "rb")


def read_10x_dir(d: Path):
    files = {f.name: f for f in d.iterdir()}
    pick = lambda pat: next((f for n, f in files.items() if re.search(pat, n, re.I)), None)
    mtx, feat, bc = pick(r"matrix\.mtx"), pick(r"(features|genes)\.tsv"), pick(r"barcodes\.tsv")
    if not (mtx and feat and bc):
        raise FileNotFoundError(f"{d}: need matrix.mtx, features/genes.tsv, barcodes.tsv; have {sorted(files)}")
    with _open(mtx) as fh:
        X = scipy.io.mmread(io.BytesIO(fh.read())).tocsr().T.tocsr()  # genes x cells -> cells x genes
    with _open(feat) as fh:
        ft = pd.read_csv(fh, sep="\t", header=None)
    with _open(bc) as fh:
        cells = pd.read_csv(fh, sep="\t", header=None)[0].astype(str).to_numpy()
    symbols = ft[1].astype(str).to_numpy() if ft.shape[1] >= 2 else ft[0].astype(str).to_numpy()
    if ft.shape[1] >= 3:  # keep only Gene Expression features (drop antibody / CRISPR rows)
        keep = ft[2].astype(str).str.contains("Gene Expression", case=False) | ~ft[2].astype(str).str.contains("Antibody|CRISPR|Custom", case=False)
        X, symbols = X[:, keep.to_numpy()], symbols[keep.to_numpy()]
    if X.shape[0] != len(cells) or X.shape[1] != len(symbols):
        raise ValueError(f"{d}: matrix {X.shape} vs {len(cells)} barcodes / {len(symbols)} features")
    return X, symbols, cells


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cells", type=int, default=200, help="skip a library with fewer barcodes than this")
    a = ap.parse_args(argv)
    sheet = pd.read_csv(a.sheet, dtype=str).fillna("")
    for c in ("sample", "condition", "dir"):
        if c not in sheet.columns:
            print(f"sheet needs a '{c}' column", file=sys.stderr)
            return 2
    if (sheet["condition"] == "").any():
        print("some rows have an empty condition; fill the sheet first:\n" + sheet[sheet["condition"] == ""][["sample", "title"]].to_string(index=False) if "title" in sheet else "", file=sys.stderr)
        return 2
    if sheet.empty:
        print("sample sheet is empty: nothing to build (check the --keep/--drop terms of fetch_geo_10x.py)", file=sys.stderr)
        return 2
    if sheet.empty:
        print("sample sheet is empty: nothing to build (check the --keep/--drop terms of fetch_geo_10x.py)", file=sys.stderr)
        return 2
    mats, obs_parts, gene_sets = [], [], []
    for _, r in sheet.iterrows():
        X, genes, cells = read_10x_dir(Path(r["dir"]))
        if X.shape[0] < a.min_cells:
            print(f"  {r['sample']}: {X.shape[0]} barcodes < {a.min_cells}, skipped")
            continue
        # collapse duplicated symbols by summing
        g = pd.Index(genes)
        if g.has_duplicates:
            codes, uniq = pd.factorize(g)
            ind = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
            X, g = (X @ ind).tocsr(), pd.Index(uniq)
        mats.append((X, g))
        gene_sets.append(set(g))
        ob = pd.DataFrame({k: [r[k]] * X.shape[0] for k in sheet.columns if k != "dir"}, index=[f"{r['sample']}_{c}" for c in cells])
        obs_parts.append(ob)
        print(f"  {r['sample']} [{r['condition']}] {X.shape[0]:,} barcodes x {X.shape[1]:,} genes")
    common = pd.Index(sorted(set.intersection(*gene_sets)))
    print(f"\ngenes shared by all libraries: {len(common):,}")
    Xs = []
    for X, g in mats:
        Xs.append(X[:, g.get_indexer(common)])
    X = sp.vstack(Xs).tocsr().astype(np.float32)
    obs = pd.concat(obs_parts)
    frac_int = float(np.mean(np.mod(X.data[:200000], 1) == 0))
    print(f"integer fraction of sampled values: {frac_int:.3f} (should be 1.0)")
    import anndata as ad

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=common.rename("gene")))
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(a.out)
    print(f"\nwrote {a.out}: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print("cells per sample x condition:\n" + pd.crosstab(adata.obs["sample"], adata.obs["condition"]).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
