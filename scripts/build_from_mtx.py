#!/usr/bin/env python3
"""Build an .h5ad from a GEO "one MatrixMarket matrix + barcodes + genes + per-cell metadata" deposit.

The third common GEO layout (after per-sample 10x directories and one dense matrix): a single
`matrix.mtx(.gz)` of genes x cells, `barcodes.tsv(.gz)`, `genes.tsv(.gz)` / `features.tsv(.gz)`, and a
metadata CSV indexed by barcode that carries the sample, condition and (often) the authors'
annotation. Cells absent from the metadata (unfiltered barcodes) are dropped. Generic: no dataset-
specific logic; every choice is an argument and is printed.

    python scripts/build_from_mtx.py --mtx GSE135893_matrix.mtx.gz --barcodes GSE135893_barcodes.tsv.gz \
        --genes GSE135893_genes.tsv.gz --meta GSE135893_IPF_metadata.csv.gz --out data/ipf/ipf.h5ad \
        --sample-col Sample_Name --condition-col Diagnosis --keep IPF --keep Control \
        --library-col orig.ident --celltype-col celltype --extra-col population --extra-col Sample_Source

Run with --dry-run first: it reads only the metadata and prints the per-sample design.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _open(p: Path):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p, "rt")


def _read_tsv_col(p: Path, col: int = 0) -> list[str]:
    out = []
    with _open(p) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            out.append(parts[min(col, len(parts) - 1)])
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mtx", required=True)
    ap.add_argument("--barcodes", required=True)
    ap.add_argument("--genes", required=True, help="genes.tsv / features.tsv; --symbol-col picks the symbol column")
    ap.add_argument("--symbol-col", type=int, default=None, help="0-based column of the gene symbol (default: last column)")
    ap.add_argument("--meta", required=True, help="per-cell metadata CSV; first column (or --cell-col) is the barcode")
    ap.add_argument("--cell-col", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample-col", required=True, help="metadata column that becomes obs['sample'] (the biological replicate)")
    ap.add_argument("--condition-col", required=True, help="metadata column that becomes obs['condition']")
    ap.add_argument("--keep", action="append", default=[], help="condition values to keep (repeatable); default all")
    ap.add_argument("--library-col", default=None, help="metadata column for the library when samples hold several (kept as obs['library'])")
    ap.add_argument("--celltype-col", default=None, help="authors' annotation -> obs['author_celltype']")
    ap.add_argument("--extra-col", action="append", default=[], help="other metadata columns to carry into obs (repeatable)")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    meta = pd.read_csv(a.meta, index_col=0 if a.cell_col is None else None)
    if a.cell_col:
        meta = meta.set_index(a.cell_col)
    meta.index = meta.index.astype(str)
    for c in [a.sample_col, a.condition_col, a.library_col, a.celltype_col] + a.extra_col:
        if c and c not in meta.columns:
            print(f"metadata has no column {c!r}; columns: {list(meta.columns)}", file=sys.stderr)
            return 2
    print(f"metadata: {meta.shape[0]} cells, columns {list(meta.columns)}")
    keep = meta if not a.keep else meta[meta[a.condition_col].astype(str).isin(a.keep)]
    dropped = meta.shape[0] - keep.shape[0]
    design = keep.groupby([a.condition_col, a.sample_col], observed=True).size()
    print(f"\nkeeping conditions {a.keep or 'all'}: {keep.shape[0]} cells ({dropped} cells of other conditions dropped)")
    print("cells per sample:")
    print(design.to_string())
    print("\nsamples per condition:", keep.groupby(a.condition_col, observed=True)[a.sample_col].nunique().to_dict())
    if a.library_col:
        libs = keep.groupby(a.sample_col, observed=True)[a.library_col].nunique()
        print(f"libraries per sample: max {int(libs.max())}; samples with >1 library: {int((libs > 1).sum())}")
    if a.dry_run:
        return 0

    import scipy.io
    import scipy.sparse as sp
    import anndata as ad

    barcodes = _read_tsv_col(Path(a.barcodes), 0)
    with _open(Path(a.genes)) as fh:
        first = fh.readline().rstrip("\n").split("\t")
    sym_col = a.symbol_col if a.symbol_col is not None else len(first) - 1
    symbols = _read_tsv_col(Path(a.genes), sym_col)
    ids = _read_tsv_col(Path(a.genes), 0) if len(first) > 1 else None
    print(f"\nbarcodes: {len(barcodes)}; genes: {len(symbols)} (symbol column {sym_col})")

    print("reading matrix (this is the slow step) ...", flush=True)
    with _open(Path(a.mtx)) as fh:
        m = scipy.io.mmread(fh)
    print(f"matrix shape {m.shape}, nnz {m.nnz}, dtype {m.dtype}", flush=True)
    if m.shape[0] == len(symbols) and m.shape[1] == len(barcodes):
        genes_by_cells = True
    elif m.shape[1] == len(symbols) and m.shape[0] == len(barcodes):
        genes_by_cells = False
    else:
        print("matrix shape matches neither genes x cells nor cells x genes", file=sys.stderr)
        return 2
    m = sp.csc_matrix(m) if genes_by_cells else sp.csr_matrix(m)

    bc_index = pd.Index(barcodes)
    present = keep.index[keep.index.isin(bc_index)]
    missing = keep.shape[0] - len(present)
    if missing:
        print(f"WARNING: {missing} metadata cells are not in the barcode file and are dropped")
    pos = bc_index.get_indexer(present)
    X = (m[:, pos].T if genes_by_cells else m[pos, :]).tocsr()
    if not np.allclose(X.data, np.round(X.data)):
        print("WARNING: matrix values are not integers; this is not a raw count matrix")
    X.data = X.data.astype(np.int32) if X.data.max() < 2**31 else X.data
    X = sp.csr_matrix(X, dtype=np.int32)

    obs = pd.DataFrame(index=present.astype(str))
    obs["sample"] = keep.loc[present, a.sample_col].astype(str).to_numpy()
    obs["condition"] = keep.loc[present, a.condition_col].astype(str).to_numpy()
    if a.library_col:
        obs["library"] = keep.loc[present, a.library_col].astype(str).to_numpy()
    if a.celltype_col:
        obs["author_celltype"] = keep.loc[present, a.celltype_col].astype(str).to_numpy()
    for c in a.extra_col:
        obs[c] = keep.loc[present, c].to_numpy()
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = obs[c].astype("category")
    var = pd.DataFrame(index=pd.Index(symbols).astype(str))
    if ids is not None:
        var["gene_id"] = ids
    var.index = ad.utils.make_index_unique(var.index) if hasattr(ad.utils, "make_index_unique") else var.index
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    adata.uns["build"] = {"builder": "build_from_mtx.py", "mtx": str(a.mtx), "meta": str(a.meta), "keep": a.keep,
                          "sample_col": a.sample_col, "condition_col": a.condition_col, "library_col": a.library_col,
                          "celltype_col": a.celltype_col, "cells_dropped_other_conditions": int(dropped),
                          "cells_missing_from_barcodes": int(missing)}
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"\nwrote {out}: {adata.n_obs} cells x {adata.n_vars} genes; obs columns {list(adata.obs.columns)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
