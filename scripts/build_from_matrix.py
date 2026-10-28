"""Build an .h5ad from a dense/TSV count matrix plus a per-cell metadata table (the common GEO "one big
matrix" deposit). Generic: discovers orientation, joins on the cell id, keeps chosen sites.

    python build_from_matrix.py --counts GSE149614_HCC.scRNAseq.S71915.count.txt.gz \
        --meta GSE149614_HCC.metadata.updated.txt.gz --out data/hcc/hcc.h5ad \
        --cell-col Cell --sample-col sample --site-col site --keep-site Tumor --keep-site Normal \
        --celltype-col celltype

Prints what it found (shapes, orientation, sites, cells per sample) so the labels can be checked
before anything else. Run with --dry-run first to see the metadata columns and site values.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cell-col", default=None, help="metadata column with the cell id (default: first column / index)")
    ap.add_argument("--sample-col", required=True)
    ap.add_argument("--site-col", required=True, help="metadata column with tissue/condition")
    ap.add_argument("--keep-site", action="append", default=[], help="site values to keep (repeatable); default all")
    ap.add_argument("--celltype-col", default=None)
    ap.add_argument("--dry-run", action="store_true", help="read metadata only and print columns/values")
    a = ap.parse_args(argv)

    meta = pd.read_csv(a.meta, sep=None, engine="python", index_col=None)
    if a.cell_col and a.cell_col in meta.columns:
        meta = meta.set_index(a.cell_col)
    else:
        meta = meta.set_index(meta.columns[0])
    meta.index = meta.index.astype(str)
    print(f"metadata: {meta.shape[0]} cells x {meta.shape[1]} columns")
    print("columns:", list(meta.columns))
    for c in (a.sample_col, a.site_col, a.celltype_col):
        if c and c in meta.columns:
            print(f"\n{c}:\n{meta[c].value_counts().to_string()}")
        elif c:
            print(f"\nWARNING: column {c!r} not in metadata")
    if a.sample_col in meta.columns and a.site_col in meta.columns:
        print("\ncells per sample x site:\n" + pd.crosstab(meta[a.sample_col], meta[a.site_col]).to_string())
    if a.dry_run:
        return 0

    if a.keep_site:
        keep = meta[a.site_col].astype(str).isin(a.keep_site)
        print(f"\nkeeping sites {a.keep_site}: {int(keep.sum())} of {len(meta)} cells")
        meta = meta[keep]

    print("\nreading counts (this can take several minutes for a dense text matrix)...")
    counts = pd.read_csv(a.counts, sep=None, engine="python", index_col=0)
    print(f"counts table: {counts.shape[0]} rows x {counts.shape[1]} columns")
    rows, cols = counts.index.astype(str), counts.columns.astype(str)
    if len(set(cols) & set(meta.index)) >= len(set(rows) & set(meta.index)):
        print("orientation: genes x cells")
        genes, cells, X = rows, cols, counts.values.T
    else:
        print("orientation: cells x genes")
        genes, cells, X = cols, rows, counts.values
    cells = pd.Index(cells)
    common = cells.intersection(meta.index)
    print(f"cells in both counts and metadata: {len(common)} (counts {len(cells)}, metadata {len(meta)})")
    if len(common) == 0:
        print("no overlap between count-matrix cell ids and metadata cell ids: check --cell-col", file=sys.stderr)
        return 2
    pos = cells.get_indexer(common)
    X = X[pos]
    frac_int = float(np.mean(np.mod(X[:2000].ravel(), 1) == 0))
    print(f"integer fraction of a sample of values: {frac_int:.3f} (should be 1.0 for raw counts)")

    import anndata as ad

    adata = ad.AnnData(X=sp.csr_matrix(X.astype(np.float32)), obs=meta.loc[common].copy(), var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.var_names_make_unique()
    adata.obs["sample"] = adata.obs[a.sample_col].astype(str)
    adata.obs["condition"] = adata.obs[a.site_col].astype(str)
    if a.celltype_col and a.celltype_col in adata.obs.columns:
        adata.obs["author_celltype"] = adata.obs[a.celltype_col].astype(str)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(a.out)
    print(f"\nwrote {a.out}: {adata.n_obs} cells x {adata.n_vars} genes")
    print("samples x condition:\n" + pd.crosstab(adata.obs["sample"], adata.obs["condition"]).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
