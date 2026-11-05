#!/usr/bin/env python3
"""Per-library marker check of the condition labels, before any pipeline run.

For each library: the fraction of cells with at least --min-count counts of each gene, grouped
by the deposited condition label. Then, per gene, each library is compared with the median of
its own condition and of the opposite condition; a library closer to the opposite condition is
listed as "resembles the opposite condition". This is the check that found swapped GEO labels
in ccRCC (CA9), a "normal" liver library expressing the tumour programme (GPC3), and tumour
cells in "distant omentum" libraries (EPCAM, PAX8).

A second rule, independent of the gap between condition medians, lists a library that carries the
other condition's marker far above the rest of its own condition (default: >= 3x its condition's
median and >= 2 points above it). It is the rule that matters when the case cells are a minority
of every library (tumour nests in omentum; a disease state inside a tissue whose composition is
otherwise kept), where the first rule reports every marker as uninformative.

The genes are the researcher's choice (a marker they trust for the case condition, ideally not
the expected target itself, so the check stays independent of the test). The matrix is read in
row chunks; the whole matrix is never loaded.

Usage:
  python label_check.py data/my.h5ad --sample-key library --condition-key tissue \
      --genes CA9 EPCAM [--layer counts] [--var-column gene_symbols] [--min-count 1] \
      [--celltype-key author_celltype --celltype Hepatocyte] [--csv label_check.csv]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _read_obs_var(path: Path):
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    obs, var = a.obs.copy(), a.var.copy()
    try:
        a.file.close()
    except Exception:
        pass
    return obs, var


def _column_values(h5: h5py.File, key: str, cols: np.ndarray, n_rows: int, chunk: int = 20_000) -> np.ndarray:
    """Return a dense (n_rows x len(cols)) array of the requested columns, reading in chunks."""
    node = h5[key]
    out = np.zeros((n_rows, len(cols)), dtype=np.float32)
    order = np.argsort(cols)
    cols_sorted = cols[order]
    if isinstance(node, h5py.Dataset):  # dense
        for s in range(0, n_rows, chunk):
            e = min(n_rows, s + chunk)
            block = node[s:e][:, cols_sorted]
            out[s:e, order] = block
        return out
    enc = node.attrs.get("encoding-type", b"csr_matrix")
    enc = enc.decode() if isinstance(enc, bytes) else str(enc)
    shape = tuple(int(x) for x in node.attrs.get("shape", (n_rows, 0)))
    data, indices, indptr = node["data"], node["indices"], node["indptr"][:]
    if enc.startswith("csr"):
        for s in range(0, n_rows, chunk):
            e = min(n_rows, s + chunk)
            lo, hi = int(indptr[s]), int(indptr[e])
            if hi == lo:
                continue
            m = sp.csr_matrix((data[lo:hi], indices[lo:hi], indptr[s : e + 1] - lo), shape=(e - s, shape[1]))
            out[s:e] = m[:, cols].toarray()
        return out
    if enc.startswith("csc"):
        for j, c in enumerate(cols):
            lo, hi = int(indptr[c]), int(indptr[c + 1])
            if hi > lo:
                out[indices[lo:hi], j] = data[lo:hi]
        return out
    raise SystemExit(f"unsupported matrix encoding {enc!r} at {key}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad")
    ap.add_argument("--sample-key", required=True)
    ap.add_argument("--condition-key", required=True)
    ap.add_argument("--genes", nargs="+", required=True, help="gene symbols to check")
    ap.add_argument("--layer", default=None, help="layer holding counts (default: X)")
    ap.add_argument("--var-column", default=None, help="var column with gene symbols, if var_names are IDs")
    ap.add_argument("--min-count", type=float, default=1.0, help="a cell is positive at >= this many counts")
    ap.add_argument("--celltype-key", default=None, help="restrict to cells of one annotated type (optional)")
    ap.add_argument("--celltype", default=None)
    ap.add_argument("--min-gap", type=float, default=0.20,
                    help="only compare libraries for a gene whose condition medians differ by at least this fraction of cells (default 0.20); a weak marker says nothing about labels")
    ap.add_argument("--outlier-fold", type=float, default=3.0,
                    help="second rule: list a library whose fraction is >= this many times its own condition's median for a marker that is higher in the other condition (default 3)")
    ap.add_argument("--outlier-gap", type=float, default=0.02,
                    help="second rule: and at least this fraction of cells above its own condition's median (default 0.02)")
    ap.add_argument("--csv", default=None, help="write the per-library table here")
    a = ap.parse_args(argv)

    path = Path(a.h5ad)
    obs, var = _read_obs_var(path)
    for k in (a.sample_key, a.condition_key):
        if k not in obs.columns:
            print(f"obs has no column {k!r}. Columns: {list(obs.columns)}", file=sys.stderr)
            return 2
    symbols = var[a.var_column].astype(str) if a.var_column else pd.Series(var.index.astype(str), index=var.index)
    sym_to_idx = {}
    for i, s in enumerate(symbols.to_numpy()):
        sym_to_idx.setdefault(s, i)
    missing = [g for g in a.genes if g not in sym_to_idx]
    genes = [g for g in a.genes if g in sym_to_idx]
    if missing:
        print(f"not in var: {missing} (check --var-column, species, or symbol case)")
    if not genes:
        return 2
    cols = np.array([sym_to_idx[g] for g in genes], dtype=int)

    key = f"layers/{a.layer}" if a.layer else "X"
    with h5py.File(path, "r") as h5:
        if key not in h5:
            print(f"{key} not in file; layers: {list(h5['layers'].keys()) if 'layers' in h5 else []}", file=sys.stderr)
            return 2
        vals = _column_values(h5, key, cols, obs.shape[0])
    sample_vals = vals[: min(vals.shape[0], 50_000)]
    if sample_vals.size and not np.all(np.abs(sample_vals - np.round(sample_vals)) < 1e-6):
        print("note: values are not integers; this looks normalised, not raw counts. Fractions are still comparable within the file, but use --layer <counts> if one exists.")

    mask = np.ones(obs.shape[0], dtype=bool)
    if a.celltype_key:
        if a.celltype_key not in obs.columns or a.celltype is None:
            print("--celltype-key needs --celltype and the column must exist", file=sys.stderr)
            return 2
        mask = obs[a.celltype_key].astype(str).to_numpy() == str(a.celltype)
        print(f"restricted to {int(mask.sum())} cells with {a.celltype_key} == {a.celltype!r}")
    pos = (vals >= a.min_count)
    df = pd.DataFrame(pos[mask], columns=genes)
    df[a.sample_key] = obs[a.sample_key].astype(str).to_numpy()[mask]
    df[a.condition_key] = obs[a.condition_key].astype(str).to_numpy()[mask]
    table = df.groupby([a.condition_key, a.sample_key], observed=True).agg(**{g: (g, "mean") for g in genes}, n_cells=(genes[0], "size"))
    table = table[table["n_cells"] > 0].reset_index()

    pd.set_option("display.width", 200)
    print(f"\nfraction of cells with >= {a.min_count:g} count(s), per library:")
    show = table.copy()
    for g in genes:
        show[g] = (show[g] * 100).round(1)
    print(show.to_string(index=False))

    conds = sorted(table[a.condition_key].unique())
    if len(conds) < 2:
        print("\nonly one condition present in the selected cells; no comparison possible")
    else:
        print(f"\nlibraries that resemble the opposite condition (per gene; only genes whose condition medians differ by >= {100*a.min_gap:.0f} points are informative):")
        any_flag = False
        for g in genes:
            med = table.groupby(a.condition_key, observed=True)[g].median()
            gaps = [abs(med[x] - med[y]) for x in conds for y in conds if x < y]
            if max(gaps) < a.min_gap:
                print(f"  {g}: condition medians differ by only {100*max(gaps):.1f} points; not an informative marker for the labels")
                continue
            for _, row in table.iterrows():
                own = row[a.condition_key]
                others = [c for c in conds if c != own]
                for oc in others:
                    if abs(med[own] - med[oc]) < a.min_gap:
                        continue
                    if abs(row[g] - med[oc]) < abs(row[g] - med[own]):
                        any_flag = True
                        print(f"  {g}: library {row[a.sample_key]} labelled {own} has {100*row[g]:.1f}% positive "
                              f"(median {own} {100*med[own]:.1f}%, median {oc} {100*med[oc]:.1f}%) -> closer to {oc}")
        if not any_flag:
            print("  none")

        # Second rule, independent of the median gap. When the case cells are a minority of every
        # library (tumour nests in omentum, a disease state inside a tissue that keeps its composition),
        # the condition medians never differ by 20 points and the first rule is blind; the informative
        # signal is a library that carries far more of the other condition's marker than the rest of
        # its own condition. Rule: for a gene whose median is higher in the other condition, a library
        # is listed when its fraction is >= --outlier-fold times its condition's median and at least
        # --outlier-gap above it. Found the tumour-bearing "distant omentum" libraries (EPCAM 10%
        # against a 0.1% median) and a benign epithelial inclusion in a control library (PAX8 5%).
        print(f"\nlibraries carrying the other condition's marker far above their own condition "
              f"(>= {a.outlier_fold:g}x the median of their condition and >= {100*a.outlier_gap:.0f} points above it; "
              f"independent of the gap between condition medians):")
        any_out = False
        for g in genes:
            med = table.groupby(a.condition_key, observed=True)[g].median()
            for _, row in table.iterrows():
                own = row[a.condition_key]
                higher = [c for c in conds if c != own and med[c] > med[own]]
                if not higher:
                    continue
                oc = max(higher, key=lambda c: med[c])
                if row[g] >= a.outlier_fold * med[own] and row[g] - med[own] >= a.outlier_gap:
                    any_out = True
                    print(f"  {g}: library {row[a.sample_key]} labelled {own} has {100*row[g]:.1f}% positive "
                          f"(median {own} {100*med[own]:.1f}%, median {oc} {100*med[oc]:.1f}%; n = {int(row['n_cells'])}) "
                          f"-> possible {oc} cells in a {own} library")
        if not any_out:
            print("  none")
        print("\nA flagged library is a prompt to look (swapped label, contaminated sample, low purity), not a verdict. "
              "Decide before any run how it will be handled; exclusion is a config decision with a written reason.")
    if a.csv:
        table.to_csv(a.csv, index=False)
        print(f"\nwrote {a.csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
