"""Build data/ccrcc/ccrcc.h5ad from GEO GSE159115 (Zhang et al. 2021, PNAS): clear cell RCC vs benign kidney.

    python scripts/build_ccrcc.py            # expects data/ccrcc/GSE159115_RAW.tar (+ optional *_anno.csv.gz)

Handles the two formats GEO RAW tarballs use for 10x data (per-sample filtered_feature_bc_matrix.h5,
or barcodes/features/matrix triplets). Writes one AnnData with raw counts in X and in layers['counts'],
obs: sample (GSM library id), patient, tissue, condition in {ccRCC, Benign}, sex; the chromophobe tumour
is excluded. If the authors' annotation CSVs are present, their cell-type column is attached as
obs['author_celltype'] by barcode (best effort; unmatched cells get 'unannotated').
"""
from __future__ import annotations

import gzip
import io
import re
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "ccrcc"
TAR = D / "GSE159115_RAW.tar"
OUT = D / "ccrcc.h5ad"

# From the GEO sample table (GSM -> library, patient, tissue, subtype, sex)
SAMPLES = {
    "GSM4819725": ("SI_18854", "SS_2005", "Tumor", "ccRCC", "M"),
    "GSM4819726": ("SI_18856", "SS_2006", "Tumor", "ccRCC", "M"),
    "GSM4819727": ("SI_18855", "SS_2006", "Benign", None, "M"),
    "GSM4819728": ("SI_19704", "SS_2007", "Tumor", "ccRCC", "M"),
    "GSM4819729": ("SI_19703", "SS_2007", "Benign", None, "M"),
    "GSM4819730": ("SI_21255", "SS_2014", "Benign", None, "F"),  # cortex
    "GSM4819731": ("SI_21256", "SS_2014", "Benign", None, "F"),  # medulla
    "GSM4819732": ("SI_21561", "SS_2016", "Tumor", "chRCC", "F"),  # excluded (chromophobe)
    "GSM4819733": ("SI_22369", "SS_2017", "Tumor", "ccRCC", "M"),
    "GSM4819734": ("SI_22368", "SS_2017", "Benign", None, "M"),
    "GSM4819735": ("SI_22605", "SS_2022", "Tumor", "ccRCC", "M"),
    "GSM4819736": ("SI_22604", "SS_2022", "Benign", None, "M"),
    "GSM4819737": ("SI_23459", "SS_2023", "Tumor", "ccRCC", "M"),
    "GSM4819738": ("SI_23843", "SS_2026", "Tumor", "ccRCC", "M"),
}


def _read_sample(files: list[Path], gsm: str) -> ad.AnnData:
    """files: all extracted files whose name starts with this GSM."""
    h5 = [f for f in files if f.suffix == ".h5" or f.name.endswith(".h5.gz")]
    if h5:
        f = h5[0]
        if f.name.endswith(".gz"):
            tmp = f.with_suffix("")
            with gzip.open(f, "rb") as src, open(tmp, "wb") as dst:
                shutil.copyfileobj(src, dst)
            f = tmp
        a = sc.read_10x_h5(f)
        a.var_names_make_unique()
        return a
    # triplet: barcodes / features|genes / matrix
    def pick(pat):
        c = [f for f in files if re.search(pat, f.name, re.I)]
        return c[0] if c else None

    bc, ft, mt = pick(r"barcodes"), pick(r"features|genes"), pick(r"matrix\.mtx")
    if not (bc and ft and mt):
        raise RuntimeError(f"{gsm}: cannot find 10x files among {[f.name for f in files]}")
    tmpdir = Path(tempfile.mkdtemp())
    for src, name in ((bc, "barcodes.tsv.gz"), (ft, "features.tsv.gz"), (mt, "matrix.mtx.gz")):
        if src.name.endswith(".gz"):
            shutil.copy(src, tmpdir / name)
        else:
            with open(src, "rb") as s_, gzip.open(tmpdir / name, "wb") as d_:
                shutil.copyfileobj(s_, d_)
    try:
        a = sc.read_10x_mtx(tmpdir, var_names="gene_symbols", cache=False)
    except Exception:
        # older 'genes.tsv' has 2 columns; features has 3. Rename if needed.
        a = sc.read_10x_mtx(tmpdir, var_names="gene_symbols", cache=False, gex_only=False)
    a.var_names_make_unique()
    shutil.rmtree(tmpdir, ignore_errors=True)
    return a


def _attach_annotation(adata: ad.AnnData) -> None:
    annos = sorted(D.glob("GSE159115_*anno*.csv*"))
    if not annos:
        print("no annotation CSVs found; obs['author_celltype'] not set")
        return
    frames = []
    for p in annos:
        df = pd.read_csv(p, index_col=0)
        df["_source"] = p.name
        frames.append(df)
    anno = pd.concat(frames)
    print(f"annotation rows: {len(anno)}; columns: {list(anno.columns)[:12]}")
    # find a cell-type column
    ct_col = next((c for c in anno.columns if re.search(r"cell.?type|celltype|annotation|cluster_name|label", c, re.I)), None)
    if ct_col is None:
        print("no cell-type column recognised; obs['author_celltype'] not set")
        return
    # barcode matching: try index as-is, then strip/add '-1', then prefix-strip
    idx = pd.Index(anno.index.astype(str))
    bc_obs = adata.obs_names.to_series()
    raw_bc = bc_obs.str.split("_").str[-1]  # our obs_names are '<GSM>_<barcode>'
    lab = pd.Series("unannotated", index=adata.obs_names)
    hit = 0
    for variant in (idx, idx.str.replace(r"-\d+$", "", regex=True), idx.str.split("_").str[-1], idx.str.split(":").str[-1]):
        m = pd.Series(anno[ct_col].astype(str).to_numpy(), index=variant)
        m = m[~m.index.duplicated()]
        matched = raw_bc.map(m)
        matched2 = raw_bc.str.replace(r"-\d+$", "", regex=True).map(m)
        got = matched.fillna(matched2)
        n = int(got.notna().sum())
        if n > hit:
            hit = n
            lab = got.fillna("unannotated")
    print(f"author annotation matched for {hit}/{adata.n_obs} cells (barcodes are not sample-unique across libraries; use as a coarse hint only)")
    adata.obs["author_celltype"] = lab.astype(str).to_numpy()


def main() -> None:
    if not TAR.exists():
        sys.exit(f"missing {TAR}: download GSE159115_RAW.tar into data/ccrcc first")
    work = D / "raw"
    work.mkdir(exist_ok=True)
    with tarfile.open(TAR) as t:
        t.extractall(work)
    files = sorted(p for p in work.rglob("*") if p.is_file())
    print(f"extracted {len(files)} files, e.g. {[f.name for f in files[:6]]}")
    parts = []
    for gsm, (lib, patient, tissue, subtype, sex) in SAMPLES.items():
        if subtype == "chRCC":
            print(f"skip {gsm} (chromophobe)")
            continue
        fs = [f for f in files if f.name.startswith(gsm)]
        if not fs:
            print(f"WARNING: no files for {gsm} ({lib})")
            continue
        a = _read_sample(fs, gsm)
        a.obs["sample"] = gsm
        a.obs["library"] = lib
        a.obs["patient"] = patient
        a.obs["tissue"] = tissue
        a.obs["condition"] = "ccRCC" if tissue == "Tumor" else "Benign"
        a.obs["sex"] = sex
        a.obs_names = [f"{gsm}_{b}" for b in a.obs_names]
        print(f"{gsm} {lib} {patient} {tissue}: {a.n_obs} cells x {a.n_vars} genes")
        parts.append(a)
    adata = ad.concat(parts, join="outer", fill_value=0, merge="first")
    adata.var_names_make_unique()
    X = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    adata.X = X.astype(np.float32)
    adata.layers["counts"] = adata.X.copy()
    for c in ("sample", "patient", "tissue", "condition", "sex", "library"):
        adata.obs[c] = adata.obs[c].astype("category")
    _attach_annotation(adata)
    print(adata)
    print(adata.obs.groupby("condition")["sample"].nunique().rename("libraries per condition"))
    print(adata.obs.groupby("condition")["patient"].nunique().rename("patients per condition"))
    adata.write_h5ad(OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
