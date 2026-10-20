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
import warnings

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
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

# GSM -> library, patient, GEO-stated tissue, subtype, sex. NOTE: for the four matched-pair patients
# (SS_2006, SS_2007, SS_2017, SS_2022) the GEO 'tissue' characteristics are SWAPPED relative to the
# authors' annotation files and to the data (CA9/NDUFA4L2 in 19-58% of cells in the libraries the
# authors call ccRCC, < 0.3% in the ones they call normal). Tissue is therefore taken from the
# annotation files (which library appears in GSE159115_ccRCC_anno vs GSE159115_normal_anno); the GEO
# value is kept in obs['geo_tissue'] for the record.
SAMPLES = {
    "GSM4819725": ("SI_18854", "SS_2005", "Tumor", "ccRCC", "M"),
    "GSM4819726": ("SI_18856", "SS_2006", "Tumor", "ccRCC", "M"),   # data: benign
    "GSM4819727": ("SI_18855", "SS_2006", "Benign", None, "M"),     # data: tumour
    "GSM4819728": ("SI_19704", "SS_2007", "Tumor", "ccRCC", "M"),   # data: benign
    "GSM4819729": ("SI_19703", "SS_2007", "Benign", None, "M"),     # data: tumour
    "GSM4819730": ("SI_21255", "SS_2014", "Benign", None, "F"),     # cortex
    "GSM4819731": ("SI_21256", "SS_2014", "Benign", None, "F"),     # medulla
    "GSM4819732": ("SI_21561", "SS_2016", "Tumor", "chRCC", "F"),   # excluded (chromophobe)
    "GSM4819733": ("SI_22369", "SS_2017", "Tumor", "ccRCC", "M"),   # data: benign
    "GSM4819734": ("SI_22368", "SS_2017", "Benign", None, "M"),     # data: tumour
    "GSM4819735": ("SI_22605", "SS_2022", "Tumor", "ccRCC", "M"),   # data: benign
    "GSM4819736": ("SI_22604", "SS_2022", "Benign", None, "M"),     # data: tumour
    "GSM4819737": ("SI_23459", "SS_2023", "Tumor", "ccRCC", "M"),
    "GSM4819738": ("SI_23843", "SS_2026", "Tumor", "ccRCC", "M"),
}

# authors' fine annotation -> broad lineage (what celltype_key should hold: tumour cells and normal
# tubule share the epithelial lineage, so normal proximal tubule is the tumour state's sibling)
BROAD = {
    "Tumor": "Epithelial",
    **{k: "Epithelial" for k in ("PT-A", "PT-B", "PT-C", "DL", "tAL", "TAL", "DCT", "CNT", "PC", "IC-A", "IC-B", "IC-PC", "UC", "Podo")},
    **{k: "Endothelial" for k in ("Endo_ACKR1", "Endo_PLVAP", "AEA-DVR", "AVR", "GC")},
    **{k: "Mesenchymal" for k in ("Peri", "vSMC", "Mesangial")},
    **{k: "Immune" for k in ("Bcell", "Plasma", "Tcell", "Tcell_CD8", "NKcell", "Macro", "Macro_MKI67", "Mono", "Mast")},
}


def _load_annotations() -> pd.DataFrame | None:
    frames = []
    for name, tissue in (("GSE159115_ccRCC_anno.csv.gz", "Tumor"), ("GSE159115_normal_anno.csv.gz", "Benign")):
        p = D / name
        if p.exists():
            df = pd.read_csv(p, index_col=0)
            df["_tissue"] = tissue
            frames.append(df)
    return pd.concat(frames) if frames else None


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


def _attach_annotation(adata: ad.AnnData, anno: pd.DataFrame | None) -> None:
    if anno is None:
        print("no annotation CSVs found; obs['author_celltype'] not set")
        adata.obs["author_celltype"] = "unannotated"
        adata.obs["author_anno_fine"] = "unannotated"
        return
    # annotation index is '<library>_<barcode>'; our obs_names are '<GSM>_<barcode>' -> rebuild as library_barcode
    lib_bc = adata.obs["library"].astype(str).to_numpy() + "_" + adata.obs_names.str.split("_").str[-1].to_numpy()
    fine = pd.Series(anno["anno"].astype(str).to_numpy(), index=anno.index.astype(str))
    fine = fine[~fine.index.duplicated()]
    got = pd.Series(lib_bc, index=adata.obs_names).map(fine)
    n = int(got.notna().sum())
    print(f"author annotation matched on library+barcode for {n}/{adata.n_obs} cells")
    adata.obs["author_anno_fine"] = got.fillna("unannotated").astype(str).to_numpy()
    adata.obs["author_celltype"] = got.map(BROAD).fillna("unannotated").astype(str).to_numpy()
    print(adata.obs["author_celltype"].value_counts().to_string())


def main() -> None:
    if not TAR.exists():
        sys.exit(f"missing {TAR}: download GSE159115_RAW.tar into data/ccrcc first")
    work = D / "raw"
    work.mkdir(exist_ok=True)
    with tarfile.open(TAR) as t:
        t.extractall(work)
    files = sorted(p for p in work.rglob("*") if p.is_file())
    print(f"extracted {len(files)} files, e.g. {[f.name for f in files[:6]]}")
    anno = _load_annotations()
    lib_tissue = anno.groupby("sample")["_tissue"].first().to_dict() if anno is not None else {}
    parts = []
    for gsm, (lib, patient, geo_tissue, subtype, sex) in SAMPLES.items():
        if subtype == "chRCC":
            print(f"skip {gsm} (chromophobe)")
            continue
        fs = [f for f in files if f.name.startswith(gsm)]
        if not fs:
            print(f"WARNING: no files for {gsm} ({lib})")
            continue
        tissue = lib_tissue.get(lib, geo_tissue)
        if tissue != geo_tissue:
            print(f"NOTE {gsm} {lib}: GEO says {geo_tissue}, authors' annotation says {tissue} -> using the annotation (confirmed by CA9/NDUFA4L2)")
        a = _read_sample(fs, gsm)
        a.obs["sample"] = gsm
        a.obs["library"] = lib
        a.obs["patient"] = patient
        a.obs["tissue"] = tissue
        a.obs["geo_tissue"] = geo_tissue
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
    for c in ("sample", "patient", "tissue", "geo_tissue", "condition", "sex", "library"):
        adata.obs[c] = adata.obs[c].astype("category")
    _attach_annotation(adata, anno)
    print(adata)
    print(adata.obs.groupby("condition", observed=True)["sample"].nunique().rename("libraries per condition"))
    print(adata.obs.groupby("condition", observed=True)["patient"].nunique().rename("patients per condition"))
    # sanity: tumour libraries must be CA9-high
    if "CA9" in adata.var_names:
        j = list(adata.var_names).index("CA9")
        ca9 = np.asarray((adata.layers["counts"][:, j] > 0).todense()).ravel() if sp.issparse(adata.layers["counts"]) else adata.layers["counts"][:, j] > 0
        chk = pd.Series(ca9, index=adata.obs_names).groupby(adata.obs["library"].astype(str)).mean().round(3)
        chk = pd.DataFrame({"CA9_frac": chk, "condition": adata.obs.groupby("library", observed=True)["condition"].first()})
        print(chk.to_string())
        bad = chk[(chk["condition"] == "ccRCC") & (chk["CA9_frac"] < 0.1)]
        if len(bad):
            print("WARNING: ccRCC libraries with low CA9:", list(bad.index))
    adata.write_h5ad(OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
