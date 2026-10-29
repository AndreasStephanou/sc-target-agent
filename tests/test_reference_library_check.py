"""check_reference_libraries: marker-free detection of a library that resembles the opposite condition.
Motivated by HCC GSE149614 library HCC04N (labelled non-tumour, tumour-like in the hepatocyte lineage):
with it in the reference set GPC3 was log2 2.6 at padj 0.20; without it, log2 6.5 at padj 8e-12."""
from __future__ import annotations

import copy

import anndata as ad

from sctarget.config import load_config
from sctarget.session import AnalysisSession
from sctarget.tools.analysis import check_reference_libraries, load_dataset
from tests.conftest import ROOT


def _session(cfg, run_dir):
    s = AnalysisSession.create(cfg, run_dir)
    load_dataset(s)
    return s


def test_clean_synthetic_has_no_flags(cfg, tmp_path):
    s = _session(cfg, tmp_path / "clean")
    r = check_reference_libraries(s, min_cells_per_library=20)
    assert r["flagged"] == [], r["note"]
    # the synthetic disease signal is a small sub-state: lineage pseudobulks carry no condition
    # structure, and the tool must say so rather than judge libraries on noise
    assert all("libraries" in v or "note" in v for v in r["lineages"].values())
    assert "could not be tested" in r["note"] or "No library" in r["note"]


def _with_condition_programme(a, dis, strength=4.0, n_genes=60, seed=0):
    """Add a condition-wide programme to every disease cell so lineage pseudobulks separate by condition."""
    import numpy as np
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    X = a.layers["counts"] if "counts" in a.layers else a.X
    X = sp.csr_matrix(X).tolil()
    genes = rng.choice(a.n_vars, n_genes, replace=False)
    cells = np.where(a.obs["condition"].astype(str).to_numpy() == dis)[0]
    for g in genes:
        col = X[:, g].toarray().ravel()
        col[cells] = np.round(col[cells] * strength + rng.poisson(2.0, len(cells)))
        X[:, g] = col.reshape(-1, 1)
    X = X.tocsr()
    if "counts" in a.layers:
        a.layers["counts"] = X
    a.X = X.copy()
    return a


def test_swapped_library_is_flagged(synthetic_dir, cfg, tmp_path):
    a = ad.read_h5ad(synthetic_dir / "synthetic.h5ad")
    dis = cfg.dataset["disease"]
    a = _with_condition_programme(a, dis)
    # relabel one DISEASE library as reference: its cells still carry the disease programme
    dis_libs = sorted(a.obs.loc[a.obs["condition"].astype(str) == dis, "sample"].astype(str).unique())
    swapped = dis_libs[0]
    ref_label = [c for c in a.obs["condition"].astype(str).unique() if c != dis][0]
    a.obs["condition"] = a.obs["condition"].astype(str)
    a.obs.loc[a.obs["sample"].astype(str) == swapped, "condition"] = ref_label
    p = tmp_path / "swapped.h5ad"
    a.write_h5ad(p)
    cfg2 = copy.deepcopy(cfg)
    cfg2.raw["dataset"]["path"] = str(p)
    s = _session(cfg2, tmp_path / "sw")
    r = check_reference_libraries(s, min_cells_per_library=20)
    assert any(v.get("libraries") for v in r["lineages"].values()), r["lineages"]  # structure now exists
    pseudo = str(s.pseudo.to_pseudo(swapped))
    flagged_libs = {f["library"] for f in r["flagged"]}
    assert pseudo in flagged_libs, (r["flagged"], r["lineages"])
    assert any(f["label"] == ref_label for f in r["flagged"] if f["library"] == pseudo)
    assert "OPPOSITE" in r["note"]


def test_scripted_critic_records_reference_library_flag(synthetic_dir, cfg, tmp_path):
    from sctarget.agents.orchestrator import run_pipeline

    a = ad.read_h5ad(synthetic_dir / "synthetic.h5ad")
    dis = cfg.dataset["disease"]
    a = _with_condition_programme(a, dis)
    dis_libs = sorted(a.obs.loc[a.obs["condition"].astype(str) == dis, "sample"].astype(str).unique())
    ref_label = [c for c in a.obs["condition"].astype(str).unique() if c != dis][0]
    a.obs["condition"] = a.obs["condition"].astype(str)
    a.obs.loc[a.obs["sample"].astype(str) == dis_libs[0], "condition"] = ref_label
    p = tmp_path / "swapped2.h5ad"
    a.write_h5ad(p)
    cfg2 = copy.deepcopy(cfg)
    cfg2.raw["dataset"]["path"] = str(p)
    res = run_pipeline(cfg2, model_name="scripted", run_dir=tmp_path / "crit")
    flags = res.facts.get("findings", {}).get("critic_flag", [])
    assert any("resembles the opposite condition" in fl.get("title", "") for fl in flags), [fl.get("title") for fl in flags]
