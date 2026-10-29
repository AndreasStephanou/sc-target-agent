"""Counts in X with no counts_layer configured (the common freshly built .h5ad) must survive normalisation.

HCC review run (2026-10-28): the first clustering call log-normalised X in place and every pseudobulk,
compartment-contrast and scoring call then failed with "needs raw counts".
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import scipy.sparse as sp

from sctarget.config import load_config
from sctarget.session import AnalysisSession
from sctarget.tools.analysis import _ensure_lognorm, _looks_like_counts, load_dataset
from tests.conftest import ROOT


def test_counts_in_x_are_stashed_before_normalisation(synthetic_dir, tmp_path_factory):
    src = ad.read_h5ad(synthetic_dir / "synthetic.h5ad")
    a = ad.AnnData(X=sp.csr_matrix(src.layers["counts"]), obs=src.obs[["sample", "condition"]].copy(), var=src.var[[]].copy())
    p = tmp_path_factory.mktemp("cx") / "counts_in_x.h5ad"
    a.write_h5ad(p)
    cfg = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(p), "counts_layer": None}, "output": {"run_dir": str(tmp_path_factory.mktemp("runs_cx"))}})
    assert cfg.dataset.get("counts_layer") is None
    s = AnalysisSession.create(cfg, tmp_path_factory.mktemp("run_cx"))
    summary = load_dataset(s)
    assert summary["counts_layer_used"] == "counts" and "copied to layer" in summary["counts_note"]
    assert s.facts["counts_layer_used"] == "counts"
    _ensure_lognorm(s)
    assert not _looks_like_counts(s.adata.X), "X should now be log-normalised"
    assert _looks_like_counts(s.counts()), "session.counts() must still return raw integer counts"


def test_missing_configured_layer_is_reported(synthetic_dir, tmp_path_factory):
    cfg = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(synthetic_dir / "synthetic.h5ad"), "counts_layer": "does_not_exist"}, "output": {"run_dir": str(tmp_path_factory.mktemp("runs_ml"))}})
    s = AnalysisSession.create(cfg, tmp_path_factory.mktemp("run_ml"))
    summary = load_dataset(s)
    assert "not in the object" in summary["counts_note"]
    assert summary["counts_layer_used"] == "counts"  # the real counts layer was found and used


def test_exclude_samples_is_applied_and_recorded(synthetic_dir, tmp_path_factory):
    import anndata as ad

    a = ad.read_h5ad(synthetic_dir / "synthetic.h5ad")
    drop = sorted(a.obs["sample"].astype(str).unique())[:1]
    cfg = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(synthetic_dir / "synthetic.h5ad"), "exclude_samples": drop + ["not_a_sample"], "exclude_samples_reason": "test"}, "output": {"run_dir": str(tmp_path_factory.mktemp("runs_ex"))}})
    s = AnalysisSession.create(cfg, tmp_path_factory.mktemp("run_ex"))
    summary = load_dataset(s)
    ex = summary["excluded_samples"]
    assert ex["n_samples_excluded"] == 1 and ex["not_found"] == ["not_a_sample"] and ex["reason"] == "test"
    assert drop[0] not in set(s.adata.obs["sample"].astype(str))
    assert s.facts["excluded_samples"]["n_samples"] == 1
    assert sum(summary["samples_per_condition"].values()) == a.obs["sample"].nunique() - 1


def test_pseudonyms_stable_under_exclusion(synthetic_dir, tmp_path_factory):
    import anndata as ad

    a = ad.read_h5ad(synthetic_dir / "synthetic.h5ad")
    libs = sorted(a.obs["sample"].astype(str).unique())
    base = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(synthetic_dir / "synthetic.h5ad")}, "output": {"run_dir": str(tmp_path_factory.mktemp("r1"))}})
    s1 = AnalysisSession.create(base, tmp_path_factory.mktemp("p1")); load_dataset(s1)
    excl = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(synthetic_dir / "synthetic.h5ad"), "exclude_samples": [libs[0]]}, "output": {"run_dir": str(tmp_path_factory.mktemp("r2"))}})
    s2 = AnalysisSession.create(excl, tmp_path_factory.mktemp("p2")); load_dataset(s2)
    for l in libs[1:]:
        assert s1.pseudo.to_pseudo(l) == s2.pseudo.to_pseudo(l)
