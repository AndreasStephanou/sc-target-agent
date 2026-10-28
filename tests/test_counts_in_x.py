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
