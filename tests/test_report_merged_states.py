"""Protocol 23q: report figures draw a state recorded under a label that was later merged.

Motivating case (CAA de novo agent run 1, round 0.4c): the analyst recorded states 16, 17 and 13;
the critic then ran merge_clusters(16, 17, 13), which made '13+16+17' the active label. Every
cell-level report figure looked the recorded labels up in the active column and drew "state 16:
0 cells". The same happened in FLC agent run 1 (state 8, merged into 8+19).
"""
from __future__ import annotations

import types
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from sctarget.report import figures as F
from sctarget.session import AnalysisSession


class _Cfg:
    dataset = {"sample_key": "sample", "condition_key": "condition", "disease": "D", "reference": ["R"]}
    analysis: dict = {}

    @staticmethod
    def get(key, default=None):
        return default


def _session(tmp_path: Path, merged: bool = True, subclustered: bool = False):
    rng = np.random.default_rng(0)
    n = 600
    base = rng.choice(["0", "1", "2", "3"], n)
    cond = np.where(rng.random(n) < 0.5, "D", "R")
    cond[np.isin(base, ["1", "2"])] = "D"  # the two pieces are disease-only, like real merged states
    sample = np.array([f"{c}{i}" for c, i in zip(cond, rng.integers(1, 4, n))])
    obs = pd.DataFrame({"sample": sample, "condition": cond, "sct_leiden_r1": base}, index=[f"c{i}" for i in range(n)])
    a = ad.AnnData(X=sp.csr_matrix(rng.poisson(1.0, (n, 12)).astype(np.float32)), obs=obs)
    a.var_names = [f"G{i}" for i in range(12)]
    a.obsm["X_umap"] = rng.normal(size=(n, 2))
    key = "sct_leiden_r1"
    if merged:
        key = "sct_leiden_r1_merged"
        a.obs[key] = pd.Series(base, index=a.obs_names).where(~np.isin(base, ["1", "2"]), "1+2").astype(str)
    if subclustered:
        key = "sct_leiden_r1_sub"
        a.obs[key] = pd.Series(base, index=a.obs_names).where(base != "3", pd.Series(np.where(rng.random(n) < 0.5, "3_0", "3_1"), index=a.obs_names)).astype(str)
    s = types.SimpleNamespace(
        adata=a, cluster_key=key, cfg=_Cfg(), run_dir=tmp_path, figures=[], pseudo=None,
        sample_key="sample", condition_key="condition", disease="D", reference=["R"], lineage_key=None, celltype_key=None,
        facts={"findings": {"disease_state": [
            {"cluster": "1", "basis": "composition", "genes": ["G1"]},
            {"cluster": "2", "basis": "composition", "genes": ["G2"]},
            {"cluster": "3", "basis": "program", "genes": ["G3"]},
        ]}},
    )
    s.require_adata = lambda: a
    s.counts = lambda: a.X
    s.save_figure = types.MethodType(AnalysisSession.save_figure, s)
    s.log = types.SimpleNamespace(note=lambda *a_, **k: None)
    (tmp_path / "figures").mkdir(exist_ok=True)
    return s, a


def test_pieces_of_a_merge_resolve_to_one_drawn_state(tmp_path):
    s, a = _session(tmp_path, merged=True)
    ents = F.resolve_states(s)
    labels = [e["label"] for e in ents]
    assert labels == ["1+2", "3"], labels
    merged = ents[0]
    assert merged["recorded"] == ["1", "2"]
    assert int(merged["cells"].sum()) == int((a.obs["sct_leiden_r1_merged"] == "1+2").sum()) > 0
    assert "merged by merge_clusters" in merged["note"]
    assert ents[1]["note"] is None
    assert F.resolved_state_label(s, "1") == "1+2" and F.resolved_state_label(s, "3") == "3"


def test_groups_and_state_shares_draw_the_merged_state(tmp_path):
    s, a = _session(tmp_path, merged=True)
    groups, order = F._groups(s)
    assert "state 1+2" in order and (groups == "state 1+2").sum() == (a.obs["sct_leiden_r1_merged"] == "1+2").sum()
    assert "state 1" not in order
    rel = F.state_library_bars(s)
    assert rel and (tmp_path / rel).exists()
    rel2 = F.state_map(s)
    assert rel2 and (tmp_path / rel2).exists()
    rel3 = F.library_expression_panel(s, ["G1", "G3"], "panel")
    assert rel3 and (tmp_path / rel3).exists()


def test_unmerged_states_resolve_as_before(tmp_path):
    s, a = _session(tmp_path, merged=False)
    ents = F.resolve_states(s)
    assert [e["label"] for e in ents] == ["1", "2", "3"]
    assert all(e["note"] is None for e in ents)


def test_derived_lineages_restore_from_the_run_record():
    from sctarget.tools.analysis import lineage_column_from_facts

    facts = {"lineages": {"key": "sct_lineage", "fine_to_lineage": {"0": "L0", "1": "L4", "2": "L4", "3": "L0", "6": "L2"}}}
    labels = np.array(["0", "1+2", "1+2", "3", "6_1", "9"])
    col = lineage_column_from_facts(labels, facts)
    assert list(col) == ["L0", "L4", "L4", "L0", "L2", "unassigned"]
    assert lineage_column_from_facts(labels, {}) is None


def test_rebuilt_session_matches_the_run(cfg, tmp_path):
    """A report rebuilt from a run directory sees the run's cells and the run's lineages, not a
    reloaded file with the QC-removed cells back and lineages re-derived under new numbers."""
    import json

    from sctarget.agents.orchestrator import STAGES, _resume_session, run_pipeline
    from sctarget.session import AnalysisSession

    res = run_pipeline(cfg, model_name="scripted", run_dir=tmp_path / "orig")
    state = json.loads((res.run_dir / "state.json").read_text())
    assign = pd.read_csv(res.run_dir / "cluster_assignments.csv", index_col=0)["cluster"].astype(str)
    s = AnalysisSession.create(cfg, res.run_dir)
    _resume_session(s, res.run_dir, keep_findings_from=set(STAGES))
    a = s.adata
    assert a.n_obs == len(assign) and set(a.obs_names) == set(assign.index)
    assert list(a.obs[state["cluster_key"]].astype(str)) == list(assign.reindex(a.obs_names))
    lin = state["facts"].get("lineages") or {}
    if lin.get("key") == AnalysisSession.DERIVED_LINEAGE_KEY:
        f2l = lin["fine_to_lineage"]
        expect = [f2l.get(c) or f2l.get(c.split("+")[0]) for c in assign.reindex(a.obs_names)]
        assert list(a.obs[AnalysisSession.DERIVED_LINEAGE_KEY].astype(str)) == expect
    if (res.run_dir / "umap_coordinates.csv").exists():
        assert "X_umap" in a.obsm and a.obsm["X_umap"].shape == (a.n_obs, 2)


def test_state_recorded_on_the_base_clustering_after_a_subcluster_step(tmp_path):
    s, a = _session(tmp_path, merged=False, subclustered=True)
    ents = F.resolve_states(s)
    by = {e["label"]: e for e in ents}
    assert set(by) == {"1", "2", "3"}
    assert by["3"]["note"] and "sct_leiden_r1" in by["3"]["note"]
    assert int(by["3"]["cells"].sum()) == int((a.obs["sct_leiden_r1"] == "3").sum())
