"""Annotation-free structure: lineages derived from the data drive sibling sets / compartments."""
import copy

from sctarget.agents.orchestrator import run_pipeline
from sctarget.eval import evaluate_run
from sctarget.session import AnalysisSession
from sctarget.tools import analysis as A


def _cfg_without_annotation(cfg):
    c = copy.deepcopy(cfg)
    c.raw["dataset"].pop("celltype_key", None)
    return c


def test_lineage_key_prefers_annotation_then_derived(cfg, tmp_path):
    s = AnalysisSession.create(cfg, tmp_path / "lk")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing=cfg.dataset.get("reuse_clustering"))
    assert s.lineage_key == cfg.dataset["celltype_key"]
    s2 = AnalysisSession.create(_cfg_without_annotation(cfg), tmp_path / "lk2")
    A.load_dataset(s2)
    assert s2.lineage_key is None
    A.integrate_and_cluster(s2, reuse_existing=None, resolution=2.0, use_harmony=False)
    r = A.derive_lineages(s2)
    assert "error" not in r and r["n_lineages"] >= 2 and s2.lineage_key == s2.DERIVED_LINEAGE_KEY
    assert set(r["fine_cluster_purity"]) == set(s2.adata.obs[s2.cluster_key].astype(str).unique())
    comp = A.cluster_composition(s2)
    assert "dominant_celltype" in comp["table"].columns  # structural logic now has lineages
    lab = A.label_lineages(s2)
    assert set(lab["labels"]) == set(r["groups"]) and all("confidence" in v for v in lab["labels"].values())


def test_tools_derive_lineages_lazily_without_annotation(cfg, tmp_path):
    s = AnalysisSession.create(_cfg_without_annotation(cfg), tmp_path / "lazy")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing=None, resolution=2.0, use_harmony=False)
    A.cluster_composition(s)  # should trigger derive_lineages internally
    assert s.lineage_key == s.DERIVED_LINEAGE_KEY and "lineages" in s.facts


def test_pipeline_passes_without_cell_type_annotation(cfg, ground_truth_path, tmp_path):
    res = run_pipeline(_cfg_without_annotation(cfg), model_name="scripted", run_dir=tmp_path / "noanno")
    assert res.facts.get("lineages"), "lineages were never derived"
    ev = evaluate_run(res.run_dir, ground_truth_path)
    assert ev.state_recovery["neuro_keratinocyte"]["found"]
    assert ev.target_rank["ERBB4"] == 1
    assert "derived from the data" in (res.run_dir / "report.md").read_text()
