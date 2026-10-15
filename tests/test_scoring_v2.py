"""Scoring v2: a second effect axis (disease vs reference within the state's compartment) for
composition states, recorded as target_class; v1 behaviour unchanged by default; rescore works."""
from sctarget.agents.orchestrator import run_pipeline
from sctarget.eval import evaluate_run
from sctarget.tools.analysis import compartment_mask, pseudobulk_de
from sctarget.tools.targets import score_targets, scoring_version


def test_default_scoring_version_is_1(cfg):
    class S:  # minimal stand-in
        pass

    s = S()
    s.cfg = cfg
    assert scoring_version(s) == 1


def test_compartment_contrast_and_mask(cfg, tmp_path):
    from sctarget.session import AnalysisSession
    from sctarget.tools.analysis import integrate_and_cluster, load_dataset

    s = AnalysisSession.create(cfg, tmp_path / "cp")
    load_dataset(s)
    integrate_and_cluster(s, reuse_existing=cfg.dataset.get("reuse_clustering"))
    a = s.adata
    key = s.cluster_key
    cl = str(a.obs[key].astype(str).value_counts().idxmax())
    in_cl = (a.obs[key].astype(str) == cl).to_numpy()
    m_all = compartment_mask(s, cl, key, shared_only=False)
    assert m_all[in_cl].all()  # the full compartment contains the cluster
    m = compartment_mask(s, cl, key)  # shared clusters only: the cluster itself is excluded
    assert not m[in_cl].any() and m.sum() > 0
    # the contrast needs reference cells in the shared clusters of the lineage; in a tiny synthetic
    # not every cluster has them, so it must either succeed or say exactly why not
    ok = []
    for c in a.obs[key].astype(str).unique():
        r = pseudobulk_de(s, cluster=str(c), contrast="disease_vs_reference_within_compartment")
        if "error" in r:
            assert "insufficient pseudobulk replicates" in r["error"], r
        else:
            ok.append(str(c))
            assert f"{c}:disease_vs_reference_within_compartment" in s.facts["pseudobulk_de"]
    assert ok, "compartment contrast never computable"


def test_scoring_v2_adds_lineage_axis_and_keeps_planted_target_on_top(synthetic_dir, cfg, tmp_path):
    import copy

    cfg2 = copy.deepcopy(cfg)
    cfg2.raw.setdefault("targets", {})["scoring_version"] = 2
    res = run_pipeline(cfg2, model_name="scripted", run_dir=tmp_path / "s2")
    assert res.facts["scoring_version"] == 2
    rows = next(iter(res.facts["ranked_targets"].values()))
    assert all("target_class" in r for r in rows)
    assert any(r["lineage_log2fc_disease_vs_reference"] is not None for r in rows), "lineage axis never computed"
    assert rows[0]["gene"] == "ERBB4"  # the planted state-specific target still wins
    md = (res.run_dir / "report.md").read_text()
    assert "scoring v2" in md and "Class" in md


def test_rescore_script_reproduces_ranking(cfg, tmp_path):
    from scripts.rescore import rescore

    base = run_pipeline(cfg, model_name="scripted", run_dir=tmp_path / "base")
    import copy

    cfg1 = copy.deepcopy(cfg)
    out = tmp_path / "base_s1"
    ch = rescore(cfg1, base.run_dir, out)
    for cl, c in ch.items():
        assert all(r["rank_before"] == r["rank_after"] for r in c["rows"]), c  # same version -> same ranking
    import yaml
    from tests.conftest import ROOT

    gt = yaml.safe_load(open(ROOT / "data" / "ground_truth" / "synthetic.yaml"))
    gt["per_cell_truth"] = None
    p = tmp_path / "gt.yaml"
    yaml.safe_dump(gt, open(p, "w"))
    ev = evaluate_run(out, p)
    assert ev.target_rank["ERBB4"] == 1
