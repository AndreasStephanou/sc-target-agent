"""Expression route: disease signal is a program shift inside a shared cell type (no composition change)."""
from pathlib import Path

import pytest

from sctarget.agents.orchestrator import run_pipeline
from sctarget.eval import evaluate_run
from sctarget.session import AnalysisSession
from sctarget.tools import analysis as A
from sctarget.tools import targets as T


@pytest.fixture(scope="module")
def program_cfg(tmp_path_factory):
    import yaml
    from scripts.make_synthetic import make_synthetic
    from sctarget.config import load_config
    from tests.conftest import ROOT

    out = tmp_path_factory.mktemp("synth_prog")
    make_synthetic(out, cells_per_sample=200, seed=2, mode="program")
    cfg = load_config(ROOT / "config" / "synthetic_program.yaml", overrides={"dataset": {"path": str(out / "synthetic.h5ad")}, "output": {"run_dir": str(tmp_path_factory.mktemp("runs_prog"))}})
    with open(ROOT / "data" / "ground_truth" / "synthetic_program.yaml") as fh:
        gt = yaml.safe_load(fh)
    gt["per_cell_truth"] = str(out / "true_state.csv")
    gtp = out / "gt.yaml"
    with open(gtp, "w") as fh:
        yaml.safe_dump(gt, fh)
    return cfg, gtp


def test_disease_programs_flags_the_shared_cluster(program_cfg, tmp_path):
    cfg, _ = program_cfg
    s = AnalysisSession.create(cfg, tmp_path / "tool")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing="leiden")
    da = A.differential_abundance(s)
    assert da["disease_enriched_clusters"] == [], "composition should be flat in program mode"
    dp = A.disease_programs(s)
    assert len(dp["disease_altered_clusters"]) == 1
    cl = dp["disease_altered_clusters"][0]
    row = dp["table"].loc[cl]
    assert row["dominant_celltype"] == "keratinocyte"
    assert "ERBB4" in row["top_up_genes"]
    assert T.state_basis(s, cl) == "program"


def test_program_route_end_to_end(program_cfg):
    cfg, gtp = program_cfg
    res = run_pipeline(cfg, model_name="scripted", run_dir=Path(cfg.run_root) / "e2e")
    for rr in res.roles.values():
        assert not [c for c in rr.tool_calls if c["error"]], rr.tool_calls
    states = res.facts["findings"]["disease_state"]
    assert any(st.get("basis") == "program" for st in states)
    ev = evaluate_run(res.run_dir, gtp)
    assert ev.state_recovery["program_keratinocyte"]["found"]
    assert ev.target_rank["ERBB4"] == 1
    assert not ev.decoys_nominated_unflagged
    assert ev.passed
