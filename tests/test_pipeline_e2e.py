import json
from pathlib import Path

import pytest

from sctarget.agents.orchestrator import run_pipeline
from sctarget.eval import evaluate_run, stability


@pytest.fixture(scope="module")
def agent_run(cfg, tmp_path_factory):
    rd = tmp_path_factory.mktemp("e2e") / "agent"
    return run_pipeline(cfg, model_name="scripted", run_dir=rd)


def test_all_roles_ran_without_tool_errors(agent_run):
    assert set(agent_run.roles) == {"planner", "analyst", "target", "critic"}
    for rr in agent_run.roles.values():
        assert rr.stopped_reason == "end_turn"
        assert not [c for c in rr.tool_calls if c["error"]], rr.tool_calls


def test_outputs_exist(agent_run):
    rd = agent_run.run_dir
    for f in ["run.jsonl", "state.json", "report.md", "report.html", "config_used.json", "cluster_assignments.csv"]:
        assert (rd / f).exists(), f
    assert list((rd / "figures").glob("*.png"))
    assert list((rd / "tables").glob("*.csv"))
    md = (rd / "report.md").read_text()
    assert "ERBB4" in md and "## 5. Critic review" in md and "## 7. Provenance" in md


def test_nothing_per_cell_reached_the_model(agent_run):
    """Every payload the model saw is small and contains no per-cell-sized arrays."""
    n_cells = agent_run.facts["dataset"]["n_cells"]
    with open(agent_run.run_dir / "run.jsonl") as fh:
        events = [json.loads(l) for l in fh if l.strip()]
    tcs = [e for e in events if e["type"] == "tool_call"]
    assert tcs

    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                yield from walk(v)
        elif isinstance(o, list):
            yield len(o)
            for v in o:
                yield from walk(v)

    for e in tcs:
        payload = json.dumps(e["result_sent"])
        assert len(payload) < 60_000, (e["tool"], len(payload))
        assert max(walk(e["result_sent"]), default=0) < n_cells / 2
        # real sample ids never appear
        assert "CAA_S" not in payload and "GINGIVA_S" not in payload


def test_evaluation_passes_and_recovers_planted_answer(agent_run, ground_truth_path):
    res = evaluate_run(agent_run.run_dir, ground_truth_path)
    assert res.state_recovery["neuro_keratinocyte"]["found"]
    assert res.target_rank["ERBB4"] == 1
    assert res.drug_recall == 1.0
    assert not res.decoys_nominated_unflagged
    assert res.passed


def test_stability_of_two_scripted_runs(cfg, agent_run, tmp_path_factory):
    rd2 = tmp_path_factory.mktemp("e2e2") / "agent2"
    run_pipeline(cfg, model_name="scripted", run_dir=rd2)
    st = stability([agent_run.run_dir, rd2], k=10)
    assert st["mean_jaccard_topk"] == 1.0 and st["rank1_consistent"]
