"""The loop must survive bad tool calls from the model and report them back as tool errors."""
from sctarget.agents.loop import run_role
from sctarget.agents.model import ModelResponse, ToolUse
from sctarget.session import AnalysisSession
from sctarget.tools import REGISTRY


class MisbehavingModel:
    name = "misbehaving"

    def __init__(self):
        self.calls = 0

    def create(self, *, system, messages, tools, role):
        self.calls += 1
        if self.calls == 1:
            tus = [ToolUse("load_dataset", {}), ToolUse("find_markers", {"cluster": "1", "bogus_arg": 2}), ToolUse("druggability_lookup", {"genes": ["X"]})]
            return ModelResponse("", tus, "tool_use", {}, [{"type": "tool_use", "id": t.id, "name": t.name, "input": t.input} for t in tus])
        # second turn: inspect the tool results and stop
        last = messages[-1]["content"]
        assert last[1]["is_error"] and "unknown arguments" in last[1]["content"]
        assert last[2]["is_error"] and "not available to role" in last[2]["content"]
        return ModelResponse("done", [], "end_turn", {}, [{"type": "text", "text": "done"}])


def test_loop_reports_errors_and_continues(cfg, tmp_path):
    s = AnalysisSession.create(cfg, tmp_path / "loop")
    m = MisbehavingModel()
    rr = run_role(role="analyst", system="x", user_prompt="go", session=s, model=m, registry=REGISTRY, max_turns=3)
    assert rr.stopped_reason == "end_turn" and rr.n_turns == 2
    errs = [c for c in rr.tool_calls if c["error"]]
    assert len(errs) == 2
    assert s.adata is not None  # load_dataset succeeded despite the others failing


def test_turn_budget(cfg, tmp_path):
    class Forever:
        name = "forever"

        def create(self, *, system, messages, tools, role):
            t = ToolUse("load_dataset", {})
            return ModelResponse("", [t], "tool_use", {}, [{"type": "tool_use", "id": t.id, "name": t.name, "input": {}}])

    s = AnalysisSession.create(cfg, tmp_path / "budget")
    rr = run_role(role="planner", system="x", user_prompt="go", session=s, model=Forever(), registry=REGISTRY, max_turns=2)
    assert rr.stopped_reason == "max_turns" and len(rr.tool_calls) == 2
