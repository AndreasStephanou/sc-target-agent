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


def test_analyst_finishing_turn_when_nothing_recorded(cfg, tmp_path):
    """An analyst that computes DE and then stops without recording gets one finishing turn."""
    from sctarget.agents.model import ModelResponse, ScriptedModel, ToolUse
    from sctarget.agents.orchestrator import run_pipeline

    class Forgetful(ScriptedModel):
        name = "forgetful"

        def __init__(self, session):
            super().__init__(session)
            self.finishing_seen = False

        def create(self, *, system, messages, tools, role):
            if role == "analyst":
                first = messages[0]["content"] if isinstance(messages[0]["content"], str) else ""
                if "FINISHING TURN" in first:
                    self.finishing_seen = True
                    if len(messages) == 1:
                        t = ToolUse("record_finding", {"kind": "summary", "title": "No disease-specific state found", "detail": "finishing turn: nothing qualifies at these thresholds", "confidence": "low"})
                        return ModelResponse("", [t], "tool_use", {}, [{"type": "tool_use", "id": t.id, "name": t.name, "input": t.input}])
                    return ModelResponse("done", [], "end_turn", {}, [{"type": "text", "text": "done"}])
                # normal pass: do the analysis but drop every record_finding, then stop silently
                resp = super().create(system=system, messages=messages, tools=tools, role=role)
                if resp.tool_uses and resp.tool_uses[0].name == "record_finding" and resp.tool_uses[0].input.get("kind") == "disease_state":
                    return ModelResponse("", [], "end_turn", {}, [{"type": "text", "text": ""}])  # DE is done by now; stop without recording
                return resp
            return super().create(system=system, messages=messages, tools=tools, role=role)

    import sctarget.agents.orchestrator as orch

    made = {}

    def fake_make_model(name, session, temperature=0.2):
        made["m"] = Forgetful(session)  # bound to the pipeline's own session
        return made["m"]

    import pytest

    mp = pytest.MonkeyPatch()
    mp.setattr(orch, "make_model", fake_make_model)
    try:
        res = run_pipeline(cfg, model_name="forgetful", run_dir=tmp_path / "fin", stages=["planner", "analyst"], write_report=False)
    finally:
        mp.undo()
    m = made["m"]
    assert m.finishing_seen, "finishing turn was not triggered"
    titles = [x["title"] for x in res.facts["findings"]["summary"]]
    assert any("No disease-specific state" in t for t in titles)
