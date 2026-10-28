"""The two gaps the HCC review run exposed (2026-10-28).

1. All 66 tests and the dense synthetic demo passed on an install that could not cluster real data
   (pynndescent 0.7 -> scipy csr_array -> scanpy .getnnz() AttributeError, only on the >8k-cell
   approximate-neighbours path and inside Scrublet). ``sctarget.envcheck`` exercises those paths.
2. The scripted analyst treated an errored tool like an empty result and reported "no disease-specific
   state found". It must abort and say the tool failed.
"""
from __future__ import annotations

import pytest

from sctarget.agents.model import ScriptedModel, ToolUse


@pytest.mark.slow
def test_envcheck_sparse_large_paths():
    from sctarget import envcheck

    failures = envcheck.run(n_cells=8500, n_genes=400, verbose=False)
    assert failures == [], "environment cannot run the sparse / approximate-neighbours tool chain: " + "; ".join(failures)


class _Session:
    facts: dict = {}
    cfg = None


def _messages_with_error(tool: str, err: str):
    return [
        {"role": "assistant", "content": [{"type": "tool_use", "id": "t1", "name": tool, "input": {}}]},
        {"role": "user", "content": [{"type": "tool_result", "tool_use_id": "t1", "content": err, "is_error": True}]},
    ]


def test_scripted_model_aborts_on_critical_tool_failure():
    m = ScriptedModel(_Session())
    m._queue["analyst"] = [ToolUse("cluster_composition", {}), ToolUse("differential_abundance", {}), "No disease-specific states found."]
    tools = [{"name": n} for n in ("record_finding", "cluster_composition", "differential_abundance", "integrate_and_cluster")]
    r = m.create(system="", messages=_messages_with_error("integrate_and_cluster", "AttributeError: 'csr_array' object has no attribute 'getnnz'"), tools=tools, role="analyst")
    assert r.tool_uses and r.tool_uses[0].name == "record_finding"
    assert "Run aborted" in r.tool_uses[0].input["title"]
    assert "integrate_and_cluster" in r.tool_uses[0].input["title"]
    assert "csr_array" in r.tool_uses[0].input["detail"]
    # the next turn ends the role with an explicit abort message, not "no state found"
    r2 = m.create(system="", messages=[], tools=tools, role="analyst")
    assert r2.stop_reason == "end_turn" and "aborted" in r2.text.lower()
    assert "No disease-specific states found" not in r2.text


def test_scripted_model_ignores_non_critical_tool_failure():
    m = ScriptedModel(_Session())
    m._queue["analyst"] = [ToolUse("cluster_composition", {}), "done"]
    tools = [{"name": n} for n in ("record_finding", "cluster_composition", "lookup_druggability")]
    r = m.create(system="", messages=_messages_with_error("lookup_druggability", "HTTP 429"), tools=tools, role="analyst")
    assert r.tool_uses and r.tool_uses[0].name == "cluster_composition"
