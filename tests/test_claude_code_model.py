"""ClaudeCodeModel is tested against a fake `claude` executable: no login, no network."""
import json
import os
import stat
import sys

import pytest

from sctarget.agents.claude_code import ClaudeCodeModel, extract_json, render_transcript
from sctarget.agents.loop import run_role
from sctarget.session import AnalysisSession
from sctarget.tools import REGISTRY

TOOLS = [{"name": "load_dataset", "description": "load", "input_schema": {"type": "object", "properties": {}}}]


def make_fake(tmp_path, replies, record=None):
    """Executable that prints the next canned reply (as a claude --output-format json envelope)."""
    q = tmp_path / "queue.json"
    q.write_text(json.dumps(replies))
    rec = tmp_path / "calls.jsonl"
    script = tmp_path / "claude"
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys, json, os\n"
        f"q = {str(q)!r}; rec = {str(rec)!r}\n"
        "prompt = sys.stdin.read()\n"
        "replies = json.load(open(q)); r = replies.pop(0); json.dump(replies, open(q, 'w'))\n"
        "open(rec, 'a').write(json.dumps({'argv': sys.argv[1:], 'prompt': prompt, 'has_key': 'ANTHROPIC_API_KEY' in os.environ, 'cwd': os.getcwd()}) + '\\n')\n"
        "print(json.dumps({'type': 'result', 'is_error': False, 'result': r, 'usage': {'input_tokens': 10, 'cache_read_input_tokens': 90, 'output_tokens': 7}}))\n"
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return str(script), rec


def calls(rec):
    return [json.loads(l) for l in rec.read_text().splitlines()]


def test_extract_json_variants():
    assert extract_json('{"a": 1}') == {"a": 1}
    assert extract_json('Sure!\n```json\n{"a": {"b": "}"}}\n```') == {"a": {"b": "}"}}
    assert extract_json('noise {not json} then {"ok": true} tail') == {"ok": True}
    assert extract_json("no json here") is None


def test_parses_tool_calls_and_final(tmp_path, monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "sk-should-not-leak")
    cmd, rec = make_fake(tmp_path, [json.dumps({"text": "loading", "tool_calls": [{"name": "load_dataset", "input": {}}]}), json.dumps({"text": "done", "tool_calls": []})])
    m = ClaudeCodeModel(model="sonnet", command=cmd)
    r1 = m.create(system="SYS", messages=[{"role": "user", "content": "go"}], tools=TOOLS, role="planner")
    assert r1.stop_reason == "tool_use" and r1.tool_uses[0].name == "load_dataset" and r1.text == "loading"
    assert r1.raw_content[1]["type"] == "tool_use" and r1.raw_content[1]["id"] == r1.tool_uses[0].id
    assert r1.usage == {"input_tokens": 100, "output_tokens": 7}  # cache reads counted as input
    r2 = m.create(system="SYS", messages=[{"role": "user", "content": "go"}], tools=TOOLS, role="planner")
    assert r2.stop_reason == "end_turn" and not r2.tool_uses and r2.text == "done"
    c = calls(rec)[0]
    assert "--tools" in c["argv"] and c["argv"][c["argv"].index("--tools") + 1] == ""
    assert "--strict-mcp-config" in c["argv"] and "--model" in c["argv"] and "sonnet" in c["argv"]
    assert not c["has_key"], "API key must not reach the CLI by default"
    assert c["cwd"] != os.getcwd()
    assert '"name":"load_dataset"' in c["argv"][c["argv"].index("--system-prompt") + 1]


def test_retries_on_garbage_then_succeeds(tmp_path):
    cmd, rec = make_fake(tmp_path, ["I will now think about it.", json.dumps({"text": "ok", "tool_calls": []})])
    m = ClaudeCodeModel(command=cmd)
    r = m.create(system="S", messages=[{"role": "user", "content": "go"}], tools=TOOLS, role="critic")
    assert r.text == "ok" and len(calls(rec)) == 2
    assert "not a single JSON object" in calls(rec)[1]["prompt"]
    assert r.usage["output_tokens"] == 14  # both attempts counted


def test_gives_up_loudly(tmp_path):
    cmd, _ = make_fake(tmp_path, ["nope"] * 3)
    with pytest.raises(RuntimeError, match="parseable"):
        ClaudeCodeModel(command=cmd).create(system="S", messages=[{"role": "user", "content": "go"}], tools=TOOLS, role="critic")


def test_transcript_rendering():
    msgs = [
        {"role": "user", "content": "start"},
        {"role": "assistant", "content": [{"type": "text", "text": "hi"}, {"type": "tool_use", "id": "t1", "name": "load_dataset", "input": {}}]},
        {"role": "user", "content": [{"type": "tool_result", "tool_use_id": "t1", "content": "{\"n\": 3}", "is_error": True}]},
    ]
    t = render_transcript(msgs)
    assert "TOOL_CALL id=t1" in t and "TOOL_RESULT id=t1 ERROR" in t and t.rstrip().endswith("tool protocol.")


def test_runs_inside_the_real_loop(cfg, tmp_path):
    """Emulated tool calls go through the registry, validation and privacy guard like native ones."""
    replies = [
        json.dumps({"text": "", "tool_calls": [{"name": "load_dataset", "input": {}}, {"name": "druggability_lookup", "input": {"genes": ["X"]}}]}),
        json.dumps({"text": "finished", "tool_calls": []}),
    ]
    cmd, rec = make_fake(tmp_path, replies)
    s = AnalysisSession.create(cfg, tmp_path / "run")
    rr = run_role(role="analyst", system="x", user_prompt="go", session=s, model=ClaudeCodeModel(command=cmd), registry=REGISTRY, max_turns=5)
    assert rr.stopped_reason == "end_turn" and rr.final_text == "finished" and rr.n_turns == 2
    assert [c["error"] is None for c in rr.tool_calls] == [True, False]  # analyst may not call druggability_lookup
    assert "not available to role" in calls(rec)[1]["prompt"]  # the model was told, as with native tool use
