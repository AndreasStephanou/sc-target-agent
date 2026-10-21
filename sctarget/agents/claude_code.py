"""Model backend that drives the Claude Code CLI (``claude -p``) instead of the Messages API.

Why it exists: it lets the pipeline run on a Claude subscription login, with no API key.
It is a DIFFERENT EXPERIMENTAL ARM from ``AnthropicModel``: Claude Code adds its own system
prompt, there is no temperature/seed control, and tool calling is emulated through a strict
JSON reply format. Report it as "Claude Code backend over the same tool registry".

How it works (stateless, one subprocess per turn):
  * the role's system prompt plus a tool-protocol block (tool names, descriptions, JSON
    schemas) is passed with ``--system-prompt``;
  * the whole transcript so far is rendered to text and sent on stdin;
  * Claude Code's own built-in tools are disabled (``--tools ""``), MCP servers are blocked,
    and the subprocess runs in an empty temp directory, so the model can only act by
    returning JSON that this loop executes through the registry (and its privacy guard);
  * the reply must be one JSON object: {"text": str, "tool_calls": [{"name": str, "input": {...}}]}.
    No tool calls means the role is finished.

ANTHROPIC_API_KEY is removed from the child's environment by default so that the CLI uses the
subscription login rather than silently billing an API key (set ``keep_api_key=True`` to keep it).
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
import uuid
from typing import Any

from .model import ModelResponse, ToolUse

PROTOCOL = """

=== TOOL PROTOCOL (read carefully) ===
You cannot call tools natively in this session. You act ONLY by replying with ONE JSON object
and nothing else (no prose outside it, no markdown fence needed):

{"text": "<short reasoning or the final summary>", "tool_calls": [{"name": "<tool name>", "input": {<arguments>}}]}

- To use tools, put 1-4 calls in "tool_calls". Calls in one reply must be independent of each
  other. You will receive their results in the next message as TOOL_RESULT blocks.
- When you are finished with your role, reply with an EMPTY "tool_calls" list and put your
  final summary in "text".
- Use only the tools listed below, with arguments that match their schemas exactly.
- Do not use any built-in tools you may believe you have (file, shell, web): they are disabled.

Available tools (name, description, JSON schema of the input):
"""


def render_tools(tools: list[dict[str, Any]]) -> str:
    return "\n".join(
        json.dumps({"name": t["name"], "description": t.get("description", ""), "input_schema": t.get("input_schema", {})}, separators=(",", ":"))
        for t in tools
    )


def render_transcript(messages: list[dict[str, Any]]) -> str:
    parts: list[str] = []
    for m in messages:
        c = m["content"]
        if isinstance(c, str):
            parts.append(f"[{m['role'].upper()}]\n{c}")
            continue
        for b in c:
            t = b.get("type")
            if t == "text":
                parts.append(f"[{m['role'].upper()}]\n{b['text']}")
            elif t == "tool_use":
                parts.append(f"[ASSISTANT TOOL_CALL id={b['id']}] {b['name']} {json.dumps(b.get('input', {}))}")
            elif t == "tool_result":
                flag = " ERROR" if b.get("is_error") else ""
                parts.append(f"[TOOL_RESULT id={b['tool_use_id']}{flag}]\n{b['content']}")
    parts.append("[YOUR TURN] Reply with the single JSON object described in the tool protocol.")
    return "\n\n".join(parts)


def extract_json(text: str) -> dict[str, Any] | None:
    """First balanced top-level JSON object in ``text`` (tolerates a ```json fence or stray prose)."""
    text = text.strip()
    fence = re.search(r"```(?:json)?\s*(\{.*\})\s*```", text, re.S)
    if fence:
        text = fence.group(1)
    start = text.find("{")
    while start != -1:
        depth, in_str, esc = 0, False, False
        for i in range(start, len(text)):
            ch = text[i]
            if in_str:
                if esc:
                    esc = False
                elif ch == "\\":
                    esc = True
                elif ch == '"':
                    in_str = False
            elif ch == '"':
                in_str = True
            elif ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    try:
                        obj = json.loads(text[start : i + 1])
                        if isinstance(obj, dict):
                            return obj
                    except json.JSONDecodeError:
                        pass
                    break
        start = text.find("{", start + 1)
    return None


class ClaudeCodeModel:
    def __init__(
        self,
        model: str | None = None,
        *,
        command: str = "claude",
        timeout_s: int = 900,
        max_parse_retries: int = 2,
        keep_api_key: bool = False,
        extra_args: list[str] | None = None,
    ) -> None:
        self.model = model
        self.name = f"claude-code:{model}" if model else "claude-code"
        self.command = command
        self.timeout_s = timeout_s
        self.max_parse_retries = max_parse_retries
        self.keep_api_key = keep_api_key
        self.extra_args = extra_args or []
        self._cwd = tempfile.mkdtemp(prefix="sctarget_cc_")

    # ---- subprocess
    def _argv(self, system: str) -> list[str]:
        argv = [self.command, "-p", "--output-format", "json", "--tools", "", "--strict-mcp-config", "--disable-slash-commands", "--no-session-persistence", "--system-prompt", system]
        if self.model:
            argv += ["--model", self.model]
        return argv + self.extra_args

    def _env(self) -> dict[str, str]:
        env = dict(os.environ)
        if not self.keep_api_key:
            env.pop("ANTHROPIC_API_KEY", None)
            env.pop("ANTHROPIC_AUTH_TOKEN", None)
        return env

    def _call(self, system: str, prompt: str) -> tuple[str, dict[str, int]]:
        proc = subprocess.run(self._argv(system), input=prompt, capture_output=True, text=True, timeout=self.timeout_s, cwd=self._cwd, env=self._env())
        if proc.returncode != 0:
            raise RuntimeError(f"claude CLI exited {proc.returncode}: {(proc.stderr or proc.stdout)[-800:]}")
        out = proc.stdout.strip()
        try:
            env = json.loads(out)
        except json.JSONDecodeError:
            return out, {"input_tokens": 0, "output_tokens": 0}
        if isinstance(env, list):  # some versions emit the event list
            env = next((e for e in reversed(env) if isinstance(e, dict) and e.get("type") == "result"), {})
        if env.get("is_error"):
            raise RuntimeError(f"claude CLI reported an error: {str(env.get('result'))[:800]}")
        u = env.get("usage") or {}
        usage = {
            "input_tokens": int(u.get("input_tokens", 0)) + int(u.get("cache_creation_input_tokens", 0)) + int(u.get("cache_read_input_tokens", 0)),
            "output_tokens": int(u.get("output_tokens", 0)),
        }
        return str(env.get("result", "")), usage

    # ---- Model interface
    def create(self, *, system: str, messages: list[dict[str, Any]], tools: list[dict[str, Any]], role: str) -> ModelResponse:
        full_system = system + PROTOCOL + render_tools(tools)
        prompt = render_transcript(messages)
        total = {"input_tokens": 0, "output_tokens": 0}
        last_text = ""
        for attempt in range(self.max_parse_retries + 1):
            text, usage = self._call(full_system, prompt)
            for k in total:
                total[k] += usage.get(k, 0)
            last_text = text
            obj = extract_json(text)
            problem = None
            if obj is None:
                problem = "your reply was not a single JSON object"
            else:
                calls = obj.get("tool_calls") or []
                if not isinstance(calls, list) or any(not isinstance(c, dict) or "name" not in c for c in calls):
                    problem = "tool_calls must be a list of {\"name\":..., \"input\":{...}} objects"
                elif not calls and not str(obj.get("text") or "").strip():
                    # an empty reply is not a decision: the model must either act or write its summary
                    problem = "your reply had no tool_calls and no text; either call the next tools or write your final summary (and record your findings first)"
            if problem is None:
                break
            prompt += f"\n\n[SYSTEM NOTICE] {problem}. Reply again with ONLY the JSON object described in the tool protocol."
        else:
            raise RuntimeError(f"claude CLI did not return a parseable tool-protocol reply after {self.max_parse_retries + 1} attempts: {last_text[:400]!r}")

        say = str(obj.get("text") or "")
        tool_uses = [
            ToolUse(name=str(c["name"]), input=dict(c.get("input") or {}), id=f"toolu_{uuid.uuid4().hex[:16]}")
            for c in (obj.get("tool_calls") or [])
        ]
        # unknown tool names are passed through: the loop reports them back to the model as tool errors
        raw: list[dict[str, Any]] = [{"type": "text", "text": say}] if say else []
        raw += [{"type": "tool_use", "id": t.id, "name": t.name, "input": t.input} for t in tool_uses]
        if not raw:
            raw = [{"type": "text", "text": ""}]
        return ModelResponse(text=say, tool_uses=tool_uses, stop_reason="tool_use" if tool_uses else "end_turn", usage=total, raw_content=raw)
