"""The tool-calling loop. One function, no framework.

For each role: send system prompt + conversation to the model; execute every tool call
it returns against the session (validated, sanitized, logged); append results; repeat
until the model stops calling tools or the turn budget is spent.
"""
from __future__ import annotations

import json
import time
import traceback
from dataclasses import dataclass, field
from typing import Any

from ..privacy import PrivacyViolation, sanitize
from ..session import AnalysisSession
from ..tools.registry import ToolRegistry, validate_args
from .model import Model


@dataclass
class RoleResult:
    role: str
    final_text: str
    n_turns: int
    tool_calls: list[dict[str, Any]] = field(default_factory=list)
    stopped_reason: str = "end_turn"


def _result_payload(obj: Any) -> str:
    return json.dumps(obj, default=str)


def run_role(
    *,
    role: str,
    system: str,
    user_prompt: str,
    session: AnalysisSession,
    model: Model,
    registry: ToolRegistry,
    max_turns: int = 25,
) -> RoleResult:
    tools = registry.schemas_for_role(role)
    messages: list[dict[str, Any]] = [{"role": "user", "content": user_prompt}]
    calls: list[dict[str, Any]] = []
    session.log.stage(role, "start", model=model.name, n_tools=len(tools))

    for turn in range(1, max_turns + 1):
        resp = model.create(system=system, messages=messages, tools=tools, role=role)
        session.log.model_call(role, model.name, len(messages), resp.usage, resp.stop_reason, resp.text or None)
        messages.append({"role": "assistant", "content": resp.raw_content})

        if not resp.tool_uses:
            session.log.stage(role, "end", n_turns=turn, reason=resp.stop_reason)
            return RoleResult(role=role, final_text=resp.text, n_turns=turn, tool_calls=calls, stopped_reason=resp.stop_reason)

        results_content: list[dict[str, Any]] = []
        for tu in resp.tool_uses:
            t0 = time.time()
            error: str | None = None
            sent: Any
            try:
                if tu.name not in registry or not registry.get(tu.name).allowed_for(role):
                    raise PermissionError(f"tool '{tu.name}' is not available to role '{role}'")
                spec = registry.get(tu.name)
                args = validate_args(spec, dict(tu.input))
                raw = spec.fn(session, **args)
                sent = sanitize(raw, session.policy)
            except PrivacyViolation as e:
                error = f"PrivacyViolation: {e}"
                sent = {"error": error}
            except Exception as e:  # tool failure is information for the agent, not a crash
                error = f"{type(e).__name__}: {e}"
                sent = {"error": error, "hint": "Check the arguments against the tool schema, or choose a different tool."}
                session.log.error(role, traceback.format_exc(limit=3))
            dur = time.time() - t0
            session.log.tool_call(role, tu.name, dict(tu.input), sent, dur, error)
            calls.append({"tool": tu.name, "args": dict(tu.input), "error": error, "duration_s": round(dur, 3)})
            results_content.append({"type": "tool_result", "tool_use_id": tu.id, "content": _result_payload(sent), "is_error": error is not None})
        messages.append({"role": "user", "content": results_content})

    session.log.stage(role, "end", n_turns=max_turns, reason="max_turns")
    return RoleResult(role=role, final_text="(turn budget exhausted)", n_turns=max_turns, tool_calls=calls, stopped_reason="max_turns")
