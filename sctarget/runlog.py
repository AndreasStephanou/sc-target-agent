"""Append-only JSONL run log. Every model call and every tool call is recorded.

Each run directory has::

    runs/<run_id>/
        run.jsonl          # ordered events: model_call, tool_call, stage, note, error
        tables/            # full result tables written by tools
        figures/           # png figures written by tools
        state.json         # final structured RunState
        report.md / .html

The log is the audit trail: a reviewer can reconstruct how any conclusion in the
report was reached by following the tool calls in order.
"""
from __future__ import annotations

import datetime as dt
import json
import time
import uuid
from pathlib import Path
from typing import Any


def _now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="milliseconds")


def new_run_id(prefix: str = "run") -> str:
    return f"{prefix}_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"


class RunLog:
    def __init__(self, run_dir: Path) -> None:
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        (self.run_dir / "tables").mkdir(exist_ok=True)
        (self.run_dir / "figures").mkdir(exist_ok=True)
        self.path = self.run_dir / "run.jsonl"
        self._seq = 0
        self._t0 = time.time()

    # ------------------------------------------------------------------ writing
    def _write(self, event: dict[str, Any]) -> int:
        self._seq += 1
        event = {"seq": self._seq, "t": _now(), "elapsed_s": round(time.time() - self._t0, 3), **event}
        with open(self.path, "a") as fh:
            fh.write(json.dumps(event, default=str) + "\n")
        return self._seq

    def stage(self, name: str, status: str, **extra: Any) -> int:
        return self._write({"type": "stage", "stage": name, "status": status, **extra})

    def note(self, role: str, text: str) -> int:
        return self._write({"type": "note", "role": role, "text": text})

    def model_call(self, role: str, model: str, n_input_msgs: int, usage: dict | None, stop_reason: str | None, text: str | None) -> int:
        return self._write(
            {
                "type": "model_call",
                "role": role,
                "model": model,
                "n_input_msgs": n_input_msgs,
                "usage": usage or {},
                "stop_reason": stop_reason,
                "text": text,
            }
        )

    def tool_call(self, role: str, tool: str, args: dict[str, Any], result_sent: Any, duration_s: float, error: str | None = None) -> int:
        return self._write(
            {
                "type": "tool_call",
                "role": role,
                "tool": tool,
                "args": args,
                "result_sent": result_sent,  # exactly what the model saw (post-sanitize)
                "duration_s": round(duration_s, 3),
                "error": error,
            }
        )

    def error(self, role: str, text: str) -> int:
        return self._write({"type": "error", "role": role, "text": text})

    # ------------------------------------------------------------------ reading
    def events(self) -> list[dict[str, Any]]:
        if not self.path.exists():
            return []
        with open(self.path) as fh:
            return [json.loads(line) for line in fh if line.strip()]

    def tool_calls(self) -> list[dict[str, Any]]:
        return [e for e in self.events() if e["type"] == "tool_call"]

    def usage_totals(self) -> dict[str, int]:
        tot: dict[str, int] = {}
        for e in self.events():
            if e["type"] == "model_call":
                for k, v in (e.get("usage") or {}).items():
                    if isinstance(v, (int, float)):
                        tot[k] = tot.get(k, 0) + int(v)
        return tot
