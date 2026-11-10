"""Token totals (report header, run_eval) sum only the numeric usage fields: since protocol 23r a
model_call's usage also carries the model id the backend reported, which must not break scoring."""
from __future__ import annotations

from sctarget.runlog import usage_totals_from_events


def test_usage_totals_skip_non_numeric_fields():
    events = [
        {"type": "stage", "stage": "planner", "status": "start"},
        {"type": "model_call", "role": "planner", "usage": {"input_tokens": 100, "output_tokens": 7, "model": "claude-fable-5-1"}},
        {"type": "tool_call", "role": "planner", "tool": "load_dataset"},
        {"type": "model_call", "role": "analyst", "usage": {"input_tokens": 250, "output_tokens": 3.0, "model": "claude-fable-5-1", "cached": True}},
        {"type": "model_call", "role": "critic", "usage": None},
    ]
    assert usage_totals_from_events(events) == {"input_tokens": 350, "output_tokens": 10}
