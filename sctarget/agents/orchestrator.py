"""Runs the roles in order and produces the report.

planner -> analyst -> target -> critic -> report

Each role gets its own conversation (fresh context) and a stage prompt built from the
structured facts the previous roles recorded. State flows through ``session.facts`` and
the run log, never through free text.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..config import Config
from ..runlog import new_run_id
from ..session import AnalysisSession
from ..tools import REGISTRY
from .loop import RoleResult, run_role
from .model import Model, make_model
from .prompts import ROLE_SYSTEM, stage_prompt

STAGES = ["planner", "analyst", "target", "critic"]


@dataclass
class PipelineResult:
    run_dir: Path
    roles: dict[str, RoleResult] = field(default_factory=dict)
    report_md: Path | None = None
    report_html: Path | None = None
    facts: dict[str, Any] = field(default_factory=dict)
    elapsed_s: float = 0.0


def run_pipeline(
    cfg: Config,
    *,
    model_name: str | None = None,
    run_dir: Path | None = None,
    stages: list[str] | None = None,
    model: Model | None = None,
    write_report: bool = True,
) -> PipelineResult:
    t0 = time.time()
    run_dir = Path(run_dir) if run_dir else cfg.run_root / new_run_id()
    session = AnalysisSession.create(cfg, run_dir)
    model = model or make_model(model_name or cfg.agent.get("model", "scripted"), session, temperature=float(cfg.agent.get("temperature", 0.2)))
    session.log.stage("pipeline", "start", model=model.name, config=cfg.raw)
    with open(run_dir / "config_used.json", "w") as fh:
        json.dump(cfg.raw, fh, indent=2, default=str)

    result = PipelineResult(run_dir=run_dir)
    for role in stages or STAGES:
        rr = run_role(
            role=role,
            system=ROLE_SYSTEM[role],
            user_prompt=stage_prompt(role, session.facts, cfg.raw),
            session=session,
            model=model,
            registry=REGISTRY,
            max_turns=int(cfg.agent.get("max_turns_per_role", 25)),
        )
        result.roles[role] = rr
        session.facts.setdefault("role_summaries", {})[role] = rr.final_text
        session.dump_state()

    # knowledge client isn't JSON-serialisable; drop before persisting
    session.facts.pop("knowledge_client", None)
    session.dump_state()

    if write_report:
        from ..report.builder import build_report

        md, html = build_report(session)
        result.report_md, result.report_html = md, html
    result.facts = session.facts
    result.elapsed_s = time.time() - t0
    session.log.stage("pipeline", "end", elapsed_s=round(result.elapsed_s, 1), usage=session.log.usage_totals())
    return result
