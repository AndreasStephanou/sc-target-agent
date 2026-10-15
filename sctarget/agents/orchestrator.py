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
from .prompts import PROMPT_VERSION, ROLE_SYSTEM, stage_prompt

STAGES = ["planner", "analyst", "target", "critic"]


@dataclass
class PipelineResult:
    run_dir: Path
    roles: dict[str, RoleResult] = field(default_factory=dict)
    report_md: Path | None = None
    report_html: Path | None = None
    facts: dict[str, Any] = field(default_factory=dict)
    elapsed_s: float = 0.0


def _resume_session(session: AnalysisSession, source: Path, keep_findings_from: set[str]) -> None:
    """Load facts from an earlier run's state.json, reload the data, and restore the cluster key,
    so later stages can run without repeating the analysis. Findings produced by stages that will
    be re-run are dropped; the rest (e.g. the analyst's disease states) are kept."""
    import shutil

    with open(source / "state.json") as fh:
        state = json.load(fh)
    facts = state["facts"]
    drop = {"planner": {"summary"}, "target": {"candidate_target"}, "critic": {"critic_flag"}}
    for role in set(drop) - keep_findings_from:
        for kind in drop[role]:
            facts.get("findings", {}).pop(kind, None)
        facts.get("role_summaries", {}).pop(role, None)
    if "target" not in keep_findings_from:  # the nominator's derived facts are recomputed
        for k in ("ranked_targets", "gene_specificity", "condition_attribution", "druggability", "literature", "orthologs"):
            facts.pop(k, None)
    if "critic" not in keep_findings_from:
        facts.get("findings", {})["summary"] = [x for x in facts.get("findings", {}).get("summary", []) if x.get("title") != "Critic summary"]
        facts.pop("differential_abundance_within_celltype", None)
    for sub in ("tables", "figures"):
        if (source / sub).exists():
            for f in (source / sub).iterdir():
                shutil.copy(f, session.run_dir / sub / f.name)
    for v in facts.get("pseudobulk_de", {}).values():
        if v.get("table_path"):
            v["table_path"] = str(session.run_dir / "tables" / Path(v["table_path"]).name)
    session.facts.update(facts)
    from ..tools.analysis import load_dataset

    load_dataset(session)
    key = state.get("cluster_key")
    a = session.adata
    if key and key not in a.obs and (source / "cluster_assignments.csv").exists():
        import pandas as pd

        assign = pd.read_csv(source / "cluster_assignments.csv", index_col=0)["cluster"].astype(str)
        a.obs[key] = assign.reindex(a.obs_names).astype(str).to_numpy()
    session.cluster_key = key
    session.figures = state.get("figures", [])
    session.log.stage("resume", "loaded", source=str(source), cluster_key=key, kept_findings=sorted(keep_findings_from))


def run_pipeline(
    cfg: Config,
    *,
    model_name: str | None = None,
    run_dir: Path | None = None,
    stages: list[str] | None = None,
    model: Model | None = None,
    write_report: bool = True,
    resume_from: Path | None = None,
) -> PipelineResult:
    t0 = time.time()
    run_dir = Path(run_dir) if run_dir else cfg.run_root / new_run_id()
    session = AnalysisSession.create(cfg, run_dir)
    if resume_from is not None:
        stages = stages or ["target", "critic"]
        _resume_session(session, Path(resume_from), keep_findings_from=set(STAGES) - set(stages))
    model = model or make_model(model_name or cfg.agent.get("model", "scripted"), session, temperature=float(cfg.agent.get("temperature", 0.2)))
    session.facts["prompt_version"] = PROMPT_VERSION
    session.log.stage("pipeline", "start", model=model.name, prompt_version=PROMPT_VERSION, config=cfg.raw)
    with open(run_dir / "config_used.json", "w") as fh:
        json.dump({**cfg.raw, "prompt_version": PROMPT_VERSION}, fh, indent=2, default=str)

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
