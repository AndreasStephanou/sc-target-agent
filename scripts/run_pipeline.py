"""Run the full agentic pipeline.

Examples
    # offline, deterministic, no API key (CI / smoke test)
    SCTARGET_OFFLINE=1 python scripts/run_pipeline.py --config config/synthetic.yaml --model scripted

    # real agent run on your data
    python scripts/run_pipeline.py --config config/caa.yaml --model claude-sonnet-4-5

    # no API key: drive the pipeline through a Claude Code subscription login (separate experimental arm)
    python scripts/run_pipeline.py --config config/caa.yaml --model claude-code            # or claude-code:sonnet

    # re-run only the nominator + critic on an existing run's analysis (no re-clustering / DE)
    python scripts/run_pipeline.py --config config/caa.yaml --model scripted --resume-from runs/caa_dry6 --run-dir runs/caa_dry6b
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sctarget.agents.orchestrator import STAGES, run_pipeline  # noqa: E402
from sctarget.config import load_config  # noqa: E402


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, help="dataset config yaml")
    ap.add_argument("--model", default=None, help="'scripted' or a Claude model id (default: agent.model from config)")
    ap.add_argument("--run-dir", default=None, help="output directory (default: runs/<auto>)")
    ap.add_argument("--stages", nargs="+", default=None, choices=STAGES)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--scoring-version", type=int, default=None, choices=[1, 2, 3], help="override targets.scoring_version (2 = shared-cluster lineage axis; 3 = all-cells lineage contrast when the shared one is unrepresentative, plus proliferation flag)")
    ap.add_argument("--no-report", action="store_true")
    ap.add_argument("--resume-from", default=None, help="reuse an earlier run directory's analysis (state.json, tables) and run only --stages (default: target critic)")
    args = ap.parse_args(argv)

    overrides: dict = {}
    if args.seed is not None:
        overrides["agent"] = {"seed": args.seed}
    if args.scoring_version is not None:
        overrides["targets"] = {"scoring_version": args.scoring_version}
    cfg = load_config(args.config, overrides)
    res = run_pipeline(cfg, model_name=args.model, run_dir=Path(args.run_dir) if args.run_dir else None, stages=args.stages, write_report=not args.no_report, resume_from=Path(args.resume_from) if args.resume_from else None)
    print(f"run dir : {res.run_dir}")
    for role, rr in res.roles.items():
        print(f"  {role:8s} turns={rr.n_turns:2d} tools={len(rr.tool_calls):2d} errors={sum(1 for c in rr.tool_calls if c['error'])}  stop={rr.stopped_reason}")
    if res.report_md:
        print(f"report  : {res.report_md}")
    print(f"elapsed : {res.elapsed_s:.1f}s")


if __name__ == "__main__":
    main()
