"""Run the full agentic pipeline.

Examples
    # offline, deterministic, no API key (CI / smoke test)
    SCTARGET_OFFLINE=1 python scripts/run_pipeline.py --config config/synthetic.yaml --model scripted

    # real agent run on your data
    python scripts/run_pipeline.py --config config/caa.yaml --model claude-sonnet-4-5

    # only some stages (e.g. re-run the critic on an existing run is NOT supported; stages
    # always start fresh)
    python scripts/run_pipeline.py --config config/caa.yaml --stages planner analyst
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
    ap.add_argument("--no-report", action="store_true")
    args = ap.parse_args(argv)

    overrides = {}
    if args.seed is not None:
        overrides = {"agent": {"seed": args.seed}}
    cfg = load_config(args.config, overrides)
    res = run_pipeline(cfg, model_name=args.model, run_dir=Path(args.run_dir) if args.run_dir else None, stages=args.stages, write_report=not args.no_report)
    print(f"run dir : {res.run_dir}")
    for role, rr in res.roles.items():
        print(f"  {role:8s} turns={rr.n_turns:2d} tools={len(rr.tool_calls):2d} errors={sum(1 for c in rr.tool_calls if c['error'])}  stop={rr.stopped_reason}")
    if res.report_md:
        print(f"report  : {res.report_md}")
    print(f"elapsed : {res.elapsed_s:.1f}s")


if __name__ == "__main__":
    main()
