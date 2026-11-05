#!/usr/bin/env python3
"""Regenerate report.md / report.html (and the report-time figures) for an existing run directory.

Reloads the dataset and the run's facts, tables and cluster labels exactly as scripts/rescore.py does,
then builds the report with the current report code. No analysis is repeated and no fact changes;
figures are added (PNG + PDF + SVG) and the two report files are overwritten in place.

    python scripts/rebuild_report.py --config config/ipf.yaml runs/ipf_calibrate_1 [runs/ipf_agent_1 ...]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sctarget.agents.orchestrator import STAGES, _resume_session  # noqa: E402
from sctarget.config import load_config  # noqa: E402
from sctarget.report.builder import build_report  # noqa: E402
from sctarget.session import AnalysisSession  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("run_dirs", nargs="+")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    for rd in a.run_dirs:
        rd = Path(rd)
        session = AnalysisSession.create(cfg, rd)
        _resume_session(session, rd, keep_findings_from=set(STAGES))
        md, html = build_report(session)
        print(f"{rd}: wrote {md.name}, {html.name}; {len(session.figures)} figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
