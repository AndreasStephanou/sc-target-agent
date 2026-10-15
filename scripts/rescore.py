"""Re-score an existing run's candidate sets under a different scoring version, without re-running
the agent. Scoring is deterministic given the stored facts, so this isolates the scoring change
from model variance. It can only re-rank genes the nominator already sent to score_targets.

    python scripts/rescore.py --config config/caa.yaml --scoring-version 2 runs/caa_cc_v2_1
    # -> runs/caa_cc_v2_1_s2/ (state.json + tables, evaluable with run_eval.py) and a rank-change table
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sctarget.agents.orchestrator import STAGES, _resume_session  # noqa: E402
from sctarget.config import load_config  # noqa: E402
from sctarget.session import AnalysisSession  # noqa: E402
from sctarget.tools.targets import score_targets  # noqa: E402


def rescore(cfg, source: Path, out: Path) -> dict:
    session = AnalysisSession.create(cfg, out)
    _resume_session(session, source, keep_findings_from=set(STAGES))  # keep everything
    before = {cl: [r["gene"] for r in rows] for cl, rows in session.facts.get("ranked_targets", {}).items()}
    human = {g: v["human_symbol"] for g, v in session.facts.get("orthologs", {}).items() if v.get("human_symbol")}
    changes = {}
    for cl, genes in before.items():
        res = score_targets(session, genes=genes, cluster=cl, human_symbols=human)
        after = [r["gene"] for r in session.facts["ranked_targets"][cl]]
        cls = {r["gene"]: r.get("target_class") for r in session.facts["ranked_targets"][cl]}
        changes[cl] = {"axes": res.get("effect_axes"), "rows": [{"gene": g, "rank_before": genes.index(g) + 1, "rank_after": after.index(g) + 1, "target_class": cls.get(g)} for g in genes]}
    session.facts.pop("knowledge_client", None)
    session.dump_state()
    with open(out / "rescore.json", "w") as fh:
        json.dump({"source": str(source), "scoring_version": session.facts.get("scoring_version"), "changes": changes}, fh, indent=2)
    return changes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--scoring-version", type=int, required=True, choices=[1, 2])
    ap.add_argument("--suffix", default=None, help="output dir suffix (default _s<version>)")
    ap.add_argument("run_dirs", nargs="+")
    a = ap.parse_args(argv)
    cfg = load_config(a.config, {"targets": {"scoring_version": a.scoring_version}})
    for d in a.run_dirs:
        src = Path(d)
        if not src.is_dir():
            continue
        out = src.with_name(src.name + (a.suffix or f"_s{a.scoring_version}"))
        ch = rescore(cfg, src, out)
        print(f"\n{src} -> {out}")
        for cl, c in ch.items():
            print(f"  cluster {cl}  axes={c['axes']}")
            print(f"  {'gene':14s} {'before':>6s} {'after':>6s}  class")
            for r in sorted(c["rows"], key=lambda r: r["rank_after"]):
                mv = r["rank_before"] - r["rank_after"]
                arrow = f"(+{mv})" if mv > 0 else (f"({mv})" if mv < 0 else "")
                print(f"  {r['gene']:14s} {r['rank_before']:6d} {r['rank_after']:6d}  {r['target_class'] or '':14s} {arrow}")


if __name__ == "__main__":
    main()
