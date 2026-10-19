"""Re-score an existing run's candidate sets under a different scoring version, without re-running
the agent. Scoring is deterministic given the stored facts, so this isolates the scoring change
from model variance. It can only re-rank genes the nominator already sent to score_targets.

    python scripts/rescore.py --config config/caa.yaml --scoring-version 2 runs/caa_cc_v2_1
    # -> runs/caa_cc_v2_1_s2/ (state.json + tables, evaluable with run_eval.py) and a rank-change table

    # what would the run have nominated with current attribution rules and no premature drop?
    python scripts/rescore.py --config config/caa_denovo.yaml --scoring-version 2 --recompute-attribution --full-candidates --suffix _reattr runs/caa_denovo_1
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


def _examined_candidates(source: Path, cluster: str) -> list[str]:
    """Every gene the nominator ran gene_specificity on for this cluster, in first-seen order:
    the candidate set it examined, before it decided which genes to score."""
    seen: list[str] = []
    with open(source / "run.jsonl") as fh:
        for line in fh:
            e = json.loads(line)
            if e.get("type") == "tool_call" and e.get("role") == "target" and e.get("tool") == "gene_specificity" and str(e["args"].get("cluster")) == str(cluster):
                for g in e["args"].get("genes") or []:
                    if g not in seen:
                        seen.append(g)
    return seen


def rescore(cfg, source: Path, out: Path, recompute_attribution: bool = False, full_candidates: bool = False) -> dict:
    session = AnalysisSession.create(cfg, out)
    _resume_session(session, source, keep_findings_from=set(STAGES))  # keep everything
    before = {cl: [r["gene"] for r in rows] for cl, rows in session.facts.get("ranked_targets", {}).items()}
    human = {g: v["human_symbol"] for g, v in session.facts.get("orthologs", {}).items() if v.get("human_symbol")}
    changes = {}
    for cl, scored in before.items():
        genes = list(scored)
        if full_candidates:
            genes = list(dict.fromkeys(_examined_candidates(source, cl) + scored)) or scored
        attr_before = dict(session.facts.get("condition_attribution", {}).get(str(cl), {}))
        if recompute_attribution:
            from sctarget.tools.analysis import condition_attribution, gene_specificity

            missing_spec = [g for g in genes if g not in {r["gene"] for r in session.facts.get("gene_specificity", {}).get(str(cl), [])}]
            if missing_spec:
                gene_specificity(session, genes=missing_spec, cluster=cl)
            condition_attribution(session, genes=genes, cluster=cl)
        attr_after = session.facts.get("condition_attribution", {}).get(str(cl), {})
        res = score_targets(session, genes=genes, cluster=cl, human_symbols=human)
        after = [r["gene"] for r in session.facts["ranked_targets"][cl]]
        cls = {r["gene"]: r.get("target_class") for r in session.facts["ranked_targets"][cl]}
        rows = []
        for g in genes:
            if g not in after:  # excluded (unnamed) genes
                continue
            rows.append({
                "gene": g,
                "rank_before": (scored.index(g) + 1) if g in scored else None,
                "rank_after": after.index(g) + 1,
                "target_class": cls.get(g),
                "attribution_before": attr_before.get(g),
                "attribution_after": attr_after.get(g) if recompute_attribution else attr_before.get(g),
            })
        changes[cl] = {"axes": res.get("effect_axes"), "rows": rows}
    session.facts.pop("knowledge_client", None)
    session.dump_state()
    with open(out / "rescore.json", "w") as fh:
        json.dump({"source": str(source), "scoring_version": session.facts.get("scoring_version"), "recompute_attribution": recompute_attribution, "full_candidates": full_candidates, "changes": changes}, fh, indent=2)
    return changes


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--scoring-version", type=int, required=True, choices=[1, 2])
    ap.add_argument("--suffix", default=None, help="output dir suffix (default _s<version>)")
    ap.add_argument("--recompute-attribution", action="store_true", help="re-run condition_attribution (current rules, e.g. the thin reference-pool guard) before scoring; needs the dataset")
    ap.add_argument("--full-candidates", action="store_true", help="score every gene the nominator examined (its gene_specificity calls), not only the ones it chose to score")
    ap.add_argument("run_dirs", nargs="+")
    a = ap.parse_args(argv)
    cfg = load_config(a.config, {"targets": {"scoring_version": a.scoring_version}})
    for d in a.run_dirs:
        src = Path(d)
        if not src.is_dir():
            continue
        out = src.with_name(src.name + (a.suffix or f"_s{a.scoring_version}"))
        ch = rescore(cfg, src, out, recompute_attribution=a.recompute_attribution, full_candidates=a.full_candidates)
        print(f"\n{src} -> {out}")
        for cl, c in ch.items():
            print(f"  cluster {cl}  axes={c['axes']}")
            print(f"  {'gene':14s} {'before':>6s} {'after':>6s}  {'class':15s} {'attribution before -> after':32s}")
            for r in sorted(c["rows"], key=lambda r: r["rank_after"]):
                rb = r["rank_before"]
                mv = (rb - r["rank_after"]) if rb is not None else None
                arrow = "(new)" if rb is None else (f"(+{mv})" if mv > 0 else (f"({mv})" if mv < 0 else ""))
                attr = f"{r['attribution_before'] or '-'} -> {r['attribution_after'] or '-'}" if r["attribution_before"] != r["attribution_after"] else (r["attribution_after"] or "-")
                print(f"  {r['gene']:14s} {rb if rb is not None else '-':>6} {r['rank_after']:6d}  {r['target_class'] or '':15s} {attr:32s} {arrow}")


if __name__ == "__main__":
    main()
