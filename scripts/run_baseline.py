"""Non-agentic baseline: the same tools, a fixed naive recipe, no model and no critic.

Recipe (deliberately what a quick first-pass script would do):
  load -> (reuse) cluster -> differential abundance -> take the most disease-enriched
  cluster -> markers vs ALL other cells -> top 30 by score -> druggability lookup ->
  rank by log2FC * druggable, nominate top 5.

It writes the same state.json layout as the agent so `scripts/run_eval.py` can score both.
The comparison agent-vs-baseline is the headline result of the project: does the agent's
sibling-cluster contrast, specificity check and critic reduce false leads?
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sctarget.config import load_config  # noqa: E402
from sctarget.runlog import new_run_id  # noqa: E402
from sctarget.session import AnalysisSession  # noqa: E402
from sctarget.tools import analysis as A  # noqa: E402
from sctarget.tools import knowledge as K  # noqa: E402


def run_baseline(cfg, run_dir: Path) -> Path:
    s = AnalysisSession.create(cfg, run_dir)
    s.log.stage("baseline", "start")
    ds = A.load_dataset(s)
    if not ds["already_processed_flag"]:
        A.run_qc(s, **cfg.analysis.get("qc", {}))
    reuse = ds["existing_clusterings"][0] if ds["already_processed_flag"] and ds["existing_clusterings"] else None
    integ = cfg.analysis.get("integration", {})
    A.integrate_and_cluster(s, reuse_existing=reuse, use_harmony=bool(integ.get("use_harmony", True)), resolution=float(integ.get("resolution", 1.0)))
    A.cluster_composition(s)
    da = A.differential_abundance(s)
    table = da["table"]
    top_cluster = str(table.index[0])  # most enriched by log2 ratio
    mk = A.find_markers(s, cluster=top_cluster, reference_clusters=None, n_top=30)
    genes = list(mk["top_markers"].index)
    species = cfg.dataset.get("species", "homo_sapiens")
    human = genes
    if species not in ("homo_sapiens", "human"):
        orth = K.map_orthologs(s, genes=genes)["mapping"]
        human = [orth[g]["human_symbol"] or g for g in genes]
    dr = K.druggability_lookup(s, genes=human)["table"]
    rows = []
    for g, h in zip(genes, human):
        lfc = float(mk["top_markers"].loc[g, "log2fc"])
        d = dr.loc[h] if h in dr.index else None
        druggable = bool(d is not None and d["is_druggable_class"])
        score = max(lfc, 0.0) * (1.0 if druggable else 0.2)
        rows.append({"gene": g, "human_symbol": h, "score": round(score, 4), "effect_size_log2fc": round(lfc, 3), "effect_size_source": "find_markers", "log2_specificity": None, "sample_consistency": None, "druggability_score": 1.0 if druggable else 0.0, "known_compound_score": 0.0, "top_inhibitors": list(d["top_inhibitors"]) if d is not None else [], "flags": []})
    rows.sort(key=lambda r: -r["score"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    s.facts["ranked_targets"] = {top_cluster: rows}
    s.facts.setdefault("findings", {})["disease_state"] = [{"title": f"Cluster {top_cluster} (baseline)", "detail": "most disease-enriched cluster by log2 ratio", "cluster": top_cluster, "genes": genes[:15], "confidence": None, "severity": None}]
    s.facts["findings"]["candidate_target"] = [{"title": f"{r['gene']} (baseline rank {r['rank']})", "detail": "baseline: log2FC x druggable", "cluster": top_cluster, "genes": [r["gene"]], "confidence": None, "severity": None} for r in rows[:5] if r["druggability_score"] > 0]
    s.facts.pop("knowledge_client", None)
    s.dump_state()
    s.log.stage("baseline", "end")
    from sctarget.report.builder import build_report

    build_report(s)
    return run_dir


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--run-dir", default=None)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    run_dir = Path(args.run_dir) if args.run_dir else cfg.run_root / new_run_id("baseline")
    out = run_baseline(cfg, run_dir)
    print(f"baseline run dir: {out}")


if __name__ == "__main__":
    main()
