"""Stability across repeated runs: how similar are the top-k rankings and the nominated sets?"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any


def _top_k(run_dir: Path, k: int) -> list[str]:
    with open(run_dir / "state.json") as fh:
        facts = json.load(fh)["facts"]
    ranked = sorted(((r["gene"], float(r["score"])) for rows in facts.get("ranked_targets", {}).values() for r in rows), key=lambda x: -x[1])
    out: list[str] = []
    for g, _ in ranked:
        if g not in out:
            out.append(g)
    return out[:k]


def _nominated(run_dir: Path) -> set[str]:
    with open(run_dir / "state.json") as fh:
        facts = json.load(fh)["facts"]
    return {c["genes"][0] for c in facts.get("findings", {}).get("candidate_target", []) if c.get("genes")}


def jaccard(a: set[str], b: set[str]) -> float:
    return len(a & b) / len(a | b) if (a | b) else 1.0


def stability(run_dirs: list[str | Path], k: int = 10) -> dict[str, Any]:
    dirs = [Path(d) for d in run_dirs]
    tops = {d.name: _top_k(d, k) for d in dirs}
    noms = {d.name: _nominated(d) for d in dirs}
    pair_j_top, pair_j_nom = [], []
    for a, b in itertools.combinations(dirs, 2):
        pair_j_top.append(jaccard(set(tops[a.name]), set(tops[b.name])))
        pair_j_nom.append(jaccard(noms[a.name], noms[b.name]))
    # which genes appear at rank 1 across runs
    rank1 = [tops[d.name][0] if tops[d.name] else None for d in dirs]
    return {
        "n_runs": len(dirs),
        "k": k,
        "mean_jaccard_topk": round(sum(pair_j_top) / len(pair_j_top), 3) if pair_j_top else None,
        "mean_jaccard_nominated": round(sum(pair_j_nom) / len(pair_j_nom), 3) if pair_j_nom else None,
        "rank1_genes": rank1,
        "rank1_consistent": len(set(rank1)) == 1,
        "top_k_per_run": tops,
        "nominated_per_run": {k_: sorted(v) for k_, v in noms.items()},
    }
