"""Score a finished run against a ground-truth file.

Metrics
  state_recovery      For each expected state: was a recorded disease_state found whose
                      top genes overlap the expected markers by >= min_marker_overlap?
                      Lenient: recorded genes + the cluster's DE/marker tables (the tools reached
                      the right place). Strict (found_strict): only genes the agent wrote into the
                      finding count (the agent noticed). The gap between them is a result.
                      (Optionally, with per_cell_truth: best-matching cluster's F1.)
  target_rank         Rank of each expected target in the ranking (None if absent).
  target_recall_at_k  Fraction of expected targets within the top k (k = targets.top_k_report
                      or 5 by default).
  drug_recall         Fraction of expected drugs that appear among the inhibitors of any
                      nominated candidate.
  false_leads         Genes ranked ABOVE the first expected target that are not expected/
                      acceptable, split into decoys_above_target (known-wrong, the error metric)
                      and novel_above_target (absent from the answer key: new hypotheses, not errors).
  nominated_precision Fraction of recorded candidate_target genes that are expected or
                      acceptable.
  critic_flagged_decoys  Decoys that were nominated but flagged by the critic (good
                      behaviour) vs nominated and unflagged (bad).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from ..runlog import usage_totals_from_events


@dataclass
class EvalResult:
    run_dir: str
    ground_truth: str
    prompt_version: int = 1
    scoring_version: int = 1
    state_recovery: dict[str, Any] = field(default_factory=dict)
    target_rank: dict[str, int | None] = field(default_factory=dict)         # best rank WITHIN a recorded state's own ranking (the key's "ranked first in its state")
    target_rank_global: dict[str, int | None] = field(default_factory=dict)  # rank in the merged ranking of all states (penalises lineage coverage; reported, not scored)
    target_state: dict[str, str | None] = field(default_factory=dict)        # the state whose ranking gives target_rank
    target_state_lineage: dict[str, str | None] = field(default_factory=dict)  # dominant lineage of that state (protocol 23i: the reader sees where the credit came from)
    target_rank_off_lineage: dict[str, Any] = field(default_factory=dict)    # when the key names the expected state's lineage: best rank in a state of ANOTHER lineage (reported, not credited)
    acceptable_rank: dict[str, int | None] = field(default_factory=dict)  # acceptable (non-expected) targets, within-state best: reported, never scored
    target_recall_at_k: float | None = None
    k: int = 5
    drug_recall: float | None = None
    drugs_found: list[str] = field(default_factory=list)
    false_leads: list[str] = field(default_factory=list)          # above the first true target and NOT expected/acceptable (= novel + decoy)
    decoys_above_target: list[str] = field(default_factory=list)  # ... of which known-wrong (decoy_targets)
    novel_above_target: list[str] = field(default_factory=list)   # ... of which not in the answer key at all: new hypotheses, not errors
    nominated: list[str] = field(default_factory=list)
    nominated_precision: float | None = None
    decoys_nominated_unflagged: list[str] = field(default_factory=list)
    decoys_nominated_flagged: list[str] = field(default_factory=list)
    n_tool_calls: int = 0
    n_tool_errors: int = 0
    tokens: dict[str, int] = field(default_factory=dict)
    passed: bool = False

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def load_ground_truth(path: str | Path) -> dict[str, Any]:
    with open(path) as fh:
        gt = yaml.safe_load(fh)
    gt.setdefault("acceptable_targets", [])
    gt.setdefault("decoy_targets", [])
    gt.setdefault("expected_drugs", [])
    gt.setdefault("expected_states", [])
    return gt


def _load_state(run_dir: Path) -> dict[str, Any]:
    with open(run_dir / "state.json") as fh:
        return json.load(fh)


def _state_recovery(gt: dict[str, Any], facts: dict[str, Any], run_dir: Path) -> dict[str, Any]:
    out: dict[str, Any] = {}
    recorded = facts.get("findings", {}).get("disease_state", [])
    # union of top genes we know about per recorded state (recorded genes + DE top_up + markers)
    def genes_for(cl: str, rec: dict[str, Any]) -> set[str]:
        g = set(rec.get("genes") or [])
        pb = facts.get("pseudobulk_de", {}).get(f"{cl}:cluster_vs_reference_clusters", {})
        g |= {r["gene"] for r in pb.get("top_up", [])[:30]}
        for k, rows in facts.get("markers", {}).items():
            if k.startswith(f"{cl}_vs_"):
                g |= {r["gene"] for r in rows[:30]}
        return g

    truth_df = None
    if gt.get("per_cell_truth"):
        p = Path(gt["per_cell_truth"])
        if not p.is_absolute():
            from ..config import REPO_ROOT

            p = REPO_ROOT / p
        if p.exists():
            truth_df = pd.read_csv(p, index_col=0)

    for exp in gt["expected_states"]:
        want = set(exp.get("marker_genes", []))
        need = int(exp.get("min_marker_overlap", max(1, len(want) // 2)))
        best = None
        for rec in recorded:
            cl = str(rec.get("cluster"))
            ov = want & genes_for(cl, rec)
            if best is None or len(ov) > best["overlap"]:
                best = {"cluster": cl, "overlap": len(ov), "overlap_genes": sorted(ov)}
        found = bool(best and best["overlap"] >= need)
        # strict variant: only genes the agent itself wrote into the disease_state finding count.
        # lenient says "the tools reached the right place"; strict says "the agent noticed".
        strict_best = None
        for rec in recorded:
            ov = want & set(rec.get("genes") or [])
            if strict_best is None or len(ov) > strict_best["overlap"]:
                strict_best = {"cluster": str(rec.get("cluster")), "overlap": len(ov), "overlap_genes": sorted(ov), "n_recorded_genes": len(rec.get("genes") or [])}
        found_strict = bool(strict_best and strict_best["overlap"] >= need)
        entry: dict[str, Any] = {"found": found, "found_strict": found_strict, "required_overlap": need, "best_match": best, "best_match_strict": strict_best}
        matched_basis = next((r.get("basis") for r in recorded if best and str(r.get("cluster")) == best["cluster"]), "composition")
        if truth_df is not None and best and matched_basis != "program" and (run_dir / "cluster_assignments.csv").exists():
            # optional: compute F1 of the matched cluster vs the true state label
            assign = pd.read_csv(run_dir / "cluster_assignments.csv", index_col=0).iloc[:, 0].astype(str)
            truth = truth_df.iloc[:, 0].reindex(assign.index)
            pred = assign == best["cluster"]
            actual = truth == exp["name"]
            tp = int((pred & actual).sum())
            prec = tp / max(1, int(pred.sum()))
            rec_ = tp / max(1, int(actual.sum()))
            entry["cell_level"] = {"precision": round(prec, 3), "recall": round(rec_, 3), "f1": round(2 * prec * rec_ / max(1e-9, prec + rec_), 3)}
        if best:
            entry["basis"] = matched_basis
        out[exp["name"]] = entry
    return out


def evaluate_run(run_dir: str | Path, ground_truth: str | Path, k: int | None = None) -> EvalResult:
    run_dir = Path(run_dir)
    gt = load_ground_truth(ground_truth)
    state = _load_state(run_dir)
    facts = state["facts"]
    res = EvalResult(run_dir=str(run_dir), ground_truth=str(ground_truth))

    # ranking: merge all clusters' rankings, keep best rank per gene
    ranked: list[tuple[str, float]] = []
    for rows in facts.get("ranked_targets", {}).values():
        for r in rows:
            ranked.append((r["gene"], float(r["score"])))
    ranked.sort(key=lambda x: -x[1])
    order: list[str] = []
    for g, _ in ranked:
        if g not in order:
            order.append(g)
    res.k = k or 5
    expected = list(gt.get("expected_targets", []))
    acceptable = set(gt["acceptable_targets"]) | set(expected)
    decoys = set(gt["decoy_targets"])

    # within-state ranks (protocol 23g): a target's rank is judged in the ranking of the state that
    # carries it, as every key is worded ("ranked first in its state"). The merged ranking is kept as
    # target_rank_global; it falls as more lineages are characterised, which is not a defect.
    per_state: dict[str, list[str]] = {str(cl): [r["gene"] for r in sorted(rows, key=lambda r: int(r.get("rank") or 10**6))] for cl, rows in facts.get("ranked_targets", {}).items()}

    # protocol 23i: the lineage of the state that carries a target is reported, and when the key
    # names the expected state's lineage (optional `lineage:` on an expected state) a rank in a state
    # of another lineage is not credited (HCC: GPC3 2nd in a CAF-like stromal state, unscored in the
    # tumour hepatocyte states the key means)
    comp = {str(r.get("cluster")): r for r in facts.get("cluster_composition", []) if isinstance(r, dict)}
    lineage_of = {cl: (comp.get(cl, {}).get("dominant_celltype") or None) for cl in per_state}
    expected_lineages = {str(e["lineage"]).lower() for e in gt.get("expected_states", []) if e.get("lineage")}

    def _on_lineage(cl: str) -> bool:
        if not expected_lineages:
            return True
        lin = (lineage_of.get(cl) or "").lower()
        return any(x == lin or x in lin or lin in x for x in expected_lineages) if lin else False

    def _best(g: str, on_lineage: bool | None = True) -> tuple[int | None, str | None]:
        hits = [(genes.index(g) + 1, cl) for cl, genes in per_state.items() if g in genes and (on_lineage is None or _on_lineage(cl) == on_lineage)]
        return min(hits) if hits else (None, None)

    res.target_rank_global = {g: (order.index(g) + 1 if g in order else None) for g in expected}
    res.target_rank, res.target_state, res.target_state_lineage, res.target_rank_off_lineage = {}, {}, {}, {}
    for g in expected:
        rk, cl = _best(g)
        res.target_rank[g], res.target_state[g] = rk, cl
        res.target_state_lineage[g] = lineage_of.get(cl) if cl else None
        if expected_lineages:
            ork, ocl = _best(g, on_lineage=False)
            if ork is not None:
                res.target_rank_off_lineage[g] = {"rank": ork, "state": ocl, "lineage": lineage_of.get(ocl)}
    acc_only = [g for g in gt.get("acceptable_targets", []) if g not in expected]
    res.acceptable_rank = {g: _best(g)[0] for g in acc_only}
    if expected:
        res.target_recall_at_k = sum(1 for g in expected if res.target_rank[g] is not None and res.target_rank[g] <= res.k) / len(expected)
    # genes above the first expected target, in the state where it ranks best; with no hit, every ranked gene
    hits = [(r, g) for g, r in res.target_rank.items() if r is not None]
    if hits:
        first_hit, first_gene = min(hits)
        above = per_state[res.target_state[first_gene]][: first_hit - 1]
    else:
        above = order
    res.false_leads = [g for g in above if g not in acceptable]
    res.decoys_above_target = [g for g in above if g in decoys]
    res.novel_above_target = [g for g in res.false_leads if g not in decoys]

    # nominated candidates
    cands = facts.get("findings", {}).get("candidate_target", [])
    res.nominated = list(dict.fromkeys(c["genes"][0] for c in cands if c.get("genes")))  # unique, order kept
    if res.nominated:
        res.nominated_precision = sum(1 for g in res.nominated if g in acceptable) / len(res.nominated)
    flagged_genes = {g for fl in facts.get("findings", {}).get("critic_flag", []) for g in (fl.get("genes") or [])}
    res.decoys_nominated_flagged = [g for g in res.nominated if g in decoys and g in flagged_genes]
    res.decoys_nominated_unflagged = [g for g in res.nominated if g in decoys and g not in flagged_genes]

    # drugs
    inhibitors: set[str] = set()
    for rows in facts.get("ranked_targets", {}).values():
        for r in rows:
            if r["gene"] in res.nominated:
                inhibitors |= {str(d).lower() for d in r.get("top_inhibitors", [])}
    for g, d in facts.get("druggability", {}).items():
        if g in res.nominated:
            inhibitors |= {str(x).lower() for x in d.get("top_inhibitors", [])}
    want_drugs = [d.lower() for d in gt["expected_drugs"]]
    # match on the compound stem so 'neratinib' matches 'neratinib maleate'
    res.drugs_found = sorted(d for d in want_drugs if any(inh == d or inh.startswith(d + " ") or d in inh.split() for inh in inhibitors))
    if want_drugs:
        res.drug_recall = len(res.drugs_found) / len(want_drugs)

    res.prompt_version = int(facts.get("prompt_version", 1))
    res.scoring_version = int(facts.get("scoring_version", 1))
    res.state_recovery = _state_recovery(gt, facts, run_dir)

    # provenance stats
    events = []
    if (run_dir / "run.jsonl").exists():
        with open(run_dir / "run.jsonl") as fh:
            events = [json.loads(l) for l in fh if l.strip()]
    tcs = [e for e in events if e["type"] == "tool_call"]
    res.n_tool_calls = len(tcs)
    res.n_tool_errors = sum(1 for e in tcs if e.get("error"))
    for kk, v in usage_totals_from_events(events).items():  # numeric fields only; usage also carries the model id (23r)
        res.tokens[kk] = res.tokens.get(kk, 0) + v

    res.passed = bool(
        all(v["found"] for v in res.state_recovery.values())
        and (res.target_recall_at_k is None or res.target_recall_at_k == 1.0)
        and not res.decoys_nominated_unflagged
    )
    return res


def format_eval(res: EvalResult) -> str:
    lines = [f"Evaluation of {res.run_dir} vs {Path(res.ground_truth).name}  (prompts v{res.prompt_version}, scoring v{res.scoring_version})", "=" * 72]
    for name, v in res.state_recovery.items():
        bm = v.get("best_match") or {}
        lines.append(f"state '{name}': {'FOUND' if v['found'] else 'missed'}  (cluster {bm.get('cluster')}, marker overlap {bm.get('overlap')}/{v['required_overlap']} needed: {', '.join(bm.get('overlap_genes', []))})")
        sb = v.get("best_match_strict") or {}
        lines.append(f"    strict (agent-recorded genes only): {'FOUND' if v.get('found_strict') else 'missed'}  (cluster {sb.get('cluster')}, overlap {sb.get('overlap')}/{v['required_overlap']}, {sb.get('n_recorded_genes', 0)} genes recorded)")
        if "cell_level" in v:
            lines.append(f"    cell-level P/R/F1: {v['cell_level']['precision']}/{v['cell_level']['recall']}/{v['cell_level']['f1']}")
    state_str = {g: (f"{cl} ({res.target_state_lineage.get(g)})" if cl and res.target_state_lineage.get(g) else cl) for g, cl in res.target_state.items()}
    lines.append(f"target ranks (within its state): {res.target_rank}  state: {state_str}  global: {res.target_rank_global}   recall@{res.k}: {res.target_recall_at_k}")
    if res.target_rank_off_lineage:
        lines.append("    off-lineage (reported, not credited; the key names the expected state's lineage): " + "; ".join(f"{g}: rank {v['rank']} in state {v['state']} ({v['lineage']})" for g, v in res.target_rank_off_lineage.items()))
    if res.acceptable_rank:
        lines.append(f"acceptable-target ranks (reported, not scored): {res.acceptable_rank}")
    lines.append(f"above first true target: {len(res.false_leads)} not in answer key -> known-wrong decoys: {res.decoys_above_target or 'none'}; novel candidates: {res.novel_above_target or 'none'}")
    lines.append(f"nominated: {res.nominated}  precision: {res.nominated_precision}")
    lines.append(f"decoys nominated - flagged by critic: {res.decoys_nominated_flagged or 'none'}; UNFLAGGED: {res.decoys_nominated_unflagged or 'none'}")
    lines.append(f"drug recall: {res.drug_recall}  found: {res.drugs_found}")
    lines.append(f"tool calls: {res.n_tool_calls} ({res.n_tool_errors} errors)   tokens: {res.tokens or 'n/a'}")
    lines.append(f"PASS: {res.passed}")
    return "\n".join(lines)
