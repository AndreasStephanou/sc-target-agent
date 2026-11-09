"""Protocol 23m (prompts v12): a stage briefing may shorten what it says about an item, never drop one.

Motivating case (CAA agent run, round 0.4c): the critic wrote "the candidate list I received was
truncated after ERBB4". Its briefing was one JSON dump cut at 9,000 characters; the claims ran
31,715. Re-measured on every agent run of rounds 0.4 and 0.4b: the cap was exceeded in all of them.
"""
from __future__ import annotations

from sctarget.agents.prompts import PROMPT_VERSION, stage_prompt


def _facts(n_states: int = 8, n_cands: int = 20, detail_chars: int = 2500) -> dict:
    long = "x" * detail_chars
    states = [{"cluster": str(i), "title": f"state {i}", "detail": long, "genes": [f"SG{i}_{j}" for j in range(15)], "basis": "composition"} for i in range(n_states)]
    cands = [{"cluster": str(i % n_states), "title": f"CAND{i} (cluster {i % n_states}): a long title", "detail": long, "genes": [f"CAND{i}"], "confidence": "medium"} for i in range(n_cands)]
    pb = {f"{i}:cluster_vs_reference_clusters": {"method": "m", "n_sig_up": 500, "top_up": [{"gene": f"PB{i}_{j}"} for j in range(60)]} for i in range(n_states)}
    pb.update({f"{i}:disease_vs_reference_within_lineage": {"n_sig_up": 400, "top_up": [{"gene": f"LN{i}_{j}"} for j in range(60)], "lineage_coverage": 0.9} for i in range(n_states)})
    mk = {f"{i}_vs_rest": [{"gene": f"MK{i}_{j}"} for j in range(40)] for i in range(n_states)}
    ranked = {str(i): [{"gene": f"CAND{i}", "score": 0.9, "rank": 1, "attribution": "state_specific", "flags": [], "extra": long} for _ in range(8)] for i in range(n_states)}
    return {
        "findings": {"disease_state": states, "candidate_target": cands},
        "pseudobulk_de": pb,
        "markers": mk,
        "ranked_targets": ranked,
        "differential_abundance": {"table": [{"cluster": str(i), "log2_ratio": 1.0} for i in range(n_states)]},
    }


def test_version_bumped():
    assert PROMPT_VERSION >= 12


def test_critic_sees_every_nomination_and_state():
    facts = _facts()
    p = stage_prompt("critic", facts, {"dataset": {}})
    assert "<truncated>" not in p.split("Differential abundance table")[0]
    for i in range(20):
        assert f"CAND{i}" in p, f"nomination CAND{i} missing from the critic briefing"
    for i in range(8):
        assert f"SG{i}_0" in p, f"state {i} missing from the critic briefing"
    # the index line comes first and names every nomination
    first = p.split("\n")[1] if p.startswith("Dataset config") else p.split("\n")[0]
    idx = [line for line in p.split("\n") if line.startswith("Nominated candidates, all of them")][0]
    assert all(f"CAND{i} [cluster" in idx for i in range(20))
    assert "every one" in p.lower() or "every nominated gene" in p.lower()
    del first


def test_nominator_sees_every_state_and_every_contrast():
    facts = _facts()
    p = stage_prompt("target", facts, {"dataset": {}})
    assert "<truncated>" not in p
    for i in range(8):
        assert f"SG{i}_0" in p and f"PB{i}_0" in p and f"LN{i}_0" in p and f"MK{i}_0" in p, f"state {i} or one of its gene lists missing"
    idx = [line for line in p.split("\n") if line.startswith("Recorded disease states, all of them")][0]
    assert idx.count("[cluster") == 8


def test_clipping_shortens_detail_but_keeps_items():
    from sctarget.agents.prompts import _fit_items

    items = [{"genes": [f"G{i}"], "detail": "y" * 5000} for i in range(30)]
    s = _fit_items(items, limit=12000)
    assert len(s) <= 12000 or "...<clipped>" in s
    assert all(f'"G{i}"' in s for i in range(30))
    assert "...<clipped>" in s


def test_small_briefing_unchanged_in_substance():
    facts = _facts(n_states=1, n_cands=2, detail_chars=50)
    p = stage_prompt("critic", facts, {"dataset": {}})
    assert "...<clipped>" not in p and "CAND0" in p and "CAND1" in p
