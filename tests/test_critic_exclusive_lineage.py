"""Protocol 23p: a condition-exclusive lineage makes the within-lineage abundance test void, not failed.

Motivating case (CAA de novo floor, round 0.4c): Harmony + derived lineages put the 59k CAA
keratinocytes into their own lineage L4 (99% disease cells); the scripted critic flagged clusters 16
and 13 "enrichment does not hold within its parent cell type" because within L4 the per-sample
proportions of one CAA sub-state against the others show no shift.
"""
from __future__ import annotations

from sctarget.agents.model import ScriptedModel, ToolUse


class _Session:
    disease = "CAA"
    facts: dict = {}

    class cfg:  # minimal stand-in; the critic review reads nothing from the config
        targets: dict = {}
        dataset: dict = {}
        analysis: dict = {}

        @staticmethod
        def get(key, default=None):
            return default


def _facts(exclusive: bool):
    # lineage L4 = clusters 16, 17, 13; lineage L0 = clusters 0, 3
    comp = [
        {"cluster": "16", "n_cells": 31000, "frac_CAA": 0.995 if exclusive else 0.40, "dominant_celltype": "L4"},
        {"cluster": "17", "n_cells": 25600, "frac_CAA": 0.993 if exclusive else 0.45, "dominant_celltype": "L4"},
        {"cluster": "13", "n_cells": 3000, "frac_CAA": 0.98 if exclusive else 0.50, "dominant_celltype": "L4"},
        {"cluster": "0", "n_cells": 23000, "frac_CAA": 0.14, "dominant_celltype": "L0"},
        {"cluster": "3", "n_cells": 3000, "frac_CAA": 0.85, "dominant_celltype": "L0"},
    ]
    return {
        "findings": {
            "disease_state": [
                {"cluster": "16", "basis": "composition", "genes": ["A"]},
                {"cluster": "13", "basis": "composition", "genes": ["B"]},
                {"cluster": "3", "basis": "composition", "genes": ["C"]},
            ],
            "candidate_target": [],
        },
        "cluster_composition": comp,
        "differential_abundance": {"table": [{"cluster": c, "log2_ratio": 4.0, "n_disease_samples_present": 6, "n_cells": 5000, "disease_enriched": True} for c in ("16", "13", "3")]},
        "differential_abundance_within_celltype": {"table": [
            {"cluster": "16", "log2_ratio": 0.07, "disease_enriched": False},
            {"cluster": "13", "log2_ratio": -0.8, "disease_enriched": False},
            {"cluster": "3", "log2_ratio": 5.4, "disease_enriched": True},
        ]},
        "reference_library_check": {"flagged": []},
    }


def _flags(facts):
    m = ScriptedModel(_Session())
    m._facts = lambda: facts
    return [s.input for s in m._critic_review() if isinstance(s, ToolUse) and s.name == "record_finding" and s.input.get("kind") == "critic_flag"]


def test_exclusive_lineage_is_reported_as_void_not_failed():
    flags = _flags(_facts(exclusive=True))
    titles = [f["title"] for f in flags]
    assert not any("does not hold" in t for t in titles), titles
    void = [f for f in flags if "condition-exclusive" in f["title"]]
    assert len(void) == 1 and void[0]["severity"] == "minor" and "L4" in void[0]["title"] and "99%" in void[0]["title"]
    assert "pan-state genes cancel" in void[0]["detail"]


def test_mixed_lineage_keeps_the_failed_enrichment_flag():
    flags = _flags(_facts(exclusive=False))
    titles = [f["title"] for f in flags]
    assert sum("does not hold" in t for t in titles) == 2  # clusters 16 and 13
    assert not any("condition-exclusive" in t for t in titles)


def test_lineage_share_computation():
    m = ScriptedModel(_Session())
    share = m._lineage_disease_share(_facts(exclusive=True))
    lin, s = share["16"]
    assert lin == "L4" and 0.99 <= s <= 1.0
    lin0, s0 = share["0"]
    assert lin0 == "L0" and 0.2 <= s0 <= 0.3
