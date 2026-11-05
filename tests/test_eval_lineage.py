"""Protocol 23i: the evaluator names the lineage of the state that carries a target, and when the
key names the expected state's lineage, a rank in a state of another lineage is reported but not
credited.

Motivating case (HCC agent run, round 0.4b): GPC3 ranked 2nd in a CAF-like stromal state and was
unscored in the tumour hepatocyte states the key means; the within-state rule reported recall 1.0.
"""
from __future__ import annotations

import json

import yaml

from sctarget.eval import evaluate_run


def _facts():
    return {
        "prompt_version": 8,
        "scoring_version": 3,
        "findings": {
            "disease_state": [
                {"cluster": "14", "basis": "composition", "genes": ["NDUFA4L2", "SQLE", "MAGEA3"]},
                {"cluster": "12", "basis": "composition", "genes": ["THY1", "FN1", "GPC3"]},
            ],
            "candidate_target": [{"genes": ["FN1"]}, {"genes": ["GPC3"]}],
            "critic_flag": [],
        },
        "cluster_composition": [
            {"cluster": "14", "dominant_celltype": "Hepatocyte"},
            {"cluster": "12", "dominant_celltype": "Fibroblast"},
        ],
        "ranked_targets": {
            "14": [{"gene": "NDUFA4L2", "score": 0.9, "rank": 1}, {"gene": "SQLE", "score": 0.8, "rank": 2}],
            "12": [{"gene": "FN1", "score": 0.95, "rank": 1}, {"gene": "GPC3", "score": 0.82, "rank": 2}],
        },
        "druggability": {},
    }


def _write(tmp_path, facts, key):
    run = tmp_path / "run"
    run.mkdir(parents=True)
    (run / "state.json").write_text(json.dumps({"facts": facts}))
    p = tmp_path / "gt.yaml"
    p.write_text(yaml.safe_dump(key))
    return run, p


def _key(lineage=None):
    st = {"name": "HCC_malignant_hepatocyte", "marker_genes": ["GPC3", "AFP", "MDK", "SPINK1", "AKR1B10"], "min_marker_overlap": 5}
    if lineage:
        st["lineage"] = lineage
    return {"expected_states": [st], "expected_targets": ["GPC3"], "acceptable_targets": ["AFP"], "expected_drugs": [], "decoy_targets": ["MKI67"], "per_cell_truth": None}


def test_lineage_of_carrying_state_is_reported_without_a_key_lineage(tmp_path):
    run, p = _write(tmp_path, _facts(), _key())
    res = evaluate_run(run, p)
    assert res.target_rank == {"GPC3": 2} and res.target_state == {"GPC3": "12"}
    assert res.target_state_lineage == {"GPC3": "Fibroblast"}
    assert res.target_rank_off_lineage == {}
    assert res.target_recall_at_k == 1.0  # frozen keys without a lineage keep their numbers


def test_off_lineage_rank_is_reported_but_not_credited_when_the_key_names_the_lineage(tmp_path):
    run, p = _write(tmp_path, _facts(), _key(lineage="Hepatocyte"))
    res = evaluate_run(run, p)
    assert res.target_rank == {"GPC3": None} and res.target_state == {"GPC3": None}
    assert res.target_recall_at_k == 0.0
    assert res.target_rank_off_lineage == {"GPC3": {"rank": 2, "state": "12", "lineage": "Fibroblast"}}
    # a target in the right lineage is credited as before
    facts = _facts()
    facts["ranked_targets"]["14"].insert(0, {"gene": "GPC3", "score": 0.99, "rank": 1})
    run2, p2 = _write(tmp_path / "b", facts, _key(lineage="hepatocyte"))
    res2 = evaluate_run(run2, p2)
    assert res2.target_rank == {"GPC3": 1} and res2.target_state_lineage == {"GPC3": "Hepatocyte"}
