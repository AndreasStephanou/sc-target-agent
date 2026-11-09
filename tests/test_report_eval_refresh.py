"""The report is built before `run_eval` scores the run; `run_eval` now rewrites section 6 from eval.json."""
from __future__ import annotations

import json

from sctarget.report.builder import _NO_KEY_TEXT, refresh_eval_section


def test_refresh_replaces_placeholder(tmp_path):
    rd = tmp_path / "run"
    rd.mkdir()
    md = "# Report\n\n## 5. Agent notes\n\nnotes\n\n" + _NO_KEY_TEXT + "\n## 7. Provenance\n\nprov\n"
    (rd / "report.md").write_text(md)
    assert refresh_eval_section(rd) is True  # no eval.json: placeholder stays, html written
    assert "No frozen answer key" in (rd / "report.md").read_text() and (rd / "report.html").exists()
    ev = {
        "ground_truth": "data/ground_truth/x.yaml",
        "prompt_version": 12,
        "scoring_version": 3,
        "passed": True,
        "state_recovery": {"S": {"found": True, "found_strict": False, "best_match": "12", "required_overlap": 5}},
        "target_rank": {"ERBB4": 4},
        "acceptable_rank": {"PEG3": 11, "BCL2": None},
        "decoys_above_target": [],
        "novel_above_target": ["GABRA4"],
        "decoys_nominated_unflagged": [],
        "nominated_precision": 0.2,
        "drug_recall": 0.75,
    }
    (rd / "eval.json").write_text(json.dumps(ev))
    assert refresh_eval_section(rd) is True
    out = (rd / "report.md").read_text()
    assert "No frozen answer key" not in out
    assert "## 6. Evaluation against the frozen answer key" in out and "ERBB4: #4" in out and "PASS" in out
    assert out.count("## 7. Provenance") == 1 and "## 5. Agent notes" in out  # neighbours untouched
    assert "ERBB4: #4" in (rd / "report.html").read_text()


def test_refresh_without_report(tmp_path):
    assert refresh_eval_section(tmp_path) is False
