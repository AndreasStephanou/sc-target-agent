import pytest

from sctarget.tools import REGISTRY, validate_args


def test_all_tools_have_valid_schemas():
    assert len(REGISTRY) >= 14
    for name in REGISTRY.names():
        spec = REGISTRY.get(name)
        s = spec.anthropic_schema()
        assert s["name"] == name
        assert len(s["description"]) > 40
        assert s["input_schema"]["type"] == "object"
        assert "properties" in s["input_schema"]


def test_role_filtering():
    analyst = {t.name for t in REGISTRY.for_role("analyst")}
    target = {t.name for t in REGISTRY.for_role("target")}
    assert "run_qc" in analyst and "run_qc" not in target
    assert "druggability_lookup" in target and "druggability_lookup" not in analyst
    assert "record_finding" in analyst and "record_finding" in target


def test_validate_args_rejects_unknown_and_missing():
    spec = REGISTRY.get("find_markers")
    with pytest.raises(ValueError):
        validate_args(spec, {"cluster": "1", "bogus": 1})
    with pytest.raises(ValueError):
        validate_args(spec, {"n_top": 5})
    with pytest.raises(ValueError):
        validate_args(spec, {"cluster": "1", "n_top": "five"})
    assert validate_args(spec, {"cluster": "1", "n_top": 5, "reference_clusters": None})["n_top"] == 5


def test_disease_state_without_genes_is_rejected(cfg, tmp_path):
    """v2 hand-off rule: a disease_state finding must carry its signature genes."""
    from sctarget.session import AnalysisSession
    from sctarget.tools.targets import record_finding

    s = AnalysisSession.create(cfg, tmp_path / "rf")
    r = record_finding(s, kind="disease_state", title="t", detail="d", cluster="3")
    assert "error" in r and "rejected" in r["error"] and r["recorded"] is None
    assert not s.facts.get("findings", {}).get("disease_state")
    r = record_finding(s, kind="disease_state", title="t", detail="d", cluster="3", genes=list("ABCDE"))
    assert r["recorded"] == "disease_state" and len(s.facts["findings"]["disease_state"]) == 1
    # other kinds are unaffected
    assert record_finding(s, kind="summary", title="t", detail="d")["recorded"] == "summary"


def test_nominator_prompt_carries_markers():
    from sctarget.agents.prompts import stage_prompt

    facts = {
        "findings": {"disease_state": [{"cluster": "3", "genes": ["G1"]}]},
        "pseudobulk_de": {"3:cluster_vs_reference_clusters": {"method": "m", "n_sig_up": 1, "top_up": [{"gene": "PB1"}]}},
        "markers": {"3_vs_1,2": [{"gene": "MK1", "log2fc": 2.0}]},
    }
    p = stage_prompt("target", facts, {"dataset": {}})
    assert "PB1" in p and "MK1" in p and "union" in p
