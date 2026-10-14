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
