import numpy as np
import pandas as pd
import pytest

from sctarget.privacy import PrivacyPolicy, PrivacyViolation, SamplePseudonymizer, sanitize


def test_truncates_long_lists_and_tables():
    pol = PrivacyPolicy(max_list_items=5, max_table_rows=3)
    out = sanitize({"genes": [f"G{i}" for i in range(20)]}, pol)
    assert len(out["genes"]) == 6 and "truncated" in out["genes"][-1]
    df = pd.DataFrame({"a": range(10)})
    out = sanitize(df, pol)
    assert out["n_rows_total"] == 10 and len(out["rows"]) == 3 and out["truncated"]


def test_rejects_per_cell_arrays():
    pol = PrivacyPolicy(per_cell_threshold=100)
    with pytest.raises(PrivacyViolation):
        sanitize(np.random.rand(500), pol)
    with pytest.raises(PrivacyViolation):
        sanitize(np.random.rand(10, 10), pol)
    with pytest.raises(PrivacyViolation):
        sanitize(list(range(1000)), pol)


def test_drops_forbidden_columns():
    pol = PrivacyPolicy(forbidden_obs_columns={"patient_id"})
    out = sanitize({"patient_id": "X123", "n_cells": 10, "table": pd.DataFrame({"patient_id": ["a"], "n": [1]})}, pol)
    assert "patient_id" not in out
    assert "patient_id" not in out["table"]["rows"][0]


def test_nan_and_numpy_scalars_become_json_safe():
    pol = PrivacyPolicy()
    out = sanitize({"x": np.float64(1.23456789), "y": float("nan"), "z": np.int64(3), "b": np.bool_(True)}, pol)
    assert out == {"x": 1.234568, "y": None, "z": 3, "b": True}


def test_pseudonymizer_roundtrip():
    p = SamplePseudonymizer(["dogB", "dogA", "dogC"])
    assert p.to_pseudo("dogA") == "S01" and p.to_pseudo("dogC") == "S03"
    assert p.to_real("S02") == "dogB"
    assert p.to_pseudo("unknown") == "unknown"
