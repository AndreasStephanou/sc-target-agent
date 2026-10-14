from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ["SCTARGET_OFFLINE"] = "1"
warnings.filterwarnings("ignore")


@pytest.fixture(scope="session")
def synthetic_dir(tmp_path_factory) -> Path:
    from scripts.make_synthetic import make_synthetic

    out = tmp_path_factory.mktemp("synth")
    make_synthetic(out, cells_per_sample=150, seed=1)  # small on purpose: the shipped leiden is too coarse, forcing the adaptive path
    return out


@pytest.fixture(scope="session")
def cfg(synthetic_dir, tmp_path_factory):
    from sctarget.config import load_config

    os.environ["SCTARGET_CACHE_DIR"] = str(tmp_path_factory.mktemp("cache"))
    return load_config(
        ROOT / "config" / "synthetic.yaml",
        overrides={"dataset": {"path": str(synthetic_dir / "synthetic.h5ad")}, "output": {"run_dir": str(tmp_path_factory.mktemp("runs"))}},
    )


@pytest.fixture(scope="session")
def ground_truth_path(synthetic_dir, tmp_path_factory) -> Path:
    import yaml

    with open(ROOT / "data" / "ground_truth" / "synthetic.yaml") as fh:
        gt = yaml.safe_load(fh)
    gt["per_cell_truth"] = str(synthetic_dir / "true_state.csv")
    p = tmp_path_factory.mktemp("gt") / "synthetic.yaml"
    with open(p, "w") as fh:
        yaml.safe_dump(gt, fh)
    return p
