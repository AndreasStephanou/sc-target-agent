"""Configuration loading and validation."""
from __future__ import annotations

import copy
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG_PATH = REPO_ROOT / "config" / "default.yaml"


def _deep_update(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_update(out[k], v)
        else:
            out[k] = v
    return out


@dataclass
class Config:
    raw: dict[str, Any] = field(default_factory=dict)

    # ---- convenience accessors -------------------------------------------------
    @property
    def dataset(self) -> dict:
        return self.raw["dataset"]

    @property
    def analysis(self) -> dict:
        return self.raw["analysis"]

    @property
    def targets(self) -> dict:
        return self.raw["targets"]

    @property
    def agent(self) -> dict:
        return self.raw["agent"]

    @property
    def privacy(self) -> dict:
        return self.raw["privacy"]

    @property
    def run_root(self) -> Path:
        return Path(self.raw.get("output", {}).get("run_dir", "runs"))

    def get(self, dotted: str, default: Any = None) -> Any:
        node: Any = self.raw
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def resolve_path(self, p: str | Path) -> Path:
        p = Path(p)
        return p if p.is_absolute() else (REPO_ROOT / p)

    def validate(self) -> None:
        ds = self.dataset
        for key in ("path", "sample_key", "condition_key", "disease", "reference"):
            if key not in ds or ds[key] in (None, ""):
                raise ValueError(f"dataset.{key} is required in config")
        if isinstance(ds["reference"], str):
            ds["reference"] = [ds["reference"]]
        w = self.targets.get("weights", {})
        if w and abs(sum(w.values()) - 1.0) > 0.05:
            raise ValueError(f"targets.weights should sum to ~1, got {sum(w.values()):.2f}")


def load_config(path: str | Path | None = None, overrides: dict | None = None) -> Config:
    with open(DEFAULT_CONFIG_PATH) as fh:
        base = yaml.safe_load(fh)
    if path is not None:
        with open(path) as fh:
            user = yaml.safe_load(fh) or {}
        base = _deep_update(base, user)
    if overrides:
        base = _deep_update(base, overrides)
    ds = base.get("dataset", {})
    if "case" in ds:  # alias: any two-condition design (the default config carries a placeholder 'disease')
        ds["disease"] = ds.pop("case")
    if "control" in ds:
        ds["reference"] = ds.pop("control")
    if isinstance(ds.get("reference"), str):
        ds["reference"] = [ds["reference"]]
    # env overrides
    if os.getenv("SCTARGET_MODEL"):
        base["agent"]["model"] = os.environ["SCTARGET_MODEL"]
    cfg = Config(base)
    cfg.validate()
    return cfg


def offline_mode() -> bool:
    return os.getenv("SCTARGET_OFFLINE", "0").strip() in ("1", "true", "True", "yes")


def cache_dir() -> Path:
    d = Path(os.getenv("SCTARGET_CACHE_DIR", REPO_ROOT / "cache"))
    d.mkdir(parents=True, exist_ok=True)
    return d
