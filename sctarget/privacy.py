"""Privacy guard: controls what may leave the machine and reach the model API.

Policy (configurable in ``privacy:`` of the config):
  * Only aggregated results leave: gene lists, fold changes, p-values, cluster
    proportions, enrichment terms. Never per-cell values, never raw count matrices.
  * Sample identifiers are pseudonymized (S01, S02, ...) before they reach the model;
    the mapping stays in the run directory on disk.
  * Columns in ``forbidden_obs_columns`` are never summarized, listed, or sent.
  * Lists and tables are truncated to a maximum size; the full tables stay on disk
    and the tool result carries their path instead.

Every tool result passes through :func:`sanitize` before being added to the model
context, and the sanitized payload is what gets logged as "sent".
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd


class PrivacyViolation(RuntimeError):
    pass


@dataclass
class PrivacyPolicy:
    max_list_items: int = 100
    max_table_rows: int = 60
    pseudonymize_samples: bool = True
    forbidden_obs_columns: set[str] = field(default_factory=lambda: {"patient_id", "donor_id", "mrn", "name", "dob", "date_of_birth"})
    # a numeric array longer than this is assumed to be per-cell data and is rejected outright
    per_cell_threshold: int = 2000

    @classmethod
    def from_config(cls, cfg: dict) -> "PrivacyPolicy":
        return cls(
            max_list_items=int(cfg.get("max_list_items", 100)),
            max_table_rows=int(cfg.get("max_table_rows", 60)),
            pseudonymize_samples=bool(cfg.get("pseudonymize_samples", True)),
            forbidden_obs_columns={c.lower() for c in cfg.get("forbidden_obs_columns", [])},
        )


class SamplePseudonymizer:
    """Stable mapping real sample id -> S01, S02, ... (order of first appearance, sorted)."""

    def __init__(self, sample_ids: list[str], enabled: bool = True) -> None:
        self.enabled = enabled
        ids = sorted({str(s) for s in sample_ids})
        self.forward = {s: f"S{i + 1:02d}" for i, s in enumerate(ids)}
        self.reverse = {v: k for k, v in self.forward.items()}

    def to_pseudo(self, s: Any) -> str:
        s = str(s)
        return self.forward.get(s, s) if self.enabled else s

    def to_real(self, s: Any) -> str:
        s = str(s)
        return self.reverse.get(s, s) if self.enabled else s

    def map_index(self, idx: pd.Index) -> pd.Index:
        return pd.Index([self.to_pseudo(s) for s in idx])


def _is_number(x: Any) -> bool:
    return isinstance(x, (int, float, np.integer, np.floating)) and not isinstance(x, bool)


def _clean_scalar(x: Any) -> Any:
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        xf = float(x)
        if math.isnan(xf) or math.isinf(xf):
            return None
        return round(xf, 6)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (pd.Timestamp,)):
        return x.isoformat()
    return x


def sanitize(obj: Any, policy: PrivacyPolicy, _depth: int = 0) -> Any:
    """Recursively convert a tool result into a JSON-safe, size-bounded payload.

    Raises :class:`PrivacyViolation` if a per-cell-sized numeric array is found.
    """
    if _depth > 12:
        return "<truncated: too deep>"

    if isinstance(obj, pd.DataFrame):
        if any(str(c).lower() in policy.forbidden_obs_columns for c in obj.columns):
            obj = obj[[c for c in obj.columns if str(c).lower() not in policy.forbidden_obs_columns]]
        n = len(obj)
        head = obj.head(policy.max_table_rows)
        rows = [{str(k): _clean_scalar(v) for k, v in r.items()} for r in head.reset_index().to_dict(orient="records")]
        out: dict[str, Any] = {"rows": rows, "n_rows_total": int(n)}
        if n > policy.max_table_rows:
            out["truncated"] = True
        return out

    if isinstance(obj, pd.Series):
        return sanitize(obj.to_frame(), policy, _depth + 1)

    if isinstance(obj, np.ndarray):
        if obj.ndim >= 2 or (obj.dtype.kind in "fiu" and obj.size > policy.per_cell_threshold):
            raise PrivacyViolation(f"refusing to send array of shape {obj.shape} to the model (looks like per-cell data)")
        obj = obj.tolist()

    if isinstance(obj, dict):
        clean: dict[str, Any] = {}
        for k, v in obj.items():
            if str(k).lower() in policy.forbidden_obs_columns:
                continue
            clean[str(k)] = sanitize(v, policy, _depth + 1)
        return clean

    if isinstance(obj, (list, tuple, set)):
        items = list(obj)
        if len(items) > policy.per_cell_threshold and all(_is_number(i) for i in items[:50]):
            raise PrivacyViolation(f"refusing to send numeric list of length {len(items)} to the model")
        n = len(items)
        items = [sanitize(i, policy, _depth + 1) for i in items[: policy.max_list_items]]
        if n > policy.max_list_items:
            items.append(f"<... {n - policy.max_list_items} more items truncated; full table on disk>")
        return items

    return _clean_scalar(obj)


def check_obs_columns(columns: list[str], policy: PrivacyPolicy) -> list[str]:
    """Return the obs columns that are safe to describe to the model."""
    return [c for c in columns if str(c).lower() not in policy.forbidden_obs_columns]
