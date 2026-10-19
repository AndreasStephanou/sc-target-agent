"""AnalysisSession: the local, stateful container that tools operate on.

The AnnData object and every full-resolution result live here, on disk and in local
memory. Tools return small JSON summaries; the privacy layer guarantees nothing
per-cell leaves. The agent only ever sees those summaries.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd

from .config import Config
from .privacy import PrivacyPolicy, SamplePseudonymizer
from .runlog import RunLog


@dataclass
class AnalysisSession:
    cfg: Config
    run_dir: Path
    log: RunLog
    policy: PrivacyPolicy
    adata: ad.AnnData | None = None
    pseudo: SamplePseudonymizer | None = None
    cluster_key: str | None = None
    tables: dict[str, Path] = field(default_factory=dict)
    figures: list[dict[str, str]] = field(default_factory=list)
    # structured facts accumulated across tools (used by the report + evaluation)
    facts: dict[str, Any] = field(default_factory=dict)
    _table_counter: int = 0

    # ------------------------------------------------------------------ factory
    @classmethod
    def create(cls, cfg: Config, run_dir: Path) -> "AnalysisSession":
        log = RunLog(run_dir)
        policy = PrivacyPolicy.from_config(cfg.privacy)
        return cls(cfg=cfg, run_dir=Path(run_dir), log=log, policy=policy)

    # ------------------------------------------------------------------ helpers
    @property
    def sample_key(self) -> str:
        return self.cfg.dataset["sample_key"]

    @property
    def condition_key(self) -> str:
        return self.cfg.dataset["condition_key"]

    @property
    def celltype_key(self) -> str | None:
        """The configured cell-type annotation column, if any (used for naming / override)."""
        return self.cfg.dataset.get("celltype_key")

    DERIVED_LINEAGE_KEY = "sct_lineage"

    @property
    def lineage_key(self) -> str | None:
        """The column that defines lineages (sibling sets, compartments, ambient controls):
        the configured annotation when present in obs, else the data-derived coarse clustering
        written by derive_lineages (None until it has been run)."""
        a = self.adata
        ck = self.celltype_key
        if ck and a is not None and ck in a.obs and not self.cfg.analysis.get("lineages", {}).get("ignore_annotation", False):
            return ck
        if a is not None and self.DERIVED_LINEAGE_KEY in a.obs:
            return self.DERIVED_LINEAGE_KEY
        return None

    @property
    def disease(self) -> str:
        return str(self.cfg.dataset["disease"])

    @property
    def reference(self) -> list[str]:
        ref = self.cfg.dataset["reference"]
        return [str(r) for r in (ref if isinstance(ref, list) else [ref])]

    def require_adata(self) -> ad.AnnData:
        if self.adata is None:
            raise RuntimeError("no dataset loaded; call load_dataset first")
        return self.adata

    def counts(self) -> Any:
        """Raw count matrix (cells x genes), from the configured layer or .X."""
        a = self.require_adata()
        layer = self.cfg.dataset.get("counts_layer")
        if layer and layer in a.layers:
            return a.layers[layer]
        return a.X

    def save_table(self, name: str, df: pd.DataFrame) -> Path:
        self._table_counter += 1
        safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)
        path = self.run_dir / "tables" / f"{self._table_counter:03d}_{safe}.csv"
        df.to_csv(path)
        self.tables[name] = path
        return path

    def save_figure(self, name: str, fig: Any, caption: str) -> Path:
        safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in name)
        path = self.run_dir / "figures" / f"{len(self.figures) + 1:02d}_{safe}.png"
        fig.savefig(path, dpi=110, bbox_inches="tight")
        import matplotlib.pyplot as plt

        plt.close(fig)
        self.figures.append({"path": str(path.relative_to(self.run_dir)), "caption": caption})
        return path

    def pseudonymize_series(self, s: pd.Series) -> pd.Series:
        if self.pseudo is None:
            return s.astype(str)
        return s.astype(str).map(self.pseudo.to_pseudo)

    def dump_state(self) -> None:
        state = {
            "cluster_key": self.cluster_key,
            "tables": {k: str(v.relative_to(self.run_dir)) for k, v in self.tables.items()},
            "figures": self.figures,
            "facts": {k: v for k, v in self.facts.items() if k != "knowledge_client"},
            "sample_pseudonyms": (self.pseudo.forward if self.pseudo else {}),
        }
        with open(self.run_dir / "state.json", "w") as fh:
            json.dump(state, fh, indent=2, default=_json_default)


def _json_default(o: Any) -> Any:
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.ndarray,)):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)
