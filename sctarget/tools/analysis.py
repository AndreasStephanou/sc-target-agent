"""Deterministic single-cell analysis tools.

Each function is registered in the tool registry with a JSON schema. They operate on
the :class:`~sctarget.session.AnalysisSession`, write full tables/figures to the run
directory, and return compact JSON summaries. Nothing per-cell is returned.
"""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

from ..session import AnalysisSession
from .registry import tool

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

ANALYST = {"planner", "analyst", "critic"}
ALL = {"planner", "analyst", "target", "critic"}


# ----------------------------------------------------------------------------- utils
def _dense_col_means(X: Any, idx: np.ndarray) -> np.ndarray:
    sub = X[idx]
    m = np.asarray(sub.mean(axis=0)).ravel()
    return m


def _pct_expressing(X: Any, idx: np.ndarray) -> np.ndarray:
    sub = X[idx]
    if sp.issparse(sub):
        return np.asarray((sub > 0).mean(axis=0)).ravel()
    return (sub > 0).mean(axis=0)


def _looks_like_counts(X: Any) -> bool:
    sample = X[:200]
    if sp.issparse(sample):
        sample = sample.toarray()
    sample = np.asarray(sample)
    if sample.size == 0:
        return True
    return bool(np.allclose(sample, np.round(sample)) and sample.max() > 30)


def _ensure_lognorm(session: AnalysisSession) -> None:
    """Make sure adata.X is log-normalized (idempotent)."""
    a = session.require_adata()
    if a.uns.get("sctarget_lognorm"):
        return
    if _looks_like_counts(a.X):
        layer = session.cfg.dataset.get("counts_layer")
        if layer and layer not in a.layers:
            a.layers[layer] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    a.uns["sctarget_lognorm"] = True


def _batches_large_enough(a: Any, key: str, min_cells: int = 100) -> bool:
    """Per-batch HVG selection is only stable when every batch has a reasonable number of cells."""
    if key not in a.obs:
        return False
    counts = a.obs[key].value_counts()
    return bool(len(counts) > 1 and counts.min() >= min_cells)


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = np.empty(n)
    ranked[order] = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity
    ranked_sorted = np.minimum.accumulate(ranked[order][::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked_sorted, 0, 1)
    return out


def _resolve_key(session: AnalysisSession, cluster: str | None, cluster_key: str | None) -> str:
    """Pick the obs column a cluster label belongs to. Explicit key wins; otherwise the active key
    if it contains the label; otherwise any other clustering column that does (e.g. the base
    clustering after a subcluster step changed the active key)."""
    a = session.require_adata()
    if cluster_key:
        if cluster_key not in a.obs:
            raise ValueError(f"cluster_key '{cluster_key}' not in obs")
        return cluster_key
    key = session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    if cluster is None or str(cluster) in set(a.obs[key].astype(str)):
        return key
    candidates = [key.split("_sub")[0]] + [c for c in a.obs.columns if c.lower().startswith(("leiden", "louvain", "sct_leiden")) or c.endswith("_sub")]
    for c in candidates:
        if c in a.obs and str(cluster) in set(a.obs[c].astype(str)):
            return c
    raise ValueError(f"cluster '{cluster}' not found in '{key}' or any other clustering column; pass cluster_key explicitly")


def _cluster_mask(a: Any, key: str, clusters: list[str] | str) -> np.ndarray:
    labels = a.obs[key].astype(str)
    if clusters == "rest":
        return np.ones(a.n_obs, dtype=bool)
    return labels.isin([str(c) for c in clusters]).to_numpy()


def _parse_cluster_list(x: Any) -> list[str] | str:
    if x in (None, "rest", ["rest"]):
        return "rest"
    if isinstance(x, (str, int)):
        return [str(x)]
    return [str(i) for i in x]


# ----------------------------------------------------------------------------- tools
@tool(
    name="load_dataset",
    description="""
    Load the configured AnnData (.h5ad) into the local session and describe it:
    cell/gene counts, available metadata columns, conditions, samples per condition,
    existing cell-type annotations and clusterings, embeddings and layers. Call this first.
    Sample names are pseudonymized (S01, S02, ...).""",
    input_schema={"properties": {}, "required": []},
    roles=ALL,
    mutates_state=True,
)
def load_dataset(session: AnalysisSession) -> dict[str, Any]:
    import anndata as ad

    ds = session.cfg.dataset
    path = session.cfg.resolve_path(ds["path"])
    a = ad.read_h5ad(path)
    for key in (session.sample_key, session.condition_key):
        if key not in a.obs:
            raise ValueError(f"obs column '{key}' not found. Available: {list(a.obs.columns)}")
    a.obs[session.sample_key] = a.obs[session.sample_key].astype(str)
    a.obs[session.condition_key] = a.obs[session.condition_key].astype(str)
    from ..privacy import SamplePseudonymizer, check_obs_columns

    session.pseudo = SamplePseudonymizer(a.obs[session.sample_key].unique().tolist(), enabled=session.policy.pseudonymize_samples)
    session.adata = a
    if "mt" not in a.var:
        a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")

    conds = a.obs[session.condition_key].value_counts()
    samples_per_cond = (
        a.obs.groupby(session.condition_key, observed=True)[session.sample_key].nunique().to_dict()
    )
    cells_per_sample = a.obs.groupby(session.sample_key, observed=True).size()
    cells_per_sample.index = session.pseudo.map_index(cells_per_sample.index)
    cond_of_sample = a.obs.drop_duplicates(session.sample_key).set_index(session.sample_key)[session.condition_key]
    cond_of_sample.index = session.pseudo.map_index(cond_of_sample.index)

    existing_clusterings = [c for c in a.obs.columns if c.lower().startswith(("leiden", "louvain")) or c.lower() in ("clusters", "cluster", "seurat_clusters")]
    existing_clusterings = [c for c in existing_clusterings if c.lower() not in session.policy.forbidden_obs_columns]
    preferred = ds.get("reuse_clustering")
    if preferred and preferred in a.obs:
        existing_clusterings = [preferred] + [c for c in existing_clusterings if c != preferred]
    celltypes: dict[str, int] = {}
    ck = session.celltype_key
    if ck and ck in a.obs:
        celltypes = a.obs[ck].astype(str).value_counts().to_dict()

    summary = {
        "n_cells": int(a.n_obs),
        "n_genes": int(a.n_vars),
        "species": ds.get("species"),
        "disease": session.disease,
        "reference_conditions": session.reference,
        "obs_columns": check_obs_columns(list(a.obs.columns), session.policy),
        "condition_counts": {str(k): int(v) for k, v in conds.items()},
        "samples_per_condition": {str(k): int(v) for k, v in samples_per_cond.items()},
        "sample_condition": {str(k): str(v) for k, v in cond_of_sample.items()},
        "cells_per_sample": {str(k): int(v) for k, v in cells_per_sample.items()},
        "celltype_key": ck if celltypes else None,
        "lineage_source": "annotation" if celltypes else "will be derived from the data (derive_lineages) after clustering",
        "celltype_counts": {str(k): int(v) for k, v in celltypes.items()},
        "existing_clusterings": existing_clusterings,
        "preferred_clustering": preferred if preferred and preferred in a.obs else None,
        "condition_descriptions": ds.get("condition_descriptions") or {},
        "embeddings": list(a.obsm.keys()),
        "layers": list(a.layers.keys()),
        "x_looks_like_raw_counts": _looks_like_counts(a.X),
        "already_processed_flag": bool(ds.get("already_processed", False)),
        "n_mt_genes": int(a.var["mt"].sum()),
    }
    missing_ref = [r for r in session.reference if r not in summary["condition_counts"]]
    if session.disease not in summary["condition_counts"]:
        raise ValueError(f"disease condition '{session.disease}' not in {list(summary['condition_counts'])}")
    if missing_ref:
        raise ValueError(f"reference conditions {missing_ref} not found in data")
    session.facts["dataset"] = summary
    return summary


@tool(
    name="run_qc",
    description="""
    Cell/gene quality control on raw counts. Call with apply=false first to PREVIEW what the
    thresholds would remove; then apply. mode='adaptive' (default) sets the mitochondrial and
    gene-count cutoffs PER SAMPLE from the data (median + n_mads median absolute deviations, in
    log space for counts), bounded by the hard caps given; this is what transfers across tissues
    - a fixed 10% mitochondrial cutoff removes most healthy kidney tubule or cardiomyocytes, whose
    normal mitochondrial fraction is 20-45%. mode='fixed' applies the caps as the thresholds.
    Returns before/after counts overall and per sample, the thresholds actually used, and, when a
    cell-type annotation exists, the fraction of each annotated type removed (a lineage losing
    > 50% of its cells is a warning that the thresholds do not fit the tissue). Skip if the dataset
    is already QC'd (already_processed_flag).""",
    input_schema={
        "properties": {
            "mode": {"type": "string", "enum": ["adaptive", "fixed"], "default": "adaptive"},
            "n_mads": {"type": "number", "default": 3.0, "description": "adaptive: cutoff = median + n_mads * MAD (per sample)"},
            "min_genes": {"type": "integer", "description": "hard floor on detected genes per cell", "default": 200},
            "max_genes": {"type": "integer", "description": "hard cap on detected genes per cell (doublet guard)", "default": 8000},
            "max_mt_pct": {"type": "number", "description": "hard cap on % mitochondrial counts (fixed mode: the threshold; adaptive mode: the ceiling the per-sample threshold may not exceed)", "default": 50.0},
            "mt_floor_pct": {"type": "number", "default": 20.0, "description": "adaptive mode: the per-sample mito threshold never falls below this, so a mitochondria-rich lineage that is a minority in a sample is not cut wholesale"},
            "min_cells_per_gene": {"type": "integer", "default": 3},
            "apply": {"type": "boolean", "default": True, "description": "false = PREVIEW: compute thresholds and what they would remove (overall, per sample, per annotated type) without changing the data. Preview first, adjust, then apply."},
        },
        "required": [],
    },
    roles=ANALYST,
    mutates_state=True,
)
def run_qc(session: AnalysisSession, mode: str = "adaptive", n_mads: float = 3.0, min_genes: int = 200, max_genes: int = 8000, max_mt_pct: float = 50.0, mt_floor_pct: float = 20.0, min_cells_per_gene: int = 3, apply: bool = True) -> dict[str, Any]:
    a = session.require_adata()
    layer = session.cfg.dataset.get("counts_layer")
    layer = layer if layer and layer in a.layers else None
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, layer=layer, percent_top=None)
    before = a.n_obs
    before_per_sample = a.obs.groupby(session.sample_key, observed=True).size()
    samples = a.obs[session.sample_key].astype(str).to_numpy()
    ng = a.obs["n_genes_by_counts"].to_numpy().astype(float)
    mt = a.obs["pct_counts_mt"].to_numpy().astype(float)
    thresholds: dict[str, dict[str, float]] = {}
    if mode == "adaptive":
        keep = np.zeros(a.n_obs, dtype=bool)
        for smp in np.unique(samples):
            m = samples == smp
            lg = np.log1p(ng[m])
            med_g, mad_g = np.median(lg), np.median(np.abs(lg - np.median(lg))) * 1.4826
            # lower gene bound stays at the fixed floor: an adaptive lower bound cuts low-complexity
            # lineages (neutrophils, erythroid, some epithelia) that are a minority in a sample.
            # The adaptive part is the upper bound (doublet guard) and the mito threshold.
            lo_g = float(min_genes)
            hi_g = min(max_genes, float(np.expm1(med_g + n_mads * mad_g)))
            med_m, mad_m = np.median(mt[m]), np.median(np.abs(mt[m] - np.median(mt[m]))) * 1.4826
            # floor: a mitochondria-rich lineage that is a MINORITY in a sample does not move the sample's
            # median/MAD, so without a floor it would be cut wholesale; dying cells sit far above 20%
            hi_m = float(min(max_mt_pct, max(mt_floor_pct, med_m + n_mads * mad_m)))
            keep[m] = (ng[m] >= lo_g) & (ng[m] <= hi_g) & (mt[m] <= hi_m)
            thresholds[str(session.pseudo.to_pseudo(smp))] = {"min_genes": round(lo_g), "max_genes": round(hi_g), "max_mt_pct": round(hi_m, 1)}
        keep = pd.Series(keep, index=a.obs_names)
    else:
        keep = (a.obs["n_genes_by_counts"] >= min_genes) & (a.obs["n_genes_by_counts"] <= max_genes) & (a.obs["pct_counts_mt"] <= max_mt_pct)
        thresholds = {"all": {"min_genes": min_genes, "max_genes": max_genes, "max_mt_pct": max_mt_pct}}
    removed_by_type = None
    ck = session.celltype_key
    if ck and ck in a.obs:
        ct = a.obs[ck].astype(str)
        rem = pd.DataFrame({"type": ct.to_numpy(), "removed": (~keep).to_numpy()}).groupby("type")["removed"].agg(["mean", "size"])
        removed_by_type = {t: {"frac_removed": round(float(r["mean"]), 3), "n": int(r["size"])} for t, r in rem.iterrows()}
    heavy = [t for t, v in (removed_by_type or {}).items() if v["frac_removed"] > 0.5 and v["n"] >= 50 and t.lower() not in ("unannotated", "unknown", "ua", "nan")]
    if not apply:
        per_sample_prev = pd.DataFrame({"before": before_per_sample, "would_keep": pd.Series(keep.to_numpy(), index=a.obs_names).groupby(samples).sum()}).fillna(0).astype(int)
        per_sample_prev["retained_frac"] = (per_sample_prev["would_keep"] / per_sample_prev["before"]).round(3)
        per_sample_prev.index = session.pseudo.map_index(per_sample_prev.index)
        dist = pd.DataFrame({"pct_mt": mt, "n_genes": ng}).groupby(samples).agg(["median", lambda v: float(np.percentile(v, 90))])
        dist.columns = ["pct_mt_median", "pct_mt_p90", "n_genes_median", "n_genes_p90"]
        dist.index = session.pseudo.map_index(dist.index)
        out = {
            "preview": True,
            "cells_before": int(before),
            "cells_would_keep": int(keep.sum()),
            "mode": mode,
            "thresholds_proposed": thresholds,
            "per_sample": per_sample_prev,
            "per_sample_distributions": dist.round(1),
            "would_remove_by_annotated_type": removed_by_type,
            "params": {"mode": mode, "n_mads": n_mads, "min_genes": min_genes, "max_genes": max_genes, "max_mt_pct": max_mt_pct, "mt_floor_pct": mt_floor_pct},
            "note": "Nothing was changed. Adjust mode / caps if a lineage would lose most of its cells, then call again with apply=true.",
        }
        if heavy:
            out["warning"] = f"These thresholds would remove > 50% of annotated types {heavy}: likely a tissue whose normal mitochondrial fraction is high. Raise mt_floor_pct / max_mt_pct or keep adaptive mode."
        return out
    session.adata = a[keep.to_numpy()].copy()
    a = session.adata
    gene_counts = np.asarray((session.counts() > 0).sum(axis=0)).ravel()
    gkeep = gene_counts >= min_cells_per_gene
    n_genes_before = a.n_vars
    session.adata = a[:, gkeep].copy()
    a = session.adata
    after_per_sample = a.obs.groupby(session.sample_key, observed=True).size()
    per_sample = pd.DataFrame({"before": before_per_sample, "after": after_per_sample}).fillna(0).astype(int)
    per_sample["retained_frac"] = (per_sample["after"] / per_sample["before"]).round(3)
    per_sample.index = session.pseudo.map_index(per_sample.index)
    out = {
        "cells_before": int(before),
        "cells_after": int(a.n_obs),
        "genes_before": int(n_genes_before),
        "genes_after": int(a.n_vars),
        "median_genes_per_cell": float(a.obs["n_genes_by_counts"].median()),
        "median_pct_mt": float(a.obs["pct_counts_mt"].median()),
        "per_sample": per_sample,
        "mode": mode,
        "thresholds_used": thresholds,
        "removed_by_annotated_type": removed_by_type,
        "params": {"mode": mode, "n_mads": n_mads, "min_genes": min_genes, "max_genes": max_genes, "max_mt_pct": max_mt_pct, "mt_floor_pct": mt_floor_pct, "min_cells_per_gene": min_cells_per_gene},
    }
    if removed_by_type:
        if heavy:
            out["warning"] = f"QC removed > 50% of annotated types {heavy}: the thresholds may not fit this tissue (e.g. mitochondria-rich epithelia). Consider mode='adaptive' with a higher max_mt_pct cap."
    session.facts["qc"] = {k: v for k, v in out.items() if k != "per_sample"}
    return out


@tool(
    name="merge_clusters",
    description="""
    The mirror of subcluster: register the UNION of several clusters as one cluster, for a disease
    state that Leiden split into pieces (typically malignant cells partitioning by patient: each
    piece is disease-dominated and present in every disease sample, but no single piece passes the
    per-sample abundance test because each sample's cells fall mostly into one piece). Give the
    clusters to merge; they should share a lineage and each be >= 80% disease cells unless you say
    otherwise. Writes a new cluster key where the merged clusters carry one label, re-runs
    differential_abundance on it and reports the merged cluster's row, so you can see whether the
    union is enriched where the pieces were not. Reversible: the previous key stays in obs.""",
    input_schema={
        "properties": {
            "clusters": {"type": "array", "items": {"type": "string"}, "description": "cluster labels to merge (>= 2)"},
            "cluster_key": {"type": ["string", "null"]},
            "min_disease_frac": {"type": "number", "default": 0.8, "description": "each piece must have at least this fraction of disease cells (set 0 to override)"},
        },
        "required": ["clusters"],
    },
    roles=ANALYST,
    mutates_state=True,
)
def merge_clusters(session: AnalysisSession, clusters: list[str], cluster_key: str | None = None, min_disease_frac: float = 0.8) -> dict[str, Any]:
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    clusters = [str(c) for c in clusters]
    if len(clusters) < 2:
        raise ValueError("give at least two clusters to merge")
    labels = a.obs[key].astype(str)
    conds = a.obs[session.condition_key].astype(str)
    missing = [c for c in clusters if c not in set(labels)]
    if missing:
        raise ValueError(f"clusters not found in {key}: {missing}")
    fracs = {c: float((conds[labels == c] == session.disease).mean()) for c in clusters}
    low = [c for c, f in fracs.items() if f < min_disease_frac]
    if low:
        return {"error": f"clusters {low} have < {min_disease_frac:.0%} disease cells ({ {c: round(fracs[c], 2) for c in low} }); merging mixed clusters would hide a composition difference. Lower min_disease_frac only with a reason.", "disease_fracs": fracs}
    merged = "+".join(sorted(clusters, key=lambda x: (len(x), x)))
    new_key = f"{key}_merged"
    new_labels = labels.where(~labels.isin(clusters), merged)
    a.obs[new_key] = new_labels.astype(str)
    session.cluster_key = new_key
    a.obs[[new_key]].rename(columns={new_key: "cluster"}).to_csv(session.run_dir / "cluster_assignments.csv")
    session.facts["clustering"] = {**session.facts.get("clustering", {}), "cluster_key": new_key, "method": f"merged {clusters} into '{merged}' on top of {key}", "n_clusters": int(a.obs[new_key].nunique())}
    da = differential_abundance(session, cluster_key=new_key)
    tab = da.get("table")
    tab = pd.DataFrame(tab) if not isinstance(tab, pd.DataFrame) else tab
    if "cluster" in tab.columns:
        row = tab[tab["cluster"].astype(str) == merged]
    else:
        row = tab.loc[[merged]] if merged in tab.index else tab.iloc[0:0]
    rec = row.iloc[0].to_dict() if len(row) else {}
    return {"cluster_key": new_key, "merged_cluster": merged, "pieces": clusters, "piece_disease_fracs": {c: round(f, 3) for c, f in fracs.items()}, "n_cells": int((a.obs[new_key] == merged).sum()), "merged_abundance": {k: rec.get(k) for k in ("log2_ratio", "mwu_p", "n_disease_samples_present", "enriched") if k in rec}, "note": "The merged label is now the active cluster; characterise it (find_markers / pseudobulk_de against its lineage siblings) like any state. Record the merge reason in your finding."}


@tool(
    name="detect_doublets",
    description="""
    Per-sample doublet detection (Scrublet, via scanpy): writes doublet_score and predicted_doublet
    for every cell and reports the doublet rate per sample. By default doublets are FLAGGED, not
    removed: cluster_composition then reports each cluster's doublet fraction and marks
    doublet_suspect clusters (a tumour-immune or tumour-endothelial doublet cluster looks exactly
    like a disease-enriched, state-specific population). Set remove=true to drop predicted doublets
    (an explicit, recorded choice). Run after run_qc and before integrate_and_cluster.""",
    input_schema={
        "properties": {
            "expected_doublet_rate": {"type": "number", "default": 0.06},
            "remove": {"type": "boolean", "default": False},
            "threshold": {"type": ["number", "null"], "default": None, "description": "fixed score threshold instead of Scrublet's automatic one"},
        },
        "required": [],
    },
    roles=ANALYST,
    mutates_state=True,
)
def detect_doublets(session: AnalysisSession, expected_doublet_rate: float = 0.06, remove: bool = False, threshold: float | None = None) -> dict[str, Any]:
    a = session.require_adata()
    layer = session.cfg.dataset.get("counts_layer")
    seed = int(session.cfg.agent.get("seed", 0))
    scores = np.full(a.n_obs, np.nan)
    pred = np.zeros(a.n_obs, dtype=bool)
    samples = a.obs[session.sample_key].astype(str).to_numpy()
    per_sample = {}
    for smp in np.unique(samples):
        m = np.where(samples == smp)[0]
        if len(m) < 100:
            per_sample[str(session.pseudo.to_pseudo(smp))] = {"n_cells": int(len(m)), "doublet_rate": None, "note": "too few cells"}
            continue
        sub = a[m].copy()
        if layer and layer in sub.layers:
            sub.X = sub.layers[layer].copy()
        try:
            sc.pp.scrublet(sub, expected_doublet_rate=expected_doublet_rate, threshold=threshold, random_state=seed, verbose=False)
            sc_ = sub.obs["doublet_score"].to_numpy().astype(float)
            pred_col = sub.obs["predicted_doublet"]
            pd_ = pred_col.fillna(False).astype(bool).to_numpy() if pred_col.dtype == object or pred_col.isna().any() else pred_col.astype(bool).to_numpy()
            fallback = None
            if threshold is None and pd_.sum() == 0 and np.isfinite(sc_).any():
                # Scrublet's automatic threshold (bimodality of simulated doublet scores) often fails;
                # fall back to calling the top expected_doublet_rate of scores doublets
                fallback = float(np.nanquantile(sc_, 1.0 - expected_doublet_rate))
                pd_ = sc_ >= fallback
        except Exception as e:  # scrublet can fail on tiny / degenerate samples
            per_sample[str(session.pseudo.to_pseudo(smp))] = {"n_cells": int(len(m)), "doublet_rate": None, "error": f"{type(e).__name__}: {str(e)[:120]}"}
            continue
        scores[m] = sc_
        pred[m] = pd_
        thr = sub.uns.get("scrublet", {}).get("threshold") if isinstance(sub.uns.get("scrublet"), dict) else None
        per_sample[str(session.pseudo.to_pseudo(smp))] = {"n_cells": int(len(m)), "doublet_rate": round(float(pd_.mean()), 4), "score_threshold": None if (thr is None or not np.isfinite(float(thr))) else round(float(thr), 3), **({"threshold_source": f"fallback quantile at expected rate ({fallback:.3f}); Scrublet's automatic threshold found no doublets"} if fallback is not None else {})}
    a.obs["doublet_score"] = scores
    a.obs["predicted_doublet"] = pred
    n_doublets = int(pred.sum())
    out: dict[str, Any] = {"n_cells": int(a.n_obs), "n_predicted_doublets": n_doublets, "overall_rate": round(n_doublets / max(1, a.n_obs), 4), "per_sample": per_sample, "removed": False}
    if remove and n_doublets:
        session.adata = a[~pred].copy()
        out["removed"] = True
        out["n_cells_after"] = int(session.adata.n_obs)
    session.facts["doublets"] = {k: v for k, v in out.items() if k != "per_sample"}
    out["note"] = "Scores are stored per cell; cluster_composition now reports doublet_frac and doublet_suspect per cluster. Flagged by default; pass remove=true to drop predicted doublets."
    return out


@tool(
    name="integrate_and_cluster",
    description="""
    Normalize, select highly variable genes, run PCA, (optionally) Harmony batch correction
    across samples, build a kNN graph, Leiden-cluster, and compute UMAP. If reuse_existing is
    given and that obs column exists, reuse that clustering instead of recomputing.
    Returns the cluster key to use downstream, cluster sizes, and a UMAP figure path.""",
    input_schema={
        "properties": {
            "reuse_existing": {"type": ["string", "null"], "description": "name of an existing obs clustering column to reuse (e.g. 'leiden'), or null to recompute"},
            "use_harmony": {"type": "boolean", "default": True},
            "n_pcs": {"type": "integer", "default": 30},
            "n_neighbors": {"type": "integer", "default": 15},
            "resolution": {"type": "number", "default": 1.0},
            "n_top_genes": {"type": "integer", "default": 2000},
        },
        "required": [],
    },
    roles=ANALYST,
    mutates_state=True,
)
def integrate_and_cluster(
    session: AnalysisSession,
    reuse_existing: str | None = None,
    use_harmony: bool = True,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolution: float = 1.0,
    n_top_genes: int = 2000,
) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    seed = int(session.cfg.agent.get("seed", 0))
    harmony_used = False

    if reuse_existing and reuse_existing in a.obs:
        key = reuse_existing
        method = f"reused existing '{reuse_existing}'"
    else:
        key = f"sct_leiden_r{resolution:g}"
        # never select (almost) every gene as "highly variable" (cap ~20% of genes, floor 300),
        # and only compute HVGs per batch when every batch is big enough for stable estimates
        n_top_genes = int(min(n_top_genes, max(300, a.n_vars // 5), a.n_vars - 1))
        batch_key = session.sample_key if _batches_large_enough(a, session.sample_key) else None
        sc.pp.highly_variable_genes(a, n_top_genes=n_top_genes, batch_key=batch_key)
        n_pcs = min(n_pcs, a.n_vars - 1, a.n_obs - 1)
        sc.tl.pca(a, n_comps=n_pcs, mask_var="highly_variable", random_state=seed)
        rep = "X_pca"
        if use_harmony:
            try:
                a.obsm["X_pca_harmony"] = _run_harmony(a, session.sample_key, seed)
                rep = "X_pca_harmony"
                harmony_used = True
            except Exception as e:
                # never swap the method silently: an auditable pipeline makes the choice explicit
                session.log.note("analyst", f"harmony failed ({e})")
                return {
                    "error": f"Harmony integration failed: {type(e).__name__}: {str(e)[:300]}",
                    "hint": "Call integrate_and_cluster again with use_harmony=false to cluster on uncorrected PCA (clusters may then split by sample; check per-sample proportions), or fix the integration environment. The choice is recorded either way.",
                }
        sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep=rep, random_state=seed)
        sc.tl.leiden(a, resolution=resolution, key_added=key, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
        sc.tl.umap(a, random_state=seed)
        method = "recomputed"

    a.obs[key] = a.obs[key].astype(str)
    session.cluster_key = key
    sizes = a.obs[key].value_counts()
    # per-cell assignments stay on disk in the run dir (never sent to the model); used by the evaluator
    a.obs[[key]].rename(columns={key: "cluster"}).to_csv(session.run_dir / "cluster_assignments.csv")

    fig_path = None
    if "X_umap" in a.obsm:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        sc.pl.umap(a, color=key, ax=axes[0], show=False, legend_loc="on data", title=f"clusters ({key})", frameon=False)
        sc.pl.umap(a, color=session.condition_key, ax=axes[1], show=False, title="condition", frameon=False)
        p = session.save_figure(f"umap_{key}", fig, f"UMAP coloured by cluster ({key}) and by condition.")
        fig_path = str(p.relative_to(session.run_dir))

    out = {
        "cluster_key": key,
        "method": method,
        "harmony_used": harmony_used,
        "n_clusters": int(len(sizes)),
        "cluster_sizes": {str(k): int(v) for k, v in sizes.items()},
        "params": {"n_pcs": n_pcs, "n_neighbors": n_neighbors, "resolution": resolution, "n_top_genes": n_top_genes},
        "figure": fig_path,
    }
    session.facts["clustering"] = out
    return out


@tool(
    name="subcluster",
    description="""
    Re-cluster only the cells of the given clusters (e.g. all keratinocyte clusters) at the
    chosen resolution, producing labels like '0.1', '0.2' in a new cluster key that keeps all
    other cells' labels unchanged. Use when a disease-specific state may be hidden inside a
    coarse cluster: subclustering one cell type isolates rare states that global clustering
    merges. Sets the new key as the session cluster key; re-run differential_abundance after.""",
    input_schema={
        "properties": {
            "clusters": {"type": "array", "items": {"type": "string"}, "description": "clusters to subcluster together"},
            "cluster_key": {"type": ["string", "null"]},
            "resolution": {"type": "number", "default": 1.0},
            "n_pcs": {"type": "integer", "default": 20},
            "n_neighbors": {"type": "integer", "default": 15},
        },
        "required": ["clusters"],
    },
    roles=ANALYST,
    mutates_state=True,
)
def subcluster(session: AnalysisSession, clusters: list[str], cluster_key: str | None = None, resolution: float = 1.0, n_pcs: int = 20, n_neighbors: int = 15) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    clusters = [str(c) for c in clusters]
    labels = a.obs[key].astype(str)
    mask = labels.isin(clusters).to_numpy()
    if mask.sum() < 50:
        raise ValueError(f"only {int(mask.sum())} cells in {clusters}; too few to subcluster")
    seed = int(session.cfg.agent.get("seed", 0))
    sub = a[mask].copy()
    n_top = int(min(2000, max(300, sub.n_vars // 5), sub.n_vars - 1))
    batch_key = session.sample_key if _batches_large_enough(sub, session.sample_key) else None
    sc.pp.highly_variable_genes(sub, n_top_genes=n_top, batch_key=batch_key)
    n_pcs = int(min(n_pcs, sub.n_vars - 1, sub.n_obs - 1))
    sc.tl.pca(sub, n_comps=n_pcs, mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(sub, n_neighbors=min(n_neighbors, sub.n_obs - 1), use_rep="X_pca", random_state=seed)
    sc.tl.leiden(sub, resolution=resolution, key_added="_sub", random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    prefix = "+".join(clusters)
    new_key = f"{key}_sub"
    new_labels = labels.copy()
    new_labels.loc[sub.obs_names] = [f"{prefix}.{s}" for s in sub.obs["_sub"].astype(str)]
    a.obs[new_key] = new_labels.astype(str)
    session.cluster_key = new_key
    a.obs[[new_key]].rename(columns={new_key: "cluster"}).to_csv(session.run_dir / "cluster_assignments.csv")
    sizes = a.obs[new_key].value_counts()
    new_sizes = {k: int(v) for k, v in sizes.items() if k.startswith(prefix + ".")}
    session.facts["clustering"] = {**session.facts.get("clustering", {}), "cluster_key": new_key, "method": f"subclustered {clusters} at resolution {resolution} on top of {key}", "n_clusters": int(len(sizes))}
    return {"cluster_key": new_key, "subclustered": clusters, "n_subclusters": len(new_sizes), "subcluster_sizes": new_sizes, "params": {"resolution": resolution, "n_pcs": n_pcs, "n_neighbors": n_neighbors, "n_top_genes": n_top}}


@tool(
    name="cluster_composition",
    description="""
    For each cluster: fraction of cells from each condition, number of samples contributing
    (>= min_cells), and the dominant existing cell-type annotation if one exists. Use it to
    see which clusters are condition-restricted and which cell type they belong to.""",
    input_schema={
        "properties": {
            "cluster_key": {"type": ["string", "null"], "description": "obs column; defaults to the session cluster key"},
            "min_cells": {"type": "integer", "default": 10},
        },
        "required": [],
    },
    roles=ALL,
)
def cluster_composition(session: AnalysisSession, cluster_key: str | None = None, min_cells: int = 10) -> dict[str, Any]:
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    ct = pd.crosstab(a.obs[key].astype(str), a.obs[session.condition_key])
    frac = (ct.T / ct.sum(axis=1)).T.round(3)
    nsamp = (
        a.obs.groupby([key, session.sample_key], observed=True).size().unstack(fill_value=0) >= min_cells
    ).sum(axis=1)
    df = frac.copy()
    df.columns = [f"frac_{c}" for c in df.columns]
    df["n_cells"] = ct.sum(axis=1)
    df["n_samples_contributing"] = nsamp.reindex(df.index).fillna(0).astype(int)
    ck = _lineage_key(session)
    if ck and ck in a.obs:
        dom = pd.crosstab(a.obs[key].astype(str), a.obs[ck].astype(str))
        df["dominant_celltype"] = dom.idxmax(axis=1)
        df["dominant_celltype_frac"] = (dom.max(axis=1) / dom.sum(axis=1)).round(3)
    if "predicted_doublet" in a.obs:
        dfrac = a.obs.groupby(a.obs[key].astype(str), observed=True)["predicted_doublet"].mean()
        df["doublet_frac"] = dfrac.reindex(df.index).round(3)
        overall = float(a.obs["predicted_doublet"].mean())
        # suspect: well above the dataset rate and above 1 in 5 cells
        df["doublet_suspect"] = (df["doublet_frac"] >= max(0.2, 3 * overall))
    df.index.name = "cluster"
    session.save_table(f"cluster_composition_{key}", df)
    session.facts["cluster_composition"] = df.reset_index().to_dict(orient="records")
    return {"cluster_key": key, "table": df}


@tool(
    name="differential_abundance",
    description="""
    Test which clusters are enriched in the disease condition versus the reference
    conditions using per-SAMPLE proportions (so that one large sample cannot drive the
    result). Reports mean proportion per group, log2 ratio, Mann-Whitney p-value, number of
    disease samples in which the cluster is present, and a disease_enriched flag.
    This is the primary tool for finding disease-specific cell states.""",
    input_schema={
        "properties": {
            "cluster_key": {"type": ["string", "null"]},
            "min_cells_per_sample": {"type": "integer", "default": 20, "description": "samples with fewer cells are ignored"},
            "enrichment_log2fc": {"type": "number", "default": 1.0, "description": "log2 ratio threshold for the disease_enriched flag"},
            "within_celltype": {"type": "boolean", "default": False, "description": "compute proportions within the cluster's dominant cell type instead of all cells"},
            "min_cluster_cells": {"type": "integer", "default": 30, "description": "clusters smaller than this are reported but never flagged disease_enriched"},
            "max_mwu_p": {"type": ["number", "null"], "default": None, "description": "Mann-Whitney p on per-sample proportions must be below this for the disease_enriched flag (guards against chance enrichment of small subclusters)"},
        },
        "required": [],
    },
    roles=ALL,
)
def differential_abundance(
    session: AnalysisSession,
    cluster_key: str | None = None,
    min_cells_per_sample: int = 20,
    enrichment_log2fc: float = 1.0,
    within_celltype: bool = False,
    min_cluster_cells: int = 30,
    max_mwu_p: float | None = None,
) -> dict[str, Any]:
    if max_mwu_p is None:
        max_mwu_p = float(session.cfg.get("analysis.differential_abundance.max_mwu_p", 0.05))
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    obs = a.obs[[key, session.sample_key, session.condition_key]].copy()
    obs[key] = obs[key].astype(str)
    ck = _lineage_key(session)
    if within_celltype and ck and ck in a.obs:
        obs["_ct"] = a.obs[ck].astype(str).to_numpy()
        dom = pd.crosstab(obs[key], obs["_ct"]).idxmax(axis=1)

    sample_sizes = obs.groupby(session.sample_key, observed=True).size()
    good_samples = sample_sizes[sample_sizes >= min_cells_per_sample].index
    obs = obs[obs[session.sample_key].isin(good_samples)]
    counts = pd.crosstab(obs[session.sample_key], obs[key])
    if within_celltype and ck and ck in a.obs:
        # numerator: cells of the cluster that carry its dominant cell-type label;
        # denominator: all cells of that cell type in the sample  -> proportion <= 1
        prop = pd.DataFrame(index=counts.index, columns=counts.columns, dtype=float)
        for cl in counts.columns:
            in_type = obs["_ct"] == dom[cl]
            num = obs[(obs[key] == cl) & in_type].groupby(session.sample_key, observed=True).size().reindex(counts.index).fillna(0)
            denom = obs[in_type].groupby(session.sample_key, observed=True).size().reindex(counts.index).fillna(0)
            prop[cl] = num / denom.replace(0, np.nan)
        prop = prop.fillna(0.0)
    else:
        prop = counts.div(counts.sum(axis=1), axis=0)
    cond_of = obs.drop_duplicates(session.sample_key).set_index(session.sample_key)[session.condition_key]
    is_dis = cond_of.reindex(prop.index) == session.disease
    is_ref = cond_of.reindex(prop.index).isin(session.reference)

    rows = []
    eps = 1e-3
    for cl in prop.columns:
        d = prop.loc[is_dis, cl].to_numpy()
        r = prop.loc[is_ref, cl].to_numpy()
        md, mr = float(d.mean()) if len(d) else 0.0, float(r.mean()) if len(r) else 0.0
        l2 = float(np.log2((md + eps) / (mr + eps)))
        if len(d) >= 2 and len(r) >= 2:
            try:
                p = float(stats.mannwhitneyu(d, r, alternative="two-sided").pvalue)
            except ValueError:
                p = 1.0
        else:
            p = float("nan")
        n_dis_present = int((d > 0).sum())
        frac_from_disease = float((obs[obs[key] == cl][session.condition_key] == session.disease).mean())
        rows.append(
            {
                "cluster": cl,
                "n_cells": int(counts[cl].sum()),
                "mean_prop_disease": round(md, 4),
                "mean_prop_reference": round(mr, 4),
                "log2_ratio": round(l2, 3),
                "mannwhitney_p": p,
                "n_disease_samples_present": n_dis_present,
                "n_disease_samples_total": int(is_dis.sum()),
                "frac_cells_from_disease": round(frac_from_disease, 3),
                "per_sample_prop_disease": [round(float(x), 4) for x in d],
                "per_sample_prop_reference": [round(float(x), 4) for x in r],
            }
        )
    df = pd.DataFrame(rows).set_index("cluster")
    p_ok = df["mannwhitney_p"].fillna(1.0) < max_mwu_p
    df["disease_enriched"] = (df["log2_ratio"] >= enrichment_log2fc) & (df["n_disease_samples_present"] >= min(2, int(is_dis.sum()))) & (df["n_cells"] >= min_cluster_cells) & p_ok
    df = df.sort_values("log2_ratio", ascending=False)
    session.save_table(f"differential_abundance_{key}", df)

    # figure: per-sample proportions
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(max(5, 0.8 * len(prop.columns) + 2), 4))
    xs = np.arange(len(prop.columns))
    rng = np.random.default_rng(0)
    for i, cl in enumerate(prop.columns):
        ax.scatter(xs[i] - 0.15 + rng.normal(0, 0.03, is_dis.sum()), prop.loc[is_dis, cl], color="#c0392b", s=18, label=session.disease if i == 0 else None)
        ax.scatter(xs[i] + 0.15 + rng.normal(0, 0.03, is_ref.sum()), prop.loc[is_ref, cl], color="#2c3e50", s=18, label="reference" if i == 0 else None)
    ax.set_xticks(xs)
    ax.set_xticklabels(prop.columns, rotation=45, ha="right")
    ax.set_ylabel("proportion of sample")
    ax.set_title("Per-sample cluster proportions")
    ax.legend(frameon=False)
    p = session.save_figure(
        f"differential_abundance_{key}{'_within_celltype' if within_celltype else ''}",
        fig,
        "Per-sample cluster proportions" + (" within the dominant cell type" if within_celltype else "") + ", disease (red) vs reference (dark).",
    )

    enriched = df.index[df["disease_enriched"]].tolist()
    # disease-dominated clusters that are present in most disease samples yet not significantly
    # enriched: the signature of ONE state split across clusters (e.g. malignant cells partitioning
    # by patient), where each piece is bimodal across samples. Reported so the model can merge them.
    near = df[(~df["disease_enriched"]) & (df["frac_cells_from_disease"] >= 0.8) & (df["n_disease_samples_present"] >= max(2, int(0.5 * is_dis.sum()))) & (df["n_cells"] >= min_cluster_cells)]
    split_candidates = near.index.astype(str).tolist()
    facts_key = "differential_abundance_within_celltype" if within_celltype else "differential_abundance"
    session.facts[facts_key] = {
        "cluster_key": key,
        "within_celltype": bool(within_celltype),
        "disease_enriched_clusters": enriched,
        "split_state_candidates": split_candidates,
        "table": df.drop(columns=["per_sample_prop_disease", "per_sample_prop_reference"]).reset_index().to_dict(orient="records"),
    }
    return {
        "cluster_key": key,
        "within_celltype": bool(within_celltype),
        "n_samples_used": {"disease": int(is_dis.sum()), "reference": int(is_ref.sum())},
        "disease_enriched_clusters": enriched,
        "split_state_candidates": split_candidates,
        "split_state_note": (f"Clusters {split_candidates} are >= 80% disease cells and present in most disease samples but NOT significantly enriched: when several such clusters share a lineage this is usually one state split by Leiden (malignant cells partitioning by patient). Test their union with merge_clusters." if len(split_candidates) >= 2 else None),
        "table": df,
        "figure": str(p.relative_to(session.run_dir)),
        "note": f"Small n per group: treat p-values as descriptive; prioritize log2_ratio and n_disease_samples_present. Clusters with < {min_cluster_cells} cells or MWU p >= {max_mwu_p} are never flagged enriched.",
    }


@tool(
    name="find_markers",
    description="""
    Marker genes for one cluster versus a reference (all other cells, or a list of clusters)
    using a Wilcoxon rank-sum test on log-normalized expression. Returns top genes with
    log2 fold change, adjusted p-value, and percent of cells expressing in/out of the group.
    Use reference_clusters to compare a state against its parent cell type's other states.""",
    input_schema={
        "properties": {
            "cluster": {"type": "string", "description": "cluster label to characterize"},
            "cluster_key": {"type": ["string", "null"]},
            "reference_clusters": {"type": ["array", "null"], "items": {"type": "string"}, "description": "clusters to compare against; null = all other cells"},
            "n_top": {"type": "integer", "default": 30},
        },
        "required": ["cluster"],
    },
    roles=ALL,
)
def find_markers(session: AnalysisSession, cluster: str, cluster_key: str | None = None, reference_clusters: list[str] | None = None, n_top: int = 30) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    cluster = str(cluster)
    key = _resolve_key(session, cluster, cluster_key)
    labels = a.obs[key].astype(str)
    ref = _parse_cluster_list(reference_clusters)
    if ref == "rest":
        sub = a
    else:
        sub = a[labels.isin([cluster] + list(ref)).to_numpy()].copy()
    sub.obs["_grp"] = np.where(sub.obs[key].astype(str) == cluster, "target", "reference")
    sc.tl.rank_genes_groups(sub, groupby="_grp", groups=["target"], reference="reference", method="wilcoxon", pts=True)
    df = sc.get.rank_genes_groups_df(sub, group="target")
    df = df.rename(columns={"names": "gene", "logfoldchanges": "log2fc", "pvals_adj": "padj"})
    pts = sub.uns["rank_genes_groups"]["pts"]
    df["pct_in"] = df["gene"].map(pts["target"]).astype(float)
    if "reference" in pts.columns:
        df["pct_out"] = df["gene"].map(pts["reference"]).astype(float)
    else:  # fall back to computing it directly
        ref_mask = (sub.obs["_grp"] == "reference").to_numpy()
        pct_ref = pd.Series(_pct_expressing(sub.X, ref_mask), index=sub.var_names)
        df["pct_out"] = df["gene"].map(pct_ref).astype(float)
    df = df[["gene", "log2fc", "padj", "pct_in", "pct_out", "scores"]].set_index("gene")
    df = df.sort_values("scores", ascending=False)
    session.save_table(f"markers_{key}_{cluster}_vs_{'rest' if ref == 'rest' else '-'.join(ref)}", df)
    up_all = df[(df["log2fc"] > 0)]
    stress = stress_fraction(list(up_all.head(30).index))
    is_stress = up_all.index.to_series().str.upper().isin(STRESS_GENES).to_numpy()
    # the dissociation programme is often SUPERIMPOSED on a real population (and differs between
    # conditions when tissues were handled differently), so: report markers with the stress genes
    # removed, and call the cluster an artefact only when little survives that removal
    up_clean = up_all[~is_stress]
    strong_clean = int(((up_clean["log2fc"] >= 1.0) & (up_clean["padj"] < 0.05)).sum())
    suspect = bool(stress >= 0.3 and strong_clean < 10)
    up = (up_clean if stress >= 0.3 else up_all).head(n_top)
    session.facts.setdefault("markers", {})[f"{cluster}_vs_{'rest' if ref == 'rest' else ','.join(ref)}"] = up.head(50).reset_index().to_dict(orient="records")
    if suspect:
        session.facts.setdefault("stress_suspect_clusters", {})[str(cluster)] = round(stress, 2)
    else:
        session.facts.get("stress_suspect_clusters", {}).pop(str(cluster), None)
    out = {
        "cluster": cluster,
        "cluster_key": key,
        "reference": "rest" if ref == "rest" else ref,
        "n_cells_target": int((labels == cluster).sum()),
        "n_cells_reference": int(sub.n_obs - (labels == cluster).sum()),
        "dissociation_stress_fraction_top30": round(stress, 2),
        "n_strong_markers_excluding_stress_genes": strong_clean,
        "stress_artifact_suspect": suspect,
        "top_markers": up,
    }
    if stress >= 0.3 and not suspect:
        out["note"] = f"{stress:.0%} of the raw top-30 markers are dissociation-stress genes (FOS, JUN, HSPA1A, DNAJB1 ...); they were removed from top_markers and {strong_clean} strong markers remain, so this is a real population carrying a stress response (handling differs between conditions). Do not nominate the stress genes."
    if suspect:
        out["warning"] = f"{stress:.0%} of the top-30 markers are dissociation-stress genes and only {strong_clean} strong markers remain without them: this cluster is likely a dissociation-stress artefact, not a disease state."
    return out


@tool(
    name="pseudobulk_de",
    description="""
    Replicate-aware differential expression: counts are summed per sample within each
    group (pseudobulk) and tested with DESeq2 (PyDESeq2; falls back to Wilcoxon on
    log-CPM if unavailable). The cluster-vs-siblings contrast is blocked by sample
    (~ sample + label: a sample's state cells vs the same sample's sibling cells) whenever
    >= 3 samples contribute both. Two contrasts:
      - 'cluster_vs_reference_clusters': cells of `cluster` vs cells of `reference_clusters`
        (default: all other clusters). Finds state-specific genes with sample-level replication.
      - 'disease_vs_reference_within_cluster': within `cluster`, disease samples vs reference
        samples. Finds disease-driven changes inside one cell type/state.
      - 'disease_vs_reference_within_compartment': within the SHARED clusters of the cluster's
        compartment (same dominant cell type, excluding the disease-enriched clusters), disease vs
        reference samples. Finds disease-wide genes of that lineage: up in disease cells of the
        ordinary, condition-shared clusters too, not only in the disease-specific state.
    Returns top up-regulated genes ranked by test statistic among genes with
    baseMean >= min_base_mean (so rare, barely-expressed genes do not dominate), plus the number
    of pseudobulk replicates per side. Requires >=2 replicates per side.""",
    input_schema={
        "properties": {
            "cluster": {"type": "string"},
            "cluster_key": {"type": ["string", "null"]},
            "contrast": {"type": "string", "enum": ["cluster_vs_reference_clusters", "disease_vs_reference_within_cluster", "disease_vs_reference_within_compartment"], "default": "cluster_vs_reference_clusters"},
            "reference_clusters": {"type": ["array", "null"], "items": {"type": "string"}},
            "min_cells_per_pseudobulk": {"type": "integer", "default": 10},
            "n_top": {"type": "integer", "default": 40},
            "padj_threshold": {"type": "number", "default": 0.05},
            "min_base_mean": {"type": "number", "default": 10.0, "description": "genes with lower mean pseudobulk counts are excluded from the ranked lists (noise guard)"},
            "min_frac_of_cluster": {"type": "number", "default": 0.01, "description": "cluster_vs_reference_clusters only: a target pseudobulk must hold at least this fraction of the cluster's cells (drops spill-over replicates)"},
        },
        "required": ["cluster"],
    },
    roles=ALL,
)
def pseudobulk_de(
    session: AnalysisSession,
    cluster: str,
    cluster_key: str | None = None,
    contrast: str = "cluster_vs_reference_clusters",
    reference_clusters: list[str] | None = None,
    min_cells_per_pseudobulk: int = 10,
    n_top: int = 40,
    padj_threshold: float = 0.05,
    min_base_mean: float = 10.0,
    min_frac_of_cluster: float = 0.01,
) -> dict[str, Any]:
    out = _run_pseudobulk(session, cluster, cluster_key, contrast, reference_clusters, min_cells_per_pseudobulk, padj_threshold, min_base_mean, record=True, min_frac_of_cluster=min_frac_of_cluster)
    if "error" in out:
        return out
    out["top_up"] = out.pop("_sig_up").head(n_top)
    return out


def _run_harmony(a, batch_key: str, seed: int) -> np.ndarray:
    """Harmony on X_pca via harmonypy directly. harmonypy < 1 returns PCs x cells, >= 1 (incl. 2.x)
    returns cells x PCs (and may return a torch tensor); scanpy's wrapper assumes the old layout and
    produces a wrong-shaped obsm. Normalise to cells x PCs here."""
    import harmonypy as hm

    ho = hm.run_harmony(np.asarray(a.obsm["X_pca"]), a.obs[[batch_key]].copy(), batch_key, random_state=seed, verbose=False)
    Z = ho.Z_corr
    if hasattr(Z, "detach"):
        Z = Z.detach().cpu().numpy()
    Z = np.asarray(Z)
    if Z.shape[0] != a.n_obs:
        Z = Z.T
    if Z.shape[0] != a.n_obs:
        raise ValueError(f"harmony returned shape {Z.shape}, expected ({a.n_obs}, n_pcs)")
    return Z


def _lineage_key(session: AnalysisSession) -> str | None:
    """Lineage column for structural logic. Uses the configured annotation when present; otherwise
    derives coarse lineages from the data (derive_lineages) the first time they are needed."""
    key = session.lineage_key
    if key is None and session.adata is not None and session.cluster_key:
        try:
            derive_lineages(session)
            key = session.lineage_key
        except Exception as e:  # structural fallback: tools then treat all clusters as one lineage
            session.log.note("analyst", f"derive_lineages failed ({e}); proceeding without lineages")
    return key


@tool(
    name="derive_lineages",
    description="""
    Derive broad lineages (compartments) FROM THE DATA, with no cell-type annotation: a coarse
    Leiden clustering (default resolution 0.2) on the same embedding as the fine clustering, with
    every fine cluster assigned to the coarse group holding most of its cells. Lineages define
    sibling sets, compartments and the ambient-control lineage downstream. Used automatically when
    the dataset has no celltype_key; call it explicitly to choose the resolution or to override an
    annotation (analysis.lineages.ignore_annotation). Returns the fine->coarse mapping and sanity
    checks (a coarse group that holds a single fine cluster cannot provide siblings). Names are NOT
    assigned here: see label_lineages, which proposes marker-based labels for interpretation only.""",
    input_schema={
        "properties": {
            "resolution": {"type": ["number", "null"], "default": None, "description": "coarse Leiden resolution (config analysis.lineages.resolution, default 0.2)"},
            "cluster_key": {"type": ["string", "null"]},
        },
        "required": [],
    },
    roles=ALL,
)
def derive_lineages(session: AnalysisSession, resolution: float | None = None, cluster_key: str | None = None) -> dict[str, Any]:
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        return {"error": "no clustering yet: run integrate_and_cluster first"}
    lcfg = session.cfg.analysis.get("lineages", {})
    resolution = float(resolution if resolution is not None else lcfg.get("resolution", 0.2))
    seed = int(session.cfg.agent.get("seed", 0))
    if "connectivities" not in a.obsp:
        rep = "X_pca_harmony" if "X_pca_harmony" in a.obsm else ("X_pca" if "X_pca" in a.obsm else None)
        if rep is None:
            _ensure_lognorm(session)
            sc.pp.highly_variable_genes(a, n_top_genes=min(2000, max(300, a.n_vars // 5)))
            sc.tl.pca(a, n_comps=min(30, a.n_vars - 1, a.n_obs - 1), mask_var="highly_variable", random_state=seed)
            rep = "X_pca"
        sc.pp.neighbors(a, n_neighbors=15, use_rep=rep, random_state=seed)
    tmp = "_sct_coarse"
    sc.tl.leiden(a, resolution=resolution, key_added=tmp, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    fine = a.obs[key].astype(str)
    coarse = a.obs[tmp].astype(str)
    # each fine cluster belongs to ONE lineage: the coarse group holding most of its cells
    ctab = pd.crosstab(fine, coarse)
    fine_to_coarse = ctab.idxmax(axis=1)
    purity = (ctab.max(axis=1) / ctab.sum(axis=1)).round(3)
    lineage = fine.map(fine_to_coarse).astype(str)
    a.obs[session.DERIVED_LINEAGE_KEY] = ("L" + lineage).to_numpy()
    del a.obs[tmp]
    groups = {}
    for g, members in fine_to_coarse.groupby(fine_to_coarse):
        groups["L" + str(g)] = {"fine_clusters": sorted(members.index.astype(str), key=lambda x: (len(x), x)), "n_cells": int((lineage == g).sum())}
    singletons = [g for g, v in groups.items() if len(v["fine_clusters"]) < int(lcfg.get("min_fine_per_group", 2))]
    session.facts["lineages"] = {"key": session.DERIVED_LINEAGE_KEY, "resolution": resolution, "groups": groups, "fine_to_lineage": {str(k): "L" + str(v) for k, v in fine_to_coarse.items()}, "fine_cluster_purity": {str(k): float(v) for k, v in purity.items()}, "singleton_groups": singletons}
    session.log.stage("analyst", "lineages_derived", resolution=resolution, n_groups=len(groups), singletons=singletons)
    return {
        "lineage_key": session.DERIVED_LINEAGE_KEY,
        "resolution": resolution,
        "n_lineages": len(groups),
        "groups": groups,
        "fine_cluster_purity": {str(k): float(v) for k, v in purity.items()},
        "singleton_groups": singletons,
        "note": "Lineages are data-derived compartments used for sibling sets, compartments and ambient controls. A singleton group (one fine cluster) gives that cluster no siblings: consider a higher resolution for the fine clustering or a lower one here. Labels are interpretation only (label_lineages).",
    }


_LINEAGE_MARKERS: dict[str, list[str]] = {
    "epithelial": ["EPCAM", "KRT5", "KRT14", "KRT8", "KRT18", "KRT19", "CDH1", "TP63"],
    "T/NK": ["CD3E", "CD3D", "CD2", "CD247", "NKG7", "GNLY", "IL7R"],
    "B/plasma": ["MS4A1", "CD79A", "CD79B", "JCHAIN", "MZB1", "IGKC"],
    "myeloid": ["CD68", "LYZ", "CSF1R", "C1QA", "C1QB", "ITGAX", "CD14"],
    "mast": ["KIT", "CPA3", "TPSAB1", "MS4A2"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "FLT1", "KDR", "CLDN5"],
    "fibroblast/mesenchymal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"],
    "smooth muscle/pericyte": ["ACTA2", "MYH11", "RGS5", "PDGFRB", "TAGLN"],
    "neural/glial": ["PLP1", "SOX10", "S100B", "NRXN1", "CDH19"],
    "melanocyte": ["PMEL", "MLANA", "TYR", "DCT"],
    "erythroid": ["HBB", "HBA1", "HBA2", "ALAS2"],
    "hepatocyte": ["ALB", "APOA1", "TTR", "SERPINA1"],
    "cardiomyocyte": ["TNNT2", "MYH6", "MYH7", "RYR2"],
}


@tool(
    name="label_lineages",
    description="""
    Propose names for lineages / clusters from canonical marker genes, FOR INTERPRETATION ONLY:
    nothing downstream reads the labels (sibling sets and compartments come from the hierarchy).
    Scores each group against a small species-agnostic marker table (symbols matched case-
    insensitively, via the human_symbol column when present) and returns the best label with a
    confidence (margin over the runner-up) and the markers that supported it. Low-confidence
    labels should be reported as such, never asserted. Groups default to the lineage key.""",
    input_schema={
        "properties": {
            "group_key": {"type": ["string", "null"], "description": "obs column to label (default: the lineage key; use the cluster key to label fine clusters)"},
        },
        "required": [],
    },
    roles=ALL,
)
def label_lineages(session: AnalysisSession, group_key: str | None = None) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    key = group_key or _lineage_key(session) or session.cluster_key
    if not key or key not in a.obs:
        return {"error": "no lineage or cluster column to label"}
    hs_col = session.cfg.dataset.get("human_symbol_var")
    names = a.var[hs_col].astype(str) if hs_col and hs_col in a.var else pd.Series(a.var_names, index=a.var_names)
    upper = {str(v).upper(): i for i, v in zip(range(a.n_vars), names)}
    groups = a.obs[key].astype(str).to_numpy()
    uniq = sorted(np.unique(groups), key=lambda x: (len(x), x))
    out = {}
    for g in uniq:
        m = groups == g
        scores = {}
        support = {}
        for lab, markers in _LINEAGE_MARKERS.items():
            idx = [upper[mk] for mk in markers if mk in upper]
            if not idx:
                continue
            sub = a.X[m][:, idx]
            frac = np.asarray((sub > 0).mean(axis=0)).ravel()
            scores[lab] = float(frac.mean())
            support[lab] = [markers[k] for k, f in zip([markers.index(mk) for mk in markers if mk in upper], frac) if f >= 0.25]
        if not scores:
            out[g] = {"label": None, "confidence": 0.0, "note": "no marker genes found in the dataset"}
            continue
        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        best, second = ranked[0], (ranked[1] if len(ranked) > 1 else (None, 0.0))
        conf = float(best[1] - second[1]) / max(best[1], 1e-6)
        out[g] = {"label": best[0], "confidence": round(conf, 3), "score": round(best[1], 3), "runner_up": second[0], "supporting_markers": support.get(best[0], [])}
    session.facts.setdefault("lineage_labels", {})[key] = out
    return {"group_key": key, "labels": out, "note": "Interpretation only. Confidence is the margin over the runner-up label; below ~0.3 report the group by its id."}


# van den Brink et al. 2017 / Denisenko 2020: dissociation-induced immediate-early / heat-shock programme
STRESS_GENES = {"FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "EGR2", "IER2", "IER3", "IER5", "KLF2", "KLF4", "KLF6", "DUSP1", "DUSP2", "ZFP36", "NR4A1", "NR4A2", "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPB1", "HSPH1", "HSP90AA1", "HSPE1", "HSPD1", "DNAJB1", "DNAJA1", "BAG3", "PPP1R15A", "SOCS3", "CEBPD", "GADD45B", "BTG2", "MCL1", "NFKBIA", "CXCL2", "RHOB"}


def stress_fraction(genes: list[str]) -> float:
    """Fraction of a gene list that belongs to the dissociation-stress programme (case-insensitive)."""
    g = [str(x).upper() for x in genes]
    return float(sum(1 for x in g if x in STRESS_GENES) / max(1, len(g)))


def compartment_mask(session: AnalysisSession, cluster: str, key: str, shared_only: bool = True) -> np.ndarray:
    """Cells of the clusters sharing `cluster`'s dominant cell type (its compartment / lineage).
    With shared_only (the default, used for the disease-wide axis) the disease-enriched clusters
    - `cluster` itself and any other cluster flagged by differential_abundance - are EXCLUDED, so
    the disease-vs-reference contrast asks whether a gene is up in disease cells of the lineage's
    ordinary, condition-shared clusters, rather than re-measuring the composition shift (a state
    that makes up most of the disease compartment would otherwise score 'disease-wide' trivially).
    Without a cell-type annotation the compartment is all cells (minus the excluded clusters)."""
    a = session.require_adata()
    labels = a.obs[key].astype(str).to_numpy()
    in_c = labels == str(cluster)
    excluded = {str(cluster)} | set(map(str, session.facts.get("differential_abundance", {}).get("disease_enriched_clusters", []))) if shared_only else set()
    ck = _lineage_key(session)
    if not ck or ck not in a.obs:
        return ~np.isin(labels, list(excluded)) if excluded else np.ones(a.n_obs, dtype=bool)
    ct = a.obs[ck].astype(str).to_numpy()
    dom = pd.Series(ct[in_c]).value_counts().idxmax() if in_c.any() else None
    if dom is None:
        return in_c
    # a cluster belongs to the compartment when the compartment cell type is its dominant one
    cl_dom = pd.Series(ct).groupby(labels).agg(lambda v: v.value_counts().idxmax())
    members = (set(cl_dom[cl_dom == dom].index.astype(str)) | {str(cluster)}) - excluded
    if not members:  # the whole lineage is disease-enriched: nothing shared to compare in
        return np.zeros(a.n_obs, dtype=bool)
    return np.isin(labels, list(members))


def _run_pseudobulk(
    session: AnalysisSession,
    cluster: str,
    cluster_key: str | None,
    contrast: str,
    reference_clusters: list[str] | None,
    min_cells_per_pseudobulk: int,
    padj_threshold: float,
    min_base_mean: float,
    record: bool = True,
    min_frac_of_cluster: float = 0.01,
) -> dict[str, Any]:
    """Shared engine for pseudobulk_de and disease_programs. Returns the summary dict with the
    full significant-up frame under '_sig_up' (caller decides how many rows to expose)."""
    a = session.require_adata()
    cluster = str(cluster)
    key = _resolve_key(session, cluster, cluster_key)
    labels = a.obs[key].astype(str).to_numpy()
    samples = a.obs[session.sample_key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    X = session.counts()
    if not _looks_like_counts(X):
        raise ValueError("pseudobulk_de needs raw counts; set dataset.counts_layer to the layer holding them")

    if contrast == "cluster_vs_reference_clusters":
        ref = _parse_cluster_list(reference_clusters)
        in_target = labels == cluster
        in_ref = (labels != cluster) if ref == "rest" else np.isin(labels, list(ref))
        group = np.where(in_target, "target", np.where(in_ref, "reference", None))
    elif contrast == "disease_vs_reference_within_compartment":
        # the cluster plus every cluster of the same dominant cell type (its lineage / compartment):
        # disease samples' compartment cells vs reference samples' compartment cells
        in_cp = compartment_mask(session, cluster, key)
        group = np.where(in_cp & (conds == session.disease), "target", np.where(in_cp & np.isin(conds, session.reference), "reference", None))
    else:
        in_cl = labels == cluster
        group = np.where(in_cl & (conds == session.disease), "target", np.where(in_cl & np.isin(conds, session.reference), "reference", None))
    pb_id = np.array([f"{s}|{g}" for s, g in zip(samples, group)], dtype=object)

    keep = np.array([g is not None for g in group])
    pb_ids = pd.Series(pb_id[keep])
    uniq = pb_ids.unique()
    Xk = X[keep]
    mat = np.zeros((len(uniq), a.n_vars))
    ncell = np.zeros(len(uniq), dtype=int)
    for i, u in enumerate(uniq):
        idx = np.where(pb_ids.to_numpy() == u)[0]
        ncell[i] = len(idx)
        mat[i] = np.asarray(Xk[idx].sum(axis=0)).ravel()
    meta = pd.DataFrame({"label": [u.split("|")[1] for u in uniq], "sample": [u.split("|")[0] for u in uniq], "n_cells": ncell}, index=uniq)
    ok = meta["n_cells"] >= min_cells_per_pseudobulk
    n_dropped_small_share = 0
    if contrast == "cluster_vs_reference_clusters":
        # a target replicate must carry a non-trivial share of the state: samples contributing
        # < min_frac_of_cluster of the cluster's cells are spill-over/doublets, and giving them a
        # full vote makes the target group bimodal and inflates dispersion
        n_cluster = int((labels == cluster).sum())
        share_ok = (meta["label"] != "target") | (meta["n_cells"] >= max(min_cells_per_pseudobulk, min_frac_of_cluster * n_cluster))
        n_dropped_small_share = int((ok & ~share_ok & (meta["label"] == "target")).sum())
        ok = ok & share_ok
    meta, mat = meta[ok], mat[ok.to_numpy()]
    n_t, n_r = int((meta["label"] == "target").sum()), int((meta["label"] == "reference").sum())
    if n_t < 2 or n_r < 2:
        return {
            "cluster": cluster,
            "contrast": contrast,
            "error": f"insufficient pseudobulk replicates (target={n_t}, reference={n_r}); need >=2 per side. Lower min_cells_per_pseudobulk or choose another contrast.",
            "replicates": {"target": n_t, "reference": n_r},
        }
    # Blocked design for the sibling contrast: compare a sample's state cells with the SAME
    # sample's sibling cells (~ sample + label). Without blocking, reference pseudobulks from
    # disease samples (which share the disease background) and from reference samples form a
    # bimodal group and the dispersion estimate swallows real effects. Within-cluster
    # disease-vs-reference cannot be blocked (condition is a sample-level variable).
    blocked = False
    if contrast == "cluster_vs_reference_clusters":
        paired = set(meta.loc[meta["label"] == "target", "sample"]) & set(meta.loc[meta["label"] == "reference", "sample"])
        if len(paired) >= 3:
            keep_rows = meta["sample"].isin(paired).to_numpy()
            meta, mat = meta[keep_rows], mat[keep_rows]
            n_t, n_r = int((meta["label"] == "target").sum()), int((meta["label"] == "reference").sum())
            blocked = True
    counts_df = pd.DataFrame(np.round(mat).astype(int), index=meta.index, columns=a.var_names)
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

    method = "pydeseq2 (~ sample + label, paired within animal)" if blocked else "pydeseq2 (~ label)"
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        meta_d = meta[["label", "sample"]].copy()
        meta_d["sample"] = meta_d["sample"].astype(str).str.replace(r"[^A-Za-z0-9_]", "_", regex=True)
        design = "~ sample + label" if blocked else "~ label"
        try:
            dds = DeseqDataSet(counts=counts_df, metadata=meta_d, design=design, quiet=True)
        except TypeError:  # older pydeseq2
            dds = DeseqDataSet(counts=counts_df, metadata=meta_d, design_factors=(["sample", "label"] if blocked else "label"), quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["label", "target", "reference"], quiet=True)
        st.summary()
        res = st.results_df.rename(columns={"log2FoldChange": "log2fc"})
        res = res[["baseMean", "log2fc", "stat", "pvalue", "padj"]]
    except Exception as e:  # fallback
        method = f"wilcoxon_logcpm (pydeseq2 unavailable: {type(e).__name__})"
        cpm = counts_df.div(counts_df.sum(axis=1), axis=0) * 1e6
        lcpm = np.log2(cpm + 1)
        t = lcpm[meta["label"] == "target"]
        r = lcpm[meta["label"] == "reference"]
        pv = np.array([stats.mannwhitneyu(t[g], r[g], alternative="two-sided").pvalue if t[g].std() + r[g].std() > 0 else 1.0 for g in lcpm.columns])
        lfc = t.mean() - r.mean()
        res = pd.DataFrame({"baseMean": counts_df.mean(axis=0), "log2fc": lfc, "stat": lfc * -np.log10(np.clip(pv, 1e-300, 1)), "pvalue": pv, "padj": _bh(pv)}, index=lcpm.columns)

    res.index.name = "gene"
    res = res.sort_values("padj")
    table_path = session.save_table(f"pseudobulk_de_{key}_{cluster}_{contrast}", res)
    # rank by the test statistic (effect size weighted by evidence), not raw log2FC, and require
    # a minimum mean expression so a handful of counts in a few nuclei cannot top the list
    sig = res[(res["padj"] < padj_threshold) & (res["baseMean"] >= min_base_mean)]
    sig_up = sig[sig["log2fc"] > 0].sort_values("stat", ascending=False)
    sig_down = sig[sig["log2fc"] < 0]
    n_low_expr_dropped = int(((res["padj"] < padj_threshold) & (res["log2fc"] > 0) & (res["baseMean"] < min_base_mean)).sum())
    if record:
        session.facts.setdefault("pseudobulk_de", {})[f"{cluster}:{contrast}"] = {
            "method": method,
            "table_path": str(table_path),
            "replicates": {"target": n_t, "reference": n_r},
            "n_sig_up": int(len(sig_up)),
            "n_sig_down": int(len(sig_down)),
            "top_up": sig_up.head(50).reset_index().to_dict(orient="records"),
        }
    return {
        "cluster": cluster,
        "contrast": contrast,
        "method": method,
        "blocked_by_sample": blocked,
        "replicates": {"target": n_t, "reference": n_r, "pseudobulk_cells": {str(session.pseudo.to_pseudo(s)) + "|" + l: int(n) for s, l, n in zip(meta["sample"], meta["label"], meta["n_cells"])}},
        "n_sig_up": int(len(sig_up)),
        "n_sig_down": int(len(sig_down)),
        "n_low_expression_excluded": n_low_expr_dropped,
        "n_target_replicates_dropped_small_share": n_dropped_small_share,
        "ranking": f"top_up is sorted by test statistic among genes with baseMean >= {min_base_mean}",
        "_sig_up": sig_up,
    }


@tool(
    name="disease_programs",
    description="""
    The EXPRESSION route to disease states, for when composition does not change: for each
    cluster (largest first, up to max_clusters, >= min_cells), run replicate-aware pseudobulk DE
    of disease samples vs reference samples WITHIN that cluster. A cluster with >= min_sig_genes
    significantly up-regulated genes is flagged 'disease_altered' - the same cell type is present
    in all conditions but runs a different program in disease. Returns per-cluster counts and the
    top genes; results are stored so pseudobulk_de/gene_specificity/score_targets can use them.
    Includes an ambient-RNA check: programs that are mostly the dominant disease state's
    signature, or that recur across unrelated cell types, are marked ambient_suspect instead of
    disease_altered, and a clean_program_genes list (program minus suspect genes) is reported.
    Use after differential_abundance finds nothing, or alongside it to catch program shifts.""",
    input_schema={
        "properties": {
            "cluster_key": {"type": ["string", "null"]},
            "max_clusters": {"type": "integer", "default": 8},
            "min_cells": {"type": "integer", "default": 200},
            "min_cells_per_pseudobulk": {"type": "integer", "default": 10},
            "padj_threshold": {"type": "number", "default": 0.05},
            "min_log2fc": {"type": "number", "default": 1.0, "description": "genes must exceed this to count towards min_sig_genes"},
            "min_sig_genes": {"type": "integer", "default": 30, "description": "absolute criterion"},
            "small_program_min": {"type": "integer", "default": 5, "description": "relative criterion: a cluster with at least this many strong genes AND at least 3x more than any other cluster is also flagged"},
            "min_base_mean": {"type": "number", "default": 10.0},
        },
        "required": [],
    },
    roles=ALL,
)
def disease_programs(
    session: AnalysisSession,
    cluster_key: str | None = None,
    max_clusters: int = 8,
    min_cells: int = 200,
    min_cells_per_pseudobulk: int = 10,
    padj_threshold: float = 0.05,
    min_log2fc: float = 1.0,
    min_sig_genes: int = 30,
    small_program_min: int = 5,
    min_base_mean: float = 10.0,
) -> dict[str, Any]:
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    sizes = a.obs[key].astype(str).value_counts()
    # clusters already enriched in disease are composition states: "disease vs reference cells
    # within" them is disease cells vs spill-over, not a program test - skip and say so
    enriched = sorted(set(map(str, session.facts.get("differential_abundance", {}).get("disease_enriched_clusters", []))))
    skipped = [c for c in sizes.index if c in enriched]
    _ensure_lognorm(session)
    labels_all = a.obs[key].astype(str).to_numpy()
    Xlog = a.X
    clusters = [c for c, n in sizes.items() if n >= min_cells and c not in enriched][:max_clusters]
    ck = _lineage_key(session)
    dom = {}
    if ck and ck in a.obs:
        dom = pd.crosstab(a.obs[key].astype(str), a.obs[ck].astype(str)).idxmax(axis=1).to_dict()
    rows, programs = [], {}
    for cl in clusters:
        r = _run_pseudobulk(session, cl, key, "disease_vs_reference_within_cluster", None, min_cells_per_pseudobulk, padj_threshold, min_base_mean, record=True)
        if "error" in r:
            rows.append({"cluster": cl, "n_cells": int(sizes[cl]), "dominant_celltype": dom.get(cl), "error": r["error"], "disease_altered": False})
            continue
        up = r["_sig_up"]
        strong = up[up["log2fc"] >= min_log2fc]
        rows.append(
            {
                "cluster": cl,
                "n_cells": int(sizes[cl]),
                "dominant_celltype": dom.get(cl),
                "replicates_disease": r["replicates"]["target"],
                "replicates_reference": r["replicates"]["reference"],
                "n_sig_up": r["n_sig_up"],
                "n_sig_up_strong": int(len(strong)),
                "n_sig_down": r["n_sig_down"],
                "top_up_genes": list(strong.index[:12]) if len(strong) else list(up.index[:12]),
                "disease_altered": bool(len(strong) >= min_sig_genes),
            }
        )
        programs[cl] = {"n_sig_up_strong": int(len(strong)), "top_up": strong.head(50).reset_index().to_dict(orient="records")}
    df = pd.DataFrame(rows).set_index("cluster")
    ambient_note = None
    if "n_sig_up_strong" in df:
        # relative criterion: a modest program that clearly stands out from every other cluster
        counts_ = df["n_sig_up_strong"].fillna(0)
        for cl in df.index:
            others = counts_.drop(cl)
            if counts_[cl] >= small_program_min and counts_[cl] >= 3 * (others.max() if len(others) else 0):
                df.loc[cl, "disease_altered"] = True

        # ---- ambient-RNA check -------------------------------------------------------------
        # When a disease-enriched population dominates the disease samples, its transcripts sit as
        # background in EVERY nucleus of those samples, so "disease vs reference within cluster"
        # finds that population's program in immune, stromal and other epithelial clusters alike.
        # A gene is suspect if it (a) is a top gene of a composition-enriched state, or (b) recurs
        # in the programs of clusters belonging to >= 2 other broad cell types.
        strong_sets = {cl: set(r["gene"] for r in programs[cl]["top_up"]) if cl in programs else set() for cl in df.index}
        dom_signature: set[str] = set()
        for e in enriched:
            m = labels_all == e
            if m.sum() >= 50:
                mean_e = np.asarray(Xlog[m].mean(axis=0)).ravel()
                dom_signature |= set(a.var_names[np.argsort(-mean_e)[:300]])
        celltype_of = {cl: dom.get(cl, f"cluster_{cl}") for cl in df.index}
        gene_celltypes: dict[str, set[str]] = {}
        for cl, gs in strong_sets.items():
            for g in gs:
                gene_celltypes.setdefault(g, set()).add(celltype_of[cl])
        clean_programs: dict[str, list[str]] = {}
        n_susp, frac_susp, clean_cols, n_clean, amb_flag = {}, {}, {}, {}, {}
        for cl in df.index:
            gs = strong_sets[cl]
            amb_flag[cl] = False
            if not gs:
                continue
            # primary criterion: within-sample lineage control (ambient ratio) on the strong genes;
            # secondary: dominant-state signature / recurrence across unrelated cell types
            glist = sorted(gs)[:200]
            idx = [list(a.var_names).index(g) for g in glist]
            Xs = Xlog[:, idx]
            Xs = Xs.toarray() if sp.issparse(Xs) else np.asarray(Xs)
            ratios = _ambient_ratios(session, np.expm1(Xs), Xs > 0, glist, labels_all == cl, disease_only=True, min_cells=min_cells_per_pseudobulk)
            amb_like = {g for g in glist if _is_ambient_like(ratios.loc[g])}
            suspect = amb_like | {g for g in gs if g in dom_signature or len(gene_celltypes.get(g, set()) - {celltype_of[cl]}) >= 2}
            frac = len(suspect) / len(gs)
            clean = [r["gene"] for r in programs[cl]["top_up"] if r["gene"] not in suspect]
            clean_programs[cl] = clean
            n_susp[cl], frac_susp[cl], clean_cols[cl], n_clean[cl] = len(suspect), round(frac, 3), ", ".join(clean[:12]), len(clean)
            # suspect whenever a non-trivial program is mostly shared/dominant-signature genes,
            # independent of whether it passed the size criteria - the recurrence itself is the evidence
            if len(gs) >= small_program_min and frac >= 0.5:
                amb_flag[cl] = True
                df.loc[cl, "disease_altered"] = False
            elif bool(df.loc[cl, "disease_altered"]) and len(clean) < small_program_min:
                df.loc[cl, "disease_altered"] = False
                amb_flag[cl] = True
        df["n_strong_suspect_ambient"] = pd.Series(n_susp)
        df["frac_suspect_ambient"] = pd.Series(frac_susp)
        df["n_clean_program_genes"] = pd.Series(n_clean)
        df["clean_program_genes"] = pd.Series(clean_cols)
        df["ambient_suspect"] = pd.Series(amb_flag).astype(bool)
        if True:
            n_amb = int(df["ambient_suspect"].sum())
            amb_types = {celltype_of[c] for c in df.index[df["ambient_suspect"]]}
            if n_amb >= 2 and len(amb_types) >= 2:
                ambient_note = (
                    f"{n_amb} cluster(s) carry a 'disease program' that is mostly the dominant disease state's signature and/or recurs across unrelated cell types - "
                    "this is the fingerprint of ambient RNA from the disease-enriched population(s) in disease samples, not cell-intrinsic biology. "
                    "Only the clean_program_genes (program minus suspect genes) should be interpreted for those clusters."
                )
        df = df.sort_values("n_sig_up_strong", ascending=False)
    session.save_table(f"disease_programs_{key}", df)
    altered = [c for c in df.index if bool(df.loc[c, "disease_altered"])]
    suspects = [c for c in df.index if "ambient_suspect" in df and bool(df.loc[c, "ambient_suspect"])]
    # for altered clusters, downstream nomination should use the clean program
    for cl in altered:
        if cl in clean_programs and f"{cl}:disease_vs_reference_within_cluster" in session.facts.get("pseudobulk_de", {}):
            pb = session.facts["pseudobulk_de"][f"{cl}:disease_vs_reference_within_cluster"]
            pb["top_up"] = [r for r in pb["top_up"] if r["gene"] in set(clean_programs[cl])] or pb["top_up"]
    session.facts["disease_programs"] = {"cluster_key": key, "disease_altered_clusters": altered, "ambient_suspect_clusters": suspects, "table": df.reset_index().to_dict(orient="records"), "ambient_note": ambient_note}
    return {
        "cluster_key": key,
        "contrast": "disease_vs_reference_within_cluster (pseudobulk per sample)",
        "skipped_composition_enriched_clusters": skipped,
        "disease_altered_clusters": altered,
        "ambient_suspect_clusters": suspects,
        "ambient_warning": ambient_note,
        "table": df,
        "note": "disease_altered = same cell type in all conditions but >= min_sig_genes genes up (padj<thr, log2FC>=min_log2fc) in disease samples AND the program is not mostly ambient-suspect. Characterise these with gene_specificity and nominate from their clean_program_genes; specificity should be judged by the within-cluster disease fold change, not cluster-vs-cluster.",
    }


def _ambient_ratios(
    session: AnalysisSession,
    X_lin: np.ndarray,
    expr: np.ndarray,
    gene_names: list[str],
    interest: np.ndarray,
    disease_only: bool = True,
    min_cells: int = 20,
) -> pd.DataFrame:
    """Within-sample lineage control for ambient RNA.

    For each gene: per disease sample, mean expression in the cells of interest vs in every other
    broad lineage; the control lineage is the one with the lowest median expression (least related
    to the gene). ambient_ratio = median over samples of interest/control. Ambient RNA reads the
    same background in both (ratio ~1); a cell-intrinsic program is several-fold higher.
    """
    a = session.require_adata()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    samples = a.obs[session.sample_key].astype(str).to_numpy()
    ck = _lineage_key(session)
    if ck and ck in a.obs:
        lineage = a.obs[ck].astype(str).to_numpy()
    else:  # no broad annotation: use clusters as lineages
        lineage = a.obs[session.cluster_key].astype(str).to_numpy() if session.cluster_key else np.array(["all"] * a.n_obs)
    scope = (conds == session.disease) if disease_only else np.ones(a.n_obs, dtype=bool)
    ref_scope = np.isin(conds, session.reference)
    dom = pd.Series(lineage[interest]).value_counts().idxmax() if interest.any() else None
    others = [l for l in np.unique(lineage[scope]) if l != dom]
    smps = np.unique(samples[scope & interest])
    ref_smps = np.unique(samples[ref_scope])
    eps = 1e-3
    out = []
    for j, g in enumerate(gene_names):
        ratios, pct_i, pct_c, ctrl_name = [], [], [], None
        # pick the control lineage for this gene: lowest median mean-expression across samples
        med = {}
        for l in others:
            vals = []
            for smp in smps:
                m = scope & (lineage == l) & (samples == smp)
                if m.sum() >= min_cells:
                    vals.append(float(X_lin[m, j].mean()))
            if len(vals) >= max(2, len(smps) // 2):
                med[l] = float(np.median(vals))
        if not med:
            out.append({"gene": g, "ambient_ratio": None, "control_lineage": None, "pct_interest": None, "pct_control": None, "control_disease_vs_reference": None})
            continue
        ctrl_name = min(med, key=med.get)
        # is the control lineage itself elevated in disease samples vs reference samples? Ambient
        # background must raise it; a gene that is flat in the control lineage cannot be ambient.
        ctrl_ref = [float(X_lin[ref_scope & (lineage == ctrl_name) & (samples == smp), j].mean()) for smp in ref_smps if (ref_scope & (lineage == ctrl_name) & (samples == smp)).sum() >= min_cells]
        ctrl_dvr = ((med[ctrl_name] + eps) / (float(np.median(ctrl_ref)) + eps)) if ctrl_ref else None
        for smp in smps:
            mi = scope & interest & (samples == smp)
            mc = scope & (lineage == ctrl_name) & (samples == smp)
            if mi.sum() >= min_cells and mc.sum() >= min_cells:
                ratios.append((float(X_lin[mi, j].mean()) + eps) / (float(X_lin[mc, j].mean()) + eps))
                pct_i.append(float(expr[mi, j].mean()))
                pct_c.append(float(expr[mc, j].mean()))
        out.append(
            {
                "gene": g,
                "ambient_ratio": round(float(np.median(ratios)), 2) if ratios else None,
                "control_lineage": ctrl_name,
                "pct_interest": round(float(np.median(pct_i)), 3) if pct_i else None,
                "pct_control": round(float(np.median(pct_c)), 3) if pct_c else None,
                "control_disease_vs_reference": None if ctrl_dvr is None else round(ctrl_dvr, 2),
                "n_samples_evaluated": len(ratios),
            }
        )
    return pd.DataFrame(out).set_index("gene")


def _is_ambient_like(row: pd.Series, ratio_threshold: float = 3.0, control_elevation: float = 1.5) -> bool:
    r, c = row.get("ambient_ratio"), row.get("control_disease_vs_reference")
    if r is None or pd.isna(r):
        return False
    if c is None or pd.isna(c):
        return bool(r < ratio_threshold)  # no reference samples for the control lineage: fall back to the ratio alone
    return bool(r < ratio_threshold and c >= control_elevation)


@tool(
    name="condition_attribution",
    description="""
    For candidate genes in one cluster, ask WHERE the expression actually sits. Within the
    cluster, compares disease-condition cells with reference-condition cells (per-sample medians
    of % expressing and mean expression); then compares with disease and reference cells in the
    sibling clusters of the same cell type. Classifies each gene:
      state_specific    high in disease cells of this cluster, low in reference cells here, low in
                        disease cells of sibling clusters
      disease_wide      also high in disease cells of the lineage's SHARED (condition-mixed) sibling
                        clusters - other disease-enriched states are excluded from this test and
                        reported separately in pct_disease_in_other_enriched_states (a disease marker, not a
                        state marker - still a legitimate target)
      shared_with_reference  also high in reference-condition cells of the same cell type (here or in
                        sibling clusters): a cell-type marker or a pan-condition (e.g. pan-tumour) gene
      reference_driven  the signal comes from the REFERENCE-condition cells inside this cluster
                        (spill-over / misassigned cells) - a false positive of cluster-vs-cluster DE
      ambient_like      within the disease samples the gene is no more than ~3x higher in these cells
                        than in the least-related lineage (immune / stromal / epithelial) of the SAME
                        samples, AND that control lineage is itself raised in disease vs reference
                        samples: the signal is sample background (ambient RNA), not cell-intrinsic
      weak              expressed in <10% of the cluster's disease cells
    The ambient_ratio column (interest / control lineage, median over disease samples) and the
    control lineage are reported for every gene so the judgement is auditable. Works symmetrically
    for tumour-cell and TME (immune, stromal) targets: the control is whichever lineage is least
    related to the gene. Run this on every candidate before nomination; score_targets penalises
    reference_driven, ambient_like and shared_with_reference genes.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
            "cluster": {"type": "string"},
            "cluster_key": {"type": ["string", "null"]},
            "min_cells": {"type": "integer", "default": 20, "description": "a condition subset smaller than this is reported as not evaluable"},
            "ambient_ratio_threshold": {"type": "number", "default": 3.0, "description": "below this interest/control-lineage ratio a gene is classed ambient_like"},
            "min_reference_frac": {"type": ["number", "null"], "default": None, "description": "reference cells inside the cluster must be at least this fraction of it for the in-cluster reference comparison to decide the class (default from config analysis.attribution.min_reference_frac_in_cluster, 0.02); below it the pool is reported as too thin and the sibling comparison carries the judgement"},
        },
        "required": ["genes", "cluster"],
    },
    roles=ALL,
)
def condition_attribution(session: AnalysisSession, genes: list[str], cluster: str, cluster_key: str | None = None, min_cells: int = 20, ambient_ratio_threshold: float = 3.0, min_reference_frac: float | None = None) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    cluster = str(cluster)
    key = _resolve_key(session, cluster, cluster_key)
    labels = a.obs[key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    samples = a.obs[session.sample_key].astype(str).to_numpy()
    is_d = conds == session.disease
    is_r = np.isin(conds, session.reference)
    in_c = labels == cluster
    ck = _lineage_key(session)
    enriched = set(map(str, session.facts.get("differential_abundance", {}).get("disease_enriched_clusters", []))) - {cluster}
    in_enriched_other = np.isin(labels, list(enriched)) if enriched else np.zeros(a.n_obs, dtype=bool)
    if ck and ck in a.obs:
        ct = a.obs[ck].astype(str).to_numpy()
        dom = pd.Series(ct[in_c]).value_counts().idxmax()
        sib = (~in_c) & (ct == dom)
        sib_desc = f"other clusters of {dom}"
    else:
        sib = ~in_c
        sib_desc = "all other clusters"
    # 'disease_wide' is judged on the lineage's SHARED (condition-mixed) clusters, excluding other
    # disease-enriched states, so the label means the same thing as score_targets' target_class:
    # a gene also up in disease cells of ordinary clusters, not merely present in a second disease state
    sib_shared = sib & ~in_enriched_other
    if enriched:
        sib_desc += f" (disease_wide judged on shared clusters, excluding enriched {sorted(enriched)})"
    # Thin reference pool guard: the in-cluster reference comparison (reference_driven /
    # shared_with_reference via reference cells INSIDE the cluster) only decides the class when
    # reference cells are a non-trivial share of the cluster. A disease state always attracts a few
    # reference cells that resemble it (doublets, misassigned nuclei, a rare real population), and
    # letting a 0.5% pool overrule the other 99.5% is the wrong weighting. Below the floor the pool
    # is reported, not consulted, and the sibling-cluster comparison carries the judgement.
    if min_reference_frac is None:
        min_reference_frac = float(session.cfg.analysis.get("attribution", {}).get("min_reference_frac_in_cluster", 0.02))
    n_in_c = int(in_c.sum())
    n_ref_in_c = int((in_c & is_r).sum())
    ref_frac_in_c = n_ref_in_c / max(1, n_in_c)
    ref_pool_ok = ref_frac_in_c >= min_reference_frac
    ref_pool_status = "ok" if ref_pool_ok else f"too_thin ({n_ref_in_c} cells, {100*ref_frac_in_c:.2f}% of cluster < {100*min_reference_frac:g}%): in-cluster reference comparison not used for the class"
    present = [g for g in genes if g in a.var_names]
    missing = [g for g in genes if g not in a.var_names]
    if not present:
        return {"error": "none of the genes are in the dataset", "missing": missing}
    X = a[:, present].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    expr = X > 0
    X = np.expm1(X)  # linear-scale means so that ratios between subsets are meaningful

    def per_sample_median(mask: np.ndarray, j: int) -> tuple[float | None, float | None, int]:
        """median over samples of (% expressing, mean normalised expr), requiring >= min_cells per sample."""
        if mask.sum() < min_cells:
            return None, None, int(mask.sum())
        pcts, means = [], []
        for smp in np.unique(samples[mask]):
            m = mask & (samples == smp)
            if m.sum() >= min_cells:
                pcts.append(float(expr[m, j].mean()))
                means.append(float(X[m, j].mean()))
        if not pcts:  # subset is big enough overall but spread thinly across samples: pool
            return float(expr[mask, j].mean()), float(X[mask, j].mean()), int(mask.sum())
        return float(np.median(pcts)), float(np.median(means)), int(mask.sum())

    amb = _ambient_ratios(session, X, expr, present, in_c, disease_only=True, min_cells=min_cells)
    rows = []
    for j, g in enumerate(present):
        a_pct, a_mean, n_a = per_sample_median(in_c & is_d, j)       # disease cells, this cluster
        b_pct, b_mean, n_b = per_sample_median(in_c & is_r, j)       # reference cells, this cluster
        d_pct, d_mean, n_d = per_sample_median(sib_shared & is_d, j)  # disease cells, SHARED sibling clusters
        r_pct, r_mean, n_r = per_sample_median(sib & is_r, j)        # reference cells, sibling clusters
        e_pct, e_mean, n_e = per_sample_median(sib & is_d & in_enriched_other, j) if enriched else (None, None, 0)  # disease cells in other enriched states
        # classify on mean expression ratios (detection % alone saturates for broadly expressed genes)
        eps = 1e-3
        if a_pct is None or a_mean is None:
            cls = "not_evaluable"
        elif ref_pool_ok and b_mean is not None and b_pct is not None and b_pct > 0.2 and b_mean >= 2 * (a_mean + eps):
            cls = "reference_driven"  # checked before 'weak': weak-in-disease but strong-in-reference IS the artefact
        elif a_pct < 0.10:
            cls = "weak"
        elif _is_ambient_like(amb.loc[g], ambient_ratio_threshold):
            cls = "ambient_like"
        elif (r_mean is not None and r_mean >= 0.5 * (a_mean + eps)) or (ref_pool_ok and b_mean is not None and b_mean >= 0.5 * (a_mean + eps)):
            cls = "shared_with_reference"
        elif d_mean is not None and d_mean >= 0.5 * (a_mean + eps):
            cls = "disease_wide"
        else:
            cls = "state_specific"
        rows.append(
            {
                "gene": g,
                "class": cls,
                "pct_disease_in_cluster": None if a_pct is None else round(a_pct, 3),
                "pct_reference_in_cluster": None if b_pct is None else round(b_pct, 3),
                "pct_disease_in_siblings": None if d_pct is None else round(d_pct, 3),
                "pct_disease_in_other_enriched_states": None if e_pct is None else round(e_pct, 3),
                "pct_reference_in_siblings": None if r_pct is None else round(r_pct, 3),
                "mean_disease_in_cluster": None if a_mean is None else round(a_mean, 3),
                "mean_reference_in_cluster": None if b_mean is None else round(b_mean, 3),
                "reference_pool_in_cluster": ref_pool_status,
                "ambient_ratio": amb.loc[g, "ambient_ratio"],
                "ambient_control_lineage": amb.loc[g, "control_lineage"],
                "pct_in_control_lineage": amb.loc[g, "pct_control"],
                "control_lineage_disease_vs_reference": amb.loc[g, "control_disease_vs_reference"],
                "n_cells": {"disease_in_cluster": n_a, "reference_in_cluster": n_b, "disease_in_siblings": n_d, "reference_in_siblings": n_r},
            }
        )
    df = pd.DataFrame(rows).set_index("gene")
    session.save_table(f"condition_attribution_{key}_{cluster}", df.drop(columns=["n_cells"]))
    session.facts.setdefault("condition_attribution", {}).setdefault(cluster, {}).update({r["gene"]: r["class"] for r in rows})
    counts = df["class"].value_counts().to_dict()
    return {
        "cluster": cluster,
        "cluster_key": key,
        "siblings": sib_desc,
        "reference_pool_in_cluster": ref_pool_status,
        "class_counts": counts,
        "table": df,
        "missing_genes": missing,
        "note": "Values are medians over samples (each sample needs >= min_cells cells in the subset; thin subsets are pooled). Classes use mean-expression ratios: reference_driven = reference cells in cluster >= 2x disease cells; shared_with_reference = reference cells (here or in siblings) >= 0.5x; disease_wide = disease cells in the SHARED sibling clusters (other disease-enriched states excluded; see pct_disease_in_other_enriched_states for those) >= 0.5x; ambient_like = within disease samples the gene is < threshold-fold higher in these cells than in the least-related lineage of the same samples AND that control lineage is itself elevated (>=1.5x) in disease vs reference samples - the two signatures of sample background rather than cell-intrinsic expression. reference_driven genes are false positives of cluster-vs-cluster DE: their expression sits in the reference-condition cells that were assigned to this cluster. The in-cluster reference comparison is only used when reference cells are >= min_reference_frac of the cluster (reference_pool_in_cluster column); a thinner pool is reported but cannot decide the class.",
    }


@tool(
    name="gene_specificity",
    description="""
    For a list of candidate genes, quantify how specific their expression is to one cluster:
    mean expression in the target cluster vs the maximum mean in any other NON-disease cluster
    (log2 specificity ratio; clusters that are flagged disease states OR are >= 50% disease cells
    are excluded, so a gene shared by two disease populations is not penalised - the all-clusters
    value is also reported), % of target cells expressing, expression in disease vs
    reference cells of the target cluster, and sample consistency (fraction of disease
    samples in which the gene is detected in >=10% of the cluster's cells). Also reports
    whether the gene is up in other condition-shared clusters (a pan-tumour or cell-type
    marker rather than a disease-state marker). Use before nominating targets.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
            "cluster": {"type": "string"},
            "cluster_key": {"type": ["string", "null"]},
        },
        "required": ["genes", "cluster"],
    },
    roles=ALL,
)
def gene_specificity(session: AnalysisSession, genes: list[str], cluster: str, cluster_key: str | None = None) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    cluster = str(cluster)
    key = _resolve_key(session, cluster, cluster_key)
    labels = a.obs[key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    samples = a.obs[session.sample_key].astype(str).to_numpy()
    present = [g for g in genes if g in a.var_names]
    missing = [g for g in genes if g not in a.var_names]
    if not present:
        return {"error": "none of the genes are in the dataset", "missing": missing}
    X = a[:, present].X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    expm = np.expm1(X)
    clusters = sorted(set(labels))
    t_mask = labels == cluster
    means = {cl: expm[labels == cl].mean(axis=0) for cl in clusters}
    # other disease states (composition-enriched or program-altered) are not "other clusters" for
    # the purpose of disease specificity: a gene shared by two disease states is still disease-specific
    disease_clusters = set(map(str, session.facts.get("differential_abundance", {}).get("disease_enriched_clusters", []))) | set(map(str, session.facts.get("disease_programs", {}).get("disease_altered_clusters", [])))
    # ... and any cluster whose cells are mostly from the disease condition, whether or not it
    # passed the enrichment test (a 70%-disease cluster that missed alpha is not a "normal" cluster)
    is_dis_cell = conds == session.disease
    for cl in clusters:
        m = labels == cl
        if m.sum() and float(is_dis_cell[m].mean()) >= 0.5:
            disease_clusters.add(cl)
    disease_clusters.discard(cluster)
    rows = []
    dis_samples = sorted(set(samples[conds == session.disease]))
    for j, g in enumerate(present):
        tm = float(means[cluster][j])
        others_all = {cl: float(means[cl][j]) for cl in clusters if cl != cluster}
        others = {cl: v for cl, v in others_all.items() if cl not in disease_clusters} or others_all
        max_other_cl = max(others, key=others.get) if others else None
        max_other = others[max_other_cl] if others else 0.0
        spec = float(np.log2((tm + 0.01) / (max_other + 0.01)))
        max_all_cl = max(others_all, key=others_all.get) if others_all else None
        spec_all = float(np.log2((tm + 0.01) / ((others_all[max_all_cl] if max_all_cl else 0.0) + 0.01)))
        pct_in = float((X[t_mask, j] > 0).mean())
        d_in = expm[t_mask & (conds == session.disease), j]
        r_in = expm[t_mask & np.isin(conds, session.reference), j]
        dvr = float(np.log2((d_in.mean() + 0.01) / (r_in.mean() + 0.01))) if len(d_in) and len(r_in) else None
        cons = []
        for s in dis_samples:
            m = t_mask & (samples == s)
            if m.sum() >= 5:
                cons.append(float((X[m, j] > 0).mean()) >= 0.10)
        consistency = float(np.mean(cons)) if cons else None
        rows.append(
            {
                "gene": g,
                "mean_in_cluster": round(tm, 3),
                "max_mean_other_cluster": round(max_other, 3),
                "max_other_cluster": max_other_cl,
                "comparison_set": "clusters < 50% disease cells and not flagged as disease states",
                "log2_specificity": round(spec, 3),
                "log2_specificity_incl_other_disease_states": round(spec_all, 3),
                "max_other_cluster_incl_disease_states": max_all_cl,
                "pct_expressing_in_cluster": round(pct_in, 3),
                "log2_disease_vs_reference_within_cluster": None if dvr is None else round(dvr, 3),
                "disease_sample_consistency": None if consistency is None else round(consistency, 3),
                "n_disease_samples_evaluated": len(cons),
            }
        )
    df = pd.DataFrame(rows).set_index("gene").sort_values("log2_specificity", ascending=False)
    session.save_table(f"gene_specificity_{key}_{cluster}", df)
    session.facts.setdefault("gene_specificity", {})[cluster] = df.reset_index().to_dict(orient="records")
    return {"cluster": cluster, "cluster_key": key, "table": df, "missing_genes": missing}


@tool(
    name="enrichment",
    description="""
    Over-representation analysis of a gene list against gene sets in a .gmt file
    (default: the configured MSigDB Hallmark-style file) with a hypergeometric test and
    BH correction. Background = all genes in the dataset. Returns top terms with overlap genes.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
            "gmt_path": {"type": ["string", "null"], "description": "override the configured gene-set file"},
            "n_top": {"type": "integer", "default": 15},
        },
        "required": ["genes"],
    },
    roles=ALL,
)
def enrichment(session: AnalysisSession, genes: list[str], gmt_path: str | None = None, n_top: int = 15) -> dict[str, Any]:
    a = session.require_adata()
    path = session.cfg.resolve_path(gmt_path or session.cfg.get("analysis.enrichment.gene_sets"))
    sets: dict[str, set[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = set(parts[2:])
    background = set(a.var_names)
    query = set(g for g in genes if g in background)
    N, n = len(background), len(query)
    rows = []
    for term, members in sets.items():
        members_bg = members & background
        K = len(members_bg)
        if K == 0:
            continue
        overlap = sorted(query & members_bg)
        k = len(overlap)
        if k == 0:
            continue
        p = float(stats.hypergeom.sf(k - 1, N, K, n))
        rows.append({"term": term, "overlap": k, "set_size_in_background": K, "p": p, "genes": overlap})
    if not rows:
        return {"n_query_genes_in_background": n, "terms": [], "note": "no overlapping terms"}
    df = pd.DataFrame(rows)
    df["padj"] = _bh(df["p"].to_numpy())
    df = df.sort_values("p").set_index("term")
    session.save_table("enrichment", df)
    session.facts.setdefault("enrichment", []).append({"query_size": n, "top": df.head(10).reset_index().to_dict(orient="records")})
    return {"n_query_genes_in_background": n, "gene_set_file": str(path.name), "terms": df.head(n_top)}


@tool(
    name="plot_genes",
    description="""
    Save a figure for a list of genes: UMAP feature plots (if an embedding exists) and a
    dot plot across clusters split by condition. Returns the figure path. Use to sanity-check
    candidate genes visually; the figure goes into the report.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}, "description": "up to 8 genes"},
            "cluster_key": {"type": ["string", "null"]},
            "title": {"type": ["string", "null"]},
        },
        "required": ["genes"],
    },
    roles=ALL,
)
def plot_genes(session: AnalysisSession, genes: list[str], cluster_key: str | None = None, title: str | None = None) -> dict[str, Any]:
    a = session.require_adata()
    _ensure_lognorm(session)
    key = cluster_key or session.cluster_key
    genes = [g for g in genes if g in a.var_names][:8]
    if not genes:
        return {"error": "no requested genes present in dataset"}
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = {}
    if "X_umap" in a.obsm:
        ncol = min(4, len(genes))
        nrow = int(np.ceil(len(genes) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3 * nrow), squeeze=False)
        for ax, g in zip(axes.ravel(), genes):
            sc.pl.umap(a, color=g, ax=ax, show=False, frameon=False, colorbar_loc=None, title=g, cmap="viridis")
        for ax in axes.ravel()[len(genes):]:
            ax.axis("off")
        fig.suptitle(title or "Candidate genes (UMAP)", y=1.02)
        p = session.save_figure("umap_genes_" + "_".join(genes[:3]), fig, f"UMAP feature plots: {', '.join(genes)}.")
        out["umap_figure"] = str(p.relative_to(session.run_dir))
    if key and key in a.obs:
        a.obs["_cl_cond"] = a.obs[key].astype(str) + " | " + a.obs[session.condition_key].astype(str)
        dp = sc.pl.dotplot(a, var_names=genes, groupby="_cl_cond", show=False, return_fig=True)
        if hasattr(dp, "make_figure"):
            dp.make_figure()
        fig = getattr(dp, "fig", None) or plt.gcf()
        p = session.save_figure("dotplot_" + "_".join(genes[:3]), fig, f"Dot plot of {', '.join(genes)} by cluster | condition.")
        out["dotplot_figure"] = str(p.relative_to(session.run_dir))
        del a.obs["_cl_cond"]
    return out
