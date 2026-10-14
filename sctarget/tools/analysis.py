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

    existing_clusterings = [c for c in a.obs.columns if c.lower() in ("leiden", "louvain", "clusters", "cluster", "seurat_clusters")]
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
        "celltype_counts": {str(k): int(v) for k, v in celltypes.items()},
        "existing_clusterings": existing_clusterings,
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
    Standard cell/gene quality control on raw counts: filters cells by number of detected
    genes and mitochondrial fraction, and genes by minimum number of cells. Returns before/after
    counts overall and per sample. Skip if the dataset is already QC'd (already_processed_flag).""",
    input_schema={
        "properties": {
            "min_genes": {"type": "integer", "description": "min detected genes per cell", "default": 200},
            "max_genes": {"type": "integer", "description": "max detected genes per cell (doublet guard)", "default": 8000},
            "max_mt_pct": {"type": "number", "description": "max % mitochondrial counts", "default": 10.0},
            "min_cells_per_gene": {"type": "integer", "default": 3},
        },
        "required": [],
    },
    roles=ANALYST,
    mutates_state=True,
)
def run_qc(session: AnalysisSession, min_genes: int = 200, max_genes: int = 8000, max_mt_pct: float = 10.0, min_cells_per_gene: int = 3) -> dict[str, Any]:
    a = session.require_adata()
    layer = session.cfg.dataset.get("counts_layer")
    layer = layer if layer and layer in a.layers else None
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, layer=layer, percent_top=None)
    before = a.n_obs
    before_per_sample = a.obs.groupby(session.sample_key, observed=True).size()
    keep = (a.obs["n_genes_by_counts"] >= min_genes) & (a.obs["n_genes_by_counts"] <= max_genes) & (a.obs["pct_counts_mt"] <= max_mt_pct)
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
        "params": {"min_genes": min_genes, "max_genes": max_genes, "max_mt_pct": max_mt_pct, "min_cells_per_gene": min_cells_per_gene},
    }
    session.facts["qc"] = {k: v for k, v in out.items() if k != "per_sample"}
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
                import harmonypy  # noqa: F401

                sc.external.pp.harmony_integrate(a, key=session.sample_key, basis="X_pca", adjusted_basis="X_pca_harmony", random_state=seed)
                rep = "X_pca_harmony"
                harmony_used = True
            except Exception as e:  # pragma: no cover - depends on optional dep
                session.log.note("analyst", f"harmony unavailable or failed ({e}); using uncorrected PCA")
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
    ck = session.celltype_key
    if ck and ck in a.obs:
        dom = pd.crosstab(a.obs[key].astype(str), a.obs[ck].astype(str))
        df["dominant_celltype"] = dom.idxmax(axis=1)
        df["dominant_celltype_frac"] = (dom.max(axis=1) / dom.sum(axis=1)).round(3)
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
) -> dict[str, Any]:
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available; run integrate_and_cluster first")
    obs = a.obs[[key, session.sample_key, session.condition_key]].copy()
    obs[key] = obs[key].astype(str)
    ck = session.celltype_key
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
    df["disease_enriched"] = (df["log2_ratio"] >= enrichment_log2fc) & (df["n_disease_samples_present"] >= min(2, int(is_dis.sum()))) & (df["n_cells"] >= min_cluster_cells)
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
    facts_key = "differential_abundance_within_celltype" if within_celltype else "differential_abundance"
    session.facts[facts_key] = {
        "cluster_key": key,
        "within_celltype": bool(within_celltype),
        "disease_enriched_clusters": enriched,
        "table": df.drop(columns=["per_sample_prop_disease", "per_sample_prop_reference"]).reset_index().to_dict(orient="records"),
    }
    return {
        "cluster_key": key,
        "within_celltype": bool(within_celltype),
        "n_samples_used": {"disease": int(is_dis.sum()), "reference": int(is_ref.sum())},
        "disease_enriched_clusters": enriched,
        "table": df,
        "figure": str(p.relative_to(session.run_dir)),
        "note": f"Small n per group: treat p-values as descriptive; prioritize log2_ratio and n_disease_samples_present. Clusters with < {min_cluster_cells} cells are never flagged enriched.",
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
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available")
    cluster = str(cluster)
    labels = a.obs[key].astype(str)
    if cluster not in set(labels):
        raise ValueError(f"cluster '{cluster}' not in {sorted(set(labels))}")
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
    up = df[(df["log2fc"] > 0)].head(n_top)
    session.facts.setdefault("markers", {})[f"{cluster}_vs_{'rest' if ref == 'rest' else ','.join(ref)}"] = up.head(50).reset_index().to_dict(orient="records")
    return {
        "cluster": cluster,
        "cluster_key": key,
        "reference": "rest" if ref == "rest" else ref,
        "n_cells_target": int((labels == cluster).sum()),
        "n_cells_reference": int(sub.n_obs - (labels == cluster).sum()),
        "top_markers": up,
    }


@tool(
    name="pseudobulk_de",
    description="""
    Replicate-aware differential expression: counts are summed per sample within each
    group (pseudobulk) and tested with DESeq2 (PyDESeq2; falls back to Wilcoxon on
    log-CPM if unavailable). Two contrasts:
      - 'cluster_vs_reference_clusters': cells of `cluster` vs cells of `reference_clusters`
        (default: all other clusters). Finds state-specific genes with sample-level replication.
      - 'disease_vs_reference_within_cluster': within `cluster`, disease samples vs reference
        samples. Finds disease-driven changes inside one cell type/state.
    Returns top up-regulated genes (log2FC, padj, baseMean) and the number of pseudobulk
    replicates per side. Requires >=2 replicates per side.""",
    input_schema={
        "properties": {
            "cluster": {"type": "string"},
            "cluster_key": {"type": ["string", "null"]},
            "contrast": {"type": "string", "enum": ["cluster_vs_reference_clusters", "disease_vs_reference_within_cluster"], "default": "cluster_vs_reference_clusters"},
            "reference_clusters": {"type": ["array", "null"], "items": {"type": "string"}},
            "min_cells_per_pseudobulk": {"type": "integer", "default": 10},
            "n_top": {"type": "integer", "default": 40},
            "padj_threshold": {"type": "number", "default": 0.05},
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
) -> dict[str, Any]:
    a = session.require_adata()
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available")
    cluster = str(cluster)
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
        pb_id = np.array([f"{s}|{g}" for s, g in zip(samples, group)], dtype=object)
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
    meta, mat = meta[ok], mat[ok.to_numpy()]
    n_t, n_r = int((meta["label"] == "target").sum()), int((meta["label"] == "reference").sum())
    if n_t < 2 or n_r < 2:
        return {
            "cluster": cluster,
            "contrast": contrast,
            "error": f"insufficient pseudobulk replicates (target={n_t}, reference={n_r}); need >=2 per side. Lower min_cells_per_pseudobulk or choose another contrast.",
            "replicates": {"target": n_t, "reference": n_r},
        }
    counts_df = pd.DataFrame(np.round(mat).astype(int), index=meta.index, columns=a.var_names)
    counts_df = counts_df.loc[:, counts_df.sum(axis=0) >= 10]

    method = "pydeseq2"
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        try:
            dds = DeseqDataSet(counts=counts_df, metadata=meta[["label"]], design="~label", quiet=True)
        except TypeError:  # older pydeseq2
            dds = DeseqDataSet(counts=counts_df, metadata=meta[["label"]], design_factors="label", quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["label", "target", "reference"], quiet=True)
        st.summary()
        res = st.results_df.rename(columns={"log2FoldChange": "log2fc"})
        res = res[["baseMean", "log2fc", "pvalue", "padj"]]
    except Exception as e:  # fallback
        method = f"wilcoxon_logcpm (pydeseq2 unavailable: {type(e).__name__})"
        cpm = counts_df.div(counts_df.sum(axis=1), axis=0) * 1e6
        lcpm = np.log2(cpm + 1)
        t = lcpm[meta["label"] == "target"]
        r = lcpm[meta["label"] == "reference"]
        pv = np.array([stats.mannwhitneyu(t[g], r[g], alternative="two-sided").pvalue if t[g].std() + r[g].std() > 0 else 1.0 for g in lcpm.columns])
        res = pd.DataFrame({"baseMean": counts_df.mean(axis=0), "log2fc": (t.mean() - r.mean()), "pvalue": pv, "padj": _bh(pv)}, index=lcpm.columns)

    res.index.name = "gene"
    res = res.sort_values("padj")
    table_path = session.save_table(f"pseudobulk_de_{key}_{cluster}_{contrast}", res)
    sig_up = res[(res["padj"] < padj_threshold) & (res["log2fc"] > 0)].sort_values("log2fc", ascending=False)
    sig_down = res[(res["padj"] < padj_threshold) & (res["log2fc"] < 0)]
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
        "replicates": {"target": n_t, "reference": n_r, "pseudobulk_cells": {str(session.pseudo.to_pseudo(s)) + "|" + l: int(n) for s, l, n in zip(meta["sample"], meta["label"], meta["n_cells"])}},
        "n_sig_up": int(len(sig_up)),
        "n_sig_down": int(len(sig_down)),
        "top_up": sig_up.head(n_top),
    }


@tool(
    name="gene_specificity",
    description="""
    For a list of candidate genes, quantify how specific their expression is to one cluster:
    mean expression in the target cluster vs the maximum mean in any other cluster
    (log2 specificity ratio), % of target cells expressing, expression in disease vs
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
    key = cluster_key or session.cluster_key
    if not key or key not in a.obs:
        raise ValueError("no cluster key available")
    cluster = str(cluster)
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
    rows = []
    dis_samples = sorted(set(samples[conds == session.disease]))
    for j, g in enumerate(present):
        tm = float(means[cluster][j])
        others = {cl: float(means[cl][j]) for cl in clusters if cl != cluster}
        max_other_cl = max(others, key=others.get) if others else None
        max_other = others[max_other_cl] if others else 0.0
        spec = float(np.log2((tm + 0.01) / (max_other + 0.01)))
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
                "log2_specificity": round(spec, 3),
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
