"""Replicate-level report figures (patch 37, report-only; no analysis change).

The pipeline's claims are about libraries, not cells, so the figures show libraries:

* library_expression_panel: for each gene, one point per library (log-CPM of the library's pseudobulk
  within each lineage, plus the recorded disease states as extra groups), coloured by condition, with a
  box per condition. Shows at a glance whether a gene is up in disease in its own lineage, how many
  libraries carry it, outlier libraries, and leakage (ambient) into other lineages.
* state_library_bars: for each recorded state, the share of its cells contributed by each library,
  coloured by condition (patient splits and single-sample states become visible).
* library_correlation_heatmaps: per lineage, the library x library correlation matrix from
  check_reference_libraries, with the condition labels on the axes.

Everything is computed from raw counts and facts already in the session; nothing here feeds back into
the analysis. Each function returns the figure path (relative to the run dir) or None.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

DISEASE_COLOR = "#c0392b"
REFERENCE_COLOR = "#2c3e50"


def _groups(session) -> tuple[np.ndarray, list[str]]:
    """Group label per cell: lineage name, with each recorded disease-state cluster as its own group."""
    from ..tools.analysis import _lineage_key

    a = session.require_adata()
    try:
        lk = _lineage_key(session)
    except Exception:  # noqa: BLE001
        lk = None
    lin = a.obs[lk].astype(str).to_numpy() if lk and lk in a.obs else np.array(["all cells"] * a.n_obs)
    groups = lin.copy()
    order = sorted(set(lin))
    ck = session.cluster_key
    if ck and ck in a.obs:
        labels = a.obs[ck].astype(str).to_numpy()
        for st in session.facts.get("findings", {}).get("disease_state", []):
            cl = str(st.get("cluster"))
            m = labels == cl
            if m.any():
                name = f"state {cl}"
                groups[m] = name
                order.append(name)
    return groups, order


def _library_pseudobulk_logcpm(session, genes: list[str]) -> pd.DataFrame | None:
    """Rows: (group, library, condition); columns: genes; values: log1p CPM of the group x library pseudobulk."""
    a = session.require_adata()
    X = session.counts()
    present = [g for g in genes if g in a.var_names]
    if not present:
        return None
    groups, _ = _groups(session)
    libs = a.obs[session.sample_key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    key = pd.Series([f"{g}\t{l}" for g, l in zip(groups, libs)])
    codes, uniq = pd.factorize(key)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    totals = np.asarray(ind @ (np.asarray(X.sum(axis=1)).ravel() if sp.issparse(X) else X.sum(axis=1))).ravel()
    gi = [a.var_names.get_loc(g) for g in present]
    sub = X[:, gi]
    sums = np.asarray((ind @ sub).todense() if sp.issparse(sub) else ind @ np.asarray(sub))
    ncell = np.asarray(ind.sum(axis=1)).ravel()
    cpm = np.log1p(sums / np.maximum(totals, 1.0)[:, None] * 1e6)
    cond_of = dict(zip(libs, conds))
    rows = []
    for i, u in enumerate(uniq):
        g, l = u.split("\t")
        rows.append({"group": g, "library": str(session.pseudo.to_pseudo(l)) if session.pseudo else l, "condition": cond_of[l], "n_cells": int(ncell[i]), **{gene: float(cpm[i, j]) for j, gene in enumerate(present)}})
    return pd.DataFrame(rows)


def library_expression_panel(session, genes: list[str], name: str, title: str | None = None, min_cells: int = 20) -> str | None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    df = _library_pseudobulk_logcpm(session, genes)
    if df is None or df.empty:
        return None
    df = df[df["n_cells"] >= min_cells]
    _, order = _groups(session)
    order = [g for g in order if g in set(df["group"])]
    genes = [g for g in genes if g in df.columns]
    if not genes or not order:
        return None
    ncol = 2 if len(genes) > 1 else 1
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(max(6, 0.9 * len(order) + 2) * ncol, 2.8 * nrow), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, gene in zip(axes.ravel(), genes):
        for xi, grp in enumerate(order):
            sub = df[df["group"] == grp]
            for off, (is_dis, color) in enumerate([(False, REFERENCE_COLOR), (True, DISEASE_COLOR)]):
                part = sub[(sub["condition"] == session.disease) == is_dis]
                if part.empty:
                    continue
                x0 = xi + (-0.18 if not is_dis else 0.18)
                vals = part[gene].to_numpy()
                if len(vals) >= 3:
                    ax.boxplot([vals], positions=[x0], widths=0.28, showfliers=False, medianprops={"color": color}, boxprops={"color": color}, whiskerprops={"color": color}, capprops={"color": color})
                ax.scatter(x0 + rng.uniform(-0.07, 0.07, len(vals)), vals, s=14, color=color, alpha=0.85, zorder=3, edgecolor="none")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=35, ha="right", fontsize=8)
        ax.set_ylabel("log1p CPM / library", fontsize=8)
        ax.set_title(gene, fontsize=10, loc="left")
        ax.grid(axis="y", alpha=0.25)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")
    axes.ravel()[0].scatter([], [], color=DISEASE_COLOR, label=f"{session.disease} libraries")
    axes.ravel()[0].scatter([], [], color=REFERENCE_COLOR, label=f"{' / '.join(session.reference)} libraries")
    axes.ravel()[0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle(title or "Per-library pseudobulk expression by lineage and recorded state (one point per library)", fontsize=10, y=1.0)
    fig.tight_layout()
    p = session.save_figure(name, fig, f"Per-library pseudobulk expression of {', '.join(genes)} by lineage and recorded state; one point per library, coloured by condition, boxes per condition (libraries with < {min_cells} cells in a group omitted).")
    return str(p.relative_to(session.run_dir))


def state_library_bars(session, name: str = "state_library_shares") -> str | None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = session.require_adata()
    ck = session.cluster_key
    states = session.facts.get("findings", {}).get("disease_state", [])
    if not ck or ck not in a.obs or not states:
        return None
    labels = a.obs[ck].astype(str).to_numpy()
    libs = a.obs[session.sample_key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    cond_of = dict(zip(libs, conds))
    lib_order = sorted(set(libs), key=lambda l: (cond_of[l] != session.disease, l))
    pseudo = {l: (str(session.pseudo.to_pseudo(l)) if session.pseudo else l) for l in lib_order}
    fig, axes = plt.subplots(len(states), 1, figsize=(max(6, 0.45 * len(lib_order) + 1), 1.9 * len(states)), squeeze=False, sharex=True)
    for ax, st in zip(axes.ravel(), states):
        cl = str(st.get("cluster"))
        m = labels == cl
        n = int(m.sum())
        shares = pd.Series(libs[m]).value_counts()
        vals = [shares.get(l, 0) / max(1, n) for l in lib_order]
        ax.bar(range(len(lib_order)), vals, color=[DISEASE_COLOR if cond_of[l] == session.disease else REFERENCE_COLOR for l in lib_order])
        ax.set_ylabel("share of\nstate cells", fontsize=8)
        ax.set_title(f"state {cl}: {n:,} cells; {int((pd.Series(vals) > 0.01).sum())} libraries contribute > 1%", fontsize=9, loc="left")
        ax.grid(axis="y", alpha=0.25)
        for s_ in ("top", "right"):
            ax.spines[s_].set_visible(False)
    axes.ravel()[-1].set_xticks(range(len(lib_order)))
    axes.ravel()[-1].set_xticklabels([pseudo[l] for l in lib_order], rotation=90, fontsize=7)
    fig.tight_layout()
    p = session.save_figure(name, fig, "Share of each recorded state's cells contributed by each library (disease libraries red, reference dark). A state carried by one or two libraries is a patient-specific claim.")
    return str(p.relative_to(session.run_dir))


def library_correlation_heatmaps(session, name: str = "library_correlation") -> str | None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chk = session.facts.get("reference_library_check", {})
    mats = chk.get("matrices") or {}
    mats = {L: m for L, m in mats.items() if m.get("corr")}
    if not mats:
        return None
    n = len(mats)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.8 * nrow), squeeze=False)
    for ax, (L, m) in zip(axes.ravel(), mats.items()):
        C = np.array(m["corr"])
        libs, labs = m["libraries"], m["labels"]
        order = sorted(range(len(libs)), key=lambda i: (labs[i] != session.disease, libs[i]))
        C = C[np.ix_(order, order)]
        im = ax.imshow(C, cmap="RdBu_r", vmin=-max(0.3, np.abs(C[~np.eye(len(C), dtype=bool)]).max()), vmax=max(0.3, np.abs(C[~np.eye(len(C), dtype=bool)]).max()))
        ticks = [f"{libs[i]}" for i in order]
        ax.set_xticks(range(len(ticks)))
        ax.set_yticks(range(len(ticks)))
        ax.set_xticklabels(ticks, rotation=90, fontsize=6)
        ax.set_yticklabels(ticks, fontsize=6)
        for t, i in zip(ax.get_yticklabels(), order):
            t.set_color(DISEASE_COLOR if labs[i] == session.disease else REFERENCE_COLOR)
        for t, i in zip(ax.get_xticklabels(), order):
            t.set_color(DISEASE_COLOR if labs[i] == session.disease else REFERENCE_COLOR)
        flagged = {f["library"] for f in chk.get("flagged", []) if f.get("lineage") == L}
        ax.set_title(f"{L} (structure {m.get('structure')}){'  flagged: ' + ', '.join(sorted(flagged)) if flagged else ''}", fontsize=8, loc="left")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="pseudobulk correlation")
    p = session.save_figure(name, fig, "Library x library correlation of lineage pseudobulks (check_reference_libraries); label colours: disease red, reference dark. A library whose row matches the other colour's block has a doubtful label or a mixed sample.")
    return str(p.relative_to(session.run_dir))


def candidate_evidence_rows(session, max_rows: int = 15) -> list[dict[str, Any]]:
    """Compact evidence table for the nominated candidates (falls back to top-ranked genes when the
    nominator recorded nothing)."""
    f = session.facts
    rt = f.get("ranked_targets", {})
    cands = f.get("findings", {}).get("candidate_target", [])
    flags_by_gene: dict[str, list[str]] = {}
    for fl in f.get("findings", {}).get("critic_flag", []):
        for g in fl.get("genes", []) or []:
            flags_by_gene.setdefault(g, []).append(fl.get("title", "").split(":")[-1].strip()[:40])
    drug = f.get("druggability", {})
    rows: list[dict[str, Any]] = []
    seen = set()
    picks: list[tuple[str, str]] = [(c["genes"][0], str(c.get("cluster"))) for c in cands if c.get("genes")]
    if not picks:
        for cl, rws in rt.items():
            picks += [(r["gene"], str(cl)) for r in rws[:5]]
    for g, cl in picks:
        if (g, cl) in seen:
            continue
        seen.add((g, cl))
        rec = next((r for r in rt.get(cl, []) if r["gene"] == g), None)
        if rec is None:
            continue
        d = drug.get(rec.get("human_symbol", g)) or drug.get(g) or {}
        n_inh, n_cmp = d.get("n_inhibitors", 0) or 0, d.get("n_compounds", 0) or 0
        direction = f"{n_inh} inhibitor(s)" if n_inh else (f"{n_cmp} compound(s), none annotated as inhibitor" if n_cmp else "no compound")
        rows.append({
            "gene": g, "state": cl, "rank": rec.get("rank"), "score": rec.get("score"), "class": rec.get("target_class"),
            "state_log2fc": rec.get("effect_size_log2fc"), "lineage_log2fc": rec.get("lineage_log2fc_disease_vs_reference"),
            "axis": (rec.get("lineage_axis_source") or "").replace("disease_vs_reference_within_", ""),
            "specificity": rec.get("log2_specificity"), "consistency": rec.get("sample_consistency"), "attribution": rec.get("attribution"),
            "drugs": direction, "top": ", ".join((rec.get("top_inhibitors") or [])[:2]), "flags": ", ".join(rec.get("flags", [])), "critic": "; ".join(flags_by_gene.get(g, [])[:2]),
        })
    rows.sort(key=lambda r: -(r.get("score") or 0))
    return rows[:max_rows]
