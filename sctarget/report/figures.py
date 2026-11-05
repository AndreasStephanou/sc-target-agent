"""Report figures (report-time; no analysis change). Every figure shows libraries or states, never
single cells, because the pipeline's claims are replicate-level. Each function returns the PNG path
relative to the run directory or None, and writes PDF + SVG beside it (see style.save_all).

* design_overview: cells per library by condition, and lineage composition per library.
* abundance_overview: per-library cluster proportions for every cluster, enriched clusters marked,
  dominant lineage coloured.
* library_expression_panel: for each gene, one point per library (log-CPM of the library's
  pseudobulk within each lineage and each recorded state), coloured by condition.
* state_library_bars: share of each recorded state's cells contributed by each library.
* library_correlation_heatmaps: per lineage, library x library correlation from check_reference_libraries.
* candidate_heatmap: the nominated (or top-ranked) candidates x the evidence components that make
  up the score, with class and attribution.
* state_volcanoes: per recorded state, the replicate-aware contrast that carried its scoring, with
  the candidates labelled.
* candidate_evidence_rows: the compact evidence table used by the report.
"""
from __future__ import annotations

import glob
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .style import (
    DISEASE_COLOR,
    DIVERGING_CMAP,
    DOUBLE_COL,
    LIGHT,
    MID_COL,
    NEUTRAL,
    SEQUENTIAL_CMAP,
    SINGLE_COL,
    apply_style,
    despine,
    lineage_palette,
    panel_label,
    REFERENCE_COLOR,
)


def _plt():
    apply_style()
    import matplotlib.pyplot as plt

    return plt


def _pseudo(session, lib: str) -> str:
    return str(session.pseudo.to_pseudo(lib)) if session.pseudo else str(lib)


def _lib_order(session, libs: np.ndarray, conds: np.ndarray) -> tuple[list[str], dict[str, str]]:
    cond_of = dict(zip(libs, conds))
    order = sorted(set(libs), key=lambda l: (cond_of[l] != session.disease, _pseudo(session, l)))
    return order, cond_of


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


# --------------------------------------------------------------------------- design / composition
def design_overview(session, name: str = "design_overview") -> str | None:
    plt = _plt()
    a = session.require_adata()
    libs = a.obs[session.sample_key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    order, cond_of = _lib_order(session, libs, conds)
    from ..tools.analysis import _lineage_key

    try:
        lk = _lineage_key(session)
    except Exception:  # noqa: BLE001
        lk = None
    lin = a.obs[lk].astype(str).to_numpy() if lk and lk in a.obs else None
    n = len(order)
    fig, axes = plt.subplots(2 if lin is not None else 1, 1, figsize=(max(SINGLE_COL, min(DOUBLE_COL, 0.16 * n + 1.2)), 1.5 * (2 if lin is not None else 1) + 0.6), squeeze=False, sharex=True)
    ax = axes[0, 0]
    counts = pd.Series(libs).value_counts()
    ax.bar(range(n), [counts.get(l, 0) for l in order], color=[DISEASE_COLOR if cond_of[l] == session.disease else REFERENCE_COLOR for l in order], width=0.75)
    ax.set_ylabel("cells")
    ax.set_title(f"Libraries: {sum(cond_of[l] == session.disease for l in order)} {session.disease}, {sum(cond_of[l] != session.disease for l in order)} reference")
    panel_label(ax, "a")
    despine(ax)
    if lin is not None:
        ax2 = axes[1, 0]
        ct = pd.crosstab(pd.Series(libs, name="lib"), pd.Series(lin, name="lineage"), normalize="index").reindex(order).fillna(0)
        pal = lineage_palette(list(ct.columns))
        bottom = np.zeros(n)
        for L in ct.columns:
            vals = ct[L].to_numpy()
            ax2.bar(range(n), vals, bottom=bottom, color=pal[L], width=0.75, label=L)
            bottom += vals
        ax2.set_ylabel("fraction of cells")
        ax2.set_ylim(0, 1)
        ax2.set_title("Lineage composition per library")
        ax2.legend(ncol=min(4, len(ct.columns)), loc="upper center", bbox_to_anchor=(0.5, -0.55), handlelength=1.0, columnspacing=0.8)
        panel_label(ax2, "b")
        despine(ax2)
    axes[-1, 0].set_xticks(range(n))
    axes[-1, 0].set_xticklabels([_pseudo(session, l) for l in order], rotation=90)
    for lab, l in zip(axes[-1, 0].get_xticklabels(), order):
        lab.set_color(DISEASE_COLOR if cond_of[l] == session.disease else REFERENCE_COLOR)
    fig.tight_layout()
    p = session.save_figure(name, fig, f"Study design. (a) Cells per library after QC, coloured by condition ({session.disease} vermilion, reference blue). (b) Lineage composition of each library; a condition that shifts composition shows here before any state is called.")
    return str(p.relative_to(session.run_dir))


def abundance_overview(session, name: str = "abundance_overview") -> str | None:
    plt = _plt()
    a = session.require_adata()
    ck = session.cluster_key
    da = session.facts.get("differential_abundance", {})
    if not ck or ck not in a.obs or not da.get("table"):
        return None
    labels = a.obs[ck].astype(str).to_numpy()
    libs = a.obs[session.sample_key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    rows = sorted(da["table"], key=lambda r: -(r.get("log2_ratio") or 0))
    clusters = [str(r["cluster"]) for r in rows]
    enriched = {str(c) for c in da.get("disease_enriched_clusters", [])}
    comp = {str(r["cluster"]): r for r in session.facts.get("cluster_composition", [])}
    lineages = [str(comp.get(c, {}).get("dominant_celltype") or "?") for c in clusters]
    pal = lineage_palette(lineages)
    lib_tot = pd.Series(libs).value_counts()
    fig, ax = plt.subplots(figsize=(max(SINGLE_COL, min(DOUBLE_COL, 0.22 * len(clusters) + 1.0)), 2.6))
    rng = np.random.default_rng(0)
    for xi, c in enumerate(clusters):
        if c in enriched:
            ax.axvspan(xi - 0.5, xi + 0.5, color=LIGHT, alpha=0.5, lw=0, zorder=0)
        m = labels == c
        per_lib = pd.Series(libs[m]).value_counts()
        for lib in lib_tot.index:
            prop = per_lib.get(lib, 0) / lib_tot[lib]
            is_dis = conds[libs == lib][0] == session.disease
            ax.scatter(xi + (0.16 if is_dis else -0.16) + rng.uniform(-0.06, 0.06), max(prop, 1e-4), s=7, color=DISEASE_COLOR if is_dis else REFERENCE_COLOR, alpha=0.85, edgecolor="none", zorder=3)
    ax.set_yscale("log")
    ax.set_ylabel("fraction of library")
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels(clusters, rotation=90)
    for lab, L in zip(ax.get_xticklabels(), lineages):
        lab.set_color(pal[L])
    for xi, r in enumerate(rows):
        if str(r["cluster"]) in enriched:
            ax.text(xi, ax.get_ylim()[1], f"{r.get('log2_ratio', 0):.1f}", ha="center", va="bottom", fontsize=5, color=NEUTRAL)
    ax.set_title("Per-library cluster proportions (shaded: enriched in disease; number = log2 ratio)")
    handles = [plt.Line2D([], [], marker="o", ls="", color=pal[L], label=L) for L in dict.fromkeys(lineages)]
    handles += [plt.Line2D([], [], marker="o", ls="", color=DISEASE_COLOR, label=f"{session.disease} library"), plt.Line2D([], [], marker="o", ls="", color=REFERENCE_COLOR, label="reference library")]
    ax.legend(handles=handles, ncol=min(5, len(handles)), loc="upper center", bbox_to_anchor=(0.5, -0.28), handlelength=1.0, columnspacing=0.8)
    despine(ax)
    fig.tight_layout()
    p = session.save_figure(name, fig, "Per-library cluster proportions (log scale), one point per library, coloured by condition; cluster labels coloured by dominant lineage; shaded clusters passed the replicate-aware enrichment test (log2 ratio shown).")
    return str(p.relative_to(session.run_dir))


# --------------------------------------------------------------------------- library-level expression
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
        rows.append({"group": g, "library": _pseudo(session, l), "condition": cond_of[l], "n_cells": int(ncell[i]), **{gene: float(cpm[i, j]) for j, gene in enumerate(present)}})
    return pd.DataFrame(rows)


def library_expression_panel(session, genes: list[str], name: str, title: str | None = None, min_cells: int = 20) -> str | None:
    plt = _plt()
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
    fig, axes = plt.subplots(nrow, ncol, figsize=(DOUBLE_COL if ncol == 2 else SINGLE_COL, 1.6 * nrow + 0.5), squeeze=False)
    rng = np.random.default_rng(0)
    for k, (ax, gene) in enumerate(zip(axes.ravel(), genes)):
        for xi, grp in enumerate(order):
            sub = df[df["group"] == grp]
            if grp.startswith("state "):
                ax.axvspan(xi - 0.5, xi + 0.5, color=LIGHT, alpha=0.5, lw=0, zorder=0)
            for is_dis, color in ((False, REFERENCE_COLOR), (True, DISEASE_COLOR)):
                part = sub[(sub["condition"] == session.disease) == is_dis]
                if part.empty:
                    continue
                x0 = xi + (0.18 if is_dis else -0.18)
                vals = part[gene].to_numpy()
                if len(vals) >= 3:
                    ax.boxplot([vals], positions=[x0], widths=0.28, showfliers=False, medianprops={"color": color, "lw": 0.8}, boxprops={"color": color, "lw": 0.6}, whiskerprops={"color": color, "lw": 0.6}, capprops={"color": color, "lw": 0.6})
                ax.scatter(x0 + rng.uniform(-0.07, 0.07, len(vals)), vals, s=7, color=color, alpha=0.85, zorder=3, edgecolor="none")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=35, ha="right")
        ax.set_ylabel("log1p CPM per library")
        ax.set_title(gene)
        panel_label(ax, chr(ord("a") + k))
        despine(ax)
    for ax in axes.ravel()[len(genes):]:
        ax.axis("off")
    axes.ravel()[0].scatter([], [], color=DISEASE_COLOR, label=f"{session.disease} libraries")
    axes.ravel()[0].scatter([], [], color=REFERENCE_COLOR, label=f"{' / '.join(session.reference)} libraries")
    axes.ravel()[0].legend(loc="upper left")
    if title:
        fig.suptitle(title, fontsize=8, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout()
    p = session.save_figure(name, fig, f"Per-library pseudobulk expression of {', '.join(genes)} by lineage and recorded state (shaded); one point per library, coloured by condition, boxes per condition (libraries with < {min_cells} cells in a group omitted). A real candidate has most disease libraries above the reference ones in its own lineage.")
    return str(p.relative_to(session.run_dir))


def state_library_bars(session, name: str = "state_library_shares") -> str | None:
    plt = _plt()
    a = session.require_adata()
    ck = session.cluster_key
    states = session.facts.get("findings", {}).get("disease_state", [])
    if not ck or ck not in a.obs or not states:
        return None
    labels = a.obs[ck].astype(str).to_numpy()
    libs = a.obs[session.sample_key].astype(str).to_numpy()
    conds = a.obs[session.condition_key].astype(str).to_numpy()
    lib_order, cond_of = _lib_order(session, libs, conds)
    fig, axes = plt.subplots(len(states), 1, figsize=(max(SINGLE_COL, min(DOUBLE_COL, 0.16 * len(lib_order) + 1.2)), 1.1 * len(states) + 0.6), squeeze=False, sharex=True)
    for k, (ax, st) in enumerate(zip(axes.ravel(), states)):
        cl = str(st.get("cluster"))
        m = labels == cl
        n = int(m.sum())
        shares = pd.Series(libs[m]).value_counts()
        vals = [shares.get(l, 0) / max(1, n) for l in lib_order]
        ax.bar(range(len(lib_order)), vals, color=[DISEASE_COLOR if cond_of[l] == session.disease else REFERENCE_COLOR for l in lib_order], width=0.75)
        ax.set_ylabel("share of\nstate cells")
        ax.set_title(f"state {cl} ({st.get('basis', 'composition')}): {n:,} cells; {int((pd.Series(vals) > 0.01).sum())} libraries contribute > 1%")
        panel_label(ax, chr(ord("a") + k))
        despine(ax)
    axes.ravel()[-1].set_xticks(range(len(lib_order)))
    axes.ravel()[-1].set_xticklabels([_pseudo(session, l) for l in lib_order], rotation=90)
    for lab, l in zip(axes.ravel()[-1].get_xticklabels(), lib_order):
        lab.set_color(DISEASE_COLOR if cond_of[l] == session.disease else REFERENCE_COLOR)
    fig.tight_layout()
    p = session.save_figure(name, fig, "Share of each recorded state's cells contributed by each library (disease libraries vermilion, reference blue). A state carried by one or two libraries is a patient-specific claim.")
    return str(p.relative_to(session.run_dir))


def library_correlation_heatmaps(session, name: str = "library_correlation") -> str | None:
    plt = _plt()
    chk = session.facts.get("reference_library_check", {})
    mats = chk.get("matrices") or {}
    mats = {L: m for L, m in mats.items() if m.get("corr")}
    if not mats:
        return None
    n = len(mats)
    ncol = min(3, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(DOUBLE_COL if ncol > 1 else SINGLE_COL, 2.3 * nrow + 0.4), squeeze=False)
    im = None
    for k, (ax, (L, m)) in enumerate(zip(axes.ravel(), mats.items())):
        C = np.array(m["corr"])
        libs, labs = m["libraries"], m["labels"]
        order = sorted(range(len(libs)), key=lambda i: (labs[i] != session.disease, libs[i]))
        C = C[np.ix_(order, order)]
        off = C[~np.eye(len(C), dtype=bool)]
        v = max(0.3, float(np.abs(off).max())) if off.size else 1.0
        im = ax.imshow(C, cmap=DIVERGING_CMAP, vmin=-v, vmax=v)
        ticks = [str(libs[i]) for i in order]
        ax.set_xticks(range(len(ticks)))
        ax.set_yticks(range(len(ticks)))
        ax.set_xticklabels(ticks, rotation=90, fontsize=5)
        ax.set_yticklabels(ticks, fontsize=5)
        for t, i in zip(ax.get_yticklabels(), order):
            t.set_color(DISEASE_COLOR if labs[i] == session.disease else REFERENCE_COLOR)
        for t, i in zip(ax.get_xticklabels(), order):
            t.set_color(DISEASE_COLOR if labs[i] == session.disease else REFERENCE_COLOR)
        flagged = {f["library"] for f in chk.get("flagged", []) if f.get("lineage") == L}
        ax.set_title(f"{L}{'  flagged: ' + ', '.join(sorted(flagged)) if flagged else ''}", fontsize=7)
        panel_label(ax, chr(ord("a") + k), x=-0.25)
        for s_ in ("top", "right", "left", "bottom"):
            ax.spines[s_].set_visible(False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    if im is not None:
        cb = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.5, pad=0.02)
        cb.set_label("pseudobulk correlation")
        cb.outline.set_visible(False)
    p = session.save_figure(name, fig, "Library x library correlation of lineage pseudobulks (check_reference_libraries); label colours: disease vermilion, reference blue. A library whose row matches the other colour's block has a doubtful label or a mixed sample.")
    return str(p.relative_to(session.run_dir))


# --------------------------------------------------------------------------- candidates
def _picks(session, max_rows: int) -> list[tuple[str, str]]:
    f = session.facts
    rt = f.get("ranked_targets", {})
    cands = f.get("findings", {}).get("candidate_target", [])
    picks: list[tuple[str, str]] = [(c["genes"][0], str(c.get("cluster"))) for c in cands if c.get("genes")]
    if not picks:
        for cl, rws in rt.items():
            picks += [(r["gene"], str(cl)) for r in rws[:5]]
    seen: set[tuple[str, str]] = set()
    out = []
    for g, cl in picks:
        if (g, cl) not in seen and any(r["gene"] == g for r in rt.get(cl, [])):
            seen.add((g, cl))
            out.append((g, cl))
    return out[:max_rows]


def candidate_evidence_rows(session, max_rows: int = 15) -> list[dict[str, Any]]:
    """Compact evidence table for the nominated candidates (falls back to top-ranked genes when the
    nominator recorded nothing)."""
    f = session.facts
    rt = f.get("ranked_targets", {})
    flags_by_gene: dict[str, list[str]] = {}
    for fl in f.get("findings", {}).get("critic_flag", []):
        for g in fl.get("genes", []) or []:
            flags_by_gene.setdefault(g, []).append(fl.get("title", "").split(":")[-1].strip()[:40])
    drug = f.get("druggability", {})
    rows: list[dict[str, Any]] = []
    for g, cl in _picks(session, max_rows * 2):
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
            "therapeutic_direction": "human check: do the listed compounds act in the direction this design needs?",
        })
    rows.sort(key=lambda r: -(r.get("score") or 0))
    return rows[:max_rows]


def _clip01(x: float | None, cap: float) -> float:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 0.0
    return float(min(max(x, 0.0), cap) / cap)


def candidate_heatmap(session, name: str = "candidate_evidence", max_rows: int = 15) -> str | None:
    """Candidates x evidence components (each scaled to 0-1 the way the score uses them)."""
    plt = _plt()
    rt = session.facts.get("ranked_targets", {})
    picks = _picks(session, max_rows)
    if not picks:
        return None
    comps = [("state log2FC", lambda r: _clip01(r.get("effect_size_log2fc"), 4.0)),
             ("lineage log2FC", lambda r: _clip01(r.get("lineage_log2fc_disease_vs_reference"), 4.0)),
             ("specificity", lambda r: _clip01(r.get("log2_specificity"), 4.0)),
             ("consistency", lambda r: _clip01(r.get("sample_consistency"), 1.0)),
             ("druggability", lambda r: _clip01(r.get("druggability_score"), 1.0)),
             ("compounds", lambda r: _clip01(r.get("known_compound_score"), 1.0)),
             ("flag penalty", lambda r: _clip01(r.get("flag_penalty"), 1.0))]
    recs = []
    for g, cl in picks:
        rec = next((r for r in rt.get(cl, []) if r["gene"] == g), None)
        if rec:
            recs.append((g, cl, rec))
    recs.sort(key=lambda t: -(t[2].get("score") or 0))
    M = np.array([[fn(rec) for _, fn in comps] for _, _, rec in recs])
    fig, (ax, axs) = plt.subplots(1, 2, figsize=(MID_COL, 0.22 * len(recs) + 0.9), gridspec_kw={"width_ratios": [len(comps), 1.2], "wspace": 0.05})
    im = ax.imshow(M, cmap=SEQUENTIAL_CMAP, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(comps)))
    ax.set_xticklabels([c for c, _ in comps], rotation=45, ha="right")
    ax.set_yticks(range(len(recs)))
    ax.set_yticklabels([f"{g}  (state {cl})" for g, cl, _ in recs])
    for i, (_, _, rec) in enumerate(recs):
        ax.text(len(comps) - 0.5 + 0.08, i, f"{rec.get('target_class') or ''} · {rec.get('attribution') or ''}", va="center", ha="left", fontsize=5, color=NEUTRAL, transform=ax.transData, clip_on=False)
    for s_ in ("top", "right", "left", "bottom"):
        ax.spines[s_].set_visible(False)
    ax.set_title("Evidence components of the nominated candidates (0-1, as used by the score)")
    axs.barh(range(len(recs)), [rec.get("score") or 0 for _, _, rec in recs], color=NEUTRAL, height=0.7)
    axs.set_ylim(len(recs) - 0.5, -0.5)
    axs.set_yticks([])
    axs.set_xlim(0, 1)
    axs.set_xlabel("score")
    despine(axs, left=True)
    cb = fig.colorbar(im, ax=[ax, axs], shrink=0.6, pad=0.25, location="bottom")
    cb.set_label("component value (0-1)")
    cb.outline.set_visible(False)
    p = session.save_figure(name, fig, "Evidence behind each nominated candidate: the components of the deterministic score, each scaled 0-1 as the score uses them (effect sizes capped at log2FC 4), with target class and attribution; right, the final score. Druggability counts compounds and does not know the therapeutic direction.")
    return str(p.relative_to(session.run_dir))


def _contrast_table(session, cluster: str) -> tuple[pd.DataFrame | None, str]:
    """Full pseudobulk table (all genes) for the contrast that carried the state's scoring."""
    from ..tools.targets import state_basis

    basis = state_basis(session, str(cluster))
    rt = session.facts.get("ranked_targets", {}).get(str(cluster), [])
    axis = next((r.get("lineage_axis_source") for r in rt if r.get("lineage_axis_source")), None)
    if basis == "program":
        contrast = "disease_vs_reference_within_cluster"
    elif axis:
        contrast = str(axis)
    else:
        contrast = "cluster_vs_reference_clusters"
    safe_cl = str(cluster).replace("+", "_")
    pats = [f"*pseudobulk_de_*_{safe_cl}_{contrast}.csv", f"*pseudobulk_de_*{safe_cl}*{contrast}*.csv"]
    for pat in pats:
        hits = sorted(glob.glob(str(session.run_dir / "tables" / pat)))
        if hits:
            df = pd.read_csv(hits[-1])
            gcol = "gene" if "gene" in df.columns else df.columns[0]
            df = df.rename(columns={gcol: "gene"})
            if "log2fc" in df.columns and "padj" in df.columns:
                return df, contrast
    return None, contrast


def state_volcanoes(session, name: str = "state_volcanoes", max_states: int = 6) -> str | None:
    plt = _plt()
    states = session.facts.get("findings", {}).get("disease_state", [])[:max_states]
    if not states:
        return None
    panels = []
    for st in states:
        cl = str(st.get("cluster"))
        df, contrast = _contrast_table(session, cl)
        if df is not None and len(df) > 20:
            panels.append((cl, st, df, contrast))
    if not panels:
        return None
    nominated = {(c["genes"][0], str(c.get("cluster"))) for c in session.facts.get("findings", {}).get("candidate_target", []) if c.get("genes")}
    ncol = min(3, len(panels))
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(DOUBLE_COL if ncol > 1 else SINGLE_COL, 2.2 * nrow + 0.3), squeeze=False)
    for k, (ax, (cl, st, df, contrast)) in enumerate(zip(axes.ravel(), panels)):
        x = df["log2fc"].to_numpy(dtype=float)
        y = -np.log10(np.clip(df["padj"].to_numpy(dtype=float), 1e-300, 1.0))
        sig = (df["padj"] < 0.05) & (np.abs(df["log2fc"]) >= 1)
        ax.scatter(x[~sig], y[~sig], s=2, color=LIGHT, edgecolor="none", rasterized=True)
        ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=2.5, color=DISEASE_COLOR, alpha=0.6, edgecolor="none", rasterized=True)
        ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=2.5, color=REFERENCE_COLOR, alpha=0.6, edgecolor="none", rasterized=True)
        genes = df["gene"].astype(str).to_numpy()
        lab = [i for i, g in enumerate(genes) if (g, cl) in nominated]
        top = list(np.argsort(-(y * (x > 0)))[:8])
        for i in dict.fromkeys(lab + top):
            if i in lab or sig.iloc[i]:
                ax.annotate(genes[i], (x[i], y[i]), fontsize=5, fontweight="bold" if i in lab else "normal", color="black" if i in lab else NEUTRAL, xytext=(2, 2), textcoords="offset points")
        ax.axhline(-np.log10(0.05), color=NEUTRAL, lw=0.5, ls=":")
        for v in (-1, 1):
            ax.axvline(v, color=NEUTRAL, lw=0.5, ls=":")
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 padj")
        ax.set_title(f"state {cl}: {contrast.replace('disease_vs_reference_within_', 'disease vs reference, ').replace('_', ' ')}", fontsize=6.5)
        panel_label(ax, chr(ord("a") + k))
        despine(ax)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    fig.tight_layout()
    p = session.save_figure(name, fig, "Replicate-aware pseudobulk contrast that carried each recorded state's scoring (vermilion: up in disease, blue: down; padj < 0.05 and |log2FC| >= 1). Nominated candidates in bold; the eight most significant up-regulated genes labelled.")
    return str(p.relative_to(session.run_dir))
