"""Build the Markdown + HTML report from the session facts and run log.

Layout (0.4): an executive summary a reader can act on (design statement, metrics, the top
candidates with their evidence and the direction question), then the evidence sections in the order
a reviewer reads them, then provenance and a figure index. Every figure is also written as PDF and
SVG (publication vector files) next to its PNG; the index lists them.
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path
from typing import Any

from ..session import AnalysisSession
from .figures import candidate_evidence_rows

CSS = """
body{font-family:Helvetica,Arial,'Liberation Sans',sans-serif;max-width:1040px;margin:2rem auto;padding:0 1rem;color:#1c1c1c;line-height:1.45}
h1{border-bottom:2px solid #333;padding-bottom:.3rem}h2{margin-top:2.2rem;border-bottom:1px solid #ccc}h3{margin-top:1.6rem}
table{border-collapse:collapse;font-size:.86rem;margin:.6rem 0}th,td{border:1px solid #ddd;padding:.3rem .5rem;text-align:left;vertical-align:top}
th{background:#f3f3f3}code{background:#f5f5f5;padding:0 .2rem}img{max-width:100%;border:1px solid #eee;margin:.5rem 0;background:#fff}
.flag-blocking{color:#b00020;font-weight:600}.flag-major{color:#c75000;font-weight:600}.flag-minor{color:#666}
.small{color:#666;font-size:.85rem}.kpi{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:.5rem;margin:.6rem 0}
.kpi div{border:1px solid #e3e3e3;border-radius:6px;padding:.45rem .6rem;background:#fafafa}.kpi b{display:block;font-size:1.25rem;color:#0072B2}
blockquote{border-left:3px solid #0072B2;margin:.6rem 0;padding:.2rem .8rem;color:#333;background:#f7fafc}
"""


def _fmt(x: Any, nd: int = 3) -> str:
    if x is None:
        return "–"
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, float):
        if abs(x) < 1e-3 and x != 0:
            return f"{x:.1e}"
        return f"{x:.{nd}f}"
    if isinstance(x, list):
        return ", ".join(str(i) for i in x) if x else "–"
    return str(x)


def _table(rows: list[dict[str, Any]], cols: list[tuple[str, str]], nd: int = 3) -> str:
    if not rows:
        return "_none_\n"
    head = "| " + " | ".join(h for _, h in cols) + " |\n|" + "|".join("---" for _ in cols) + "|\n"
    body = "".join("| " + " | ".join(_fmt(r.get(k), nd) for k, _ in cols) + " |\n" for r in rows)
    return head + body


def _kv(rows: list[tuple[str, Any]]) -> str:
    return "| | |\n|---|---|\n" + "".join(f"| {k} | {_fmt(v)} |\n" for k, v in rows)


def _fig_md(x: dict[str, Any] | None) -> str:
    return f"\n![{x['caption']}]({x['path']})\n" if x else ""


def _report_time_figures(session: AnalysisSession) -> dict[str, Any]:
    """Figures built when the report is rendered. Each one is optional: a failure is logged and the
    report goes on without it. Data-dependent figures need the object in memory."""
    out: dict[str, Any] = {"candidate_panels": []}
    existing = {Path(x["path"]).name for x in session.figures}

    def _fig(rel):
        return next((x for x in session.figures if x["path"] == rel), None)

    def _have(tag: str):
        return next((x for x in session.figures if tag in Path(x["path"]).name), None)

    def _try(tag: str, fn, *args, **kw):
        if any(tag in n for n in existing):
            return _have(tag)
        if fn is None:  # a process that imported an older figures module simply skips the figure
            return None
        try:
            rel = fn(session, *args, **kw)
            return _fig(rel) if rel else None
        except Exception as e:  # noqa: BLE001 - figures never break a report
            try:
                session.log.note("report", f"figure {tag} skipped: {type(e).__name__}: {e}")
            except Exception:  # noqa: BLE001
                pass
            return None

    from . import figures as F

    g = lambda n: getattr(F, n, None)  # noqa: E731 - tolerate a stale module in a long-running process
    if session.adata is not None:
        out["design"] = _try("design_overview", g("design_overview"))
        out["abundance"] = _try("abundance_overview", g("abundance_overview"))
        out["state_map"] = _try("state_map", g("state_map"))
        out["state_shares"] = _try("state_library_shares", g("state_library_bars"))
        f = session.facts
        cands = f.get("findings", {}).get("candidate_target", [])
        rt = f.get("ranked_targets", {})
        per_state: dict[str, list[str]] = {}
        for c in cands:
            if c.get("genes"):
                per_state.setdefault(str(c.get("cluster")), []).append(c["genes"][0])
        if not per_state:
            for cl, rows in rt.items():
                per_state[str(cl)] = [r["gene"] for r in rows[:6]]
        resolve = g("resolved_state_label")
        for cl, genes in per_state.items():
            genes = list(dict.fromkeys(genes))[:8]
            if not genes:
                continue
            name = f"library_expression_state{cl}_" + "_".join(genes[:3])
            shown = None
            try:
                shown = resolve(session, cl) if resolve else None
            except Exception:  # noqa: BLE001
                shown = None
            where = f"state {cl}" if not shown or shown == str(cl) else f"state {cl} (drawn as merged state {shown})"
            x = _try(name, g("library_expression_panel"), genes, name, title=f"Nominated candidates for {where}: per-library pseudobulk expression")
            if x:
                out["candidate_panels"].append(x)
    out["library_corr"] = _try("library_correlation", g("library_correlation_heatmaps"))
    out["candidate_heatmap"] = _try("candidate_evidence", g("candidate_heatmap"))
    out["volcanoes"] = _try("state_volcanoes", g("state_volcanoes"))
    return out


def _eval_section(session: AnalysisSession) -> str:
    return eval_section_from_dir(session.run_dir)


_NO_KEY_TEXT = "## 6. Evaluation\n\n_No frozen answer key was scored for this run (no `eval.json`). Without a key the sanity checks of the no-key protocol apply: lineage markers classed shared_with_reference, artefact genes flagged, stability across repeat runs, per-library panels, direction._\n"


def refresh_eval_section(run_dir: Path | str) -> bool:
    """Rewrite section 6 of an existing report.md (and report.html) from the run's eval.json.

    The report is built at the end of the run, before `run_eval` scores it, so every scored run
    used to carry the "no key" placeholder in section 6 while the result sat in eval.json and the
    log. `run_eval` calls this after writing eval.json. Returns True when the report was rewritten."""
    run_dir = Path(run_dir)
    md_path = run_dir / "report.md"
    if not md_path.exists():
        return False
    new = eval_section_from_dir(run_dir) or _NO_KEY_TEXT
    md = md_path.read_text()
    start = md.find("\n## 6. ")
    end = md.find("\n## 7. ", start + 1) if start >= 0 else -1
    if start < 0 or end < 0:
        return False
    md = md[: start + 1] + new.rstrip("\n") + "\n\n" + md[end + 1 :]
    md_path.write_text(md)
    _write_html(run_dir, md)
    return True


def eval_section_from_dir(run_dir: Path) -> str:
    p = Path(run_dir) / "eval.json"
    if not p.exists():
        return ""
    try:
        ev = json.loads(p.read_text())
    except Exception:  # noqa: BLE001
        return ""
    md = ["## 6. Evaluation against the frozen answer key\n"]  # numbering: sections 7-8 follow whether or not this one exists
    md.append(f"Key: `{ev.get('ground_truth')}` · prompts v{ev.get('prompt_version')} · scoring v{ev.get('scoring_version')} · **{'PASS' if ev.get('passed') else 'not passed'}** under the key's own criteria.\n\n")
    rows = []
    for name, r in (ev.get("state_recovery") or {}).items():
        rows.append({"state": name, "lenient": r.get("found"), "strict": r.get("found_strict"), "cluster": r.get("best_match"), "need": r.get("required_overlap")})
    md.append("**State recovery** (lenient: the matched cluster's tables carry the signature; strict: the agent itself recorded it):\n\n" + _table(rows, [("state", "Expected state"), ("lenient", "Found"), ("strict", "Strict"), ("cluster", "Best-matching cluster"), ("need", "Markers required")]))
    md.append("\n**Target ranks:** " + ", ".join(f"{g}: {('#' + str(v)) if v else 'not ranked'}" for g, v in (ev.get("target_rank") or {}).items()) + "\n")
    acc = {g: v for g, v in (ev.get("acceptable_rank") or {}).items() if v}
    md.append("\n**Acceptable targets ranked:** " + (", ".join(f"{g} #{v}" for g, v in sorted(acc.items(), key=lambda kv: kv[1])) or "none") + "\n")
    md.append(f"\n**Decoys above the first expected target:** {_fmt(ev.get('decoys_above_target'))} · **novel genes above it:** {len(ev.get('novel_above_target') or [])} · **decoys nominated unflagged:** {_fmt(ev.get('decoys_nominated_unflagged'))} · **nominated precision:** {_fmt(ev.get('nominated_precision'))} · **drug recall:** {_fmt(ev.get('drug_recall'))}\n")
    md.append("\n> A missed expected state can be a resolution problem (the state sits inside a larger cluster), a key-drafting problem (literature markers that are not consistent in this dataset) or a method miss; the matched cluster and its tables tell which.\n")
    return "\n".join(md)


def build_markdown(session: AnalysisSession) -> str:
    f = session.facts
    cfg_ds = session.cfg.dataset
    ds = f.get("dataset", {})
    findings = f.get("findings", {})
    states = findings.get("disease_state", [])
    cands = findings.get("candidate_target", [])
    flags = findings.get("critic_flag", [])
    summaries = findings.get("summary", [])
    events = session.log.events()
    # backend name as configured, plus the model ids the backend reported having used (a CLI backend
    # without --model runs on a default the pipeline does not set; facts["model_ids"] records it)
    model_names = sorted({e.get("model") for e in events if e["type"] == "model_call" and e.get("model")})
    model_ids = sorted({str(m) for m in (session.facts.get("model_ids") or []) if m})
    if model_ids and not any(any(m.endswith(":" + i) or m == i for m in model_names) for i in model_ids):
        model_names = [f"{n} ({', '.join(model_ids)})" for n in model_names] or model_ids
    usage = session.log.usage_totals()
    calls = [e for e in events if e["type"] == "tool_call"]
    n_err = len([e for e in calls if e.get("error")])
    elapsed = max((e.get("elapsed_s") or 0) for e in events) if events else None
    report_figs = _report_time_figures(session)
    fig_by_name = {Path(x["path"]).name: x for x in session.figures}
    objective = str(session.cfg.targets.get("objective", "target"))
    comp = {str(r["cluster"]): r for r in f.get("cluster_composition", [])}
    da = f.get("differential_abundance", {})
    rt = f.get("ranked_targets", {})
    uniq_cands = list(dict.fromkeys(c["genes"][0] for c in cands if c.get("genes")))
    sev = {s: len([fl for fl in flags if fl.get("severity") == s]) for s in ("blocking", "major", "minor")}
    chk = f.get("reference_library_check") or {}
    flag_by_gene: dict[str, list[str]] = {}
    for fl in flags:
        for g in fl.get("genes") or []:
            flag_by_gene.setdefault(g, []).append(f"{fl.get('severity', '')}: {fl['title']}")

    md: list[str] = []
    md.append(f"# {'Target discovery' if objective == 'target' else 'Response-gene'} report — {ds.get('disease', '?')} vs {', '.join(ds.get('reference_conditions', []))}\n")
    md.append(
        f"<span class='small'>Run `{session.run_dir.name}` · {dt.datetime.now().strftime('%Y-%m-%d %H:%M')} · model(s): {', '.join(model_names) or '–'} · prompts v{f.get('prompt_version', 1)} · scoring v{f.get('scoring_version', 1)} · objective `{objective}` · "
        f"tokens in/out: {usage.get('input_tokens', 0):,}/{usage.get('output_tokens', 0):,} · {len(calls)} tool calls</span>\n"
    )
    md.append(
        "> **Read this first.** Every number below was produced by a deterministic tool call recorded in `run.jsonl`; "
        "the model chose which tools to run and wrote the interpretation. Candidates are hypotheses for orthogonal validation, not conclusions. "
        "Druggability counts compounds and does not know the therapeutic direction. Values marked *fixture* come from bundled illustrative data, not a live database.\n"
    )

    # ---- 0. executive summary
    md.append("## Executive summary\n")
    md.append("### Design\n")
    descr = cfg_ds.get("condition_descriptions") or {}
    cond_lines = []
    for k, v in ds.get("condition_counts", {}).items():
        role = "case" if str(k) == str(ds.get("disease")) else "reference"
        cond_lines.append(f"- **{k}** ({role}): {v:,} cells, {ds.get('samples_per_condition', {}).get(k, '?')} libraries" + (f" — {descr[k]}" if k in descr else ""))
    md.append("\n".join(cond_lines) + "\n")
    lin_src = ("annotation `" + str(ds.get("celltype_key")) + "`") if ds.get("celltype_key") and not f.get("lineages") else (("derived from the data (coarse Leiden r" + str(f.get("lineages", {}).get("resolution")) + ", " + str(len(f.get("lineages", {}).get("groups", {}))) + " groups)") if f.get("lineages") else "none")
    design_rows: list[tuple[str, Any]] = [
        ("Objective", "druggable targets (`target`)" if objective == "target" else "condition-specific response genes (`marker`; druggability reported, not scored)"),
        ("Species", ds.get("species")),
        ("Replicates", "libraries are pseudonymised S01…; the per-sample tests treat each library as one replicate — if libraries are nested in donors, the config comments and the design check say so"),
        ("Lineages", lin_src),
        ("Clustering", f"`{f.get('clustering', {}).get('cluster_key', '–')}` ({f.get('clustering', {}).get('method', '–')}, harmony={f.get('clustering', {}).get('harmony_used')}); base clustering `{f.get('clustering', {}).get('base_cluster_key', '–')}`"),
    ]
    if "qc" in f:
        q = f["qc"]
        design_rows.append(("QC", f"{q['cells_before']:,} → {q['cells_after']:,} cells, {q['genes_before']:,} → {q['genes_after']:,} genes ({q.get('mode', 'fixed')} thresholds; min_genes={q['params']['min_genes']}, max_mt%={q['params']['max_mt_pct']})"))
    if f.get("excluded_samples"):
        ex = f["excluded_samples"]
        design_rows.append(("Excluded by config", f"{ex.get('n_samples')} libraries ({ex.get('n_cells', 0):,} cells). Reason: {ex.get('reason') or 'none recorded'}"))
    else:
        design_rows.append(("Excluded by config", "none (primary arm)"))
    md.append(_kv(design_rows))
    md.append(_fig_md(report_figs.get("design")))

    md.append("### Key numbers\n")
    n_states_comp = len([s for s in states if s.get("basis", "composition") != "program"])
    n_states_prog = len(states) - n_states_comp
    kpis: list[tuple[str, Any]] = [
        ("Cells after QC", f.get("qc", {}).get("cells_after", ds.get("n_cells"))),
        ("Libraries (case / reference)", f"{ds.get('samples_per_condition', {}).get(ds.get('disease'), '?')} / {sum(v for k, v in ds.get('samples_per_condition', {}).items() if str(k) != str(ds.get('disease')))}"),
        ("Clusters", f.get("clustering", {}).get("n_clusters")),
        ("Disease-enriched clusters", len(da.get("disease_enriched_clusters", []))),
        ("States recorded (composition / programme)", f"{n_states_comp} / {n_states_prog}"),
        ("Candidates nominated (unique genes)", len(uniq_cands)),
        ("Critic flags (blocking / major / minor)", f"{sev['blocking']} / {sev['major']} / {sev['minor']}"),
        ("Library/lineage pairs flagged by the marker-free label check", len(chk.get("flagged", [])) if chk else "not run"),
        ("Tool calls / errors", f"{len(calls)} / {n_err}"),
        ("Model calls · tokens in / out", f"{len([e for e in events if e['type']=='model_call'])} · {usage.get('input_tokens', 0):,} / {usage.get('output_tokens', 0):,}"),
        ("Elapsed", f"{elapsed/60:.1f} min" if elapsed else "–"),
    ]
    ev_path = session.run_dir / "eval.json"
    if ev_path.exists():
        try:
            ev = json.loads(ev_path.read_text())
            sr = ev.get("state_recovery") or {}
            kpis.append(("Answer key: states found (lenient / strict)", f"{sum(1 for r in sr.values() if r.get('found'))} / {sum(1 for r in sr.values() if r.get('found_strict'))} of {len(sr)}"))
            kpis.append(("Answer key: expected target rank", ", ".join(f"{g} {('#' + str(v)) if v else '–'}" for g, v in (ev.get('target_rank') or {}).items())))
            kpis.append(("Answer key: decoys nominated unflagged", len(ev.get("decoys_nominated_unflagged") or [])))
        except Exception:  # noqa: BLE001
            pass
    md.append("<div class='kpi'>" + "".join(f"<div><b>{_fmt(v)}</b>{k}</div>" for k, v in kpis) + "</div>\n")
    md.append("\n" + _kv(kpis) + "\n")

    md.append("### Top candidates\n")
    if not cands:
        md.append("_No candidate targets were recorded._\n")
    else:
        rows = []
        drug = f.get("druggability", {})
        for c in cands:
            g = c["genes"][0] if c.get("genes") else "?"
            rec = next((r for r in rt.get(c.get("cluster"), []) if r["gene"] == g), {})
            d = drug.get(rec.get("human_symbol", g)) or drug.get(g) or {}
            n_inh = d.get("n_inhibitors", 0) or 0
            rows.append({
                "gene": g, "cluster": c.get("cluster"), "class": rec.get("target_class"), "attribution": rec.get("attribution"), "score": rec.get("score"),
                "compounds": rec.get("top_inhibitors", [])[:3],
                "direction": (f"{n_inh} inhibitor-type interaction(s)" if n_inh else "no inhibitor-type interaction listed") + " — human check",
                "confidence": c.get("confidence"), "critic": flag_by_gene.get(g, []),
            })
        rows.sort(key=lambda r: -(r.get("score") or 0))
        md.append(_table(rows, [("gene", "Gene"), ("cluster", "State"), ("class", "Class"), ("attribution", "Attribution"), ("score", "Score"), ("compounds", "Known compounds"), ("direction", "Therapeutic direction"), ("confidence", "Nominator confidence"), ("critic", "Critic flags")]))
        if objective == "marker":
            md.append("\n> **Objective `marker`:** these are condition-specific *response genes*, ranked by effect size, within-cluster specificity and consistency across replicates. Compounds and direction are listed for information only and did not enter the score.\n")
        md.append("\n> **Direction is a human check.** For each candidate, state whether the listed compounds inhibit or activate it and whether that is the direction this design needs. A `reference_driven` or `ambient_like` attribution means the signal is not from these cells; `shared_with_reference` means a lineage or pan-condition gene; for a programme state (`basis: program`) it means *not induced*. A candidate set aside by attribution while the per-library evidence supports it should be reported as conditional on the reference cells identified, not dropped.\n")
    md.append(_fig_md(report_figs.get("candidate_heatmap")))
    for s in summaries:
        md.append(f"\n**{s['title']}.** {s['detail']}\n")

    # ---- 1. dataset
    md.append("## 1. Dataset and QC\n")
    md.append(
        f"- {ds.get('n_cells', 0):,} cells × {ds.get('n_genes', 0):,} genes · species `{ds.get('species')}` · "
        f"conditions: {', '.join(f'{k} ({v:,})' for k, v in ds.get('condition_counts', {}).items())}\n"
        f"- samples per condition: {', '.join(f'{k}: {v}' for k, v in ds.get('samples_per_condition', {}).items())} (pseudonymized S01…)\n"
        f"- lineages: {lin_src}; clustering used: `{f.get('clustering', {}).get('cluster_key', '–')}` ({f.get('clustering', {}).get('method', '–')}, harmony={f.get('clustering', {}).get('harmony_used')})\n"
    )
    if f.get("counts_layer_used"):
        md.append(f"- raw counts: layer `{f['counts_layer_used']}`\n")
    if "qc" in f:
        q = f["qc"]
        md.append(f"- QC: {q['cells_before']:,} → {q['cells_after']:,} cells, {q['genes_before']:,} → {q['genes_after']:,} genes (min_genes={q['params']['min_genes']}, max_mt%={q['params']['max_mt_pct']})\n")
        lost = q.get("lineage_loss") or q.get("per_celltype_removed")
        if isinstance(lost, dict) and lost:
            md.append("- per-lineage QC loss: " + ", ".join(f"{k}: {_fmt(v)}" for k, v in list(lost.items())[:12]) + "\n")
    if f.get("doublets"):
        d = f["doublets"]
        md.append(f"- doublets: {d.get('method', '–')}, {d.get('n_flagged', d.get('n_doublets', '–'))} cells flagged\n")
    for name, x in fig_by_name.items():
        if name.startswith("01_umap") or "umap_sct" in name or "umap_leiden" in name:
            md.append(_fig_md(x))
            break

    # ---- 2. states
    md.append("## 2. Condition-associated cell states\n")
    if da.get("table"):
        md.append(f"Per-sample differential abundance (cluster key `{da.get('cluster_key')}`); enriched clusters: **{', '.join(da.get('disease_enriched_clusters', [])) or 'none'}**.\n\n")
        md.append(
            _table(
                sorted(da["table"], key=lambda r: -(r.get("log2_ratio") or 0)),
                [("cluster", "Cluster"), ("n_cells", "Cells"), ("mean_prop_disease", "Mean prop. disease"), ("mean_prop_reference", "Mean prop. reference"), ("log2_ratio", "log2 ratio"), ("mannwhitney_p", "MWU p"), ("n_disease_samples_present", "Disease samples present"), ("disease_enriched", "Enriched")],
            )
        )
        md.append(_fig_md(report_figs.get("abundance")))
        if not report_figs.get("abundance"):
            for name, x in fig_by_name.items():
                if "differential_abundance" in name:
                    md.append(_fig_md(x))
                    break
    dp = f.get("disease_programs", {})
    if dp.get("table"):
        md.append(f"\nExpression route (`disease_programs`, on clustering `{dp.get('cluster_key')}`): disease-altered programmes in clusters **{', '.join(dp.get('disease_altered_clusters', [])) or 'none'}**" + (f"; ambient-suspect: {', '.join(dp.get('ambient_suspect_clusters', []))}" if dp.get("ambient_suspect_clusters") else "") + ".\n")
    md.append(_fig_md(report_figs.get("state_map")))
    md.append(_fig_md(report_figs.get("state_shares")))
    for st in states:
        cl = st.get("cluster")
        md.append(f"\n### {st['title']}\n\n{st['detail']}\n\n_Basis: {st.get('basis', 'composition')}; confidence: {st.get('confidence', '–')}; dominant lineage: {comp.get(str(cl), {}).get('dominant_celltype', '–')}._\n")
        pb = f.get("pseudobulk_de", {}).get(f"{cl}:cluster_vs_reference_clusters")
        if pb:
            md.append(f"\n**Pseudobulk DE vs sibling clusters** ({pb['method']}; replicates target={pb['replicates']['target']}, reference={pb['replicates']['reference']}; {pb['n_sig_up']} up / {pb['n_sig_down']} down at padj<0.05). Top up-regulated:\n\n")
            md.append(_table(pb["top_up"][:12], [("gene", "Gene"), ("log2fc", "log2FC"), ("padj", "padj"), ("baseMean", "baseMean")]))
        pbw = f.get("pseudobulk_de", {}).get(f"{cl}:disease_vs_reference_within_cluster")
        if pbw and st.get("basis") == "program":
            md.append(f"\n**Pseudobulk DE, disease vs reference within the cluster** ({pbw['method']}; replicates target={pbw['replicates']['target']}, reference={pbw['replicates']['reference']}; {pbw['n_sig_up']} up / {pbw['n_sig_down']} down). Top up-regulated:\n\n")
            md.append(_table(pbw["top_up"][:12], [("gene", "Gene"), ("log2fc", "log2FC"), ("padj", "padj"), ("baseMean", "baseMean")]))
        mk_key = next((k for k in f.get("markers", {}) if k.startswith(f"{cl}_vs_")), None)
        if mk_key:
            md.append(f"\n**Cell-level markers** ({mk_key.replace('_vs_', ' vs ')}):\n\n")
            md.append(_table(f["markers"][mk_key][:10], [("gene", "Gene"), ("log2fc", "log2FC"), ("padj", "padj"), ("pct_in", "% in"), ("pct_out", "% out")]))
        for en in f.get("enrichment", []):
            if en.get("top"):
                md.append("\n**Enrichment (ORA):**\n\n" + _table(en["top"][:6], [("term", "Term"), ("overlap", "Overlap"), ("set_size_in_background", "Set size"), ("padj", "padj"), ("genes", "Genes")]))
                break
    md.append(_fig_md(report_figs.get("volcanoes")))

    # ---- 3. ranking
    md.append("## 3. Target ranking\n")
    if not rt:
        md.append("_No ranking was produced._\n")
    ev_rows = candidate_evidence_rows(session)
    if ev_rows:
        md.append("**Evidence summary for the nominated candidates** (one row per gene and state; full per-state tables below). `axis` says which lineage contrast carried the lineage log2FC; `drugs` counts compounds but does not know the therapeutic direction - the critic and the reader judge that.\n\n")
        md.append(_table(ev_rows, [("rank", "#"), ("gene", "Gene"), ("state", "State"), ("score", "Score"), ("class", "Class"), ("state_log2fc", "State log2FC"), ("lineage_log2fc", "Lineage log2FC"), ("axis", "Axis"), ("specificity", "log2 spec."), ("consistency", "Consistency"), ("attribution", "Attribution"), ("drugs", "Drugs"), ("top", "Top compounds"), ("flags", "Flags"), ("critic", "Critic")]))
    for x in report_figs.get("candidate_panels", []):
        md.append(_fig_md(x))
    for cl, rows in rt.items():
        basis = next((st.get("basis", "composition") for st in states if str(st.get("cluster")) == str(cl)), "composition")
        md.append(f"\n**State {cl}** ({basis}) — score = weighted sum of effect size, specificity, druggability, known compounds, sample consistency (weights in `config_used.json`; objective `{objective}`).\n\n")
        md.append(
            _table(
                rows[:15],
                [("rank", "#"), ("gene", "Gene"), ("human_symbol", "Human"), ("score", "Score"), ("target_class", "Class"), ("effect_size_log2fc", "log2FC"), ("lineage_log2fc_disease_vs_reference", "lineage log2FC"), ("log2_specificity", "log2 spec."), ("sample_consistency", "Consistency"), ("druggability_score", "Druggable"), ("known_compound_score", "Compounds"), ("top_inhibitors", "Inhibitors"), ("attribution", "Attribution"), ("flags", "Flags")],
            )
        )
    drug = f.get("druggability", {})
    if drug:
        src = {r.get("source_dgidb") for r in drug.values()} | {r.get("source_opentargets") for r in drug.values()}
        md.append(f"\n<span class='small'>Druggability sources used: {', '.join(sorted(s for s in src if s))}.</span>\n")
    lit = f.get("literature", {})
    if lit:
        md.append("\n**Literature (PubMed hit counts):** " + "; ".join(f"{g}: {v.get('count')}" for g, v in list(lit.items())[:12]) + "\n")

    # ---- 4. critic
    md.append("## 4. Critic review\n")
    if chk:
        md.append(f"Marker-free library label check: {chk.get('n_lineages_judged', 0)} of {chk.get('n_lineages_total', 0)} lineages showed condition structure and were judged; {len(chk.get('flagged', []))} library/lineage pair(s) flagged.\n")
        md.append(_fig_md(report_figs.get("library_corr")))
    if not flags:
        md.append("_No flags raised._\n")
    for s in ("blocking", "major", "minor", None):
        for fl in flags:
            if fl.get("severity") == s:
                md.append(f"- <span class='flag-{s or 'minor'}'>[{s or 'note'}]</span> **{fl['title']}** — {fl['detail']}\n")

    # ---- 5. notes
    md.append("## 5. Agent notes\n")
    for role, text in f.get("role_summaries", {}).items():
        if text:
            md.append(f"**{role}:** {text}\n\n")

    # ---- 6. evaluation (only with a key; without one the heading reads 'not scored')
    ev_md = _eval_section(session)
    md.append(ev_md if ev_md else _NO_KEY_TEXT)

    # ---- 7. provenance
    md.append("## 7. Provenance\n")
    md.append(f"{len(calls)} tool calls ({n_err} errors), {len([e for e in events if e['type']=='model_call'])} model calls. Full log: `run.jsonl`; tables: `tables/`; state: `state.json`; config: `config_used.json`.\n\n")
    md.append(_table([{"seq": e["seq"], "role": e["role"], "tool": e["tool"], "args": json.dumps(e["args"])[:90], "err": "✗" if e.get("error") else ""} for e in calls], [("seq", "#"), ("role", "Role"), ("tool", "Tool"), ("args", "Arguments"), ("err", "Error")]))

    # ---- 8. figure index
    md.append("## 8. Figure index\n")
    md.append("Each figure is written as PNG (300 dpi, shown above) and, where available, as PDF and SVG with editable text for publication.\n\n")
    md.append(_table([{"png": x["path"], "vector": ", ".join(x.get("vector") or []) or "–", "caption": x["caption"][:140]} for x in session.figures], [("png", "PNG"), ("vector", "Vector"), ("caption", "Caption")]))
    return "\n".join(md)


def build_report(session: AnalysisSession) -> tuple[Path, Path]:
    md = build_markdown(session)
    md_path = session.run_dir / "report.md"
    md_path.write_text(md)
    html_path = _write_html(session.run_dir, md)
    return md_path, html_path


def _write_html(run_dir: Path, md: str) -> Path:
    html_path = Path(run_dir) / "report.html"
    try:
        import markdown as mdlib

        body = mdlib.markdown(md, extensions=["tables", "fenced_code"])
    except Exception:  # pragma: no cover
        body = "<pre>" + md.replace("<", "&lt;") + "</pre>"
    # inline figures so the HTML is a single self-contained file that can be emailed
    import base64
    import re

    def _inline(m: "re.Match[str]") -> str:
        rel = m.group(1)
        p = Path(run_dir) / rel
        if p.exists() and p.suffix.lower() == ".png":
            b64 = base64.b64encode(p.read_bytes()).decode()
            return f'src="data:image/png;base64,{b64}"'
        return m.group(0)

    body = re.sub(r'src="(figures/[^"]+)"', _inline, body)
    html_path.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>sc-target-agent report</title><style>{CSS}</style></head><body>{body}</body></html>")
    return html_path
