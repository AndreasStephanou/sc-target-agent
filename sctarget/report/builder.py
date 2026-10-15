"""Build the Markdown + HTML report from the session facts and run log."""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path
from typing import Any

from ..session import AnalysisSession

CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1000px;margin:2rem auto;padding:0 1rem;color:#1c1c1c;line-height:1.45}
h1{border-bottom:2px solid #333;padding-bottom:.3rem}h2{margin-top:2.2rem;border-bottom:1px solid #ccc}
table{border-collapse:collapse;font-size:.9rem;margin:.6rem 0}th,td{border:1px solid #ddd;padding:.3rem .5rem;text-align:left;vertical-align:top}
th{background:#f3f3f3}code{background:#f5f5f5;padding:0 .2rem}img{max-width:100%;border:1px solid #eee;margin:.5rem 0}
.flag-blocking{color:#b00020;font-weight:600}.flag-major{color:#c75000;font-weight:600}.flag-minor{color:#666}
.small{color:#666;font-size:.85rem}
"""


def _fmt(x: Any, nd: int = 3) -> str:
    if x is None:
        return "–"
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


def build_markdown(session: AnalysisSession) -> str:
    f = session.facts
    ds = f.get("dataset", {})
    findings = f.get("findings", {})
    states = findings.get("disease_state", [])
    cands = findings.get("candidate_target", [])
    flags = findings.get("critic_flag", [])
    summaries = findings.get("summary", [])
    events = session.log.events()
    model_names = sorted({e.get("model") for e in events if e["type"] == "model_call" and e.get("model")})
    usage = session.log.usage_totals()
    fig_by_name = {Path(x["path"]).name: x for x in session.figures}

    md: list[str] = []
    md.append(f"# Target discovery report — {ds.get('disease', '?')} vs {', '.join(ds.get('reference_conditions', []))}\n")
    md.append(
        f"<span class='small'>Run `{session.run_dir.name}` · {dt.datetime.now().strftime('%Y-%m-%d %H:%M')} · model(s): {', '.join(model_names) or '–'} · prompts v{session.facts.get('prompt_version', 1)} · "
        f"tokens in/out: {usage.get('input_tokens', 0):,}/{usage.get('output_tokens', 0):,} · {len([e for e in events if e['type']=='tool_call'])} tool calls</span>\n"
    )
    md.append(
        "> **Read this first.** Every number below was produced by a deterministic tool call recorded in `run.jsonl`; "
        "the model chose which tools to run and wrote the interpretation. Candidates are hypotheses for wet-lab validation, not conclusions. "
        "Druggability values marked *fixture* come from bundled illustrative data, not a live database.\n"
    )

    # ---- dataset
    md.append("## 1. Dataset\n")
    md.append(
        f"- {ds.get('n_cells', 0):,} cells × {ds.get('n_genes', 0):,} genes · species `{ds.get('species')}` · "
        f"conditions: {', '.join(f'{k} ({v:,})' for k, v in ds.get('condition_counts', {}).items())}\n"
        f"- samples per condition: {', '.join(f'{k}: {v}' for k, v in ds.get('samples_per_condition', {}).items())} (pseudonymized S01…)\n"
        f"- existing annotations: {ds.get('celltype_key') or 'none'}; clustering used: `{f.get('clustering', {}).get('cluster_key', '–')}` ({f.get('clustering', {}).get('method', '–')}, harmony={f.get('clustering', {}).get('harmony_used')})\n"
    )
    if "qc" in f:
        q = f["qc"]
        md.append(f"- QC: {q['cells_before']:,} → {q['cells_after']:,} cells, {q['genes_before']:,} → {q['genes_after']:,} genes (min_genes={q['params']['min_genes']}, max_mt%={q['params']['max_mt_pct']})\n")
    for name, x in fig_by_name.items():
        if name.startswith("01_umap") or "umap_sct" in name or "umap_leiden" in name:
            md.append(f"\n![{x['caption']}]({x['path']})\n")
            break

    # ---- executive summary
    md.append("## 2. Summary of nominations\n")
    if not cands:
        md.append("_No candidate targets were recorded._\n")
    else:
        flag_by_gene: dict[str, list[str]] = {}
        for fl in flags:
            for g in fl.get("genes") or []:
                flag_by_gene.setdefault(g, []).append(f"{fl.get('severity', '')}: {fl['title']}")
        rows = []
        for c in cands:
            g = c["genes"][0] if c.get("genes") else "?"
            rec = next((r for r in f.get("ranked_targets", {}).get(c.get("cluster"), []) if r["gene"] == g), {})
            rows.append(
                {
                    "gene": g,
                    "cluster": c.get("cluster"),
                    "score": rec.get("score"),
                    "compounds": rec.get("top_inhibitors", [])[:3],
                    "confidence": c.get("confidence"),
                    "critic": flag_by_gene.get(g, []),
                }
            )
        md.append(_table(rows, [("gene", "Gene"), ("cluster", "State (cluster)"), ("score", "Score"), ("compounds", "Known inhibitors"), ("confidence", "Nominator confidence"), ("critic", "Critic flags")]))
    for s in summaries:
        md.append(f"\n**{s['title']}.** {s['detail']}\n")

    # ---- disease states
    md.append("## 3. Disease-specific cell states\n")
    da = f.get("differential_abundance", {})
    if da.get("table"):
        md.append(f"Per-sample differential abundance (cluster key `{da.get('cluster_key')}`); enriched clusters: **{', '.join(da.get('disease_enriched_clusters', [])) or 'none'}**.\n\n")
        md.append(
            _table(
                sorted(da["table"], key=lambda r: -(r.get("log2_ratio") or 0)),
                [("cluster", "Cluster"), ("n_cells", "Cells"), ("mean_prop_disease", "Mean prop. disease"), ("mean_prop_reference", "Mean prop. reference"), ("log2_ratio", "log2 ratio"), ("mannwhitney_p", "MWU p"), ("n_disease_samples_present", "Disease samples present"), ("disease_enriched", "Enriched")],
            )
        )
        for name, x in fig_by_name.items():
            if "differential_abundance" in name:
                md.append(f"\n![{x['caption']}]({x['path']})\n")
                break
    comp = {r["cluster"]: r for r in f.get("cluster_composition", [])}
    for st in states:
        cl = st.get("cluster")
        md.append(f"\n### {st['title']}\n\n{st['detail']}\n\n_Confidence: {st.get('confidence', '–')}; dominant cell type: {comp.get(cl, {}).get('dominant_celltype', '–')}._\n")
        pb = f.get("pseudobulk_de", {}).get(f"{cl}:cluster_vs_reference_clusters")
        if pb:
            md.append(f"\n**Pseudobulk DE vs sibling clusters** ({pb['method']}; replicates target={pb['replicates']['target']}, reference={pb['replicates']['reference']}; {pb['n_sig_up']} up / {pb['n_sig_down']} down at padj<0.05). Top up-regulated:\n\n")
            md.append(_table(pb["top_up"][:12], [("gene", "Gene"), ("log2fc", "log2FC"), ("padj", "padj"), ("baseMean", "baseMean")]))
        mk_key = next((k for k in f.get("markers", {}) if k.startswith(f"{cl}_vs_")), None)
        if mk_key:
            md.append(f"\n**Cell-level markers** ({mk_key.replace('_vs_', ' vs ')}):\n\n")
            md.append(_table(f["markers"][mk_key][:10], [("gene", "Gene"), ("log2fc", "log2FC"), ("padj", "padj"), ("pct_in", "% in"), ("pct_out", "% out")]))
        for en in f.get("enrichment", []):
            if en.get("top"):
                md.append("\n**Enrichment (ORA):**\n\n" + _table(en["top"][:6], [("term", "Term"), ("overlap", "Overlap"), ("set_size_in_background", "Set size"), ("padj", "padj"), ("genes", "Genes")]))
                break
        for x in session.figures:
            n = Path(x["path"]).name
            if "umap_genes" in n or "dotplot" in n:
                md.append(f"\n![{x['caption']}]({x['path']})\n")

    # ---- ranking
    md.append("## 4. Target ranking\n")
    rt = f.get("ranked_targets", {})
    if not rt:
        md.append("_No ranking was produced._\n")
    for cl, rows in rt.items():
        md.append(f"\n**Cluster {cl}** — score = weighted sum of effect size, specificity, druggability, known compounds, sample consistency (weights in `config_used.json`).\n\n")
        md.append(
            _table(
                rows[:15],
                [("rank", "#"), ("gene", "Gene"), ("human_symbol", "Human"), ("score", "Score"), ("effect_size_log2fc", "log2FC"), ("log2_specificity", "log2 spec."), ("sample_consistency", "Consistency"), ("druggability_score", "Druggable"), ("known_compound_score", "Compounds"), ("top_inhibitors", "Inhibitors"), ("flags", "Flags")],
            )
        )
    drug = f.get("druggability", {})
    if drug:
        src = {r.get("source_dgidb") for r in drug.values()} | {r.get("source_opentargets") for r in drug.values()}
        md.append(f"\n<span class='small'>Druggability sources used: {', '.join(sorted(s for s in src if s))}.</span>\n")
    lit = f.get("literature", {})
    if lit:
        md.append("\n**Literature (PubMed hit counts):** " + "; ".join(f"{g}: {v.get('count')}" for g, v in list(lit.items())[:12]) + "\n")

    # ---- critic
    md.append("## 5. Critic review\n")
    if not flags:
        md.append("_No flags raised._\n")
    for sev in ("blocking", "major", "minor", None):
        for fl in flags:
            if fl.get("severity") == sev:
                md.append(f"- <span class='flag-{sev or 'minor'}'>[{sev or 'note'}]</span> **{fl['title']}** — {fl['detail']}\n")

    # ---- role notes
    md.append("## 6. Agent notes\n")
    for role, text in f.get("role_summaries", {}).items():
        if text:
            md.append(f"**{role}:** {text}\n\n")

    # ---- provenance
    md.append("## 7. Provenance\n")
    calls = [e for e in events if e["type"] == "tool_call"]
    md.append(f"{len(calls)} tool calls, {len([e for e in events if e['type']=='model_call'])} model calls. Full log: `run.jsonl`; tables: `tables/`; state: `state.json`.\n\n")
    md.append(_table([{"seq": e["seq"], "role": e["role"], "tool": e["tool"], "args": json.dumps(e["args"])[:90], "err": "✗" if e.get("error") else ""} for e in calls], [("seq", "#"), ("role", "Role"), ("tool", "Tool"), ("args", "Arguments"), ("err", "Error")]))
    return "\n".join(md)


def build_report(session: AnalysisSession) -> tuple[Path, Path]:
    md = build_markdown(session)
    md_path = session.run_dir / "report.md"
    md_path.write_text(md)
    html_path = session.run_dir / "report.html"
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
        p = session.run_dir / rel
        if p.exists() and p.suffix.lower() == ".png":
            b64 = base64.b64encode(p.read_bytes()).decode()
            return f'src="data:image/png;base64,{b64}"'
        return m.group(0)

    body = re.sub(r'src="(figures/[^"]+)"', _inline, body)
    html_path.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>Target discovery report</title><style>{CSS}</style></head><body>{body}</body></html>")
    return md_path, html_path
