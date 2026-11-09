"""Model backends behind one small interface.

``AnthropicModel`` wraps the Anthropic Python SDK's Messages API with tool use.
``ScriptedModel`` is a deterministic test double that issues a fixed, sensible sequence
of tool calls per role so the whole pipeline can run in CI and offline without an API
key. It exists to test the machinery (loop, tools, privacy, report), not to stand in for
the reasoning of a real model.
"""
from __future__ import annotations

import json

import numpy as np
import os
import uuid
from dataclasses import dataclass, field
from typing import Any, Protocol


@dataclass
class ToolUse:
    name: str
    input: dict[str, Any]
    id: str = field(default_factory=lambda: f"toolu_{uuid.uuid4().hex[:16]}")


@dataclass
class ModelResponse:
    text: str
    tool_uses: list[ToolUse]
    stop_reason: str
    usage: dict[str, int]
    raw_content: list[dict[str, Any]]  # assistant content blocks to append to the conversation


class Model(Protocol):
    name: str

    def create(self, *, system: str, messages: list[dict[str, Any]], tools: list[dict[str, Any]], role: str) -> ModelResponse: ...


# --------------------------------------------------------------------------- Anthropic
class AnthropicModel:
    def __init__(self, model: str, temperature: float = 0.2, max_tokens: int = 4096, api_key: str | None = None) -> None:
        import anthropic

        self.name = model
        self.temperature = temperature
        self.max_tokens = max_tokens
        self.client = anthropic.Anthropic(api_key=api_key or os.getenv("ANTHROPIC_API_KEY"))

    def create(self, *, system: str, messages: list[dict[str, Any]], tools: list[dict[str, Any]], role: str) -> ModelResponse:
        resp = self.client.messages.create(
            model=self.name,
            system=system,
            messages=messages,
            tools=tools,
            max_tokens=self.max_tokens,
            temperature=self.temperature,
        )
        text_parts, tool_uses, raw = [], [], []
        for block in resp.content:
            raw.append(block.model_dump(exclude_none=True))
            if block.type == "text":
                text_parts.append(block.text)
            elif block.type == "tool_use":
                tool_uses.append(ToolUse(name=block.name, input=dict(block.input), id=block.id))
        usage = {"input_tokens": resp.usage.input_tokens, "output_tokens": resp.usage.output_tokens}
        return ModelResponse(text="\n".join(text_parts), tool_uses=tool_uses, stop_reason=resp.stop_reason or "end_turn", usage=usage, raw_content=raw)


# --------------------------------------------------------------------------- Scripted
class ScriptedModel:
    """Deterministic stand-in. Reads prior tool results from the conversation and the
    session facts, then issues the next tool call of a fixed per-role plan."""

    name = "scripted"

    def __init__(self, session: Any) -> None:
        self.session = session
        self._step: dict[str, Any] = {}
        self._queue: dict[str, list[ToolUse | str]] = {}

    # ---- helpers -----------------------------------------------------------
    def _facts(self) -> dict[str, Any]:
        return self.session.facts

    def _plan(self, role: str) -> list[ToolUse | str]:
        f = self._facts()
        cfg = self.session.cfg
        if role == "planner":
            return [
                ToolUse("load_dataset", {}),
                ToolUse(
                    "record_finding",
                    {
                        "kind": "summary",
                        "title": "Analysis plan",
                        "detail": "Reuse the existing clustering if the object is processed, otherwise QC + integrate. Find disease-enriched clusters by per-sample differential abundance, characterise them against sibling clusters of the same cell type, then nominate druggable, state-specific genes and have the critic test for single-sample and pan-condition artefacts.",
                    },
                ),
                "Plan recorded.",
            ]
        if role == "analyst":
            return [ToolUse("load_dataset", {}), "__analyst_dynamic__"]
        if role == "target":
            return ["__target_dynamic__"]
        if role == "critic":
            return ["__critic_dynamic__"]
        return ["Nothing to do."]

    # ---- dynamic sub-plans (filled once the needed facts exist) -------------
    def _analyst_dynamic(self) -> list[ToolUse | str]:
        f = self._facts()
        ds = f.get("dataset", {})
        cfg = self.session.cfg
        steps: list[ToolUse | str] = []
        processed = bool(ds.get("already_processed_flag"))
        if not processed:
            steps.append(ToolUse("run_qc", dict(cfg.analysis.get("qc", {}))))
            steps.append(ToolUse("detect_doublets", {}))
        reuse = ds["existing_clusterings"][0] if processed and ds.get("existing_clusterings") else None
        integ = cfg.analysis.get("integration", {})
        steps.append(ToolUse("integrate_and_cluster", {"reuse_existing": reuse, "use_harmony": bool(integ.get("use_harmony", True)), "resolution": float(integ.get("resolution", 1.0)), "n_pcs": int(integ.get("n_pcs", 30))}))
        steps.append(ToolUse("cluster_composition", {}))
        steps.append(ToolUse("differential_abundance", {"enrichment_log2fc": float(cfg.get("analysis.differential_abundance.enrichment_log2fc", 1.0))}))
        steps.append("__analyst_characterise__")
        return steps

    def _analyst_characterise(self) -> list[ToolUse | str]:
        f = self._facts()
        da = f.get("differential_abundance", {})
        comp = {r["cluster"]: r for r in f.get("cluster_composition", [])}
        rows = {r["cluster"]: r for r in da.get("table", [])}
        steps: list[ToolUse | str] = []
        # split states: several disease-dominated near-miss clusters of one lineage -> test their union once
        split = [str(c) for c in da.get("split_state_candidates", [])]
        if len(split) >= 2 and not self._step.get("merged"):
            by_lin: dict[str, list[str]] = {}
            for c in split:
                by_lin.setdefault(str(comp.get(c, {}).get("dominant_celltype") or "?"), []).append(c)
            groups = [g for g in by_lin.values() if len(g) >= 2]
            if groups:
                group = max(groups, key=lambda g: sum(comp.get(c, {}).get("n_cells", 0) for c in g))
                self._step["merged"] = True
                steps.append(ToolUse("record_finding", {"kind": "summary", "title": f"Merging split clusters {', '.join(group)}", "detail": "These clusters are each >= 80% disease cells and present in most disease samples but none passes the per-sample enrichment test on its own (bimodal across samples): the pattern of one state partitioned by patient. Testing their union.", "confidence": "medium"}))
                steps.append(ToolUse("merge_clusters", {"clusters": group}))
                steps.append(ToolUse("cluster_composition", {}))
                steps.append("__analyst_characterise__")
                return steps
        # size-aware ordering: effect x log10(cells), so a large consistent state outranks a tiny extreme one;
        # dissociation-stress suspects are skipped
        stress = set(map(str, f.get("stress_suspect_clusters", {}).keys()))
        stress |= {str(c) for c, r in comp.items() if r.get("doublet_suspect")}
        def prio(c):
            r = rows.get(c, {})
            return -((r.get("log2_ratio") or 0) * np.log10(max(10, r.get("n_cells") or 10)))
        def prio2(c):
            return (0 if "+" in str(c) else 1, prio(c))  # a merged cluster was built on purpose: characterise it first
        ranked_enriched = [c for c in sorted(da.get("disease_enriched_clusters", []), key=prio2) if str(c) not in stress]
        # lineage coverage (protocol 23e): in a multi-lineage disease (IPF: epithelium, macrophages,
        # fibroblasts all shift) the top enriched cluster of EACH lineage is characterised before a
        # second cluster of the same lineage; the old "top 3 overall" filled the cap with one lineage
        # and left the macrophage states without a marker table
        seen_lin: set[str] = set()
        first_per_lin: list = []
        rest: list = []
        for c in ranked_enriched:
            lin = str(comp.get(str(c), {}).get("dominant_celltype") or "?")
            (rest if lin in seen_lin else first_per_lin).append(c)
            seen_lin.add(lin)
        enriched = (first_per_lin + rest)[:5]
        if not enriched:
            retries = self._step.get("analyst_recluster", 0)
            da_thr = float(self.session.cfg.get("analysis.differential_abundance.enrichment_log2fc", 1.0))
            integ = self.session.cfg.analysis.get("integration", {})
            if retries < 2:
                # adaptive step 1/2: a rare disease state can hide inside a large, condition-balanced
                # cell type. Subcluster one cell type at a time, largest first (most room to hide).
                tried: set[str] = self._step.setdefault("subclustered_types", set())
                base_comp = {r["cluster"]: r for r in f.get("cluster_composition", []) if "." not in str(r["cluster"])}
                by_type: dict[str, list[str]] = {}
                for c, r in base_comp.items():
                    dom = r.get("dominant_celltype") or f"cluster_{c}"
                    by_type.setdefault(dom, []).append(c)
                sizes = {t: sum(base_comp[c]["n_cells"] for c in cs) for t, cs in by_type.items()}
                order = [t for t in sorted(sizes, key=lambda t: -sizes[t]) if t not in tried]
                if order:
                    self._step["analyst_recluster"] = retries + 1
                    dom = order[0]
                    tried.add(dom)
                    group = by_type[dom]
                    steps.append(ToolUse("record_finding", {"kind": "summary", "title": f"Subclustering {dom} ({', '.join(group)})", "detail": f"No cluster was disease-enriched; subclustering the {dom} cells (largest cell type not yet examined) to look for a rarer disease-specific state hidden inside a condition-balanced cluster.", "confidence": "medium"}))
                    base_key = str(f.get("clustering", {}).get("cluster_key", "")).split("_sub")[0] or None
                    steps.append(ToolUse("subcluster", {"clusters": group, "cluster_key": base_key, "resolution": 1.0}))
                    steps.append(ToolUse("cluster_composition", {}))
                    steps.append(ToolUse("differential_abundance", {"enrichment_log2fc": da_thr}))
                    steps.append("__analyst_characterise__")
                    return steps
            if retries <= 2:
                # adaptive step 3: global re-clustering at higher resolution
                self._step["analyst_recluster"] = 3
                steps.append(ToolUse("record_finding", {"kind": "summary", "title": "Re-clustering at resolution 2.0", "detail": "Subclustering did not reveal a disease-enriched state; recomputing the global clustering at higher resolution.", "confidence": "medium"}))
                steps.append(ToolUse("integrate_and_cluster", {"reuse_existing": None, "use_harmony": bool(integ.get("use_harmony", True)), "resolution": 2.0, "n_pcs": int(integ.get("n_pcs", 30))}))
                steps.append(ToolUse("cluster_composition", {}))
                steps.append(ToolUse("differential_abundance", {"enrichment_log2fc": da_thr}))
                steps.append("__analyst_characterise__")
                return steps
            # adaptive step 4: composition is flat -> look for disease PROGRAMS inside shared clusters
            if not self._step.get("programs_tried"):
                self._step["programs_tried"] = True
                steps.append(ToolUse("record_finding", {"kind": "summary", "title": "Switching to the expression route", "detail": "No cluster is enriched in the disease condition after subclustering and re-clustering. Testing whether shared cell types run a disease-specific program instead (pseudobulk disease vs reference within each cluster).", "confidence": "medium"}))
                steps.append(ToolUse("disease_programs", {"cluster_key": None}))  # tool defaults to the base clustering (protocol 23)
                steps.append("__analyst_programs__")
                return steps
            steps.append(ToolUse("record_finding", {"kind": "summary", "title": "No disease-specific state found", "detail": "Neither the composition route (differential abundance, subclustering, re-clustering) nor the expression route (within-cluster disease programs) found a disease-specific state at the configured thresholds.", "confidence": "medium"}))
            steps.append("No disease-specific states found.")
            return steps
        # composition found something; ALSO test for disease PROGRAMS in shared cell types first
        # (on the base clustering, so subclusters do not fragment the test) - cheap, and it catches
        # the case where the composition hits are weak and the real signal is a program shift
        if not self._step.get("programs_tried"):
            self._step["programs_tried"] = True
            steps.append(ToolUse("disease_programs", {"cluster_key": None}))  # tool defaults to the base clustering (protocol 23)
            steps.append("__analyst_programs__")
        for cl in enriched:
            dom = comp.get(cl, {}).get("dominant_celltype")
            siblings = [c for c, r in comp.items() if c != cl and dom and r.get("dominant_celltype") == dom] or None
            steps.append(ToolUse("find_markers", {"cluster": cl, "reference_clusters": siblings, "n_top": 30}))
            steps.append(ToolUse("pseudobulk_de", {"cluster": cl, "reference_clusters": siblings, "n_top": 40}))
            # the disease-wide axis: disease vs reference cells within the lineage's shared clusters.
            # When a lineage holds several disease sub-states, the sibling contrast only sees what
            # separates them from each other; this contrast sees what separates disease from normal.
            steps.append(ToolUse("pseudobulk_de", {"cluster": cl, "contrast": "disease_vs_reference_within_compartment", "n_top": 40}))
            if int(self.session.cfg.targets.get("scoring_version", 1)) >= 3:
                # scoring v3: the plain lineage contrast (all disease cells of the lineage vs all reference
                # cells). In tumours most disease cells sit in the enriched clusters, so the shared-cluster
                # contrast above compares an atypical minority; this one compares tumour with normal.
                steps.append(ToolUse("pseudobulk_de", {"cluster": cl, "contrast": "disease_vs_reference_within_lineage", "n_top": 40}))
            steps.append(f"__analyst_enrich__{cl}")
        return steps

    def _analyst_enrich(self, cl: str) -> list[ToolUse | str]:
        f = self._facts()
        stress = f.get("stress_suspect_clusters", {})
        if str(cl) in stress:
            return [ToolUse("record_finding", {"kind": "summary", "title": f"Cluster {cl} not recorded: dissociation-stress artefact", "detail": f"{int(100*stress[str(cl)])}% of its top 30 markers are immediate-early / heat-shock genes (FOS, JUN, HSPA1A, DNAJB1 ...), the dissociation-induced programme; its enrichment in disease samples reflects sample handling, not biology.", "cluster": cl, "confidence": "medium"})]
        pb = f.get("pseudobulk_de", {}).get(f"{cl}:cluster_vs_reference_clusters", {})
        from ..tools.targets import excluded_gene

        up = [r["gene"] for r in pb.get("top_up", []) if not excluded_gene(r["gene"])][:60]
        if not up:
            mk_key = next((k for k in f.get("markers", {}) if k.startswith(f"{cl}_vs_")), None)
            up = [r["gene"] for r in f.get("markers", {}).get(mk_key, []) if not excluded_gene(r["gene"])][:60]
        cp = f.get("pseudobulk_de", {}).get(f"{cl}:disease_vs_reference_within_compartment", {})
        cp_genes = [r["gene"] for r in cp.get("top_up", []) if not excluded_gene(r["gene"])][:60]
        ln = f.get("pseudobulk_de", {}).get(f"{cl}:disease_vs_reference_within_lineage", {})
        ln_genes = [r["gene"] for r in ln.get("top_up", []) if not excluded_gene(r["gene"])][:60]
        if cp_genes or ln_genes:  # recorded signature: what defines the state among its siblings AND what separates it from normal
            up = list(dict.fromkeys(up[:8] + ln_genes[:8] + cp_genes[:8] + up[8:] + ln_genes[8:] + cp_genes[8:]))
        cycling = f.get("cycling_clusters", {}).get(str(cl))
        da_rows = {r["cluster"]: r for r in f.get("differential_abundance", {}).get("table", [])}
        row = da_rows.get(cl, {})
        comp = {r["cluster"]: r for r in f.get("cluster_composition", [])}
        dom = comp.get(cl, {}).get("dominant_celltype", "unknown")
        steps: list[ToolUse | str] = [
            ToolUse("enrichment", {"genes": up}),
            ToolUse("plot_genes", {"genes": up[:4], "title": f"cluster {cl} top genes"}),
            ToolUse(
                "record_finding",
                {
                    "kind": "disease_state",
                    "title": f"Cluster {cl}: disease-enriched {dom} state",
                    "detail": f"Cluster {cl} has log2 abundance ratio {row.get('log2_ratio')} (disease vs reference, per-sample), present in {row.get('n_disease_samples_present')}/{row.get('n_disease_samples_total')} disease samples; {round(100*row.get('frac_cells_from_disease',0))}% of its cells come from the disease condition. Pseudobulk DE vs sibling clusters yields {pb.get('n_sig_up')} significantly up-regulated genes (method: {pb.get('method')}){'; disease vs reference within the shared clusters of the lineage yields ' + str(cp.get('n_sig_up')) + ' (top: ' + ', '.join(cp_genes[:6]) + ')' if cp_genes else ''}{'; disease vs reference across the whole lineage yields ' + str(ln.get('n_sig_up')) + ' (top: ' + ', '.join(ln_genes[:6]) + ')' if ln_genes else ''}.{' CYCLING SUB-STATE: ' + str(int(100*cycling)) + '% of its top markers are cell-cycle genes; they describe proliferation, not disease identity.' if cycling else ''} Top genes: {', '.join(up[:8])}.",
                    "cluster": cl,
                    "genes": up[:15],
                    "confidence": "high" if (row.get("n_disease_samples_present", 0) >= 2 and (row.get("log2_ratio") or 0) >= 2) else "medium",
                },
            ),
        ]
        return steps

    def _analyst_programs(self) -> list[ToolUse | str]:
        f = self._facts()
        dp = f.get("disease_programs", {})
        altered = dp.get("disease_altered_clusters", [])[:2]
        rows = {r["cluster"]: r for r in dp.get("table", [])}
        if not altered:
            if f.get("differential_abundance", {}).get("disease_enriched_clusters"):
                return []  # composition states are being characterised; nothing extra to add
            return ["__analyst_characterise__"]  # falls through to the final "nothing found" summary
        from ..tools.targets import excluded_gene

        steps: list[ToolUse | str] = []
        for cl in altered:
            r = rows[cl]
            pb = f.get("pseudobulk_de", {}).get(f"{cl}:disease_vs_reference_within_cluster", {})
            up = [x["gene"] for x in pb.get("top_up", []) if not excluded_gene(x["gene"])][:60]
            steps.append(ToolUse("enrichment", {"genes": up}))
            steps.append(ToolUse("plot_genes", {"genes": up[:4], "title": f"cluster {cl} disease program"}))
            steps.append(
                ToolUse(
                    "record_finding",
                    {
                        "kind": "disease_state",
                        "basis": "program",
                        "title": f"Cluster {cl}: disease-altered {r.get('dominant_celltype') or 'cell'} program",
                        "detail": f"Cluster {cl} is present in all conditions (no abundance shift) but disease samples show {r.get('n_sig_up_strong')} genes up (padj<0.05, log2FC>=1) vs reference samples in a replicate-aware pseudobulk test ({r.get('replicates_disease')} disease vs {r.get('replicates_reference')} reference pseudobulks). Top genes: {', '.join(up[:8])}.",
                        "cluster": cl,
                        "genes": up[:15],
                        "confidence": "high" if (r.get("n_sig_up_strong") or 0) >= 100 else "medium",
                    },
                )
            )
        if not f.get("differential_abundance", {}).get("disease_enriched_clusters"):
            steps.append("Recorded disease-altered programs.")
        return steps

    def _target_dynamic(self) -> list[ToolUse | str]:
        f = self._facts()
        states = f.get("findings", {}).get("disease_state", [])
        if not states:
            return ["No disease states were recorded by the analyst; nothing to nominate."]
        # composition states carry the stronger evidence (a population that exists only in disease);
        # score them first so a turn budget never silences them in favour of program states
        states = sorted(states, key=lambda st: (0 if st.get("basis", "composition") == "composition" else 1))
        self._step["target_order"] = [st["cluster"] for st in states[:6]]
        steps: list[ToolUse | str] = []
        species = self.session.cfg.dataset.get("species", "homo_sapiens")
        for st in states[:6]:
            cl = st["cluster"]
            contrast = "disease_vs_reference_within_cluster" if st.get("basis") == "program" else "cluster_vs_reference_clusters"
            pb = f.get("pseudobulk_de", {}).get(f"{cl}:{contrast}", {})
            from ..tools.targets import excluded_gene

            def _pool(rows, n_stat=20, n_effect=20, min_log2fc=1.0):
                # the n most significant, plus the n_effect largest significant effects among every
                # other stored gene (protocol 23f, boundary fixed in 23h: the effect set is taken over
                # the whole stored pool, not only over the rows the tool labelled 'effect', so a gene
                # at significance rank 21-40 with a top effect is kept)
                rows = [r for r in rows if not excluded_gene(r["gene"])]
                head = [r["gene"] for r in rows if r.get("pool") != "effect"][:n_stat]
                rest = [r for r in rows if r["gene"] not in head and (r.get("log2fc") or 0) >= min_log2fc]
                rest = sorted(rest, key=lambda r: -(r.get("log2fc") or 0))[:n_effect]
                return head + [r["gene"] for r in rest]

            pb_genes = _pool(pb.get("top_up", [])) or st.get("genes", [])[:20]
            # cell-level markers are kept as well: for a pure state they often carry the biology
            # the pseudobulk list misses (e.g. receptors expressed in 70% of cells at modest counts)
            mk_key = next((k for k in f.get("markers", {}) if k.startswith(f"{cl}_vs_")), None)
            mk_genes = [r["gene"] for r in f.get("markers", {}).get(mk_key, []) if not excluded_gene(r["gene"])][:15]
            cp = f.get("pseudobulk_de", {}).get(f"{cl}:disease_vs_reference_within_compartment", {})
            cp_genes = _pool(cp.get("top_up", []))
            ln = f.get("pseudobulk_de", {}).get(f"{cl}:disease_vs_reference_within_lineage", {})
            ln_genes = _pool(ln.get("top_up", []))
            genes = list(dict.fromkeys(pb_genes + mk_genes + cp_genes + ln_genes))  # markers and disease-wide genes never truncated away
            steps.append(ToolUse("gene_specificity", {"genes": genes, "cluster": cl}))
            steps.append(ToolUse("condition_attribution", {"genes": genes, "cluster": cl}))
            if species not in ("homo_sapiens", "human"):
                steps.append(ToolUse("map_orthologs", {"genes": genes}))
            steps.append(f"__target_lookup__{cl}")
        return steps

    def _target_lookup(self, cl: str) -> list[ToolUse | str]:
        f = self._facts()
        spec = f.get("gene_specificity", {}).get(cl, [])
        genes = [r["gene"] for r in spec]
        orth = f.get("orthologs", {})
        human = [orth.get(g, {}).get("human_symbol") or g for g in genes]
        human = [h for h in dict.fromkeys(human) if h]
        disease = self.session.disease
        return [
            ToolUse("druggability_lookup", {"genes": human}),
            ToolUse("literature_count", {"genes": human[:10], "disease_terms": [disease, "ameloblastoma"] if disease.upper() == "CAA" else [disease]}),
            ToolUse("score_targets", {"genes": genes, "cluster": cl}),
            f"__target_record__{cl}",
        ]

    def _target_record(self, cl: str) -> list[ToolUse | str]:
        f = self._facts()
        ranked = f.get("ranked_targets", {}).get(cl, [])
        steps: list[ToolUse | str] = []
        # a gene already nominated from another state (two clusters of one lineage sharing a programme)
        # is nominated once; the report lists the states it came from (protocol 23d)
        already: set[str] = self._step.setdefault("nominated_genes", set())
        for r in ranked[:5]:
            if r["druggability_score"] == 0 and r["known_compound_score"] == 0:
                continue
            if r["score"] < 0.3 or r.get("attribution") in ("reference_driven", "ambient_like"):
                continue  # weak, heavily penalised, or a wrong-condition artefact: not worth a nomination
            if r["gene"] in already:
                continue
            already.add(r["gene"])
            steps.append(
                ToolUse(
                    "record_finding",
                    {
                        "kind": "candidate_target",
                        "title": f"{r['gene']} (rank {r['rank']}, score {r['score']})",
                        "detail": f"log2FC {r['effect_size_log2fc']} ({r['effect_size_source']}), log2 specificity {r['log2_specificity']}, sample consistency {r['sample_consistency']}, druggability {r['druggability_score']}, known compounds {r['known_compound_score']} ({', '.join(r['top_inhibitors'][:3]) or 'none'}). Flags: {', '.join(r['flags']) or 'none'}.",
                        "cluster": cl,
                        "genes": [r["gene"]],
                        "confidence": "high" if r["score"] >= 0.75 and not r["flags"] else ("medium" if r["score"] >= 0.5 else "low"),
                    },
                )
            )
        # a plain string ends the role, so only say "done" after the LAST state has been processed
        states = self._step.get("target_order") or [st["cluster"] for st in f.get("findings", {}).get("disease_state", [])][:6]
        if cl == states[-1]:
            n_total = len(f.get("findings", {}).get("candidate_target", [])) + len([s for s in steps if isinstance(s, ToolUse)])
            steps.append(f"Nominated {n_total} candidate targets across clusters {', '.join(states)}.")
        return steps

    def _critic_dynamic(self) -> list[ToolUse | str]:
        f = self._facts()
        steps: list[ToolUse | str] = [ToolUse("differential_abundance", {"within_celltype": True})]
        # marker-free label check: does any library resemble the opposite condition within its lineage?
        steps.append(ToolUse("check_reference_libraries", {}))
        steps.append("__critic_review__")
        return steps

    EXCLUSIVE_LINEAGE_SHARE = 0.9  # protocol 23p: a lineage with >= 90% disease cells is condition-exclusive

    def _lineage_disease_share(self, f: dict[str, Any]) -> dict[str, tuple[str, float]]:
        """cluster -> (dominant lineage, share of disease cells in that whole lineage), from cluster_composition."""
        rows = f.get("cluster_composition") or []
        if isinstance(rows, dict):
            rows = rows.get("table", [])
        disease = getattr(getattr(self, "session", None), "disease", None)
        key = f"frac_{disease}" if disease else None
        if rows and (not key or key not in rows[0]):
            cands = [k for k in rows[0] if k.startswith("frac_") and k != "frac_cells_from_disease"]
            key = key if key in (rows[0] if rows else {}) else (cands[0] if len(cands) == 1 else key)
        tot: dict[str, float] = {}
        dis: dict[str, float] = {}
        for r in rows:
            lin = str(r.get("dominant_celltype"))
            n = float(r.get("n_cells") or 0)
            fr = r.get(key) if key else None
            if fr is None:
                continue
            tot[lin] = tot.get(lin, 0.0) + n
            dis[lin] = dis.get(lin, 0.0) + n * float(fr)
        out: dict[str, tuple[str, float]] = {}
        for r in rows:
            lin = str(r.get("dominant_celltype"))
            if tot.get(lin):
                out[str(r.get("cluster"))] = (lin, dis[lin] / tot[lin])
        return out

    def _critic_review(self) -> list[ToolUse | str]:
        f = self._facts()
        steps: list[ToolUse | str] = []
        da_rows = {r["cluster"]: r for r in f.get("differential_abundance", {}).get("table", [])}
        da_within = {r["cluster"]: r for r in f.get("differential_abundance_within_celltype", {}).get("table", [])}
        for fl in f.get("reference_library_check", {}).get("flagged", []):
            steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"Library {fl['library']} (labelled {fl['label']}) resembles the opposite condition in lineage {fl['lineage']}",
                                                     "detail": f"Its {fl['n_cells']} {fl['lineage']} cells correlate more with the opposite-condition libraries (mean r {fl['mean_r_opposite_label']}) than with its own label (mean r {fl['mean_r_same_label']}); {int(100*fl['frac_nearest_opposite_label'])}% of its nearest libraries carry the other label. Either the label is wrong or the sample is contaminated with the other condition. Every disease-vs-reference contrast in this lineage is weakened by it. Verify the label; if unresolved, run a second arm with dataset.exclude_samples and a recorded reason.",
                                                     "severity": "major"}))
        lineage_share = self._lineage_disease_share(f)
        exclusive_noted: set[str] = set()
        for st in f.get("findings", {}).get("disease_state", []):
            cl = st["cluster"]
            row = da_rows.get(cl, {})
            wrow = da_within.get(cl, {})
            if st.get("basis") == "program":
                continue  # program states are judged by within-cluster DE, not by abundance
            if wrow and not wrow.get("disease_enriched"):
                lin, share = lineage_share.get(str(cl), (None, None))
                if share is not None and share >= self.EXCLUSIVE_LINEAGE_SHARE:
                    # Protocol 23p: the lineage is itself condition-exclusive (derived lineages do
                    # this when the disease cells are distinct enough to form their own coarse
                    # group). The within-lineage test then compares disease sub-states with each
                    # other and is void, not failed; the state is the lineage.
                    if lin not in exclusive_noted:
                        exclusive_noted.add(lin)
                        steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"Lineage {lin} is condition-exclusive ({int(round(100 * share))}% disease cells): the within-lineage test is void for its clusters", "detail": f"Cluster {cl} (and every other cluster of lineage {lin}) is almost entirely disease cells, so recomputing proportions within the lineage compares disease sub-states with each other (log2 ratio {wrow.get('log2_ratio')}; all-cells ratio {row.get('log2_ratio')}). That is not a failed enrichment: the lineage itself is the disease-associated population and its composition test is the all-cells one. Consequences: the sibling contrast within this lineage is blind to genes shared by all its clusters (pan-state genes cancel), so read the disease-wide axis (within-lineage and compartment contrasts) for those; the reference cells inside the lineage are few, so attribution classes are decided on a thin reference pool.", "cluster": cl, "severity": "minor"}))
                else:
                    steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"Cluster {cl} enrichment does not hold within its parent cell type", "detail": f"Recomputing proportions within the dominant cell type gives log2 ratio {wrow.get('log2_ratio')} (all-cells ratio was {row.get('log2_ratio')}); the enrichment may reflect a composition shift rather than a disease-specific state.", "cluster": cl, "severity": "major"}))
            if row and row.get("n_disease_samples_present", 0) < 2:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"Cluster {cl} is driven by a single disease sample", "detail": "Only one disease sample contributes cells; this cannot be distinguished from a sample-specific artefact.", "cluster": cl, "severity": "blocking"}))
            if row and row.get("n_cells", 0) < 50:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"Cluster {cl} is small ({row.get('n_cells')} cells)", "detail": "Marker and DE statistics on fewer than 50 cells are unstable; confirm at a different clustering resolution.", "cluster": cl, "severity": "minor"}))
        for ct in f.get("findings", {}).get("candidate_target", []):
            g = ct["genes"][0] if ct.get("genes") else None
            cl = ct.get("cluster")
            rec = next((r for r in f.get("ranked_targets", {}).get(cl, []) if r["gene"] == g), None)
            if not rec:
                continue
            if "low_specificity" in rec["flags"] or "shared_across_conditions" in rec["flags"]:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: not specific to the disease state", "detail": f"log2 specificity {rec['log2_specificity']} - expression is comparable in other clusters/conditions, so {g} looks like a cell-type or pan-tumour marker rather than a disease-state target.", "cluster": cl, "genes": [g], "severity": "major"}))
            if "rare_reference_counterpart" in rec["flags"]:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: expressed by a rare reference population inside the state", "detail": f"A small number of reference-condition cells (< 0.5% of all reference cells) sit inside cluster {cl} and express {g}. Too few to stand for the reference condition, so they do not change the class, but they say {g} is expressed by some normal cell type present at low frequency: an on-target consideration for any therapy, to be identified (which lineage, which libraries) before validation.", "cluster": cl, "genes": [g], "severity": "minor"}))
            if "proliferation_marker" in rec["flags"]:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: cell-cycle gene, not a disease-specific target", "detail": f"{g} is a proliferation (S / G2-M) programme gene. It marks that the cells of cluster {cl} divide, which is true of every tumour and of regenerating tissue; it says nothing about what makes this state disease-specific and would score in any cancer dataset. Treat the cycling sub-state's identity genes (lineage contrast) as the candidates instead.", "cluster": cl, "genes": [g], "severity": "major"}))
            if "inconsistent_across_samples" in rec["flags"]:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: inconsistent across disease samples", "detail": f"Detected in >=10% of state cells in only {rec['sample_consistency']} of disease samples.", "cluster": cl, "genes": [g], "severity": "major"}))
            if rec.get("attribution") == "reference_driven":
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: expression sits in the reference-condition cells of cluster {cl}", "detail": f"condition_attribution classifies {g} as reference_driven: within cluster {cl} it is expressed mainly by {', '.join(self.session.reference)} cells, so its cluster-vs-sibling fold change is an artefact of mixed-condition cells in the cluster, not a disease signal.", "cluster": cl, "genes": [g], "severity": "blocking"}))
            elif rec.get("attribution") == "ambient_like":
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: ambient RNA, not cell-intrinsic", "detail": f"Within the disease samples, {g} is no more than ~3x higher in cluster {cl} cells than in the least-related lineage of the same samples - the level expected from sample background, not from expression by these cells.", "cluster": cl, "genes": [g], "severity": "blocking"}))
            elif rec.get("attribution") == "shared_with_reference":
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: shared with reference-condition cells, not disease-specific", "detail": f"{g} is expressed at comparable levels in reference-condition cells of the same cell type (a cell-type marker or a pan-condition gene).", "cluster": cl, "genes": [g], "severity": "major"}))
            elif rec.get("attribution") == "disease_wide":
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: disease-wide rather than state-specific", "detail": f"{g} is also high in disease cells of sibling clusters. Fine as a disease target; do not describe it as defining this state.", "cluster": cl, "genes": [g], "severity": "minor"}))
            if rec["effect_size_source"] == "find_markers":
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: effect size lacks replicate-aware support", "detail": "Only the cell-level Wilcoxon result is available; pseudobulk DE did not confirm it (or was not run).", "cluster": cl, "genes": [g], "severity": "minor"}))
        amb = f.get("disease_programs", {})
        if amb.get("ambient_note"):
            steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": "Ambient RNA from the dominant disease population", "detail": amb["ambient_note"] + f" Clusters affected: {', '.join(amb.get('ambient_suspect_clusters', []))}. Any 'disease-wide' label on state genes should be read with this in mind.", "severity": "major"}))
        # a gene nominated from two states would otherwise get the same flag twice
        seen_flags: set[tuple[str, str]] = set()
        deduped: list[ToolUse | str] = []
        for item in steps:
            if isinstance(item, ToolUse) and item.name == "record_finding" and item.input.get("kind") == "critic_flag":
                key = (",".join(item.input.get("genes") or []), item.input.get("title", ""))
                if key in seen_flags:
                    continue
                seen_flags.add(key)
            deduped.append(item)
        steps = deduped
        n_flags = len([s for s in steps if isinstance(s, ToolUse)])
        steps.append(ToolUse("record_finding", {"kind": "summary", "title": "Critic summary", "detail": f"Raised {n_flags} flag(s). Candidates without flags are supported by replicate-aware DE, cluster specificity and consistency across disease samples; treat all nominations as hypotheses for IHC / pharmacologic validation.", "confidence": "medium"}))
        steps.append("Critic review complete.")
        return steps

    # ---- Model protocol -----------------------------------------------------
    # a failure in one of these leaves nothing for the rest of the plan to work on; the script must say so
    # rather than carry on and narrate "no state found" (HCC review run, 2026-10-28: clustering crashed on a
    # library version conflict and the calibration report read as a clean negative)
    CRITICAL_TOOLS = {"load_dataset", "run_qc", "integrate_and_cluster", "differential_abundance", "merge_clusters", "derive_lineages"}

    def _last_tool_errors(self, messages: list[dict[str, Any]]) -> list[tuple[str, str]]:
        """(tool name, error text) for errored tool_result blocks in the most recent user message."""
        if not messages or messages[-1].get("role") != "user" or not isinstance(messages[-1].get("content"), list):
            return []
        # map tool_use ids to names from the preceding assistant turn
        names: dict[str, str] = {}
        if len(messages) >= 2 and isinstance(messages[-2].get("content"), list):
            for b in messages[-2]["content"]:
                if b.get("type") == "tool_use":
                    names[b["id"]] = b["name"]
        out = []
        for b in messages[-1]["content"]:
            if b.get("type") == "tool_result" and b.get("is_error"):
                out.append((names.get(b.get("tool_use_id"), "?"), str(b.get("content"))[:300]))
        return out

    def create(self, *, system: str, messages: list[dict[str, Any]], tools: list[dict[str, Any]], role: str) -> ModelResponse:
        if role not in self._queue:
            self._queue[role] = self._plan(role)
        q = self._queue[role]
        allowed = {t["name"] for t in tools}
        critical = [(n, e) for n, e in self._last_tool_errors(messages) if n in self.CRITICAL_TOOLS]
        if critical and not self._step.get(f"{role}_aborted"):
            self._step[f"{role}_aborted"] = True
            name, err = critical[0]
            q.clear()
            q.extend([
                ToolUse("record_finding", {"kind": "summary", "title": f"Run aborted: tool '{name}' failed",
                                           "detail": f"The scripted procedure stopped because a required tool failed; downstream results would be meaningless. Error: {err}. This is an environment or data-format problem, not a biological negative: fix it (python scripts/sctarget.py check-env) and re-run.",
                                           "confidence": "high"}),
                f"Scripted {role} aborted: {name} failed ({err[:120]}). Nothing below this point is a result.",
            ])
        while q:
            item = q.pop(0)
            if isinstance(item, str) and item.startswith("__"):
                expand = {
                    "__analyst_dynamic__": self._analyst_dynamic,
                    "__analyst_characterise__": self._analyst_characterise,
                    "__analyst_programs__": self._analyst_programs,
                    "__target_dynamic__": self._target_dynamic,
                    "__critic_dynamic__": self._critic_dynamic,
                    "__critic_review__": self._critic_review,
                }
                if item in expand:
                    q[:0] = expand[item]()
                elif item.startswith("__analyst_enrich__"):
                    q[:0] = self._analyst_enrich(item.split("__analyst_enrich__", 1)[1])
                elif item.startswith("__target_lookup__"):
                    q[:0] = self._target_lookup(item.split("__target_lookup__", 1)[1])
                elif item.startswith("__target_record__"):
                    q[:0] = self._target_record(item.split("__target_record__", 1)[1])
                continue
            if isinstance(item, str):
                return ModelResponse(text=item, tool_uses=[], stop_reason="end_turn", usage={"input_tokens": 0, "output_tokens": 0}, raw_content=[{"type": "text", "text": item}])
            if item.name not in allowed:
                continue  # skip tools this role may not use
            raw = [{"type": "tool_use", "id": item.id, "name": item.name, "input": item.input}]
            return ModelResponse(text="", tool_uses=[item], stop_reason="tool_use", usage={"input_tokens": 0, "output_tokens": 0}, raw_content=raw)
        return ModelResponse(text="Done.", tool_uses=[], stop_reason="end_turn", usage={"input_tokens": 0, "output_tokens": 0}, raw_content=[{"type": "text", "text": "Done."}])


def make_model(name: str, session: Any, temperature: float = 0.2) -> Model:
    if name == "scripted":
        return ScriptedModel(session)
    if name == "claude-code" or name.startswith("claude-code:"):
        # subscription-login backend (no API key); a separate experimental arm, see claude_code.py
        from .claude_code import ClaudeCodeModel

        return ClaudeCodeModel(model=name.split(":", 1)[1] if ":" in name else None)
    return AnthropicModel(model=name, temperature=temperature)
