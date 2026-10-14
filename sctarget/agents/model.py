"""Model backends behind one small interface.

``AnthropicModel`` wraps the Anthropic Python SDK's Messages API with tool use.
``ScriptedModel`` is a deterministic test double that issues a fixed, sensible sequence
of tool calls per role so the whole pipeline can run in CI and offline without an API
key. It exists to test the machinery (loop, tools, privacy, report), not to stand in for
the reasoning of a real model.
"""
from __future__ import annotations

import json
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
        enriched = sorted(da.get("disease_enriched_clusters", []), key=lambda c: -next((r.get("log2_ratio") or 0) for r in da.get("table", []) if r["cluster"] == c))[:2]
        steps: list[ToolUse | str] = []
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
                steps.append(ToolUse("disease_programs", {"cluster_key": str(f.get("clustering", {}).get("cluster_key", "")).split("_sub")[0] or None}))
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
            steps.append(ToolUse("disease_programs", {"cluster_key": str(f.get("clustering", {}).get("cluster_key", "")).split("_sub")[0] or None}))
            steps.append("__analyst_programs__")
        for cl in enriched:
            dom = comp.get(cl, {}).get("dominant_celltype")
            siblings = [c for c, r in comp.items() if c != cl and dom and r.get("dominant_celltype") == dom] or None
            steps.append(ToolUse("find_markers", {"cluster": cl, "reference_clusters": siblings, "n_top": 30}))
            steps.append(ToolUse("pseudobulk_de", {"cluster": cl, "reference_clusters": siblings, "n_top": 40}))
            steps.append(f"__analyst_enrich__{cl}")
        return steps

    def _analyst_enrich(self, cl: str) -> list[ToolUse | str]:
        f = self._facts()
        pb = f.get("pseudobulk_de", {}).get(f"{cl}:cluster_vs_reference_clusters", {})
        from ..tools.targets import excluded_gene

        up = [r["gene"] for r in pb.get("top_up", []) if not excluded_gene(r["gene"])][:40]
        if not up:
            mk_key = next((k for k in f.get("markers", {}) if k.startswith(f"{cl}_vs_")), None)
            up = [r["gene"] for r in f.get("markers", {}).get(mk_key, []) if not excluded_gene(r["gene"])][:40]
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
                    "detail": f"Cluster {cl} has log2 abundance ratio {row.get('log2_ratio')} (disease vs reference, per-sample), present in {row.get('n_disease_samples_present')}/{row.get('n_disease_samples_total')} disease samples; {round(100*row.get('frac_cells_from_disease',0))}% of its cells come from the disease condition. Pseudobulk DE vs sibling clusters yields {pb.get('n_sig_up')} significantly up-regulated genes (method: {pb.get('method')}). Top genes: {', '.join(up[:8])}.",
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
            up = [x["gene"] for x in pb.get("top_up", []) if not excluded_gene(x["gene"])][:40]
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
        steps: list[ToolUse | str] = []
        species = self.session.cfg.dataset.get("species", "homo_sapiens")
        for st in states[:4]:
            cl = st["cluster"]
            contrast = "disease_vs_reference_within_cluster" if st.get("basis") == "program" else "cluster_vs_reference_clusters"
            pb = f.get("pseudobulk_de", {}).get(f"{cl}:{contrast}", {})
            from ..tools.targets import excluded_gene

            genes = [r["gene"] for r in pb.get("top_up", []) if not excluded_gene(r["gene"])][:25] or st.get("genes", [])[:25]
            # add the cell-level markers too (a naive analyst's first pick)
            mk_key = next((k for k in f.get("markers", {}) if k.startswith(f"{cl}_vs_")), None)
            genes += [r["gene"] for r in f.get("markers", {}).get(mk_key, []) if not excluded_gene(r["gene"])][:15]
            genes = list(dict.fromkeys(genes))[:30]
            steps.append(ToolUse("gene_specificity", {"genes": genes, "cluster": cl}))
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
        for r in ranked[:5]:
            if r["druggability_score"] == 0 and r["known_compound_score"] == 0:
                continue
            if r["score"] < 0.3:
                continue  # weak or heavily penalised: not worth a nomination
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
        states = [st["cluster"] for st in f.get("findings", {}).get("disease_state", [])][:4]
        if cl == states[-1]:
            n_total = len(f.get("findings", {}).get("candidate_target", [])) + len([s for s in steps if isinstance(s, ToolUse)])
            steps.append(f"Nominated {n_total} candidate targets across clusters {', '.join(states)}.")
        return steps

    def _critic_dynamic(self) -> list[ToolUse | str]:
        f = self._facts()
        steps: list[ToolUse | str] = [ToolUse("differential_abundance", {"within_celltype": True})]
        steps.append("__critic_review__")
        return steps

    def _critic_review(self) -> list[ToolUse | str]:
        f = self._facts()
        steps: list[ToolUse | str] = []
        da_rows = {r["cluster"]: r for r in f.get("differential_abundance", {}).get("table", [])}
        da_within = {r["cluster"]: r for r in f.get("differential_abundance_within_celltype", {}).get("table", [])}
        for st in f.get("findings", {}).get("disease_state", []):
            cl = st["cluster"]
            row = da_rows.get(cl, {})
            wrow = da_within.get(cl, {})
            if st.get("basis") == "program":
                continue  # program states are judged by within-cluster DE, not by abundance
            if wrow and not wrow.get("disease_enriched"):
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
            if "inconsistent_across_samples" in rec["flags"]:
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: inconsistent across disease samples", "detail": f"Detected in >=10% of state cells in only {rec['sample_consistency']} of disease samples.", "cluster": cl, "genes": [g], "severity": "major"}))
            if rec["effect_size_source"] == "find_markers":
                steps.append(ToolUse("record_finding", {"kind": "critic_flag", "title": f"{g}: effect size lacks replicate-aware support", "detail": "Only the cell-level Wilcoxon result is available; pseudobulk DE did not confirm it (or was not run).", "cluster": cl, "genes": [g], "severity": "minor"}))
        n_flags = len([s for s in steps if isinstance(s, ToolUse)])
        steps.append(ToolUse("record_finding", {"kind": "summary", "title": "Critic summary", "detail": f"Raised {n_flags} flag(s). Candidates without flags are supported by replicate-aware DE, cluster specificity and consistency across disease samples; treat all nominations as hypotheses for IHC / pharmacologic validation.", "confidence": "medium"}))
        steps.append("Critic review complete.")
        return steps

    # ---- Model protocol -----------------------------------------------------
    def create(self, *, system: str, messages: list[dict[str, Any]], tools: list[dict[str, Any]], role: str) -> ModelResponse:
        if role not in self._queue:
            self._queue[role] = self._plan(role)
        q = self._queue[role]
        allowed = {t["name"] for t in tools}
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
    return AnthropicModel(model=name, temperature=temperature)
