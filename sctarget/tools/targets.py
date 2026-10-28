"""Target scoring and structured findings.

``score_targets`` is deterministic and explainable: it combines evidence already gathered
by other tools (pseudobulk effect size, cluster specificity, druggability, known
compounds, cross-sample consistency) with the weights in ``targets.weights``. The agent
decides *which* genes to score and writes the narrative; the numbers come from here.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from ..session import AnalysisSession
from .registry import tool

ALL = {"planner", "analyst", "target", "critic"}


DEFAULT_EXCLUDE_PATTERNS = [r"^ENS[A-Z]*G\d+", r"^LOC\d+", r"^U\d+-\d+", r"^MT-", r"^RP[LS]\d", r"^SNOR", r"^MIR\d", r"^LINC\d", r"^[A-Z0-9]+\.\d+$"]


def excluded_gene(name: str, patterns: list[str] | None = None) -> bool:
    """Unnamed / placeholder / housekeeping-RNA genes that can never be nominated as targets."""
    import re

    return any(re.match(p, str(name)) for p in (patterns or DEFAULT_EXCLUDE_PATTERNS))


def _clip01(x: float | None, scale: float) -> float:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return 0.0
    return float(np.clip(x / scale, 0.0, 1.0))


def state_basis(session: AnalysisSession, cluster: str) -> str:
    """'composition' (disease-enriched cluster) or 'program' (disease-altered program within a shared cluster)."""
    for st in session.facts.get("findings", {}).get("disease_state", []):
        if str(st.get("cluster")) == str(cluster) and st.get("basis"):
            return str(st["basis"])
    if str(cluster) in session.facts.get("disease_programs", {}).get("disease_altered_clusters", []) and str(cluster) not in session.facts.get("differential_abundance", {}).get("disease_enriched_clusters", []):
        return "program"
    return "composition"


SCORING_VERSIONS = {
    1: "effect size = cluster vs sibling clusters (composition states) or disease vs reference within cluster (program states)",
    2: "as v1, plus a second effect axis for composition states: disease vs reference within the SHARED clusters of the state's compartment (its lineage minus the disease-enriched clusters); the larger axis is used and recorded as target_class state_specific | disease_wide (disease_wide requires lineage log2FC >= 1 and > the state log2FC)",
    3: "as v2, with the lineage axis taken from the ALL-cells lineage contrast (disease vs reference cells of the whole lineage, enriched clusters included) whenever the shared-cluster contrast is unrepresentative (it holds < 25% of the lineage's disease cells, or < 3 reference replicates, or failed); the source is recorded per gene (lineage_axis_source). Plus a proliferation_marker flag (penalty 0.5) for cell-cycle genes. Motivated by the HCC review run: 90-100% of malignant hepatocytes sat in enriched clusters, so the shared-cluster contrast compared an atypical 9% of them and GPC3 (log2FC +2.4) did not reach significance; TOP2A from a cycling sub-state was nominated unflagged",
}

LINEAGE_SHARED_MIN_FRAC = 0.25
LINEAGE_MIN_REF_REPLICATES = 3


def scoring_version(session: AnalysisSession) -> int:
    return int(session.cfg.targets.get("scoring_version", 1))


def _effect_size(session: AnalysisSession, gene: str, cluster: str, basis: str = "composition", contrast: str | None = None) -> tuple[float | None, str]:
    pb = session.facts.get("pseudobulk_de", {})
    wanted = contrast or ("disease_vs_reference_within_cluster" if basis == "program" else "cluster_vs_reference_clusters")
    for key, val in pb.items():
        cl, contrast = key.split(":", 1)
        if cl == cluster and contrast == wanted:
            for row in val.get("top_up", []):
                if row.get("gene") == gene:
                    return float(row["log2fc"]), "pseudobulk_de"
            # not among the significant hits: read the full table from disk
            path = val.get("table_path")
            if path:
                try:
                    full = pd.read_csv(path, index_col=0)
                    if gene in full.index:
                        padj = full.loc[gene, "padj"]
                        lfc = float(full.loc[gene, "log2fc"])
                        # down-weight weak evidence instead of discarding it outright
                        if pd.isna(padj) or padj >= 0.25:
                            return 0.0, "pseudobulk_de"
                        return (lfc if padj < 0.05 else 0.5 * lfc), "pseudobulk_de"
                except Exception:
                    pass
    if wanted == "cluster_vs_reference_clusters":
        for key, rows in session.facts.get("markers", {}).items():
            if key.startswith(f"{cluster}_vs_"):
                for row in rows:
                    if row.get("gene") == gene:
                        return float(row["log2fc"]), "find_markers"
    return None, "missing"


def _ensure_compartment_de(session: AnalysisSession, cluster: str, contrast: str = "disease_vs_reference_within_compartment") -> bool:
    """Scoring v2/v3: make sure the lineage-level pseudobulk exists for this cluster (computed once,
    cached in facts like any pseudobulk_de result)."""
    key = f"{cluster}:{contrast}"
    if key in session.facts.get("pseudobulk_de", {}):
        return True
    if session.adata is None:
        return False
    from .analysis import _run_pseudobulk

    try:
        r = _run_pseudobulk(session, cluster, None, contrast, None, 10, 0.05, 10.0, record=True)
    except Exception as e:  # noqa: BLE001 - a failed contrast is recorded as unavailable, not fatal
        session.facts.setdefault("pseudobulk_de_errors", {})[key] = f"{type(e).__name__}: {str(e)[:200]}"
        return False
    if "error" in r:
        session.facts.setdefault("pseudobulk_de_errors", {})[key] = str(r["error"])[:200]
        return False
    return True


def lineage_axis_source(session: AnalysisSession, cluster: str) -> tuple[str | None, dict[str, Any]]:
    """Scoring v3: which lineage-level contrast carries the disease-wide axis for this cluster, and why.
    The shared-cluster contrast is preferred (it separates 'disease-wide programme' from 'the state
    itself'), but only when it is representative: >= 25% of the lineage's disease cells and >= 3
    reference replicates. Otherwise the all-cells lineage contrast is used. Both are attempted."""
    pb = session.facts.get("pseudobulk_de", {})
    have_shared = _ensure_compartment_de(session, cluster, "disease_vs_reference_within_compartment")
    have_all = _ensure_compartment_de(session, cluster, "disease_vs_reference_within_lineage")
    shared = pb.get(f"{cluster}:disease_vs_reference_within_compartment", {}) if have_shared else {}
    cov = shared.get("lineage_coverage") or {}
    frac = cov.get("frac_of_lineage_disease_cells_in_contrast")
    n_ref = (shared.get("replicates") or {}).get("reference", 0)
    why = {"shared_available": have_shared, "all_available": have_all, "shared_frac_of_lineage_disease_cells": frac, "shared_reference_replicates": n_ref}
    if have_shared and (frac is None or frac >= LINEAGE_SHARED_MIN_FRAC) and n_ref >= LINEAGE_MIN_REF_REPLICATES:
        return "disease_vs_reference_within_compartment", why
    if have_all:
        why["reason"] = "shared-cluster contrast unrepresentative or unavailable" + (f" ({frac:.0%} of lineage disease cells, {n_ref} reference replicates)" if frac is not None else "")
        return "disease_vs_reference_within_lineage", why
    if have_shared:
        why["reason"] = "all-cells lineage contrast unavailable; using the shared-cluster contrast despite low coverage"
        return "disease_vs_reference_within_compartment", why
    return None, why


@tool(
    name="score_targets",
    description="""
    Deterministically score and rank candidate target genes for one disease-specific
    cluster using evidence already collected: effect size (pseudobulk_de or find_markers),
    specificity and sample consistency (gene_specificity), druggability and known
    compounds (druggability_lookup). Weights come from config targets.weights. Returns the
    ranked table with per-component scores and flags (low_specificity, not_druggable,
    no_known_compounds, inconsistent_across_samples, missing_evidence). Call the evidence
    tools first; missing evidence scores 0 and is flagged. Flags also penalise the score
    (low_specificity x0.6, shared_across_conditions x0.6, inconsistent_across_samples x0.8,
    reference_driven / ambient_like x0.2, shared_with_reference / weak x0.5; disease_wide is informational). Run
    condition_attribution on the candidates first or they are flagged missing_evidence:attribution.
    With targets.scoring_version 2, composition states get a second effect axis (disease vs
    reference within the state's compartment / lineage) and each gene's target_class says which
    axis carried it: state_specific or disease_wide.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
            "cluster": {"type": "string"},
            "human_symbols": {"type": ["object", "null"], "description": "optional mapping dataset symbol -> human symbol used for druggability lookups"},
        },
        "required": ["genes", "cluster"],
    },
    roles=ALL,
)
def score_targets(session: AnalysisSession, genes: list[str], cluster: str, human_symbols: dict[str, str] | None = None) -> dict[str, Any]:
    w = dict(session.cfg.targets.get("weights", {}))
    objective = str(session.cfg.targets.get("objective", "target"))
    if objective == "marker":
        # a perturbation / response study: rank the strongest, most specific, most consistent
        # condition-specific genes; druggability is reported but does not drive the ranking
        w = {"effect_size": 0.40, "specificity": 0.35, "druggability": 0.0, "known_compounds": 0.0, "consistency_across_samples": 0.25, **{k: v for k, v in w.items() if k in ("effect_size", "specificity", "consistency_across_samples") and session.cfg.targets.get("weights_override_marker")}}
    w_eff = float(w.get("effect_size", 0.25))
    w_spec = float(w.get("specificity", 0.25))
    w_drug = float(w.get("druggability", 0.25))
    w_cmp = float(w.get("known_compounds", 0.15))
    w_cons = float(w.get("consistency_across_samples", 0.10))
    spec_rows = {r["gene"]: r for r in session.facts.get("gene_specificity", {}).get(str(cluster), [])}
    drug_rows = session.facts.get("druggability", {})
    human_symbols = human_symbols or {}
    for g, v in session.facts.get("orthologs", {}).items():
        if v.get("human_symbol") and g not in human_symbols:
            human_symbols[g] = v["human_symbol"]

    patterns = session.cfg.targets.get("exclude_gene_patterns") or None
    excluded = [g for g in dict.fromkeys(genes) if excluded_gene(g, patterns)]
    basis = state_basis(session, str(cluster))
    sv = scoring_version(session)
    lin_contrast, lin_why = None, {}
    if sv >= 3 and basis == "composition":
        lin_contrast, lin_why = lineage_axis_source(session, str(cluster))
    elif sv == 2 and basis == "composition" and _ensure_compartment_de(session, str(cluster)):
        lin_contrast = "disease_vs_reference_within_compartment"
    use_lineage_axis = lin_contrast is not None
    session.facts["scoring_version"] = sv
    from .analysis import CELL_CYCLE_GENES
    rows = []
    for g in dict.fromkeys(genes):
        if g in excluded:
            continue
        hs = human_symbols.get(g, g)
        flags: list[str] = []
        eff, eff_src = _effect_size(session, g, str(cluster), basis)
        eff_lin, target_class = None, ("program" if basis == "program" else "state_specific")
        if use_lineage_axis:
            eff_lin, _ = _effect_size(session, g, str(cluster), basis, contrast=lin_contrast)
            # the larger axis carries the effect component; which one won is recorded, not hidden.
            # disease_wide needs a real lineage effect (>= 1 log2) that beats the state effect.
            if eff_lin is not None and eff_lin >= 1.0 and (eff is None or eff_lin > eff):
                target_class = "disease_wide"
        if eff is None and eff_lin is None:
            flags.append("missing_evidence:effect_size")
        s_eff = max(_clip01(eff, 4.0), _clip01(eff_lin, 4.0) if use_lineage_axis else 0.0)

        sp = spec_rows.get(g)
        if sp is None:
            flags.append("missing_evidence:specificity")
            s_spec, cons = 0.0, None
        else:
            cons = sp.get("disease_sample_consistency")
            dvr = sp.get("log2_disease_vs_reference_within_cluster")
            if basis == "program":
                # the state is a shared cell type running a disease program: what matters is how
                # disease-specific the gene is INSIDE that cluster, not whether the cluster is unique
                s_spec = _clip01(dvr, 3.0)
                if dvr is None or dvr < 0.5:
                    flags.append("shared_across_conditions")
            else:
                s_spec = _clip01(sp.get("log2_specificity"), 4.0)
                if (sp.get("log2_specificity") or 0) < 1.0:
                    flags.append("low_specificity")
                if dvr is not None and dvr < 0.5 and (sp.get("log2_specificity") or 0) < 1.0:
                    flags.append("shared_across_conditions")
        s_cons = 0.5 if cons is None else float(cons)
        if cons is not None and cons < 0.67:
            flags.append("inconsistent_across_samples")

        attr = session.facts.get("condition_attribution", {}).get(str(cluster), {}).get(g)
        if attr is None:
            flags.append("missing_evidence:attribution")
        elif attr in ("reference_driven", "ambient_like", "shared_with_reference", "disease_wide", "weak"):
            flags.append(attr)

        d = drug_rows.get(hs) or drug_rows.get(g)
        if d is None:
            flags.append("missing_evidence:druggability")
            s_drug, s_cmp = 0.0, 0.0
        else:
            tract = [t for t in (d.get("ot_sm_tractability"), d.get("ot_ab_tractability")) if t]
            if d.get("is_druggable_class") and any("Approved" in t or "Clinical" in t for t in tract):
                s_drug = 1.0
            elif d.get("is_druggable_class"):
                s_drug = 0.6
            elif tract:
                s_drug = 0.3
            else:
                s_drug = 0.0
                flags.append("not_druggable")
            # approved compounds (typed inhibitor or untyped) drive the score; typed-but-unapproved
            # inhibitors add a little. DGIdb often leaves interaction type empty for well-known drugs.
            n_appr = int(d.get("n_approved_compounds", d.get("n_approved_inhibitors")) or 0)
            n_cmp = int(d.get("n_compounds", d.get("n_inhibitors")) or 0)
            s_cmp = float(min(1.0, n_appr / 3.0 + (0.1 if n_cmp > n_appr else 0.0)))
            if n_cmp == 0:
                flags.append("no_known_compounds")

        score = w_eff * s_eff + w_spec * s_spec + w_drug * s_drug + w_cmp * s_cmp + w_cons * s_cons
        # flags are not just annotations: a gene that is equally high elsewhere, or inconsistent
        # across disease samples, is a worse target however druggable it is
        penalty = 1.0
        if "low_specificity" in flags:
            penalty *= 0.6
        if "shared_across_conditions" in flags:
            penalty *= 0.6
        if "inconsistent_across_samples" in flags:
            penalty *= 0.8
        if "reference_driven" in flags or "ambient_like" in flags:
            penalty *= 0.2  # wrong condition's cells, or sample background: not a cell-intrinsic disease signal
        if "shared_with_reference" in flags or "weak" in flags:
            penalty *= 0.5
        if sv >= 3 and hs.upper() in CELL_CYCLE_GENES:
            flags.append("proliferation_marker")
            penalty *= 0.5  # describes that the cells divide, not what the disease is; a decoy in every key
        # disease_wide is informational only: a disease-wide gene is still a legitimate target
        score *= penalty
        rows.append(
            {
                "gene": g,
                "human_symbol": hs,
                "basis": basis,
                "score": round(float(score), 4),
                "effect_size_log2fc": None if eff is None else round(eff, 3),
                "effect_size_source": eff_src,
                "lineage_log2fc_disease_vs_reference": None if eff_lin is None else round(eff_lin, 3),
                "lineage_axis_source": lin_contrast,
                "target_class": target_class,
                "log2_specificity": None if sp is None else sp.get("log2_specificity"),
                "log2_disease_vs_reference_within_cluster": None if sp is None else sp.get("log2_disease_vs_reference_within_cluster"),
                "sample_consistency": cons,
                "druggability_score": s_drug,
                "known_compound_score": round(s_cmp, 3),
                "top_inhibitors": (d or {}).get("top_inhibitors", []),
                "attribution": attr,
                "flag_penalty": round(penalty, 2),
                "flags": flags,
            }
        )
    df = pd.DataFrame(rows).set_index("gene").sort_values("score", ascending=False)
    df["rank"] = np.arange(1, len(df) + 1)
    session.save_table(f"ranked_targets_{cluster}", df)
    session.facts.setdefault("ranked_targets", {})[str(cluster)] = df.reset_index().to_dict(orient="records")
    return {"cluster": str(cluster), "state_basis": basis, "scoring_version": sv, "objective": objective, "lineage_axis": lin_why | {"used": lin_contrast}, "effect_axes": (["cluster_vs_siblings", lin_contrast] if use_lineage_axis else ["cluster_vs_siblings" if basis == "composition" else "disease_vs_reference_within_cluster"]), "specificity_component": ("within-cluster disease vs reference log2FC" if basis == "program" else "cluster vs max other cluster log2 ratio"), "weights": {"effect_size": w_eff, "specificity": w_spec, "druggability": w_drug, "known_compounds": w_cmp, "consistency": w_cons}, "excluded_unnamed_genes": excluded, "ranking": df}


MIN_STATE_GENES = 5  # a disease_state finding without its signature cannot be acted on downstream


@tool(
    name="record_finding",
    description="""
    Record a structured conclusion so it appears in the report and can be evaluated.
    kinds: 'disease_state' (a cluster you conclude is disease-specific), 'candidate_target'
    (a gene you nominate), 'critic_flag' (a concern about a claim), 'summary' (free-text
    conclusion for the report). Always include the evidence you relied on.
    A 'disease_state' finding MUST carry `genes`: the state's top up-regulated genes (>= 5,
    ideally 15, from pseudobulk_de / find_markers against sibling clusters). Record the state
    AFTER characterising it; a disease_state without genes is rejected.""",
    input_schema={
        "properties": {
            "kind": {"type": "string", "enum": ["disease_state", "candidate_target", "critic_flag", "summary"]},
            "title": {"type": "string"},
            "detail": {"type": "string", "description": "2-5 sentences of justification citing tool results"},
            "cluster": {"type": ["string", "null"]},
            "genes": {"type": ["array", "null"], "items": {"type": "string"}},
            "confidence": {"type": ["string", "null"], "enum": ["high", "medium", "low", None]},
            "severity": {"type": ["string", "null"], "enum": ["blocking", "major", "minor", None], "description": "critic_flag only"},
            "basis": {"type": ["string", "null"], "enum": ["composition", "program", None], "description": "disease_state only: 'composition' = cluster enriched in disease; 'program' = shared cluster with a disease-specific expression program"},
        },
        "required": ["kind", "title", "detail"],
    },
    roles=ALL,
)
def record_finding(
    session: AnalysisSession,
    kind: str,
    title: str,
    detail: str,
    cluster: str | None = None,
    genes: list[str] | None = None,
    confidence: str | None = None,
    severity: str | None = None,
    basis: str | None = None,
) -> dict[str, Any]:
    if kind == "disease_state" and len(genes or []) < MIN_STATE_GENES:
        return {
            "error": f"disease_state finding for cluster {cluster} rejected: it carries {len(genes or [])} genes, at least {MIN_STATE_GENES} are required.",
            "hint": "Characterise the state first (pseudobulk_de and find_markers against sibling clusters of the same cell type), then record it again with its top up-regulated genes (ideally 15) in `genes`. Downstream roles only see the genes you record.",
            "recorded": None,
        }
    entry = {"title": title, "detail": detail, "cluster": cluster, "genes": genes or [], "confidence": confidence, "severity": severity}
    if kind == "disease_state":
        entry["basis"] = basis or "composition"
    session.facts.setdefault("findings", {}).setdefault(kind, []).append(entry)
    return {"recorded": kind, "n_total_of_kind": len(session.facts["findings"][kind])}
