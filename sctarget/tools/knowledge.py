"""Knowledge tools: orthologs, druggability, known compounds, literature counts.

Live sources (all free, no key required):
  * Ensembl REST          https://rest.ensembl.org           (ortholog mapping)
  * DGIdb GraphQL v5      https://dgidb.org/api/graphql       (drug-gene interactions, gene categories)
  * Open Targets GraphQL  https://api.platform.opentargets.org/api/v4/graphql (tractability, known drugs)
  * NCBI E-utilities      https://eutils.ncbi.nlm.nih.gov     (PubMed counts)

Every response is cached on disk (``SCTARGET_CACHE_DIR``) keyed by a hash of the request,
so a re-run is deterministic and offline. With ``SCTARGET_OFFLINE=1`` only the cache and
the bundled fixture file are consulted; a gene not covered returns ``source: "none"``
rather than failing, so the pipeline keeps running on air-gapped nodes.

Only gene symbols and aggregate statistics ever leave the machine in these requests;
no expression values are sent.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from ..config import REPO_ROOT, cache_dir, offline_mode
from ..session import AnalysisSession
from .registry import tool

FIXTURE_PATH = REPO_ROOT / "data" / "fixtures" / "knowledge_fixtures.json"
TARGET_ROLES = {"planner", "target", "critic"}

_fixtures: dict[str, Any] | None = None


def fixtures() -> dict[str, Any]:
    global _fixtures
    if _fixtures is None:
        with open(FIXTURE_PATH) as fh:
            _fixtures = json.load(fh)
    return _fixtures


class KnowledgeClient:
    """Thin HTTP client with disk cache + fixture fallback."""

    def __init__(self, offline: bool | None = None, timeout: float = 30.0) -> None:
        self.offline = offline_mode() if offline is None else offline
        self.timeout = timeout
        self.cache = cache_dir()
        self.calls: list[dict[str, Any]] = []

    # ------------------------------------------------------------------ cache
    def _key(self, kind: str, payload: Any) -> Path:
        h = hashlib.sha1(json.dumps([kind, payload], sort_keys=True).encode()).hexdigest()[:20]
        return self.cache / f"{kind}_{h}.json"

    def _get_cached(self, kind: str, payload: Any) -> Any | None:
        p = self._key(kind, payload)
        if p.exists():
            with open(p) as fh:
                return json.load(fh)
        return None

    def _put_cached(self, kind: str, payload: Any, value: Any) -> None:
        with open(self._key(kind, payload), "w") as fh:
            json.dump(value, fh)

    def _request(self, kind: str, method: str, url: str, payload: Any, **kw: Any) -> tuple[Any | None, str]:
        """Returns (json, source) where source in {'cache','live','none'}."""
        cached = self._get_cached(kind, payload)
        if cached is not None:
            return cached, "cache"
        if self.offline:
            return None, "none"
        t0 = time.time()
        try:
            if method == "GET":
                r = requests.get(url, params=payload if isinstance(payload, dict) else None, timeout=self.timeout, headers={"Accept": "application/json"}, **kw)
            else:
                r = requests.post(url, json=payload, timeout=self.timeout, headers={"Accept": "application/json", "Content-Type": "application/json"}, **kw)
            r.raise_for_status()
            data = r.json()
            self._put_cached(kind, payload, data)
            self.calls.append({"kind": kind, "url": url, "ok": True, "s": round(time.time() - t0, 2)})
            return data, "live"
        except Exception as e:  # network down / API changed
            self.calls.append({"kind": kind, "url": url, "ok": False, "error": f"{type(e).__name__}: {e}", "s": round(time.time() - t0, 2)})
            return None, "none"

    # ------------------------------------------------------------------ Ensembl
    def orthologs(self, genes: list[str], from_species: str) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        if from_species in ("homo_sapiens", "human"):
            return {g: {"human_symbol": g, "type": "identity", "source": "n/a"} for g in genes}
        fx = fixtures().get("orthologs", {}).get(from_species, {})
        for g in genes:
            payload = {"species": from_species, "symbol": g, "target_species": "homo_sapiens", "type": "orthologues", "format": "condensed"}
            data, src = self._request("ensembl_homology", "GET", f"https://rest.ensembl.org/homology/symbol/{from_species}/{g}", {"target_species": "homo_sapiens", "type": "orthologues", "format": "condensed", "content-type": "application/json"})
            if data is not None:
                homs = (data.get("data") or [{}])[0].get("homologies", [])
                hs = [h for h in homs if h.get("species") == "homo_sapiens"]
                if hs:
                    ids = [h["id"] for h in hs]
                    look, _ = self._request("ensembl_lookup", "POST", "https://rest.ensembl.org/lookup/id", {"ids": ids})
                    symbols = sorted({(look or {}).get(i, {}).get("display_name") for i in ids} - {None}) if look else []
                    out[g] = {"human_symbol": symbols[0] if symbols else None, "all_human_symbols": symbols, "type": hs[0].get("type"), "source": src}
                else:
                    out[g] = {"human_symbol": None, "type": "no_ortholog", "source": src}
            elif g in fx:
                out[g] = {**fx[g], "source": "fixture"}
            else:
                # conservative default: assume same symbol (true for most dog/human genes) but flag it
                out[g] = {"human_symbol": g, "type": "assumed_same_symbol", "source": "none"}
        return out

    # ------------------------------------------------------------------ DGIdb
    DGIDB_QUERY = """
    query($names: [String!]!) {
      genes(names: $names) {
        nodes {
          name
          geneCategoriesWithSources { name }
          interactions {
            interactionScore
            interactionTypes { type directionality }
            drug { name approved }
            sources { sourceDbName }
          }
        }
      }
    }"""

    def dgidb(self, genes: list[str]) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        fx = fixtures().get("dgidb", {})
        payload = {"query": self.DGIDB_QUERY, "variables": {"names": sorted(genes)}}
        data, src = self._request("dgidb", "POST", "https://dgidb.org/api/graphql", payload)
        found: dict[str, dict[str, Any]] = {}
        if data is not None:
            for node in (data.get("data", {}).get("genes", {}) or {}).get("nodes", []):
                inters = []
                for it in node.get("interactions", []) or []:
                    inters.append(
                        {
                            "drug": (it.get("drug") or {}).get("name"),
                            "approved": bool((it.get("drug") or {}).get("approved")),
                            "types": [t.get("type") for t in it.get("interactionTypes", []) or []],
                            "score": it.get("interactionScore"),
                            "sources": [s.get("sourceDbName") for s in it.get("sources", []) or []],
                        }
                    )
                found[node["name"].upper()] = {"categories": [c["name"] for c in node.get("geneCategoriesWithSources", []) or []], "interactions": inters, "source": src}
        for g in genes:
            if g.upper() in found:
                out[g] = found[g.upper()]
            elif g in fx:
                out[g] = {**fx[g], "source": "fixture"}
            else:
                out[g] = {"categories": [], "interactions": [], "source": "none"}
        return out

    # ------------------------------------------------------------------ Open Targets
    OT_SEARCH = """
    query($q: String!) { search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) { hits { id name } } }"""
    OT_TARGET = """
    query($id: String!) {
      target(ensemblId: $id) {
        approvedSymbol
        tractability { label modality value }
        knownDrugs(size: 50) { count rows { drug { name } phase mechanismOfAction } }
      }
    }"""

    def open_targets(self, genes: list[str]) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        fx = fixtures().get("open_targets", {})
        url = "https://api.platform.opentargets.org/api/v4/graphql"
        for g in genes:
            data, src = self._request("ot_search", "POST", url, {"query": self.OT_SEARCH, "variables": {"q": g}})
            hits = (((data or {}).get("data") or {}).get("search") or {}).get("hits") or []
            hit = next((h for h in hits if h.get("name", "").upper() == g.upper()), None) if hits else None
            if hit:
                tdata, tsrc = self._request("ot_target", "POST", url, {"query": self.OT_TARGET, "variables": {"id": hit["id"]}})
                t = ((tdata or {}).get("data") or {}).get("target") or {}
                sm = [x["label"] for x in t.get("tractability", []) if x.get("modality") == "SM" and x.get("value")]
                ab = [x["label"] for x in t.get("tractability", []) if x.get("modality") == "AB" and x.get("value")]
                kd = t.get("knownDrugs") or {}
                phases = [r.get("phase") for r in kd.get("rows", []) if r.get("phase") is not None]
                out[g] = {
                    "ensembl_id": hit["id"],
                    "tractability": {"small_molecule": sm[0] if sm else None, "antibody": ab[0] if ab else None},
                    "known_drugs_count": kd.get("count", 0),
                    "max_clinical_phase": max(phases) if phases else None,
                    "known_drugs": sorted({r["drug"]["name"] for r in kd.get("rows", []) if r.get("drug")})[:25],
                    "source": tsrc,
                }
            elif g in fx:
                out[g] = {**fx[g], "source": "fixture"}
            else:
                out[g] = {"ensembl_id": None, "tractability": {"small_molecule": None, "antibody": None}, "known_drugs_count": 0, "max_clinical_phase": None, "source": "none"}
        return out

    # ------------------------------------------------------------------ PubMed
    def pubmed_count(self, gene: str, disease_terms: list[str]) -> dict[str, Any]:
        q = f"{gene} AND ({' OR '.join(disease_terms)})"
        fx = fixtures().get("pubmed", {})
        data, src = self._request("pubmed", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", {"db": "pubmed", "term": q, "retmode": "json", "retmax": 5})
        if data is not None:
            res = data.get("esearchresult", {})
            return {"query": q, "count": int(res.get("count", 0)), "pmids": res.get("idlist", []), "source": src}
        key = f"{gene} {disease_terms[0]}" if disease_terms else gene
        if key in fx:
            return {"query": q, **fx[key], "source": "fixture"}
        return {"query": q, "count": None, "pmids": [], "source": "none"}


# ----------------------------------------------------------------------------- tools
def _client(session: AnalysisSession) -> KnowledgeClient:
    if "knowledge_client" not in session.facts:
        session.facts["knowledge_client"] = KnowledgeClient()
    return session.facts["knowledge_client"]


@tool(
    name="map_orthologs",
    description="""
    Map gene symbols from the dataset's species to human symbols (Ensembl orthologs).
    Human datasets return identity. Flags one-to-many and missing orthologs so they are not
    silently dropped. Run before druggability lookups for non-human data.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
            "from_species": {"type": ["string", "null"], "description": "defaults to dataset.species"},
        },
        "required": ["genes"],
    },
    roles=TARGET_ROLES,
)
def map_orthologs(session: AnalysisSession, genes: list[str], from_species: str | None = None) -> dict[str, Any]:
    sp = from_species or session.cfg.dataset.get("species", "homo_sapiens")
    res = _client(session).orthologs(genes, sp)
    session.facts.setdefault("orthologs", {}).update(res)
    unmapped = [g for g, v in res.items() if not v.get("human_symbol")]
    assumed = [g for g, v in res.items() if v.get("type") == "assumed_same_symbol"]
    return {"from_species": sp, "mapping": res, "unmapped": unmapped, "assumed_same_symbol": assumed}


@tool(
    name="druggability_lookup",
    description="""
    For human gene symbols, query DGIdb (gene categories such as KINASE / CLINICALLY
    ACTIONABLE, drug-gene interactions with approval status and interaction type) and
    Open Targets (small-molecule / antibody tractability, number of known drugs, max
    clinical phase). Returns one row per gene plus the inhibitor list. Results are cached;
    'source' tells you whether each value came from a live call, the cache, or the bundled
    fixture (fixture values are illustrative only).""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}, "description": "human symbols (<= 40)"},
        },
        "required": ["genes"],
    },
    roles=TARGET_ROLES,
)
def druggability_lookup(session: AnalysisSession, genes: list[str]) -> dict[str, Any]:
    genes = list(dict.fromkeys(g for g in genes if g))[:40]
    client = _client(session)
    dg = client.dgidb(genes)
    ot = client.open_targets(genes)
    wanted = {c.upper() for c in session.cfg.targets.get("druggable_classes", [])}
    rows = []
    for g in genes:
        d, o = dg[g], ot[g]
        inhibitors = [i for i in d["interactions"] if any((t or "").lower() in ("inhibitor", "antagonist", "blocker", "antibody", "inverse agonist", "negative modulator") for t in i["types"])]
        approved_inh = [i for i in inhibitors if i.get("approved")]
        cats = [c.upper() for c in d["categories"]]
        rows.append(
            {
                "gene": g,
                "druggable_categories": sorted(set(cats) & wanted) or sorted(set(cats))[:4],
                "is_druggable_class": bool(set(cats) & wanted),
                "n_interactions": len(d["interactions"]),
                "n_inhibitors": len(inhibitors),
                "n_approved_inhibitors": len(approved_inh),
                "top_inhibitors": [i["drug"] for i in sorted(inhibitors, key=lambda x: -(x.get("score") or 0))[:6]],
                "ot_sm_tractability": o["tractability"].get("small_molecule"),
                "ot_ab_tractability": o["tractability"].get("antibody"),
                "ot_known_drugs": o.get("known_drugs_count", 0),
                "ot_max_phase": o.get("max_clinical_phase"),
                "source_dgidb": d.get("source"),
                "source_opentargets": o.get("source"),
            }
        )
    df = pd.DataFrame(rows).set_index("gene")
    session.save_table("druggability", df)
    session.facts.setdefault("druggability", {}).update({r["gene"]: r for r in rows})
    live_failures = [c for c in client.calls if not c["ok"]]
    return {
        "table": df,
        "offline_mode": client.offline,
        "n_live_failures": len(live_failures),
        "note": "source=fixture means illustrative bundled data, not a live database result." if any(r["source_dgidb"] == "fixture" for r in rows) else None,
    }


@tool(
    name="literature_count",
    description="""
    PubMed hit counts for each gene combined with disease terms (e.g. 'ameloblastoma'),
    with up to 5 PMIDs. A low count on a strong candidate indicates novelty; a high count
    means the target is already well studied in this disease. Cached.""",
    input_schema={
        "properties": {
            "genes": {"type": "array", "items": {"type": "string"}},
            "disease_terms": {"type": "array", "items": {"type": "string"}},
        },
        "required": ["genes", "disease_terms"],
    },
    roles=TARGET_ROLES,
)
def literature_count(session: AnalysisSession, genes: list[str], disease_terms: list[str]) -> dict[str, Any]:
    client = _client(session)
    out = {g: client.pubmed_count(g, disease_terms) for g in genes[:30]}
    session.facts.setdefault("literature", {}).update(out)
    return {"results": out}
