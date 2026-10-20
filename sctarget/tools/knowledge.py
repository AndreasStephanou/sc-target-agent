"""Knowledge tools: orthologs, druggability, known compounds, literature counts.

Live sources (all free, no key required):
  * Ensembl REST          https://rest.ensembl.org           (ortholog mapping)
  * DGIdb GraphQL v5      https://dgidb.org/api/graphql       (drug-gene interactions, gene categories)
  * Open Targets GraphQL  https://api.platform.opentargets.org/api/v4/graphql (tractability, known drugs)
  * NCBI E-utilities      https://eutils.ncbi.nlm.nih.gov     (PubMed counts)

Every response is cached on disk (``SCTARGET_CACHE_DIR``) keyed by a hash of the request,
so a re-run is deterministic and offline.

Provenance rules (the ``source`` field on every result):
  * ``live``    fetched now;   ``cache``  fetched in an earlier run
  * ``fixture`` bundled illustrative data - used ONLY when ``SCTARGET_OFFLINE=1``
  * ``none``    online mode but the call failed (after retries); the gene gets empty
                evidence and downstream scoring flags it as missing, rather than
                quietly substituting fixture values.

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


# --------------------------------------------------------------------------- drug names
_SALT_WORDS = {
    "MALEATE", "DIMALEATE", "ANHYDROUS", "HYDROCHLORIDE", "DIHYDROCHLORIDE", "MESYLATE", "TOSYLATE",
    "DITOSYLATE", "SODIUM", "POTASSIUM", "CITRATE", "SULFATE", "PHOSPHATE", "ACETATE", "TARTRATE",
    "HYDRATE", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "FUMARATE", "SUCCINATE", "BESYLATE", "MALATE",
    "CAMSYLATE", "XINAFOATE", "ETHANOLATE", "HBR", "HCL", "BROMIDE", "CHLORIDE", "LACTATE", "GLUCONATE",
}
INHIBITOR_TYPES = {"inhibitor", "antagonist", "blocker", "antibody", "inverse agonist", "negative modulator", "suppressor", "inhibitory allosteric modulator", "antisense oligonucleotide"}


def normalize_drug_name(name: str) -> str:
    """'NERATINIB MALEATE' -> 'NERATINIB'; 'DACOMITINIB ANHYDROUS' -> 'DACOMITINIB'."""
    toks = str(name).upper().replace(",", " ").split()
    while len(toks) > 1 and toks[-1] in _SALT_WORDS:
        toks.pop()
    return " ".join(toks)


def summarize_interactions(interactions: list[dict[str, Any]]) -> dict[str, Any]:
    """Collapse DGIdb interactions to one row per compound and count the useful ones.

    A compound counts if it is approved OR has an inhibitor-like interaction type.
    Untyped, unapproved entries (endogenous ligands, bare ChEMBL ids, screening hits)
    are kept in ``n_interactions`` but excluded from the compound lists.
    """
    merged: dict[str, dict[str, Any]] = {}
    for it in interactions:
        if not it.get("drug"):
            continue
        key = normalize_drug_name(it["drug"])
        m = merged.setdefault(key, {"drug": key, "approved": False, "types": set(), "score": 0.0, "sources": set(), "aliases": set()})
        m["approved"] = m["approved"] or bool(it.get("approved"))
        m["types"] |= {str(t).lower() for t in (it.get("types") or []) if t}
        m["score"] = max(m["score"], float(it.get("score") or 0.0))
        m["sources"] |= set(it.get("sources") or [])
        m["aliases"].add(str(it["drug"]).upper())
    rows = []
    for m in merged.values():
        inhib = bool(m["types"] & INHIBITOR_TYPES)
        if not (m["approved"] or inhib):
            continue
        rows.append({**m, "types": sorted(m["types"]), "sources": sorted(m["sources"]), "aliases": sorted(m["aliases"]), "inhibitor_like": inhib})
    rows.sort(key=lambda r: (not (r["approved"] and r["inhibitor_like"]), not r["approved"], not r["inhibitor_like"], -r["score"]))
    return {
        "n_interactions": len(interactions),
        "n_compounds": len(rows),
        "n_inhibitors": sum(1 for r in rows if r["inhibitor_like"]),
        "n_approved_compounds": sum(1 for r in rows if r["approved"]),
        "n_approved_inhibitors": sum(1 for r in rows if r["approved"] and r["inhibitor_like"]),
        "top_compounds": [r["drug"] for r in rows[:8]],
        "compounds": rows,
    }


class KnowledgeClient:
    """Thin HTTP client with disk cache, retries, and offline fixture fallback."""

    def __init__(self, offline: bool | None = None, timeout: float = 30.0, retries: int = 2) -> None:
        self.offline = offline_mode() if offline is None else offline
        self.timeout = timeout
        self.retries = retries
        self.cache = cache_dir()
        self.calls: list[dict[str, Any]] = []

    # ------------------------------------------------------------------ cache
    def _key(self, kind: str, url: str, payload: Any, params: Any) -> Path:
        # the URL is part of the key: for Ensembl GETs the gene symbol lives in the path
        h = hashlib.sha1(json.dumps([kind, url, payload, params], sort_keys=True).encode()).hexdigest()[:20]
        return self.cache / f"{kind}_{h}.json"

    def _request(self, kind: str, method: str, url: str, payload: Any, params: dict[str, Any] | None = None, timeout: float | None = None) -> tuple[Any | None, str]:
        """Returns (json, source) where source in {'cache','live','none'}.

        Retries on timeouts / 5xx with backoff. HTTP 4xx errors are not retried, and the
        response body (which carries GraphQL validation messages) is kept in ``calls``.
        """
        p = self._key(kind, url, payload, params)
        if p.exists():
            with open(p) as fh:
                return json.load(fh), "cache"
        if self.offline:
            return None, "none"
        timeout = timeout or self.timeout
        last_err = ""
        for attempt in range(self.retries + 1):
            t0 = time.time()
            try:
                if method == "GET":
                    r = requests.get(url, params=params, timeout=timeout, headers={"Accept": "application/json"})
                else:
                    r = requests.post(url, json=payload, params=params, timeout=timeout, headers={"Accept": "application/json", "Content-Type": "application/json"})
                if r.status_code >= 400:
                    body = r.text[:400].replace("\n", " ")
                    last_err = f"HTTP {r.status_code}: {body}"
                    if 400 <= r.status_code < 500:
                        break  # not transient
                    raise requests.HTTPError(last_err)
                data = r.json()
                with open(p, "w") as fh:
                    json.dump(data, fh)
                self.calls.append({"kind": kind, "url": url, "ok": True, "attempt": attempt + 1, "s": round(time.time() - t0, 2)})
                return data, "live"
            except Exception as e:  # timeouts, connection errors, 5xx
                last_err = last_err or f"{type(e).__name__}: {e}"
                if attempt < self.retries:
                    time.sleep(2.0 * (attempt + 1) ** 2)  # 2 s, 8 s
                    last_err = ""
                    continue
        self.calls.append({"kind": kind, "url": url, "ok": False, "error": last_err, "s": round(time.time() - t0, 2)})
        return None, "none"

    def _fixture_or_none(self, section: str, key: str, empty: dict[str, Any]) -> dict[str, Any]:
        fx = fixtures().get(section, {})
        if self.offline and key in fx:
            return {**fx[key], "source": "fixture"}
        return {**empty, "source": "none"}

    # ------------------------------------------------------------------ Ensembl
    def orthologs(self, genes: list[str], from_species: str) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        if from_species in ("homo_sapiens", "human"):
            return {g: {"human_symbol": g, "type": "identity", "source": "n/a"} for g in genes}
        fx = fixtures().get("orthologs", {}).get(from_species, {})
        params = {"target_species": "homo_sapiens", "type": "orthologues", "format": "condensed", "sequence": "none", "content-type": "application/json"}
        for g in genes:
            data, src = self._request("ensembl_homology", "GET", f"https://rest.ensembl.org/homology/symbol/{from_species}/{g}", None, params=params, timeout=60)
            if data is not None:
                homs = (data.get("data") or [{}])[0].get("homologies", [])
                hs = [h for h in homs if h.get("species") == "homo_sapiens"]
                if hs:
                    ids = sorted({h["id"] for h in hs})
                    look, _ = self._request("ensembl_lookup", "POST", "https://rest.ensembl.org/lookup/id", {"ids": ids}, timeout=60)
                    symbols = sorted({(look or {}).get(i, {}).get("display_name") for i in ids} - {None}) if look else []
                    if symbols:
                        # prefer the human symbol identical to the query when one-to-many
                        best = g if g in symbols else symbols[0]
                        out[g] = {"human_symbol": best, "all_human_symbols": symbols, "human_ensembl_ids": ids, "type": hs[0].get("type"), "source": src}
                    else:
                        # ortholog exists but the symbol lookup failed: keep the query symbol, flag it
                        out[g] = {"human_symbol": g, "all_human_symbols": [], "human_ensembl_ids": ids, "type": hs[0].get("type"), "symbol_unverified": True, "source": "partial"}
                else:
                    out[g] = {"human_symbol": None, "type": "no_ortholog", "source": src}
            elif self.offline and g in fx:
                out[g] = {**fx[g], "source": "fixture"}
            else:
                # conservative default: assume the same symbol (true for most dog/human genes) but flag it
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
            elif data is not None:
                out[g] = {"categories": [], "interactions": [], "source": src}  # queried live, gene unknown to DGIdb
            else:
                out[g] = self._fixture_or_none("dgidb", g, {"categories": [], "interactions": []})
        return out

    # ------------------------------------------------------------------ Open Targets
    OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
    OT_SEARCH = """
    query($q: String!) { search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 1}) { hits { id name } } }"""
    # progressively simpler queries: if the first is rejected (schema drift), the next is tried
    # Open Targets removed `knownDrugs` from Target (confirmed 2026-10-14); compounds come from
    # DGIdb. Progressive fallback is kept so future schema drift degrades gracefully.
    OT_TARGET_QUERIES = [
        ("tractability_only", """
        query($id: String!) { target(ensemblId: $id) { approvedSymbol tractability { label modality value } } }"""),
        ("symbol_only", """
        query($id: String!) { target(ensemblId: $id) { approvedSymbol } }"""),
    ]

    def open_targets(self, genes: list[str]) -> dict[str, dict[str, Any]]:
        out: dict[str, dict[str, Any]] = {}
        empty = {"ensembl_id": None, "tractability": {"small_molecule": None, "antibody": None}, "known_drugs_count": 0, "max_clinical_phase": None, "known_drugs": []}
        for g in genes:
            data, src = self._request("ot_search", "POST", self.OT_URL, {"query": self.OT_SEARCH, "variables": {"q": g}})
            hits = (((data or {}).get("data") or {}).get("search") or {}).get("hits") or []
            hit = next((h for h in hits if h.get("name", "").upper() == g.upper()), None) if hits else None
            if not hit:
                out[g] = self._fixture_or_none("open_targets", g, empty)
                continue
            merged: dict[str, Any] = {**empty, "ensembl_id": hit["id"], "query_variants_used": []}
            got_any = False
            for name, q in self.OT_TARGET_QUERIES:
                tdata, tsrc = self._request(f"ot_target_{name}", "POST", self.OT_URL, {"query": q, "variables": {"id": hit["id"]}})
                t = ((tdata or {}).get("data") or {}).get("target") or None
                if not t:
                    continue
                got_any = True
                merged["query_variants_used"].append(name)
                merged["source"] = tsrc
                if "tractability" in t and t["tractability"] is not None:
                    sm = [x["label"] for x in t["tractability"] if x.get("modality") == "SM" and x.get("value")]
                    ab = [x["label"] for x in t["tractability"] if x.get("modality") == "AB" and x.get("value")]
                    merged["tractability"] = {"small_molecule": sm[0] if sm else None, "antibody": ab[0] if ab else None}
                if t.get("knownDrugs"):
                    kd = t["knownDrugs"]
                    phases = [r.get("phase") for r in kd.get("rows", []) if r.get("phase") is not None]
                    merged["known_drugs_count"] = kd.get("count", 0)
                    merged["max_clinical_phase"] = max(phases) if phases else None
                    merged["known_drugs"] = sorted({normalize_drug_name(r["drug"]["name"]) for r in kd.get("rows", []) if r.get("drug")})[:25]
                if name == "tractability_only":
                    break  # we have everything the current schema offers
            if got_any:
                out[g] = merged
            else:
                out[g] = {**self._fixture_or_none("open_targets", g, empty), "ensembl_id": hit["id"]}
        return out

    # ------------------------------------------------------------------ PubMed
    def pubmed_count(self, gene: str, disease_terms: list[str]) -> dict[str, Any]:
        q = f"{gene} AND ({' OR '.join(disease_terms)})"
        data, src = self._request("pubmed", "GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", None, params={"db": "pubmed", "term": q, "retmode": "json", "retmax": 5})
        if data is not None:
            res = data.get("esearchresult", {})
            return {"query": q, "count": int(res.get("count", 0)), "pmids": res.get("idlist", []), "source": src}
        key = f"{gene} {disease_terms[0]}" if disease_terms else gene
        fx = fixtures().get("pubmed", {})
        if self.offline and key in fx:
            return {"query": q, **fx[key], "source": "fixture"}
        return {"query": q, "count": None, "pmids": [], "source": "none"}


# ----------------------------------------------------------------------------- tools
def _client(session: AnalysisSession) -> KnowledgeClient:
    if "knowledge_client" not in session.facts:
        session.facts["knowledge_client"] = KnowledgeClient()
    return session.facts["knowledge_client"]


HUMAN_ALIASES = {"homo_sapiens", "human", "hsapiens", "h_sapiens", "homo sapiens"}


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
    res: dict[str, dict[str, Any]] = {}
    remaining = list(genes)
    if str(sp).strip().lower() in HUMAN_ALIASES:
        # human data: symbols are already the knowledge-base vocabulary; no lookup, no network
        res = {g: {"human_symbol": g, "type": "identity", "source": "human dataset"} for g in genes}
        session.facts.setdefault("orthologs", {}).update(res)
        return {"from_species": sp, "mapping": res, "unmapped": [], "assumed_or_unverified": [], "n_from_var_column": 0, "note": "human dataset: identity mapping"}
    # 1) a precomputed mapping shipped inside the AnnData (dataset.human_symbol_var) wins:
    #    it is offline, reproducible, and under the analyst's control (e.g. made once with biomaRt)
    col = session.cfg.dataset.get("human_symbol_var")
    if col and session.adata is not None and col in session.adata.var:
        var = session.adata.var[col]
        for g in genes:
            if g in var.index:
                hs = var.loc[g]
                if isinstance(hs, str) and hs and hs.lower() not in ("nan", "none", ""):
                    res[g] = {"human_symbol": hs, "type": "var_column", "source": f"var[{col}]"}
                else:
                    res[g] = {"human_symbol": None, "type": "no_ortholog", "source": f"var[{col}]"}
        remaining = [g for g in genes if g not in res]
    # 2) everything else goes to Ensembl (cached)
    if remaining:
        res.update(_client(session).orthologs(remaining, sp))
    res = {g: res[g] for g in genes if g in res}
    session.facts.setdefault("orthologs", {}).update(res)
    unmapped = [g for g, v in res.items() if not v.get("human_symbol")]
    assumed = [g for g, v in res.items() if v.get("type") == "assumed_same_symbol" or v.get("symbol_unverified")]
    return {"from_species": sp, "mapping": res, "unmapped": unmapped, "assumed_or_unverified": assumed, "n_from_var_column": sum(1 for v in res.values() if v.get("type") == "var_column")}


@tool(
    name="druggability_lookup",
    description="""
    For human gene symbols, query DGIdb (gene categories such as KINASE / CLINICALLY
    ACTIONABLE, drug-gene interactions with approval status and interaction type) and
    Open Targets (small-molecule / antibody tractability, number of known drugs, max
    clinical phase). Compounds are de-duplicated across salt forms; a compound counts if it
    is approved or has an inhibitor-like interaction type. Returns one row per gene plus the
    top compounds. 'source' tells you whether each value came from a live call, the cache,
    the bundled fixture (offline mode only, illustrative), or 'none' (lookup failed).""",
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
        s = summarize_interactions(d["interactions"])
        cats = [c.upper() for c in d["categories"]]
        rows.append(
            {
                "gene": g,
                "druggable_categories": sorted(set(cats) & wanted) or sorted(set(cats))[:4],
                "is_druggable_class": bool(set(cats) & wanted),
                "n_interactions": s["n_interactions"],
                "n_compounds": s["n_compounds"],
                "n_inhibitors": s["n_inhibitors"],
                "n_approved_compounds": s["n_approved_compounds"],
                "n_approved_inhibitors": s["n_approved_inhibitors"],
                "top_inhibitors": s["top_compounds"],  # name kept for backward compatibility: approved or inhibitor-like compounds, best first
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
    notes = []
    if any(r["source_dgidb"] == "fixture" or r["source_opentargets"] == "fixture" for r in rows):
        notes.append("source=fixture means illustrative bundled data (offline mode), not a live database result.")
    if any(r["source_dgidb"] == "none" for r in rows):
        notes.append("source=none for DGIdb means the lookup failed; treat druggability for those genes as unknown.")
    if any(r["source_opentargets"] == "none" for r in rows):
        notes.append("Open Targets lookup failed for some genes; tractability is unknown for them (DGIdb evidence still applies).")
    return {
        "table": df,
        "offline_mode": client.offline,
        "n_live_failures": len(live_failures),
        "last_failure": (live_failures[-1]["error"][:300] if live_failures else None),
        "note": " ".join(notes) or None,
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
