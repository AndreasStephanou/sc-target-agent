"""Parsing of real-world knowledge API payloads (captured from a live DGIdb call, 2026-10-14)."""
from sctarget.tools.knowledge import KnowledgeClient, normalize_drug_name, summarize_interactions

# Abbreviated but faithful slice of the live DGIdb response for ERBB4: note the many
# approved drugs with EMPTY interaction types, salt-form duplicates, endogenous ligands
# and bare ChEMBL ids - all of which the summariser must handle.
LIVE_ERBB4 = [
    {"drug": "AEE788", "approved": False, "types": ["inhibitor"], "score": 0.086, "sources": ["GuideToPharmacology"]},
    {"drug": "ERLOTINIB", "approved": True, "types": [], "score": 0.0215, "sources": ["DTC", "ClinPGx"]},
    {"drug": "CHEMBL:CHEMBL578061", "approved": False, "types": [], "score": 0.029, "sources": ["DTC"]},
    {"drug": "AFATINIB", "approved": True, "types": [], "score": 0.1009, "sources": ["CIViC", "CGI"]},
    {"drug": "AFATINIB DIMALEATE", "approved": True, "types": ["inhibitor"], "score": 0.2439, "sources": ["ChEMBL"]},
    {"drug": "LAPATINIB", "approved": True, "types": [], "score": 0.1423, "sources": ["DTC", "CIViC", "DoCM"]},
    {"drug": "DACOMITINIB ANHYDROUS", "approved": True, "types": ["inhibitor"], "score": 0.2582, "sources": ["ChEMBL"]},
    {"drug": "DACOMITINIB", "approved": True, "types": ["inhibitor"], "score": 0.2439, "sources": ["ChEMBL"]},
    {"drug": "NERATINIB MALEATE", "approved": True, "types": ["inhibitor"], "score": 0.2439, "sources": ["ChEMBL"]},
    {"drug": "NERATINIB", "approved": True, "types": ["inhibitor"], "score": 0.0585, "sources": ["ChEMBL", "Drug Repurposing Hub"]},
    {"drug": "NEUREGULIN-1", "approved": False, "types": [], "score": 0.3658, "sources": ["GuideToPharmacology"]},
    {"drug": "HB-EGF", "approved": False, "types": [], "score": 0.3658, "sources": ["GuideToPharmacology"]},
    {"drug": "PALIPERIDONE", "approved": True, "types": [], "score": 0.1045, "sources": ["ClinPGx"]},
    {"drug": "IBRUTINIB", "approved": True, "types": ["inhibitor"], "score": 0.017, "sources": ["GuideToPharmacology"]},
]


def test_normalize_drug_name():
    assert normalize_drug_name("NERATINIB MALEATE") == "NERATINIB"
    assert normalize_drug_name("Afatinib Dimaleate") == "AFATINIB"
    assert normalize_drug_name("DACOMITINIB ANHYDROUS") == "DACOMITINIB"
    assert normalize_drug_name("CANERTINIB DIHYDROCHLORIDE") == "CANERTINIB"
    assert normalize_drug_name("BOTULINUM TOXIN TYPE A") == "BOTULINUM TOXIN TYPE A"
    assert normalize_drug_name("SODIUM") == "SODIUM"  # never strip to nothing


def test_summarize_live_payload_counts_untyped_approved_drugs_and_dedupes_salts():
    s = summarize_interactions(LIVE_ERBB4)
    assert s["n_interactions"] == len(LIVE_ERBB4)
    names = set(s["top_compounds"]) | {r["drug"] for r in s["compounds"]}
    # salt forms collapsed
    assert {"NERATINIB", "AFATINIB", "DACOMITINIB"} <= names
    assert not any("MALEATE" in n or "ANHYDROUS" in n for n in names)
    # approved-but-untyped drugs are kept as compounds
    assert {"ERLOTINIB", "LAPATINIB"} <= names
    # ligands / bare ids / unapproved-untyped are excluded
    assert not ({"NEUREGULIN-1", "HB-EGF", "CHEMBL:CHEMBL578061"} & names)
    assert s["n_approved_compounds"] >= 7
    assert s["n_approved_inhibitors"] >= 4
    # approved + typed inhibitors rank first
    assert s["top_compounds"][0] in {"DACOMITINIB", "NERATINIB", "AFATINIB"}


def test_online_mode_never_uses_fixtures_on_failure(tmp_path, monkeypatch):
    """When a live call fails in online mode the result is 'none', not a fixture value."""
    monkeypatch.setenv("SCTARGET_CACHE_DIR", str(tmp_path))
    c = KnowledgeClient(offline=False, timeout=0.001, retries=0)
    # unroutable address -> immediate failure
    monkeypatch.setattr(c, "_request", lambda *a, **k: (None, "none"))
    d = c.dgidb(["ERBB4"])
    assert d["ERBB4"]["source"] == "none" and d["ERBB4"]["interactions"] == []
    o = c.open_targets(["ERBB4"])
    assert o["ERBB4"]["source"] == "none"
    m = c.orthologs(["ERBB4"], "canis_lupus_familiaris")
    assert m["ERBB4"]["type"] == "assumed_same_symbol" and m["ERBB4"]["source"] == "none"


def test_offline_mode_uses_fixtures(tmp_path, monkeypatch):
    monkeypatch.setenv("SCTARGET_CACHE_DIR", str(tmp_path))
    c = KnowledgeClient(offline=True)
    assert c.dgidb(["ERBB4"])["ERBB4"]["source"] == "fixture"
    assert c.open_targets(["ERBB4"])["ERBB4"]["source"] == "fixture"
    assert c.orthologs(["ERBB4"], "canis_lupus_familiaris")["ERBB4"]["source"] == "fixture"


def test_cache_key_includes_url(tmp_path, monkeypatch):
    """Regression: two GETs that differ only in the URL path must not share a cache entry."""
    monkeypatch.setenv("SCTARGET_CACHE_DIR", str(tmp_path))
    c = KnowledgeClient(offline=True)
    k1 = c._key("ensembl_homology", "https://rest.ensembl.org/homology/symbol/dog/PEG3", None, {"a": 1})
    k2 = c._key("ensembl_homology", "https://rest.ensembl.org/homology/symbol/dog/KRT14", None, {"a": 1})
    assert k1 != k2
