from pathlib import Path

import pytest

from sctarget.session import AnalysisSession
from sctarget.tools import analysis as A
from sctarget.tools import knowledge as K
from sctarget.tools import targets as T


@pytest.fixture(scope="module")
def session(cfg, tmp_path_factory):
    s = AnalysisSession.create(cfg, tmp_path_factory.mktemp("tools_run"))
    A.load_dataset(s)
    return s


def test_load_dataset_summary(session):
    d = session.facts["dataset"]
    assert d["disease"] == "CAA" and set(d["reference_conditions"]) == {"COSCC", "GINGIVA"}
    assert d["samples_per_condition"] == {"CAA": 3, "COSCC": 3, "GINGIVA": 3}
    assert all(k.startswith("S") for k in d["cells_per_sample"])  # pseudonymized
    assert "leiden" in d["existing_clusterings"]


def test_qc_then_recompute_clustering_and_find_enriched(session):
    qc = A.run_qc(session, min_genes=100, max_mt_pct=15.0)
    assert qc["cells_after"] <= qc["cells_before"]
    r = A.integrate_and_cluster(session, reuse_existing=None, use_harmony=False, resolution=1.0, n_top_genes=500)
    assert r["method"] == "recomputed" and r["n_clusters"] >= 5
    assert (session.run_dir / "cluster_assignments.csv").exists()
    comp = A.cluster_composition(session)
    assert "dominant_celltype" in comp["table"].columns
    da = A.differential_abundance(session)
    assert len(da["disease_enriched_clusters"]) >= 1
    top = da["table"].index[0]
    assert da["table"].loc[top, "frac_cells_from_disease"] > 0.7


def test_markers_pseudobulk_specificity_scoring(session):
    da = session.facts["differential_abundance"]
    cl = da["disease_enriched_clusters"][0]
    comp = {r["cluster"]: r for r in session.facts["cluster_composition"]}
    dom = comp[cl]["dominant_celltype"]
    sibs = [c for c, r in comp.items() if c != cl and r["dominant_celltype"] == dom]
    mk = A.find_markers(session, cluster=cl, reference_clusters=sibs, n_top=20)
    top = set(mk["top_markers"].index)
    assert {"ERBB4", "PEG3"} & top
    de = A.pseudobulk_de(session, cluster=cl, reference_clusters=sibs)
    assert "error" not in de and de["replicates"]["target"] >= 2
    planted = {"ERBB4", "PEG3", "SYP", "SNAP25", "STMN2", "GAP43", "NEFL", "SYT1", "RAB3A", "DLG4"}
    assert len(planted & set(de["top_up"].index)) >= 4
    sp = A.gene_specificity(session, genes=["ERBB4", "EGFR", "PDGFRB"], cluster=cl)["table"]
    assert sp.loc["ERBB4", "log2_specificity"] > 1.5 > sp.loc["EGFR", "log2_specificity"]
    assert sp.loc["ERBB4", "disease_sample_consistency"] >= 0.66
    en = A.enrichment(session, genes=list(mk["top_markers"].index))
    assert en["terms"].index[0] == "MINI_NEURONAL_SYNAPTIC_SIGNALING"
    K.druggability_lookup(session, genes=["ERBB4", "EGFR", "PDGFRB", "PEG3"])
    rk = T.score_targets(session, genes=["ERBB4", "EGFR", "PDGFRB", "PEG3"], cluster=cl)["ranking"]
    assert rk.index[0] == "ERBB4"
    assert "low_specificity" in rk.loc["EGFR", "flags"]
    assert "not_druggable" in rk.loc["PEG3", "flags"]


def test_knowledge_offline_fixture_and_orthologs(session):
    m = K.map_orthologs(session, genes=["ERBB4", "LOC_FAKE_DOG_GENE", "NOVELGENE1"], from_species="canis_lupus_familiaris")
    assert m["mapping"]["ERBB4"]["human_symbol"] == "ERBB4"
    assert "LOC_FAKE_DOG_GENE" in m["unmapped"]
    assert "NOVELGENE1" in m["assumed_or_unverified"]
    d = K.druggability_lookup(session, genes=["ERBB4", "ZZZ_NOT_A_GENE"])["table"]
    assert d.loc["ERBB4", "n_approved_inhibitors"] >= 1 and d.loc["ERBB4", "source_dgidb"] == "fixture"
    assert d.loc["ZZZ_NOT_A_GENE", "source_dgidb"] == "none" and d.loc["ZZZ_NOT_A_GENE", "n_interactions"] == 0


def test_tool_errors_are_informative(session):
    with pytest.raises(ValueError):
        A.find_markers(session, cluster="not-a-cluster")
    out = A.gene_specificity(session, genes=["NOPE1", "NOPE2"], cluster=session.facts["differential_abundance"]["disease_enriched_clusters"][0])
    assert "error" in out
