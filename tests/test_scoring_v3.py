"""Scoring v3 (HCC review run, protocol section 19): the all-cells lineage contrast, the choice of
lineage axis recorded per gene, the proliferation flag, the cycling-sub-state note."""
import copy

from sctarget.agents.orchestrator import run_pipeline
from sctarget.tools.analysis import CELL_CYCLE_GENES, cell_cycle_fraction, pseudobulk_de
from sctarget.tools.targets import SCORING_VERSIONS, lineage_axis_source


def test_cell_cycle_fraction_and_list():
    assert {"TOP2A", "MKI67", "UBE2C", "BIRC5", "PCNA", "MCM2"} <= CELL_CYCLE_GENES
    assert "GPC3" not in CELL_CYCLE_GENES and "ERBB4" not in CELL_CYCLE_GENES and "CA9" not in CELL_CYCLE_GENES
    assert cell_cycle_fraction(["TOP2A", "MKI67", "GPC3", "AFP"]) == 0.5
    assert cell_cycle_fraction([]) == 0.0


def test_within_lineage_contrast_reports_coverage(cfg, tmp_path):
    from sctarget.session import AnalysisSession
    from sctarget.tools.analysis import integrate_and_cluster, load_dataset

    s = AnalysisSession.create(cfg, tmp_path / "ln")
    load_dataset(s)
    integrate_and_cluster(s, reuse_existing=None, resolution=2.0, use_harmony=False)
    key = s.cluster_key
    ok = 0
    for c in s.adata.obs[key].astype(str).unique():
        r = pseudobulk_de(s, cluster=str(c), contrast="disease_vs_reference_within_lineage")
        if "error" in r:
            assert "insufficient pseudobulk replicates" in r["error"], r
            continue
        ok += 1
        cov = r["lineage_coverage"]
        # the all-cells contrast uses every disease cell of the lineage, by construction
        assert cov["frac_of_lineage_disease_cells_in_contrast"] == 1.0, cov
        assert cov["disease_cells_in_contrast"] == cov["disease_cells_in_lineage"] > 0
        rs = pseudobulk_de(s, cluster=str(c), contrast="disease_vs_reference_within_compartment")
        if "error" not in rs:
            assert rs["lineage_coverage"]["frac_of_lineage_disease_cells_in_contrast"] <= 1.0
    assert ok, "lineage contrast never computable"
    # the axis choice is explicit and explained
    for c in s.adata.obs[key].astype(str).unique():
        src, why = lineage_axis_source(s, str(c))
        assert src in (None, "disease_vs_reference_within_compartment", "disease_vs_reference_within_lineage")
        assert "shared_available" in why and "all_available" in why


def test_scoring_v3_records_axis_source_flags_proliferation_and_keeps_target(cfg, tmp_path, monkeypatch):
    import sctarget.tools.analysis as an

    cfg3 = copy.deepcopy(cfg)
    cfg3.raw.setdefault("targets", {})["scoring_version"] = 3
    # the synthetic genes are not real symbols: declare one of the planted decoys a cell-cycle gene
    # so the flag path is exercised (TOP2A is the planted proliferation decoy in the synthetic key)
    monkeypatch.setattr(an, "CELL_CYCLE_GENES", set(an.CELL_CYCLE_GENES) | {"TOP2A"})
    res = run_pipeline(cfg3, model_name="scripted", run_dir=tmp_path / "s3")
    assert res.facts["scoring_version"] == 3
    rows = next(iter(res.facts["ranked_targets"].values()))
    assert all("lineage_axis_source" in r for r in rows)
    srcs = {r["lineage_axis_source"] for r in rows}
    assert srcs <= {None, "disease_vs_reference_within_compartment", "disease_vs_reference_within_lineage"}
    assert rows[0]["gene"] == "ERBB4", [r["gene"] for r in rows[:5]]  # planted target still first
    for r in rows:
        if r["gene"] == "TOP2A":
            assert "proliferation_marker" in r["flags"] and r["flag_penalty"] <= 0.5
    assert "3" in str(res.facts.get("scoring_version")) and 3 in SCORING_VERSIONS
    md = (res.run_dir / "report.md").read_text()
    assert "scoring v3" in md


def test_scoring_v2_unchanged_by_v3_code(cfg, tmp_path):
    cfg2 = copy.deepcopy(cfg)
    cfg2.raw.setdefault("targets", {})["scoring_version"] = 2
    res = run_pipeline(cfg2, model_name="scripted", run_dir=tmp_path / "s2b")
    rows = next(iter(res.facts["ranked_targets"].values()))
    assert all(r["lineage_axis_source"] in (None, "disease_vs_reference_within_compartment") for r in rows)
    assert not any("proliferation_marker" in r["flags"] for r in rows)
    assert rows[0]["gene"] == "ERBB4"
