"""Ambient RNA from a dominant disease state must not be mistaken for disease programs in other cell types."""
from pathlib import Path

import pytest

from sctarget.session import AnalysisSession
from sctarget.tools import analysis as A


@pytest.fixture(scope="module")
def ambient_sess(tmp_path_factory):
    from scripts.make_synthetic import make_synthetic
    from sctarget.config import load_config
    from tests.conftest import ROOT

    out = tmp_path_factory.mktemp("synth_amb")
    make_synthetic(out, cells_per_sample=300, seed=0, mode="state", ambient=0.15)
    cfg = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(out / "synthetic.h5ad")}})
    s = AnalysisSession.create(cfg, tmp_path_factory.mktemp("run_amb"))
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing=None, use_harmony=False, resolution=1.5, n_top_genes=500)
    A.cluster_composition(s)
    A.differential_abundance(s)
    return s


def test_ambient_programs_are_flagged_not_recorded(ambient_sess):
    s = ambient_sess
    dp = A.disease_programs(s, min_cells=100)
    t = dp["table"]
    # the non-keratinocyte clusters show the state genes "up in disease" purely from ambient RNA ...
    non_kc = t[(t["dominant_celltype"] != "keratinocyte") & t["n_sig_up_strong"].notna()]
    assert len(non_kc) >= 2
    assert (non_kc["n_sig_up_strong"] > 0).any()
    # ... but none of them is accepted as a disease-altered program, and the warning fires
    assert not any(t.loc[c, "disease_altered"] for c in non_kc.index)
    assert dp["ambient_warning"] is not None
    assert set(dp["ambient_suspect_clusters"]) >= set(non_kc.index[non_kc["n_sig_up_strong"] >= 5])


def test_ambient_ratio_separates_real_state_genes_from_background(ambient_sess):
    """In the ambient dataset the planted state genes are REAL in the state cluster (ratio >> 3 vs
    immune/stromal cells of the same CAA samples) but AMBIENT in a fibroblast cluster (ratio ~1)."""
    s = ambient_sess
    enriched = s.facts["differential_abundance"]["disease_enriched_clusters"]
    comp = {r["cluster"]: r for r in s.facts["cluster_composition"]}
    fib = next(c for c, r in comp.items() if r["dominant_celltype"] == "fibroblast")
    fib_attr = A.condition_attribution(s, genes=["ERBB4", "PEG3", "SNAP25", "COL1A1"], cluster=fib, min_cells=10)["table"]
    assert fib_attr.loc["ERBB4", "class"] in ("ambient_like", "weak", "not_evaluable")
    assert fib_attr.loc["PEG3", "class"] in ("ambient_like", "weak", "not_evaluable")
    assert fib_attr.loc["COL1A1", "class"] == "shared_with_reference"   # a real fibroblast gene, every condition
    if enriched:
        st_attr = A.condition_attribution(s, genes=["ERBB4", "PEG3", "SNAP25"], cluster=enriched[0], min_cells=10)["table"]
        assert st_attr.loc["ERBB4", "class"] in ("state_specific", "disease_wide")
        assert st_attr.loc["ERBB4", "ambient_ratio"] > 3
