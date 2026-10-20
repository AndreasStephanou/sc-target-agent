"""Adaptive QC, merge_clusters, split-state reporting, dissociation-stress flag."""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from sctarget.session import AnalysisSession
from sctarget.tools import analysis as A


def test_adaptive_qc_keeps_mito_rich_lineage(cfg, tmp_path):
    """A lineage with a high but normal mitochondrial fraction survives adaptive QC and is reported removed under fixed 10%."""
    s = AnalysisSession.create(cfg, tmp_path / "qc")
    A.load_dataset(s)
    a = s.adata
    # make one annotated type mitochondria-rich by adding counts to MT- genes (or inventing one)
    ct = a.obs[cfg.dataset["celltype_key"]].astype(str)
    target = ct.value_counts().index[-1]
    m = (ct == target).to_numpy()
    mt_genes = [g for g in a.var_names if g.startswith("MT-")]
    assert mt_genes, "synthetic data should carry MT- genes"
    layer = cfg.dataset.get("counts_layer")
    X0 = a.layers[layer] if layer and layer in a.layers else a.X
    X = sp.lil_matrix(X0) if not sp.issparse(X0) else X0.tolil()
    j = list(a.var_names).index(mt_genes[0])
    tot = np.asarray(X.sum(axis=1)).ravel()
    X[np.where(m)[0], j] = (0.09 * tot[m]).astype(int)  # lifts that lineage to ~15% mito: normal for many epithelia, above a 10% cutoff
    X = X.tocsr()
    if layer and layer in a.layers:
        a.layers[layer] = X
    a.X = X.copy()
    import copy

    def fresh(name):
        t = AnalysisSession.create(cfg, tmp_path / name)
        A.load_dataset(t)  # sets up the pseudonymiser
        t.adata = a.copy()
        return t

    s_fixed = fresh("qc_fixed")
    r_fixed = A.run_qc(s_fixed, mode="fixed", max_mt_pct=10.0)
    assert r_fixed["removed_by_annotated_type"][target]["frac_removed"] > 0.9
    assert "warning" in r_fixed
    s_ad = fresh("qc_ad")
    r_ad = A.run_qc(s_ad, mode="adaptive")
    assert r_ad["removed_by_annotated_type"][target]["frac_removed"] < 0.5
    assert all(v["max_mt_pct"] <= 50 for v in r_ad["thresholds_used"].values())


def test_stress_fraction_and_marker_flag():
    assert A.stress_fraction(["FOS", "JUN", "HSPA1A", "ERBB4"]) == 0.75
    assert A.stress_fraction(["ERBB4", "PEG3"]) == 0.0


def test_merge_clusters_reports_union_abundance(cfg, tmp_path):
    s = AnalysisSession.create(cfg, tmp_path / "mg")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing=None, resolution=2.0, use_harmony=False)
    A.cluster_composition(s)
    da = A.differential_abundance(s)
    tab = da["table"]
    lab = s.adata.obs[s.cluster_key].astype(str)
    cond = s.adata.obs[s.condition_key].astype(str)
    dis = [c for c in lab.unique() if (cond[lab == c] == s.disease).mean() >= 0.8 and (lab == c).sum() >= 20]
    if len(dis) < 2:
        import pytest

        pytest.skip("synthetic clustering has < 2 disease-dominated clusters at this resolution")
    r = A.merge_clusters(s, clusters=dis[:2])
    assert "error" not in r, r
    assert r["merged_cluster"] in set(s.adata.obs[s.cluster_key].astype(str))
    assert r["n_cells"] == sum(int((lab == c).sum()) for c in dis[:2])
    assert "log2_ratio" in r["merged_abundance"]
    # refuses to merge a mixed cluster by default
    mixed = [c for c in lab.unique() if (cond[lab == c] == s.disease).mean() < 0.5]
    if mixed:
        r2 = A.merge_clusters(s, clusters=[dis[0], mixed[0]])
        assert "error" in r2


def test_split_state_candidates_reported(cfg, tmp_path):
    s = AnalysisSession.create(cfg, tmp_path / "sp")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing=None, resolution=2.0, use_harmony=False)
    da = A.differential_abundance(s)
    assert "split_state_candidates" in da and isinstance(da["split_state_candidates"], list)
