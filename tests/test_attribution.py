"""condition_attribution: where does a gene's expression actually sit?"""
import numpy as np
import pandas as pd
import pytest
import anndata as ad

from sctarget.session import AnalysisSession
from sctarget.tools import analysis as A
from sctarget.tools import targets as T


@pytest.fixture(scope="module")
def sess(cfg, tmp_path_factory):
    s = AnalysisSession.create(cfg, tmp_path_factory.mktemp("attr"))
    A.load_dataset(s)
    A.run_qc(s, min_genes=100, max_mt_pct=15.0)
    A.integrate_and_cluster(s, reuse_existing=None, use_harmony=False, resolution=1.0, n_top_genes=500)  # the shipped leiden is too coarse on this small set
    A.cluster_composition(s)
    A.differential_abundance(s)
    assert s.facts["differential_abundance"]["disease_enriched_clusters"]
    return s


def test_classes_on_planted_synthetic(sess):
    cl = sess.facts["differential_abundance"]["disease_enriched_clusters"][0]
    out = A.condition_attribution(sess, genes=["ERBB4", "PEG3", "EGFR", "KRT14", "GENE0001"], cluster=cl)
    cls = out["table"]["class"].to_dict()
    assert cls["ERBB4"] == "state_specific" and cls["PEG3"] == "state_specific"
    assert cls["EGFR"] == "shared_with_reference"   # pan-tumour: also high in COSCC (a reference condition) keratinocytes
    assert cls["KRT14"] == "shared_with_reference"  # all keratinocytes, every condition
    assert cls["GENE0001"] in ("weak", "shared_with_reference", "not_evaluable")
    assert sess.facts["condition_attribution"][cl]["EGFR"] == "shared_with_reference"


def test_reference_driven_detection(cfg, tmp_path):
    """A gene expressed only by the reference-condition cells that sit inside a disease-enriched
    cluster must be classified reference_driven (the NETO1 case)."""
    rng = np.random.default_rng(0)
    n = 600
    cond = np.array(["CAA"] * 300 + ["COSCC"] * 200 + ["GINGIVA"] * 100)
    sample = np.array([f"{c}_S{(i % 3) + 1}" for i, c in enumerate(cond)])
    cluster = np.array(["X"] * 240 + ["Y"] * 60 + ["X"] * 60 + ["Y"] * 140 + ["Y"] * 100)  # X = disease-enriched with 60 OSCC spill-over
    X = rng.poisson(0.2, size=(n, 3)).astype(np.float32)
    X[(cluster == "X") & (cond == "COSCC"), 0] = rng.poisson(6, size=60)        # SPILL: only in OSCC cells of X
    X[(cluster == "X") & (cond == "CAA"), 1] = rng.poisson(6, size=240)        # STATE: only in CAA cells of X
    X[:, 2] = rng.poisson(5, size=n)                                            # HOUSE: everywhere
    a = ad.AnnData(X=X, obs=pd.DataFrame({"sample": sample, "condition": cond, "cell_type": "keratinocyte", "leiden": cluster}, index=[f"c{i}" for i in range(n)]), var=pd.DataFrame(index=["SPILL", "STATE", "HOUSE"]))
    a.layers["counts"] = a.X.copy()
    p = tmp_path / "mini.h5ad"
    a.write_h5ad(p)
    from sctarget.config import load_config
    from tests.conftest import ROOT

    c = load_config(ROOT / "config" / "synthetic.yaml", overrides={"dataset": {"path": str(p)}, "privacy": {"pseudonymize_samples": True}})
    s = AnalysisSession.create(c, tmp_path / "run")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing="leiden")
    out = A.condition_attribution(s, genes=["SPILL", "STATE", "HOUSE"], cluster="X", min_cells=10)
    cls = out["table"]["class"].to_dict()
    assert cls["SPILL"] == "reference_driven"
    assert cls["STATE"] == "state_specific"
    assert cls["HOUSE"] == "shared_with_reference"
    # and the scorer punishes it
    s.facts["gene_specificity"] = {"X": [{"gene": g, "log2_specificity": 2.0, "disease_sample_consistency": 1.0, "log2_disease_vs_reference_within_cluster": 1.0} for g in ["SPILL", "STATE"]]}
    s.facts["druggability"] = {g: {"gene": g, "is_druggable_class": True, "ot_sm_tractability": "Approved Drug", "n_approved_compounds": 3, "n_compounds": 3, "top_inhibitors": ["X"]} for g in ["SPILL", "STATE"]}
    rk = T.score_targets(s, genes=["SPILL", "STATE"], cluster="X")["ranking"]
    assert rk.index[0] == "STATE" and rk.loc["SPILL", "score"] < 0.3 * rk.loc["STATE", "score"]
    assert "reference_driven" in rk.loc["SPILL", "flags"]


def test_run_harmony_handles_both_orientations(monkeypatch):
    """harmonypy >= 1 returns cells x PCs (sometimes a torch tensor); older returns PCs x cells."""
    import types

    import anndata as ad
    import numpy as np
    import pandas as pd

    from sctarget.tools import analysis as A

    n, k = 40, 5
    a = ad.AnnData(X=np.zeros((n, 3)), obs=pd.DataFrame({"s": ["a", "b"] * 20}, index=[f"c{i}" for i in range(n)]))
    a.obsm["X_pca"] = np.random.default_rng(0).normal(size=(n, k))

    class Tensorish:
        def __init__(self, arr):
            self.arr = arr

        def detach(self):
            return self

        def cpu(self):
            return self

        def numpy(self):
            return self.arr

    for Z in (np.ones((n, k)), np.ones((k, n)), Tensorish(np.ones((n, k)))):
        fake = types.SimpleNamespace(run_harmony=lambda X, obs, key, random_state, verbose, Z=Z: types.SimpleNamespace(Z_corr=Z))
        monkeypatch.setitem(__import__("sys").modules, "harmonypy", fake)
        out = A._run_harmony(a, "s", 0)
        assert out.shape == (n, k)


def test_integrate_and_cluster_fails_loudly_when_harmony_breaks(cfg, tmp_path, monkeypatch):
    from sctarget.session import AnalysisSession
    from sctarget.tools import analysis as A

    s = AnalysisSession.create(cfg, tmp_path / "h")
    A.load_dataset(s)

    def boom(*a, **k):
        raise ValueError("shape mismatch")

    monkeypatch.setattr(A, "_run_harmony", boom)
    r = A.integrate_and_cluster(s, reuse_existing=None, use_harmony=True)
    assert "error" in r and "Harmony" in r["error"] and "use_harmony=false" in r["hint"]
    r2 = A.integrate_and_cluster(s, reuse_existing=None, use_harmony=False)
    assert "error" not in r2 and r2.get("harmony_used") is False



def test_thin_reference_pool_does_not_decide_class(cfg, tmp_path):
    """A handful of reference cells inside a disease cluster must not flip a gene to
    shared_with_reference / reference_driven; the pool is reported instead."""
    import numpy as np
    import scipy.sparse as sp

    from sctarget.session import AnalysisSession
    from sctarget.tools import analysis as A

    s = AnalysisSession.create(cfg, tmp_path / "thin")
    A.load_dataset(s)
    A.integrate_and_cluster(s, reuse_existing=cfg.dataset.get("reuse_clustering"))
    A._ensure_lognorm(s)
    a = s.adata
    key = s.cluster_key
    labels = a.obs[key].astype(str).to_numpy()
    conds = a.obs[s.condition_key].astype(str).to_numpy()
    is_r = np.isin(conds, s.reference)
    frac = {c: float((conds[labels == c] == s.disease).mean()) for c in np.unique(labels)}
    cl = max(frac, key=frac.get)
    # make the reference pool in this cluster thin but present: keep exactly 2 reference cells
    in_c = labels == cl
    ref_idx = np.where(in_c & is_r)[0]
    if len(ref_idx) < 2:
        extra = [i for i in np.where(is_r)[0] if i not in ref_idx][: 2 - len(ref_idx)]
        labels[extra] = cl
        ref_idx = np.array(list(ref_idx) + extra)
    else:
        for i in ref_idx[2:]:
            labels[i] = "other_ref"
        ref_idx = ref_idx[:2]
    a.obs[key] = labels
    in_c = labels == cl
    thin_frac = len(ref_idx) / in_c.sum()
    floor = thin_frac * 2  # the synthetic clusters are small; test the guard's logic at a matching scale
    g = "ERBB4"
    j = list(a.var_names).index(g)
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    X[ref_idx, j] = X[in_c, j].max() + 1.0  # the thin reference pool is strongly positive
    a.X = sp.csr_matrix(X)
    thin = A.condition_attribution(s, genes=[g], cluster=cl, min_cells=1, min_reference_frac=floor)["table"].loc[g]
    assert str(thin["reference_pool_in_cluster"]).startswith("too_thin")
    assert thin["class"] not in ("shared_with_reference", "reference_driven")
    full = A.condition_attribution(s, genes=[g], cluster=cl, min_cells=1, min_reference_frac=0.0)["table"].loc[g]
    assert full["reference_pool_in_cluster"] == "ok"
    assert full["class"] in ("shared_with_reference", "reference_driven")  # with the guard off the 5 cells decide
