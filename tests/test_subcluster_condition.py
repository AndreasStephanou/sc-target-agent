"""Protocol 23k: `subcluster` with `condition` splits only the cells of one condition inside the
chosen clusters, as a diagnostic that leaves the session cluster key unchanged. It is the settling
move when condition_attribution calls a candidate reference_driven: what are the reference cells
inside the disease state's cluster?
"""
from __future__ import annotations

from sctarget.session import AnalysisSession
from sctarget.tools.analysis import find_markers, integrate_and_cluster, load_dataset, subcluster


def test_condition_restricted_subcluster_is_diagnostic(cfg, tmp_path):
    s = AnalysisSession.create(cfg, tmp_path / "sc")
    load_dataset(s)
    integrate_and_cluster(s, reuse_existing=None, resolution=1.0, use_harmony=False)
    a = s.adata
    key = s.cluster_key
    cond = a.obs[s.condition_key].astype(str)
    # pick the cluster with the most reference-condition cells
    ref = s.reference[0] if isinstance(s.reference, (list, tuple)) else str(s.reference)
    counts = a.obs.loc[(cond == str(ref)).to_numpy(), key].astype(str).value_counts()
    cl = str(counts.index[0])
    n_ref_in_cl = int(counts.iloc[0])
    before = a.obs[key].astype(str).copy()

    r = subcluster(s, clusters=[cl], condition=str(ref), resolution=0.5)
    new_key = r["cluster_key"]
    assert new_key != key and new_key in a.obs
    assert s.cluster_key == key, "a condition-restricted split must not change the session cluster key"
    assert r["condition"] == str(ref) and "note" in r
    labels = a.obs[new_key].astype(str)
    split = labels.str.startswith(f"{cl}.{ref}.")
    # exactly the reference cells of the cluster were relabelled; everything else is unchanged
    assert int(split.sum()) == n_ref_in_cl == sum(r["subcluster_sizes"].values())
    assert (labels[~split.to_numpy()] == before[~split.to_numpy()]).all()
    assert (cond[split.to_numpy()] == str(ref)).all()
    # the new labels are usable for markers on the new key
    lab = next(iter(r["subcluster_sizes"]))
    m = find_markers(s, cluster=lab, cluster_key=new_key)
    assert "error" not in m


def test_condition_must_be_a_level_of_the_condition_column(cfg, tmp_path):
    s = AnalysisSession.create(cfg, tmp_path / "sc2")
    load_dataset(s)
    integrate_and_cluster(s, reuse_existing=None, resolution=1.0, use_harmony=False)
    cl = str(s.adata.obs[s.cluster_key].astype(str).value_counts().index[0])
    try:
        subcluster(s, clusters=[cl], condition="NoSuchCondition")
    except ValueError as e:
        assert "not a value of" in str(e)
    else:
        raise AssertionError("an unknown condition must raise")
