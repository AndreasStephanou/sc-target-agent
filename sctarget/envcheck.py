"""Environment check: exercise the code paths a real dataset hits that the dense synthetic demo does not.

Why: a fresh install in October 2026 passed all tests and the demo, then failed at clustering on the
first real dataset (pynndescent 0.7 returned a scipy ``csr_array``; scanpy 1.11 called ``.getnnz()``
on it). The approximate-neighbours path only runs above ~8k cells, and Scrublet always uses it, so
nothing small or dense can catch the break. This check builds a sparse random dataset just over that
threshold and runs the real tool chain on it: normalise, HVG, PCA, neighbours, Leiden, Scrublet.

    python -m sctarget.envcheck          # exit 0 = OK; prints each step
"""
from __future__ import annotations

import sys
import time

import numpy as np


def run(n_cells: int = 9000, n_genes: int = 600, seed: int = 0, verbose: bool = True) -> list[str]:
    """Return a list of failure messages (empty = environment OK)."""
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp

    failures: list[str] = []
    rng = np.random.default_rng(seed)
    t0 = time.time()

    def step(name, fn):
        try:
            fn()
            if verbose:
                print(f"  ok   {name:30s} ({time.time() - t0:5.1f}s)")
        except Exception as e:  # noqa: BLE001 - report everything
            msg = f"{name}: {type(e).__name__}: {str(e)[:200]}"
            failures.append(msg)
            if verbose:
                print(f"  FAIL {msg}")

    # three loose "cell types" so neighbours/Leiden have structure; sparse integer counts like 10x data
    centres = rng.gamma(2.0, 1.0, size=(3, n_genes))
    labels = rng.integers(0, 3, n_cells)
    lam = centres[labels] * rng.gamma(2.0, 0.5, size=(n_cells, 1))
    X = sp.csr_matrix(rng.poisson(lam).astype(np.float32))
    a = ad.AnnData(X=X)
    a.obs["sample"] = [f"s{i % 4}" for i in range(n_cells)]
    a.var_names = [f"G{i}" for i in range(n_genes)]

    step("normalize+log1p", lambda: (sc.pp.normalize_total(a, target_sum=1e4), sc.pp.log1p(a)))
    step("highly_variable_genes", lambda: sc.pp.highly_variable_genes(a, n_top_genes=300))
    step("pca", lambda: sc.tl.pca(a, n_comps=20, mask_var="highly_variable", random_state=seed))
    step("neighbors (approximate path)", lambda: sc.pp.neighbors(a, n_neighbors=15, use_rep="X_pca", random_state=seed))
    step("leiden (igraph)", lambda: sc.tl.leiden(a, resolution=1.0, random_state=seed, flavor="igraph", n_iterations=2, directed=False))

    def scrub():
        sub = ad.AnnData(X=X[: n_cells // 2].copy())
        sub.var_names = a.var_names
        sc.pp.scrublet(sub, expected_doublet_rate=0.06, random_state=seed, verbose=False)

    step("scrublet (per-sample)", scrub)

    def harmony():
        import harmonypy as hm

        ho = hm.run_harmony(a.obsm["X_pca"], a.obs, "sample", verbose=False, random_state=seed)
        Z = np.asarray(ho.Z_corr)
        if Z.shape != a.obsm["X_pca"].shape:
            Z = Z.T
        assert Z.shape == a.obsm["X_pca"].shape, f"harmony output shape {Z.shape}"

    step("harmonypy", harmony)

    def deseq():
        from pydeseq2.dds import DeseqDataSet

        counts = np.asarray(X[:8].sum(axis=0)).ravel()
        import pandas as pd

        pb = pd.DataFrame(rng.poisson(np.tile(counts, (6, 1)) + 1), index=[f"p{i}" for i in range(6)], columns=a.var_names)
        meta = pd.DataFrame({"cond": ["a", "a", "a", "b", "b", "b"]}, index=pb.index)
        dds = DeseqDataSet(counts=pb, metadata=meta, design="~cond", quiet=True)
        dds.deseq2()

    step("pydeseq2", deseq)
    return failures


def main(argv=None) -> int:
    print("sc-target-agent environment check (sparse, >8k cells: the paths the synthetic demo does not hit)")
    failures = run()
    if failures:
        print("\nENVIRONMENT NOT USABLE for real data:")
        for f in failures:
            print("  -", f)
        print("\nMost likely cause: package versions drifted. Rebuild from the pinned environment.yml:\n  conda env create -f environment.yml   (or: pip install 'pynndescent<0.7' 'umap-learn<0.6')")
        return 1
    print("environment OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
