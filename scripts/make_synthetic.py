"""Generate a synthetic snRNA-seq dataset shaped like a CAA / COSCC / gingiva study.

The dataset has a planted, disease-specific keratinocyte state (``neuro_keratinocyte``)
that over-expresses a neuronal program plus ERBB4 and PEG3, mirroring the published
finding. Several *decoy* druggable genes are planted so a naive ranking produces false
leads: EGFR (up in both tumours), PDGFRB (fibroblast marker), KDR/FLT1 (endothelial),
MKI67/TOP2A (proliferation in both tumours). A good pipeline should rank ERBB4 first
and down-weight the decoys for lack of disease specificity.

Usage:
    python scripts/make_synthetic.py --out data/synthetic --cells-per-sample 400 --seed 0
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

# --------------------------------------------------------------------------- design
CONDITIONS = {"CAA": 3, "COSCC": 3, "GINGIVA": 3}  # condition -> n samples

CELL_TYPES = {
    #  type: (base fraction, marker genes)
    "keratinocyte": (0.40, ["KRT5", "KRT14", "KRT17", "KRT15", "TP63", "SFN", "DSP", "PKP1"]),
    "fibroblast": (0.25, ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "FBLN1"]),
    "T_cell": (0.12, ["CD3E", "CD3D", "CD2", "PTPRC", "IL7R", "TRAC", "CD247"]),
    "macrophage": (0.13, ["CD68", "CD163", "LYZ", "C1QA", "C1QB", "MS4A7", "AIF1"]),
    "endothelial": (0.10, ["PECAM1", "VWF", "CDH5", "KDR", "FLT1", "CLDN5", "EMCN"]),
}

# planted disease-specific state (a keratinocyte sub-state)
NEURO_STATE_GENES = ["ERBB4", "PEG3", "SYP", "SNAP25", "STMN2", "GAP43", "NEFL", "SYT1", "RAB3A", "DLG4"]
NEURO_STATE_FRACTION = {"CAA": 0.30, "COSCC": 0.01, "GINGIVA": 0.005}  # of keratinocytes

# decoys: (gene, which conditions it is elevated in, which cell type)
PAN_TUMOR_GENES = ["EGFR", "MKI67", "TOP2A", "KRAS", "CCND1"]      # up in CAA and COSCC keratinocytes
MT_GENES = ["MT-CO1", "MT-CO2", "MT-ND1", "MT-ND2", "MT-ATP6", "MT-CYB"]

N_FILLER_GENES = 1400


def _gene_universe() -> list[str]:
    genes: list[str] = []
    for _, (_, markers) in CELL_TYPES.items():
        genes += markers
    genes += NEURO_STATE_GENES + PAN_TUMOR_GENES + MT_GENES
    # a few extra real druggable symbols as pure noise (should never be nominated)
    genes += ["KIT", "MET", "ALK", "BRAF", "JAK2", "ABL1", "SRC", "CDK4", "CDK6", "PIK3CA"]
    genes += [f"GENE{i:04d}" for i in range(N_FILLER_GENES)]
    seen, out = set(), []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def make_synthetic(out_dir: Path, cells_per_sample: int = 400, seed: int = 0) -> Path:
    rng = np.random.default_rng(seed)
    genes = _gene_universe()
    gidx = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)

    # gene-level baseline mean expression (log-normal), low for most
    base_mu = np.exp(rng.normal(loc=-1.5, scale=1.0, size=n_genes))
    base_mu[[gidx[g] for g in MT_GENES]] *= 25.0  # mito genes abundant

    rows_obs: list[dict] = []
    count_blocks: list[np.ndarray] = []
    truth: list[str] = []

    sample_id = 0
    for cond, n_samples in CONDITIONS.items():
        for _ in range(n_samples):
            sample_id += 1
            sname = f"{cond}_S{sample_id}"
            # per-sample batch effect on every gene (mild)
            batch = np.exp(rng.normal(0, 0.15, size=n_genes))
            # per-sample tweak of cell-type composition
            fracs = np.array([CELL_TYPES[ct][0] for ct in CELL_TYPES])
            fracs = fracs * rng.dirichlet(np.ones(len(fracs)) * 40)
            fracs = fracs / fracs.sum()
            n_per_type = rng.multinomial(cells_per_sample, fracs)

            for ct, n_ct in zip(CELL_TYPES, n_per_type):
                if n_ct == 0:
                    continue
                mu = base_mu * batch
                mu = np.tile(mu, (n_ct, 1))
                # cell-type markers
                for g in CELL_TYPES[ct][1]:
                    mu[:, gidx[g]] *= rng.uniform(8, 15)
                # tumour programs in keratinocytes
                if ct == "keratinocyte" and cond in ("CAA", "COSCC"):
                    for g in PAN_TUMOR_GENES:
                        mu[:, gidx[g]] *= rng.uniform(4, 7)
                # disease-specific neuro state within keratinocytes
                state = np.array([ct] * n_ct, dtype=object)
                if ct == "keratinocyte":
                    is_neuro = rng.random(n_ct) < NEURO_STATE_FRACTION[cond]
                    for g in NEURO_STATE_GENES:
                        fold = rng.uniform(10, 20) if g in ("ERBB4", "PEG3") else rng.uniform(6, 12)
                        mu[is_neuro, gidx[g]] *= fold
                    state[is_neuro] = "neuro_keratinocyte"
                # library size variation + NB sampling
                libsize = np.exp(rng.normal(0, 0.35, size=n_ct))[:, None]
                mu = mu * libsize
                theta = 2.0  # NB dispersion
                p = theta / (theta + mu)
                counts = rng.negative_binomial(theta, p).astype(np.int32)
                count_blocks.append(counts)
                for i in range(n_ct):
                    rows_obs.append({"sample": sname, "condition": cond, "cell_type": ct})
                truth += list(state)

    X = np.vstack(count_blocks)
    obs = pd.DataFrame(rows_obs)
    obs.index = [f"cell_{i:06d}" for i in range(len(obs))]
    var = pd.DataFrame(index=genes)
    var["mt"] = var.index.str.startswith("MT-")

    adata = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()

    # Make it look like a processed object (what a lab would hand over)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=600, batch_key="sample")
    sc.tl.pca(adata, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=seed)
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden", random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(adata, random_state=seed)
    adata.uns["synthetic"] = {"seed": seed, "cells_per_sample": cells_per_sample}

    out_dir.mkdir(parents=True, exist_ok=True)
    h5 = out_dir / "synthetic.h5ad"
    adata.write_h5ad(h5)
    # ground-truth per-cell state is kept OUTSIDE the h5ad so the agent cannot see it
    pd.Series(truth, index=adata.obs_names, name="true_state").to_csv(out_dir / "true_state.csv")
    with open(out_dir / "design.json", "w") as fh:
        json.dump(
            {
                "conditions": CONDITIONS,
                "cell_types": {k: v[1] for k, v in CELL_TYPES.items()},
                "neuro_state_genes": NEURO_STATE_GENES,
                "neuro_state_fraction": NEURO_STATE_FRACTION,
                "pan_tumor_genes": PAN_TUMOR_GENES,
            },
            fh,
            indent=2,
        )
    return h5


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/synthetic")
    ap.add_argument("--cells-per-sample", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args(argv)
    path = make_synthetic(Path(args.out), args.cells_per_sample, args.seed)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
