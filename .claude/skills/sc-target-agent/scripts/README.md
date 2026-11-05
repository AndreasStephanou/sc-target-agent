# Scripts

The pipeline's code is not duplicated here. It lives in the `sc-target-agent` repository
(`../references/repository.md`), which has its own tests, pinned environment and versioned method
history; the skill drives `scripts/sctarget.py` there, so a fix lands once and the version a run
recorded (`config_used.json`) always points at real code.

Two small, deterministic pre-run checks that the repository's `inspect` does not yet perform are
shipped with the skill. They read `obs`/`var` and the count matrix in chunks, never load the whole
matrix, modify nothing, and need only the repository's conda environment (anndata, h5py, numpy,
pandas, scipy). Run them with `python -I` from outside the data directory.

| script | what it prints | when |
|---|---|---|
| `design_check.py H5AD --sample-key S --condition-key C [--case V] [--donor-key D] [--json out.json]` | where integer counts are; samples and cells per condition; imbalance; libraries nested in donors (pseudo-replication) and donors present in both conditions (pairing); every sample-level categorical covariate cross-tabulated against the condition, with aliases of the sample or condition column skipped and perfectly confounded covariates flagged | workflow step 2, after `inspect`, before any run |
| `label_check.py H5AD --sample-key S --condition-key C --genes G1 G2 ... [--layer L] [--var-column COL] [--celltype-key K --celltype V] [--min-gap 0.2] [--csv out.csv]` | per-library fraction of cells positive for each gene, by deposited condition; libraries closer to the opposite condition's median than to their own, for genes whose condition medians differ by at least `--min-gap` | workflow step 3, with a marker the researcher trusts (not the expected target) |

Both were run on the repository's synthetic dataset and on the public ccRCC object (GSE159115,
29k cells, 13 libraries). The design check located the counts, reported 7 vs 6 libraries, skipped
the `library`/`tissue` aliases, and with `--donor-key patient` reported that 5 of 8 patients
contributed two libraries and that 4 patients appear in both conditions (the matched-pair design
the config comments describe); without the donor key it reported `patient` as partly aligned with
condition (Cramér's V 0.62). The label check with CA9 separated tumour from benign libraries in
all 13; with NDUFA4L2 it listed two tumour libraries of lower tumour content as closer to benign
(a purity prompt, not a label error); with a weak marker (PTPRC) it declined to compare because
the condition medians differed by only 12 points.

Candidates for moving into the repository's `inspect`: both checks, and a `--stability` mode of
`run_eval.py` that does not require an answer key (`../../PROPOSED_CODE_CHANGES.md`).
