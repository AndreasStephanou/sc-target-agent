# Evaluation cases

Eight cases following `docs/testing-a-skill.md` (straightforward, problem, ambiguous or
difficult, plus one out-of-scope case that must not trigger the skill). Each has an `input.md` (the researcher's request and what the data look like) and an
`expected-behavior.md` (what a scientifically appropriate response does; wording is not
prescribed). `test-cases.json` lists them for the CoPAIR structural checker; it executes nothing.

| case | type | data | what it tests |
|---|---|---|---|
| case-01 | straightforward | public ccRCC scRNA-seq, GSE159115, with a known marker (CA9) | the full workflow: five inputs, key frozen before the run, calibration read correctly, agent result read with the critic |
| case-02 | problem | the same dataset as it was actually encountered | swapped public metadata, QC that removes a lineage, a protein-level target in the key, a stress-artefact cluster at the top of the ranking |
| case-03 | ambiguous | a hypothetical mouse dataset, no annotation, 2 vs 2 samples, no known answer | refusing to guess the condition assignment, the design check and a confounded batch, stating the replicate limitation, running without a key, not over-claiming |
| case-04 | problem | public HCC scRNA-seq, GSE149614, as encountered in the author's reviewer run | a fresh install that passes its tests but cannot cluster, a reference library that expresses the disease programme, a known target that fails the replicate-aware test because of it, a cycling sub-state's cell-cycle markers |
| case-05 | difficult | public ovarian omental-metastasis scRNA-seq, GSE286057, as encountered | choosing a dataset by composition, a key for the pure composition case, a known target set aside by attribution for a real but minor reason, knowing when to stop changing rules |
| case-06 | difficult (new in 0.4.0) | a hypothetical human lung-fibrosis study, 34 libraries from 14 donors, chemistry partly aligned with condition, no known target | classifying a non-malignant design, naming the independent replicate, reading the confounding table, the no-key protocol, composition findings, stating that the design type is untested |
| case-07 | ambiguous (new in 0.4.0) | a hypothetical mouse knockout vs wild-type, 4 vs 4, `marker` mode | refusing the knocked-out gene as the expected target, the sign-convention control, sex genes as decoys, reading a `marker`-mode ranking without druggability, programme shifts vs new states |
| case-08 | out of scope (new in 0.4.0) | bulk RNA-seq, one sample per arm, a request to "validate" a gene | the skill does not trigger; it explains the input it needs and that no output validates a target |

The author ran cases 01, 02, 04 and 05 for real on the public data (results in
`../references/validation-results.md`). Cases 03, 06 and 07 are reasoning tests with fabricated
inputs. Nobody has yet run the 0.4.0 workflow end-to-end on a new dataset; the with-skill vs
without-skill comparison is the next test.

A machine-checkable floor exists in the repository itself: `python scripts/sctarget.py demo`
(offline, synthetic data) must find the planted state, rank the planted target first with its
drugs, nominate no decoys and PASS; `pytest -q` must pass. The skill's helper scripts were run on
the synthetic dataset and the ccRCC object (`../scripts/README.md`).
