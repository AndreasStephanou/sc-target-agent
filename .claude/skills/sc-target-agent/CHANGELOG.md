# Changelog

## 0.4.0 — 2026-11-05 (draft; not yet run end-to-end on a new dataset)

Rewrite of the skill text for accuracy and generality beyond tumour-vs-normal designs; no change
to the repository's code or method (still tag `results-v3`, prompts v6, scoring v3).

- **Design classification first.** New `references/design-types.md`: seven design rows (tumour vs
  tissue of origin; tumour vs a tissue lacking the lineage; tumour vs another tumour; non-malignant
  disease vs healthy; perturbation; paired-only; multi-level) with, for each, what the reference
  contains, the objective (`target` / `marker`), which flags are informative, what a fair key and
  its decoys look like, and whether it has been validated. Workflow step 0 asks for the row.
- **Fifth required input:** the design type and objective, and what the reference represents.
  `objective: marker` (response genes; druggability reported, not scored) is documented for the
  first time.
- **Two helper scripts** (`scripts/`): `design_check.py` (counts location, replicates, imbalance,
  pseudo-replication via a donor column, pairing, confounding of sample-level covariates with the
  condition) and `label_check.py` (per-library marker fractions by deposited label, libraries that
  resemble the opposite condition). Run on the synthetic dataset and on the public ccRCC object.
- **Labels verified before any run, by two routes:** marker-based with a marker that is not the
  expected target, and marker-free via the calibration run's `check_reference_libraries`.
  Exclusions pre-registered in writing.
- **No-key protocol** for the common case of no known answer: positive and negative controls,
  three-run stability, replicate panels, direction; stated as sanity checks, not validation.
- **Attribution reading rule:** a blocking `reference_driven` verdict on a candidate that the
  per-sample and compartment evidence support is a prompt to identify the reference cells; the
  candidate is reported as conditional on them, not dropped (generalised from the ovarian FOLR1 case).
- **Therapeutic direction** is a mandatory per-candidate statement in the output, not an optional ask.
- **Stop conditions** listed in one place.
- **Limitations labelled** G (general) or O (observed in tumour datasets; untested elsewhere), with
  an explicit list of designs that have never been run.
- **`references/flags.md`**: one table per flag family with the generalised and the tumour readings.
- Dataset anecdotes moved out of the workflow into the references; the workflow reads as a procedure.
- Evals: three fabricated cases added (06 non-malignant disease with nested libraries; 07 knockout
  vs wild-type in `marker` mode; 08 out-of-scope bulk request that must not trigger the skill);
  `evals/test-cases.json` for the CoPAIR checker; case 03 expectations reference the design check.
- README gains a plain `Status:` line (the bold form failed the CoPAIR validator) and the
  Anthropic skills acknowledgement the validator asks for.

## 0.3.0 — 2026-11-02

- Fourth validation dataset, reviewer-style from a fresh clone, 2025 GEO deposit chosen by search:
  ovarian cancer omental metastases (GSE286057, FOLR1). First model-driven runs since scoring v3;
  repository tag `results-v3`.
- Workflow: compartment-composition check in step 1 (a gastric atlas was dropped as immune-only);
  GEO per-sample 10x fetcher and sample-sheet builder; attribution status fields
  (`rare_in_reference`, `rare_reference_counterpart`) explained in step 6.
- Method: attribution representativeness guard (lineage-relative); prompts v6 (critic reads the
  pool status; exclusion is a human config decision). Limitation stated: attribution verdicts depend
  on derived-lineage cleanliness.
- Infrastructure: 429 backoff with Retry-After; transient-error retry in the Claude Code backend.

## 0.2.0 — 2026-10-29

- Third validation dataset, run as a reviewer from a fresh clone: human HCC (GSE149614, GPC3).
  Repository tag `results-v2`.
- Workflow step 0 (`check-env`), GEO matrix builder, config-level library exclusion with a recorded
  reason, marker-free library label check (`check_reference_libraries`).
- Scoring v3 / prompts v5: whole-lineage disease-vs-reference contrast with the axis choice
  recorded per gene; proliferation flag for cell-cycle genes; cycling sub-state note.
- Report: per-library pseudobulk panels, state library shares, library-correlation heatmaps,
  candidate evidence table (dot plots and feature UMAPs removed).

## 0.1.0 — 2026-10-27

- First incubator release, status Experimental. Procedure skill over the `sc-target-agent`
  repository (tag `results-v1`): inspect → verify labels → freeze key → calibrate → run → read with
  the critic → optional arms. Validated by the author on canine acanthomatous ameloblastoma
  snRNA-seq (ERBB4) and public human ccRCC scRNA-seq (GSE159115; CA9).
