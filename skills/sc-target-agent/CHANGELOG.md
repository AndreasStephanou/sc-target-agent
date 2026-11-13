# Changelog

## 0.4.5 — 2026-11-13 (exact versions and a container image for the repository)

No helper-script change. Documentation of repository protocol 23s: `environment/*.lock` holds the
exact versions behind the current results tag (written by `scripts/freeze_env.sh`), a `Dockerfile`
installs from those lock files only, and CI builds, tests and publishes the image to GHCR on every
push and tag. Workflow step 1 and `references/repository.md` give the lock-file and container
commands beside the conda one; the agent run stays on the host because it uses the Claude Code
CLI login. Conda remains the primary path.

## 0.4.4 — 2026-11-10 (the model behind a run is recorded; one run on a named model)

No helper-script change. Documentation of one repository fix and one run (repository tag
results-v7; prompts v13, scoring v3):

- Workflow step 6: launch the agent with an explicit model id, `--model claude-code:<model id>`
  (for example `claude-code:claude-fable-5-1`). Bare `claude-code` runs on the Claude Code CLI's
  own default model, which depends on the user's settings and changes over time; every agent run of
  rounds 0.4-0.4c was launched that way and the logs recorded only the backend name. Repository
  protocol 23r: the run log (`usage.model` on every model call), the run facts (`model_ids`) and the
  report header now carry the model id the CLI reports, whether or not one was requested.
- `validation-results.md`: the CAA primary arm on Claude Fable 5.1 (prompts v13): PASS, ERBB4 4th
  within its state, FGFR2 nominated and flagged, 2.1M input tokens, top-10 Jaccard 0.54 with the
  v12 run; same answer, more tokens, more microenvironment nominations, each flagged.
- Context: the skill package is about 20K tokens with every reference open and does not approach
  any limit. The agent roles do grow: per-role peaks across the 16 recorded runs reach 217K input
  tokens in the nominator on the largest datasets, because the stateless CLI backend re-sends every
  tool result each turn. The API backend therefore needs a model with a context window above 200K
  for datasets of 100k cells and more; the configs' `agent.model` default (claude-sonnet-4-5, 200K)
  is not that model.
- Maintainer line (`README.md`, `copair.yaml`) names both affiliations: De Vlaminck lab and
  Sethupathy lab, Cornell University.

## 0.4.3 — 2026-11-08 (round 0.4c: the development dataset, canine ameloblastoma, under skill 0.4)

No helper-script change. Documentation of the round and of four repository fixes the skill relies on
(repository tag results-v6; prompts v13, scoring v3):

- Workflow step 7 adds: read the agent notes (report section 5) for what each role says about its
  own inputs, not only its conclusions. The round's most consequential finding came from one such
  sentence: the critic reported that the candidate list it had been given was cut short. Every
  stage briefing had been truncated at a flat character cap since the method was frozen, in every
  agent run; the evaluator could not have shown it (repository protocol 23m, prompts v12: a
  briefing never drops an item).
- `validation-results.md`: round 0.4c section and the CAA table rows (three floors, five agent
  runs, including a de novo arm with Harmony on the library and derived lineages: ERBB4 1st within
  its state there); the statement that the "decoys nominated unflagged" counts of rounds 0.3-0.4b
  measured a critic working from a partial list; how to read a de novo arm's stability.
- Repository fix recorded (23p, prompts v13): when a derived lineage is itself condition-exclusive
  (>= 90% disease cells), the within-lineage abundance test is void, not failed; the scripted critic
  and the critic prompt say so instead of flagging the state. Found on the CAA de novo arm, where
  the CAA keratinocytes formed their own lineage.
- Repository fix recorded (23q, report only): a state recorded before `merge_clusters` folded its
  cluster into a merged label was drawn as empty in the report figures; the figures now draw the
  merged pieces as one state and the report gains a map of the recorded states on the embedding
  under the final cluster key. `scripts/rebuild_report.py` applies it to an existing run; no number
  changes. If a report shows a state with 0 cells, rebuild it.
- Repository fixes recorded: an empty entry in a precomputed ortholog table is "unknown", not "no
  ortholog", and a symbol Ensembl Compara cannot pair but that names a human gene is kept on a
  flagged assumption (23n); `run_eval` rewrites the report's evaluation section (23o).
- Lesson for step 3 written into the text: a whole-library marker check cannot see a reference
  library that is partly case-like (10% marker-positive against 60-87% in the case libraries); the
  floor's marker-free lineage-wise correlation check can, and the exclusion arm is pre-registered
  from that result before the agent run on the arm.

## 0.4.2 — 2026-11-05 (round 0.4b: HCC, ovarian and FLC under skill 0.4; in progress)

Helper-script fixes from the pre-run checks on the three rebuilt tumour datasets; no change to the
repository's method (prompts v8, scoring v3):

- `label_check.py` gains a second rule, independent of the gap between condition medians: a library
  carrying the other condition's marker at >= 3x its own condition's median (`--outlier-fold`,
  `--outlier-gap`) is listed. The first rule is blind when the case cells are a minority of every
  library (ovarian: EPCAM medians 20% vs 0.2%); the second names the two tumour-bearing "distant
  omentum" libraries and the healthy donor's benign epithelial inclusion before any run. Workflow
  step 3 says when to read which list.
- `design_check.py` distinguishes a covariate **nested within the condition** (several levels, each
  inside one condition: regions, sub-sites, grades) from one perfectly confounded with it, and
  counts pseudo-replication per donor-condition pair so a paired design is no longer reported as
  over-counted replication. Workflow step 2 explains the nested case.
- Attribution reading rule extended (step 6): when "the reference holds the lineage" is answered by
  a derived lineage, the settling move is a reference-only subcluster of the state's cluster
  (`subcluster(..., condition=<reference>)`, repository protocol 23k, prompts v10); the omentum case
  is given as the worked example. Limitation added: with one condition per library the integration
  on the library can absorb the disease signal (FLC tumour nuclei co-clustered with normal
  hepatocytes); the expression route is what finds the state, so it runs in every analysis.
- Repository changes this round that the skill relies on: scripted floor's effect pool taken over the
  whole stored list (23h); evaluator prints the lineage of the state carrying a target and honours an
  optional `lineage:` on an expected state (23i; new keys should name it); prompts v9/v10 report a
  reference-driven candidate with per-sample support as conditional rather than dropping it (23j).

## 0.4.1 — 2026-11-05 (validation round in progress; working copy of the repository, not pushed)

Skill text follows the method changes made during the round (repository protocol sections 23-23f,
prompts v7-v8, scoring unchanged at v3):

- Expression route runs on the base clustering; same-lineage clusters do not compete in the
  programme criterion; programme-basis attribution and consistency follow the within-cluster
  contrast; one nomination per gene (23, 23b, 23d). The synthetic `marker`-mode calibration, which
  failed before the round, now passes.
- Candidate pool keeps the largest significant effects as well as the most significant statistics
  (23f); `flags.md` documents the `pool` field.
- Scripted floor and analyst prompt cover the top enriched cluster of every lineage; the expression
  route runs once in every analysis; a programme cluster of the case lineage is subclustered when a
  rare state is plausible (23e, prompts v8).
- `inspect` matches canonical markers case-insensitively (mouse objects had an empty composition
  line) and the epithelial warning is conditional on the question (23c).
- Label checks run within the marker's compartment in non-malignant and perturbation designs;
  covariates known only from the deposit are added to the sample sheet so they reach `obs`.
- Third key-drafting failure mode recorded (time-point-dependent read-out).
- Evaluator reports a target's rank within the state that carries it (23g); results tables report every expected state.
- Round outcomes: IPF SPP1 macrophage state recovered strictly under prompts v8 (1 of 3 key states; ITGB6 unranked); doxorubicin programme found, expected read-out absent at the sampled time; ccRCC regression clean with a two-run stability pair.
- Report rewritten: executive summary (design statement, key numbers, top candidates with a
  direction column), evaluation section, figure index; publication style with PNG + PDF + SVG per
  figure; new design, abundance, candidate-evidence and volcano figures; `rebuild_report.py`.

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
