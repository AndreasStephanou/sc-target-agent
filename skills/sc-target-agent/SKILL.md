---
name: sc-target-agent
description: Run auditable, agent-driven target discovery on a two-condition single-cell or single-nucleus RNA-seq dataset with the sc-target-agent pipeline. Two objectives - druggable targets (tumour vs normal, non-malignant disease vs healthy) or condition-specific response genes (treated vs control, knockout vs wild-type, infected vs mock) - with critic flags and a frozen answer key when a target is already known. Use when a researcher has an .h5ad or per-sample 10x matrices with per-library condition labels and at least two (preferably three) biological replicates per condition. Do not use for single-replicate designs, for more than two conditions or a continuous covariate in one run, or to declare a candidate validated.
---

# Purpose

Walk a researcher through the sc-target-agent pipeline on their own dataset: classify the study
design, check the data and labels, freeze an answer key if one exists, run a free calibration pass,
run the agent, and read the ranking with the critic. The pipeline is a four-role agent (planner,
analyst, nominator, critic) that acts only through deterministic single-cell tools, logs every tool
call with the exact payload the model saw, and never sends per-cell data or real sample identifiers
to the model. The code lives in the `sc-target-agent` repository (`references/repository.md`);
this skill drives its scripts and does not re-implement the analysis. Commands written as
`python scripts/...` run inside the repository clone; the two pre-run checks that ship with this
skill are written as `<skill>/scripts/...` and run from the skill's own directory.

The pipeline is design-agnostic in its code (any two conditions, any species, `objective: target`
or `marker`) but it has only been validated on four tumour datasets. This skill therefore does two
things the previous version did not: it makes the researcher classify the design before anything
runs (`references/design-types.md`), and it labels every rule as *general* or *observed in one
dataset family* so that a non-tumour user knows what has and has not been tested.

Do not use this skill to:

- analyse a design with fewer than two independent biological replicates per condition;
- compare more than two conditions or a continuous covariate in one run (split into pairwise runs);
- declare any nominated gene a validated target or mechanism: every output is a hypothesis for
  orthogonal follow-up (IHC, functional assay, perturbation), and the report says so;
- edit an answer key after seeing results, or change code, prompts or scoring weights to make one
  dataset pass.

## Required inputs

Ask for, or locate in the object, these before anything else. If one is missing, say which and
stop; do not infer the case/reference assignment from a column name or a letter.

1. **Raw counts.** An `.h5ad` with integer counts in `X` or a layer, per-sample 10x matrices, or a
   GEO "one matrix + metadata" deposit (`scripts/build_from_matrix.py --dry-run`, then build).
   `<skill>/scripts/design_check.py` reports where integer counts are without loading the matrix.
2. **The sample column**: one value per library. State also what the **independent biological unit**
   is (donor, animal, culture). If one donor contributes several libraries, the sample column
   over-counts replication; decide with the researcher whether to pool to donor or to state the inflation.
3. **The condition column**, which value is the case, and **what the reference represents**
   (tissue of origin, a different disease, a different tissue, vehicle, wild-type). This decides the
   design type below and how attribution verdicts should be read.
4. **Species.** Human needs nothing further; other species are mapped to human symbols for drug
   lookups and unmapped genes are marked.
5. **The design type and objective** (`references/design-types.md`): which row the study falls in,
   and whether the question is *druggable targets* (`targets.objective: target`) or
   *condition-specific response genes* (`targets.objective: marker`, druggability reported but not scored).

Optional: a cell-type annotation (lineages are derived from the data otherwise); a known target
or published marker set for an answer key; the researcher's tolerance for run time and cost.

## Workflow

0. **Classify the design.** Pick the row in `references/design-types.md` and write it into the
   config as comments and `dataset.condition_descriptions`. The row sets the objective, which flags
   are informative, what the decoys should be, and which parts of the method have been tested on
   that design. Tell the researcher plainly when their design is in the *untested* set.
1. **Check the environment.** `python scripts/sctarget.py check-env` exercises the sparse,
   large-data tool chain on a synthetic object. A fresh install once passed every test and the demo
   yet could not cluster a real dataset. If it fails, rebuild from the pinned `environment.yml`.
2. **Inspect, then check the design.** `python scripts/sctarget.py inspect data/<file>.h5ad --name
   <name>` guesses the columns, prints the compartment composition per condition (immune, epithelial,
   fibroblast, endothelial marker-positive fractions), and scaffolds `config/<name>.yaml` and a key
   template. Then run `<skill>/scripts/design_check.py` with the sample, condition and, if any,
   donor columns. Act on every warning: no integer counts; fewer than three samples per condition;
   more than two levels; sample imbalance; libraries nested in donors; paired donors; and any
   sample-level covariate (batch, chemistry, sex, site) that is **perfectly confounded** with the
   condition - in that case no result can be attributed to the condition and the researcher must
   decide whether to proceed at all. A covariate whose levels are **nested** within the condition
   (tissue regions, sub-sites, grades: several levels, each inside one condition) is reported
   separately: harmless if its levels are biological sub-groups, confounding if they are processing
   batches; say which in the write-up. If the compartment the question is about is nearly absent in
   both conditions (an immune-sorted deposit, epithelium lost at dissociation), stop: a target test
   aimed at cells the data barely hold is wasted. For a GEO series with per-sample files, fetch two
   libraries first (`scripts/fetch_geo_10x.py`, `scripts/build_from_10x_samples.py`) and look.
3. **Verify the labels before any run.** Public metadata has been wrong in two of three public
   validation sets. Run `<skill>/scripts/label_check.py` with a marker the researcher trusts
   for the case condition - ideally **not the expected target**, so the check stays independent of
   the test - and read the libraries it lists as resembling the opposite condition. In a tumour the
   marker can be read over the whole library; in a non-malignant disease or a perturbation (design
   rows 4-5) the disease shifts composition rather than replacing the tissue, so whole-library
   fractions are uninformative (IPF: MMP7 medians within 13 points) and the check must run **within
   the marker's compartment** (`--celltype-key population --celltype Epithelial`: MMP7 36% vs 10%).
   When the case cells are a minority of every library and there is no annotation to restrict to
   (tumour nests in omentum: EPCAM medians 20% vs 0.2%, below the 20-point gap), the first rule
   calls every marker uninformative; read the script's **second list** instead - libraries carrying
   the other condition's marker at 3x or more their own condition's median. That list named the two
   tumour-bearing "distant omentum" libraries and a healthy donor's benign epithelial inclusion
   (PAX8 5%), the cells that later blocked FOLR1 in round 0.3. Covariates that exist only in the deposit's metadata (processing batch, chemistry, site) must be
   written into `obs` at build time (a column in the sample sheet; `build_from_10x_samples.py`
   carries every sheet column) or the design check cannot see them. The calibration
   run in step 5 adds the marker-free `check_reference_libraries` (within each lineage, does any
   library correlate with the opposite condition?). Decide **now, in writing**, how a doubtful
   library will be handled: keep it in the primary run, and pre-register a second arm with
   `dataset.exclude_samples: [...]` and `exclude_samples_reason`, both recorded in the run. Never
   relabel a library silently, and never clean the reference with the expected target and call the
   result blind.
4. **Decide about an answer key.** If a target (or a response gene) for this condition is already
   known and independently validated, write `data/ground_truth/<name>.yaml` now, commit it, tag it
   (`git tag <name>-key-frozen`), and never open it again after the first run. Rules that
   generalise across designs:
   - the expected target must be **transcriptionally distinctive in the case condition**; a gene
     regulated at the protein level (HIF-2alpha/EPAS1) belongs in `acceptable_targets`, with its
     transcriptional read-out as the expected target;
   - state markers must be **consistent across replicates in this dataset family**, not a textbook
     list (the HCC key's literature markers were patient-heterogeneous and never significant);
   - decoys are chosen **per design type** (`references/design-types.md`): proliferation and lineage
     genes for tumours; generic injury, interferon and stress programmes for non-malignant disease;
     the perturbed gene itself and the trivial downstream read-out for perturbations;
   - say in the key what the test is: "ranked first" when the reference holds the lineage,
     "found and ranked among the real targets" when it does not (pure composition case).
   With no known answer, skip the key and follow the **no-key protocol** in `references/design-types.md`:
   positive controls that must behave (lineage markers classed `shared_with_reference`, not nominated;
   a known condition-associated gene, if any, visible somewhere in the tables), three repeat runs for
   stability, and the explicit statement that the evaluator does not apply.
5. **Calibrate (no model tokens).** `python scripts/sctarget.py calibrate --config config/<name>.yaml
   [--ground-truth ...]` runs a fixed scripted procedure through the same tools, including the
   marker-free library check. Read `runs/<name>_calibrate/report.md` in this order: QC loss per
   lineage (a lineage that lost more than half its cells means the thresholds do not fit the tissue;
   adaptive per-sample QC is the default, and mitochondrial-rich tissues such as kidney tubule,
   cardiomyocytes and hepatocytes need a higher cap), enriched clusters (tiny or infiltrate-only
   means the composition route favours consistent infiltrates over case cells that split by
   replicate; look at `split_state_candidates`), doublet and dissociation-stress suspects (not
   disease states however enriched), the recorded signatures (do they look like the biology?), the
   library check, and the eval line if a key exists. Iterate on **config only**: thresholds, caps,
   pre-registered exclusions. Never change code or prompts to fit one dataset.
6. **Run the agent.** `python scripts/sctarget.py run --config config/<name>.yaml --model
   claude-code:<model id>` uses the Claude Code login on a named model (for example
   `claude-code:claude-fable-5-1`). Bare `claude-code` runs on the CLI's own default model, which
   depends on the user's Claude Code settings; whichever model answers is recorded in the run log
   and the report header (`usage.model`, `facts.model_ids`), so a report can always say which model
   produced it. `--model <claude model id>` alone uses the API with `ANTHROPIC_API_KEY` from the
   environment; on datasets of 100k cells and more the nominator's context reaches 150-220K tokens,
   so that model needs a window above 200K. Expect
   0.6-2.3M input tokens and 20-90 minutes for 30k-200k cells; run one dataset at a time (the
   druggability services rate-limit). Tell the researcher what leaves the machine: cluster-level
   tables, marker lists and gene symbols go to the model; gene symbols go to DGIdb, Open Targets and
   NCBI; no per-cell data and no real sample identifiers.
7. **Read the result with the researcher.** Start with the report's executive summary (design
   statement, key numbers, top candidates), then the critic's flags (section 4), then the ranking
   (section 3): the evidence table, then the per-library panels (one point per library
   by lineage and state, coloured by condition; a real candidate shows most case libraries above the
   reference ones in its own lineage; a single high library or a shift in every lineage is a
   replicate effect or ambient RNA). For each top candidate check, using `references/flags.md`:
   - `target_class` (`state_specific` or `disease_wide`: two different claims) and
     `lineage_axis_source`;
   - `proliferation_marker`: in a tumour it means "the cells divide"; in a non-malignant design a
     genuinely expanding population (activated T cells, regenerating epithelium) will also carry it,
     so read it as "describe the expansion, do not nominate the cycle gene", not as "artefact";
   - attribution: `reference_driven` or `ambient_like` means the signal is not from these cells;
     `shared_with_reference` means a lineage or pan-condition gene. **When attribution sets aside a
     candidate that the per-sample abundance and the compartment contrast both support, do not
     discard it: identify the reference cells** (which lineage, which libraries, how many, what
     markers) and report the candidate as *conditional on that population*. A rare benign
     population expressing the gene (`rare_reference_counterpart`) is an on-target consideration,
     not a disqualification. When the tools answer "the reference holds this lineage in bulk" with
     a **derived** lineage, that answer is only as fine as the lineage: in the omentum design the
     derived epithelial lineage absorbed mesothelium, so 151 benign epithelial cells read as a
     representative counterpart and FOLR1 was held back in five agent runs (prompts v8 and v9, the
     latter stating the conditional rule). The settling move, which the critic named and did not
     take, is a **reference-only subcluster of the state's cluster** (`subcluster` restricted to
     reference cells, then markers): mesothelium (CALB2, UPK3B, MSLN) and benign Müllerian-type
     epithelium (PAX8, KRT19, WFDC2, FOLR1) separate at once, and the candidate is then reported as
     conditional on the named population, with the histology or exclusion arm that would settle it;
   - **therapeutic direction**: druggability counts compounds and does not know whether inhibiting
     the gene helps. For every top candidate, state whether the listed compounds inhibit or
     activate and whether that is the therapeutic direction for this design (a PHD inhibitor
     stabilises HIF; an agonist antibody on a checkpoint is not an inhibitor); ask the critic
     explicitly if the flag is absent;
   - ortholog status for non-human data (`same_symbol_human_gene` and `assumed_same_symbol` are
     flagged assumptions, not mappings); "novel" genes above a known target are hypotheses to look
     up, not errors.
   The nominator over-nominates (precision 0.1-0.3): treat the top 3-5 per state as the list.
   Then read the **agent notes** (report section 5) for what each role says about its own inputs
   and limits, not only its conclusions. In round 0.4c one such sentence - the critic reporting
   that the candidate list it had been given was cut short - exposed a briefing truncated at a flat
   character cap in every agent run since the method was frozen (repository protocol 23m); no
   evaluator metric could have shown it. A role saying "I could not run X" or "the list I received
   ended at Y" is a method finding to record, whatever the scores say.
8. **Optional arms and reporting.** Annotation-free run (`analysis.lineages.ignore_annotation:
   true`) to compare with an inherited annotation; three repeat runs for a stability estimate
   (`scripts/run_eval.py --stability`, which currently needs a key; without one, compare the top-k
   of `tables/*ranked_targets*.csv` across runs by hand); `python scripts/sctarget.py report` to
   regenerate every number from the run directories.

**Stop conditions.** Stop and say so, rather than proceed, when: `check-env` fails; a condition has
fewer than two independent replicates; a covariate is perfectly confounded with the condition; the
compartment the question concerns is absent in both conditions; a run ends with tool errors in
`run.jsonl` ("no disease-specific state found" after errors is a broken chain, not a negative);
or a guard would have to be changed a second time to reach the key (record the limitation instead).

## Output format

Report in this order:

1. **Design statement** - design type and objective, what the reference represents, independent
   replicates per condition (donors, not libraries), covariates checked for confounding, libraries
   flagged by the label checks and the pre-registered decision about them, what parts of the method
   have been validated on this design type and what has not.
2. **Data summary** - cells and libraries per condition after QC, species, inherited or derived
   lineages, per-lineage QC losses.
3. **Condition-associated states** - each recorded state with lineage, enrichment across replicates,
   signature (sibling, compartment and whole-lineage contrasts), and any doublet, stress,
   replicate-split or composition-shift caveat.
4. **Candidate ranking** - per state, the top 3-5 with score, `target_class`, attribution (and, for
   a set-aside candidate, the reference population identified), druggability with **direction**,
   ortholog status, and the critic's flag. In `marker` mode: effect, specificity, consistency, and
   druggability as information only.
5. **Critic flags** - verbatim, with what each means for the ranking.
6. **Evaluation** - only with a frozen key: state recovery (lenient and strict), target rank,
   decoys above the target, decoys nominated unflagged, drug recall, method versions.
   Without a key: the positive-control checks and the stability estimate.
7. **Recommended next checks** - the orthogonal experiment each top candidate needs, and the
   config-level changes worth trying if calibration showed a problem.
8. **Human-review points** - listed explicitly.

Point to the run directory: `report.md`, `report.html`, `state.json`, `tables/`, `figures/`, and
`run.jsonl` (every tool call with the exact payload the model saw).

## Human review required

A knowledgeable researcher must decide or confirm:

- the design type, the objective, what the reference represents, and which column is the
  independent biological unit;
- that the condition labels were verified and how each doubtful library is handled (pre-registered);
- whether a confounded or imbalanced design should be run at all, and how the limitation is stated;
- whether the answer key is fair for this design type and was frozen before the first run;
- whether QC thresholds fit the tissue, and whether a lineage lost to QC matters for the question;
- whether a recorded state is biology or an artefact (doublets, stress, one replicate, a
  composition shift that is the whole story);
- for each top candidate: therapeutic direction, relevance of the compounds, and - when attribution
  set it aside - whether the reference population identified changes the claim;
- whether a "novel" candidate above the known target is an omission in the key or a false positive;
- every proposed change to code, prompts or scoring, which must become a labelled version with a
  stated failure, never a dataset-specific fix.

## Limitations

Labelled **G** (general to the method) or **O** (observed in the tumour datasets; untested elsewhere).

- **G** Expression-based: targets regulated post-translationally are invisible.
- **G** Replicate-aware tests need independent replicates; libraries nested in donors, or a
  covariate aligned with the condition, cannot be repaired downstream.
- **G** Druggability is a compound count without direction; direction is a human check.
- **G** Drug knowledge is human-centric; other species rest on ortholog mapping.
- **G** Attribution is the most conservative step and its verdicts depend on how clean the derived
  lineages are; a rare reference population expressing the candidate can set a real target aside.
- **O** The composition route favours infiltrates consistent across patients over malignant cells
  that split by patient; merging mitigates this. Expected to matter less when the condition shifts
  programmes within existing lineages (non-malignant disease), but this is untested.
- **O** The proliferation flag was motivated by tumours; its behaviour on genuinely expanding
  populations in non-malignant designs is untested.
- **O** Nomination precision 0.1-0.3; the nominator rarely subclusters unasked.
- **O** The candidate pool is finite: the most significant genes by statistic plus the largest
  significant effects (protocol 23f, after FOLR1 sat 39th by statistic and the cardiomyocyte
  atrophy genes sat outside the pool in the doxorubicin floor); a real gene can still fall outside it.
- **G** Programme-basis states (a lineage that exists in both conditions and runs a different
  programme) are judged by the within-cluster, replicate-aware contrast: attribution, consistency and
  the relative programme criterion follow that basis (protocol 23/23b/23d). Before those rules the
  expression route reported false negatives on exactly the designs it exists for.
- **G** When every library is one condition, the integration (Harmony on the library) cannot tell
  library from condition and may absorb part of the disease signal: in the fibrolamellar carcinoma
  snRNA-seq, tumour nuclei co-clustered with normal hepatocytes and the tumour state was reached
  only by the expression route (2,009 genes up within the cluster; PRKACA 2nd behind CA12). The
  route therefore runs in every analysis, also in tumour designs; a composition-only reading of
  such a dataset would have been a false negative.
- Validated only on four tumour datasets (canine ameloblastoma vs another oral carcinoma; ccRCC,
  HCC and ovarian omental metastasis vs matched tissue). `marker` mode has run only on synthetic
  data. No non-malignant disease, perturbation, blood, brain or paired-only design has been run.
- The subscription backend has no seed or temperature control; stability is behavioural.
- Needs a Claude Code login or an Anthropic API key for the agent run; calibration and the demo run
  offline.

## Final checks

Before completing the task:

- do not invent cells, genes, scores, drugs or run outcomes; every number traces to `run.jsonl` or
  a table in the run directory, and every helper-script result to its printed output;
- distinguish what the tools observed from what the agent or you interpret; name the design type
  and say which claims are untested on it;
- state the uncertainty: replicates, confounding, thin reference pools, replicate splits, unflagged
  decoys, lenient-only state recovery, set-aside candidates and why;
- confirm that no per-cell data, real sample identifiers, raw files, API keys or restricted data were
  committed, pasted or sent anywhere, and that the privacy guard was not disabled;
- confirm that the answer key, if any, was not edited after a run, and that every exclusion was
  pre-registered with a reason in the config;
- list the human-review points with the researcher's answers, or mark them open.
