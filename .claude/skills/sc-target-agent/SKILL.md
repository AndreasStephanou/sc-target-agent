---
name: sc-target-agent
description: Run auditable, agent-driven drug-target discovery on a single-cell or single-nucleus RNA-seq dataset with two conditions (disease vs reference, treated vs control, knockout vs wild-type). Use when someone has an .h5ad (or 10x matrices) with per-sample condition labels and wants disease-specific cell states and druggable targets with evidence, or wants to validate the method against a known answer.
---

# sc-target-agent

You are helping a researcher run the sc-target-agent pipeline on their own dataset. The pipeline is
a four-role agent (planner, analyst, nominator, critic) over ~20 deterministic single-cell tools,
with an audit log and an evaluation harness. Your job is to get the data into shape, run the
calibration, interpret what it says, run the agent, and read the critic's flags with the researcher.
You do not re-implement the analysis; you drive the repository's scripts.

## Before anything: the four facts you need

Ask for (or find in the object) these; nothing works without them.

1. **Where the raw counts are.** A `.h5ad` with integer counts in `X` or a layer, or per-sample
   10x matrices, or GEO's "one dense count matrix + metadata table" deposit (then
   `python scripts/build_from_matrix.py --dry-run ...` to see the columns, and build; `scripts/build_ccrcc.py`
   is the worked example for 10x directories).
2. **The sample column.** One value per library/biological replicate. At least 3 samples per
   condition, or the per-sample tests are weak and you must say so.
3. **The condition column** and which value is the case (disease/treated) and which are reference.
4. **Species.** Human needs nothing; mouse/dog/other map to human symbols for drug lookups.

A cell-type annotation is optional (lineages are derived from the data otherwise). If one exists,
note that it came from the original analysis and may inherit its choices.

## Step 0: check the environment

```bash
python scripts/sctarget.py check-env
```

This runs the sparse, large-data tool chain (approximate neighbours, Leiden, Scrublet, Harmony,
DESeq2) on a 9k-cell synthetic object. A fresh install once passed every test and the demo yet
could not cluster any real dataset (a library-version conflict); this catches that in a minute
instead of after a 10-minute run. If it fails, rebuild from the pinned `environment.yml`.

## Step 1: inspect and scaffold

```bash
python scripts/sctarget.py inspect data/<file>.h5ad --name <short-name>
```

Read its guesses aloud with the researcher and fix `config/<name>.yaml`. Check the warnings:
fewer than 3 samples in a condition, no counts found, more than two condition levels.

**Do not trust public metadata.** GEO/ArrayExpress sample labels have been wrong before (in the
ccRCC test set, tumour/benign were swapped for 4 of 13 libraries). Verify condition labels
against a marker the researcher knows (`CA9` for ccRCC, a known disease gene for theirs) per
library before any run.

## Step 2: decide whether there is an answer key

If the researcher already knows a target for this disease (published, functionally validated),
write `data/ground_truth/<name>.yaml` NOW, before any run, and commit it:
expected state markers (6-20 genes), expected target(s), acceptable targets, expected drugs,
decoys (druggable but obviously wrong: proliferation genes, lineage markers, pan-condition genes).
Then `git tag <name>-key-frozen`. Never edit it after the first run; if it was wrong, say so in
the write-up (see the EPAS1 lesson below).

**An expected target must be transcriptionally distinctive in the disease state.** The pipeline
nominates from expression. A target activated post-translationally (HIF-2alpha / EPAS1 in ccRCC)
is invisible to it and would test the key, not the method; put such a gene in `acceptable_targets`
and choose its transcriptional read-out as the expected target instead.

No known answer? Skip the key. The pipeline still runs and produces hypotheses with critic flags;
the evaluator is simply not applicable.

## Step 3: calibrate (free, no tokens)

```bash
python scripts/sctarget.py calibrate --config config/<name>.yaml [--ground-truth data/ground_truth/<name>.yaml]
```

This runs the fixed scripted procedure through the same tools. Read `runs/<name>_calibrate/report.md`
with the researcher and look for these, in order:

* **QC removed a lineage.** Section 1: if a cell type lost > 50% of its cells, the thresholds do
  not fit the tissue (kidney tubule and cardiomyocytes run 20-45% mitochondrial reads normally).
  The default is adaptive per sample; raise `max_mt_pct` / `mt_floor_pct` in the config if needed.
* **Enriched clusters are tiny or immune.** The composition route favours infiltrates that are
  consistent across patients over malignant cells that split by patient. Look at
  `split_state_candidates`; the agent has `merge_clusters` for this.
* **Doublet and stress suspects.** Clusters flagged `doublet_suspect` or `stress_artifact_suspect`
  are not disease states however enriched. Dissociated data always has some.
* **The recorded signatures.** Do they look like the biology the researcher expects? If the state
  genes are all "sub-state vs sub-state" differences, check the compartment contrast lists.
* **The eval line** (if a key exists): state found? decoys above the target? If the scripted floor
  already recovers the answer, the agent run measures what reasoning adds; if it misses, the
  calibration tells you where (and a fix should be a generic rule, never a gene).

Iterate on config only (thresholds, caps). Do not change code or prompts to fit one dataset.

## Step 4: run the agent

```bash
python scripts/sctarget.py run --config config/<name>.yaml            # Claude Code login, no API key
python scripts/sctarget.py run --config config/<name>.yaml --model claude-sonnet-4-5   # API key
```

Expect 0.6-1.7M input tokens and 20-90 minutes for 30k-200k cells on the subscription backend.
Run one at a time (the druggability databases rate-limit). The run directory has `report.md`,
`report.html`, `state.json`, `tables/`, `figures/`, and `run.jsonl` (every tool call with the
exact payload the model saw). Nothing per-cell ever reaches the model.

## Step 5: read it with the researcher

Start with **section 5, the critic flags**, then **section 4, the ranking**. For each candidate:

* `target_class`: `state_specific` (defines this state among its lineage) or `disease_wide`
  (up in disease cells across the lineage). Both can be targets; they are different claims.
* `attribution`: `reference_driven` or `ambient_like` means the signal is not from these cells;
  `shared_with_reference` means a lineage or pan-condition gene; `reference_pool_in_cluster`
  says whether the in-cluster reference comparison was even usable.
* **Direction.** Druggability counts compounds; it does not know whether inhibiting the target
  helps. The critic has caught this (EGLN3 inhibitors stabilise HIF). Ask it explicitly.
* **Orthologs.** For non-human data, `assumed_or_unverified` genes rest on a same-name assumption.
* **Novel candidates.** Genes above the known target that are not decoys are hypotheses, not
  errors; the eval reports them separately. Check the literature before dismissing them (CD70
  and FOLH1 in ccRCC were real clinical-stage targets the key had omitted).

The nominator tends to over-nominate (precision 0.1-0.3) and rarely subclusters on its own;
treat the top 3-5 per state as the real list.

## Step 6: optional arms and reporting

* **Annotation-free:** set `analysis.lineages.ignore_annotation: true` and run again; lineages are
  derived from the data. Compare the two.
* **Repeat runs:** 3 runs give a stability estimate (`run_eval.py --stability`).
* **Results table:** `python scripts/sctarget.py report --ground-truth data/ground_truth/<name>.yaml runs/<name>_*`
  regenerates every number from the run directories.

## Things this skill must not do

* Edit `data/ground_truth/*.yaml` after a run, or write a key after seeing results.
* Change prompts, scoring weights or thresholds in code to make one dataset pass. Config-level
  caps are fine; anything else becomes a new labelled version with a reason in
  `docs/EVALUATION_PROTOCOL.md`.
* Present a nominated gene as validated. Every candidate is a hypothesis for IHC / functional
  follow-up; the report says so and so should you.
* Send per-cell data or real sample identifiers anywhere; the privacy guard enforces this, do
  not disable it.

## Where things are

`docs/EVALUATION_PROTOCOL.md` is the method history (every rule and the failure that motivated
it). `docs/REPORT.md` is the write-up of the two validation datasets (canine ameloblastoma,
human ccRCC). `tests/` (`pytest -q`) must pass after any change.
