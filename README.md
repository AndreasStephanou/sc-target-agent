# sc-target-agent

Agentic target discovery from single-cell / single-nucleus RNA-seq.
Deterministic analysis tools, orchestrated by an auditable LLM tool loop, scored against
ground truth you already know.

Built as a CoPAIR fellowship project at Cornell and validated on two datasets: it rediscovered
the ERBB4 / neratinib vulnerability of canine acanthomatous ameloblastoma from our own snRNA-seq,
and the HIF programme (CA9 first) of human clear cell renal cell carcinoma from public scRNA-seq,
from counts and condition labels alone, with no decoy ever nominated unflagged.
`docs/REPORT.md` is the write-up; `docs/EVALUATION_PROTOCOL.md` is the method history.

**Works on any species and any two conditions** (disease vs reference, treated vs control,
knockout vs wild-type). Human data need no ortholog step; a cell-type annotation is optional
(lineages are derived from the data). Runs on a Claude subscription login (no API key) or the API.

```
   .h5ad  ──►  planner ──► analyst ──► target nominator ──► critic ──► report.md/.html
               (plan)      (states)    (druggable genes)    (flags)    + run.jsonl audit log
                  │            │              │                │
                  └────────────┴──────────────┴────────────────┘
                         18 deterministic tools (scanpy, PyDESeq2, DGIdb, Open Targets …)
                         only aggregated results ever reach the model
```

## Why this design

* **Tools, not free-form code.** The model chooses *which* analysis to run and with what
  parameters; the numbers come from tested Python functions. Every call is logged with its
  arguments and the exact payload the model saw (`run.jsonl`), so any sentence in the report
  can be traced to a tool result.
* **Replicate-aware by default.** Disease states are found by per-*sample* differential
  abundance; state genes by pseudobulk DESeq2 against *sibling clusters of the same cell
  type* (so you get state markers, not cell-type markers), **blocked by sample**
  (`~ sample + label`: a sample's state cells vs the same sample's sibling cells), so the
  shared background of disease samples cannot inflate the reference dispersion. Specificity is
  judged against clusters that are mostly reference cells, not against other disease populations. A `gene_specificity` tool and a
  critic role exist specifically to catch pan-tumour and cell-type decoys.
* **Where does the signal sit?** `condition_attribution` asks, for every candidate gene, whether
  its expression in a cluster comes from the disease-condition cells or from reference-condition
  cells that landed in the same cluster (`reference_driven`: a false positive of cluster-vs-cluster
  DE that mixed-condition clusters produce), whether it is shared with reference cells of the same
  cell type (`shared_with_reference`: cell-type marker or pan-condition gene), or whether it is
  disease-wide rather than state-specific. The scorer penalises the first two and the critic
  raises them as flags.
* **Ambient-RNA check.** When one disease population dominates the disease samples, its
  transcripts sit as background in every nucleus from those samples, so any "disease vs
  reference" comparison inside a shared cluster finds that population's program in immune,
  stromal and other epithelial cells alike. Every gene claim gets an *ambient ratio*: its level in
  the cells of interest vs the least-related lineage in the same disease samples. A gene is
  `ambient_like` when that ratio is < 3 **and** the control lineage is itself raised in disease vs
  reference samples (background must show up there too; a gene flat in the control lineage cannot
  be ambient). The test is symmetric - a TME target in macrophages is checked against epithelium -
  and never privileges the tumour compartment. `disease_programs` uses the same ratio to strip
  programs down to `clean_program_genes` and marks mostly-ambient programs `ambient_suspect`.
  Exercised by `make_synthetic --ambient 0.15`. Upstream decontamination (CellBender / SoupX /
  decontX) remains the principled fix; this is the safeguard for uncorrected data.
* **Two routes to a disease state.** *Composition*: a cluster enriched in disease samples
  (`differential_abundance`, with subclustering / re-clustering fallbacks). *Expression*:
  a cell type present in every condition that runs a disease-specific program
  (`disease_programs`: pseudobulk disease-vs-reference *within* each cluster). States carry a
  `basis` (composition | program) and `score_targets` switches its specificity component
  accordingly (cluster-vs-cluster for composition, within-cluster disease fold change for
  program). Both routes are exercised by the synthetic datasets (`make_synthetic --mode state|program`).
* **Privacy layer.** Raw counts and per-cell values never leave the machine. Tool outputs are
  truncated, sample IDs pseudonymized, and any per-cell-sized array is rejected before it can
  be sent (`sctarget/privacy.py`; tested in `tests/test_privacy.py` and the e2e test).
* **Evaluable.** A `ground_truth.yaml` per dataset plus `scripts/run_eval.py` give state
  recovery, target rank, drug recall, false leads, decoy handling, and run-to-run stability.
  A non-agentic baseline (`scripts/run_baseline.py`) uses the same tools with a naive recipe,
  so the agent's added value is measured, not assumed.
* **Runs without an API key.** `--model scripted` is a deterministic stand-in that exercises
  the whole machinery (used in tests/CI). Swap in any Claude model ID for real runs.

## Quick start (5 minutes, offline)

```bash
conda env create -f environment.yml && conda activate sctarget     # or: pip install -e ".[full,dev]"
python scripts/sctarget.py demo        # environment check + offline synthetic demo        # generates a synthetic dataset with a planted answer, runs the scripted
                                       # pipeline and the naive baseline offline, scores both (no key, no network)
SCTARGET_OFFLINE=1 pytest -q           # 66 tests, ~4 minutes
```

The individual steps behind `demo` are `scripts/make_synthetic.py`, `scripts/run_pipeline.py --model scripted`,
`scripts/run_baseline.py` and `scripts/run_eval.py`; see `scripts/sctarget.py --help`.

Open `runs/agent/report.html`. On the synthetic data the agent ranks ERBB4 first with
neratinib/afatinib/lapatinib, nominates no decoys, and passes; the baseline also finds the
state but nominates KRAS and TOP2A (pan-tumour decoys) unflagged and fails.

## Running on the real CAA data

1. Export your processed object as `.h5ad` with: raw UMI counts in a layer (e.g. `counts`),
   one obs column per animal/library (`sample`), one for condition (`CAA` / `COSCC` /
   `GINGIVA`), and optionally your cell-type annotation. Existing `leiden`/`X_umap` are reused.
2. `cp config/caa_template.yaml config/caa.yaml` and edit the column names / species.
3. Download MSigDB Hallmark (`h.all.*.symbols.gmt`) and point `analysis.enrichment.gene_sets`
   at it (the bundled `mini_hallmark.gmt` is a test stand-in, not real gene sets).
4. **Before running:** fill `data/ground_truth/caa_template.yaml` → `caa.yaml` from the
   published results, commit it, and do not edit it afterwards. That is what makes the blind
   run honest (see `docs/EVALUATION_PROTOCOL.md`).
5. Run with a real model:
   ```bash
   export ANTHROPIC_API_KEY=sk-ant-...
   unset SCTARGET_OFFLINE                   # allow live DGIdb / Open Targets / Ensembl / PubMed
   python scripts/run_pipeline.py --config config/caa.yaml --model claude-sonnet-4-5
   python scripts/run_baseline.py --config config/caa.yaml
   python scripts/run_eval.py --ground-truth data/ground_truth/caa.yaml runs/<agent_run> runs/<baseline_run>
   ```
   Repeat the agent run 3-5× (different `--seed`, or just rerun) and add `--stability` to
   `run_eval.py` to measure consistency.

Data handling: the model API receives gene symbols, fold changes, p-values, proportions,
enrichment terms and cluster sizes – never counts or per-cell values. Check the payloads
yourself: every `tool_call` event in `run.jsonl` has a `result_sent` field that is exactly
what the model saw. Confirm this fits Cornell's data policy for your dataset before the
first live run; `SCTARGET_OFFLINE=1` keeps even gene-symbol lookups local.

## Your own dataset in four commands

```bash
python scripts/sctarget.py inspect   data/my.h5ad --name mydata          # guesses sample/condition/species, writes config/mydata.yaml + a key template
#   -> edit config/mydata.yaml; verify condition labels against a marker you trust (public metadata is wrong more often than you think)
#   -> if you KNOW a target for this disease, fill data/ground_truth/mydata.yaml now and commit it (frozen before any run)
python scripts/sctarget.py calibrate --config config/mydata.yaml         # scripted procedure, no tokens: QC per lineage, states, doublet/stress flags
python scripts/sctarget.py run       --config config/mydata.yaml         # the agent (Claude Code login by default; --model <api model> with a key)
python scripts/sctarget.py report    --ground-truth data/ground_truth/mydata.yaml runs/mydata_*   # results table from the run directories
```

Read `runs/<run>/report.md` section 5 (critic flags) before section 4 (ranking). Every number in
it traces to a line in `runs/<run>/run.jsonl`. Nothing per-cell ever reaches the model.

Cost guide (one node, no GPU): the scripted calibration takes ~10 minutes for 30k cells and
~65 minutes for 200k cells (QC, Harmony, Leiden, pseudobulk DESeq2 per state). An agent run adds
20-90 minutes and 0.6-1.7M input tokens on the subscription backend. Run agent runs one at a time:
the druggability databases rate-limit.

**Using Claude Code?** Copy `skills/sc-target-agent/` into your project's `.claude/skills/` and
Claude will walk through the steps above with the lessons from both validation datasets built in
(what to check in QC, how to read a thin reference pool, why an expression method cannot nominate
a protein-level target, and so on).

## What is in the box

```
sctarget/
  config.py            config loading (config/default.yaml + per-dataset yaml + env)
  session.py           AnalysisSession: holds AnnData + results locally; tools operate on it
  privacy.py           sanitize(): truncation, per-cell rejection, pseudonymization
  runlog.py            append-only JSONL audit log
  tools/
    registry.py        @tool decorator -> JSON schema, role permissions, arg validation
    analysis.py        load_dataset, run_qc, integrate_and_cluster, subcluster, cluster_composition,
                       differential_abundance, disease_programs, find_markers, pseudobulk_de,
                       gene_specificity, condition_attribution, enrichment, plot_genes
    knowledge.py       map_orthologs (Ensembl), druggability_lookup (DGIdb + Open Targets),
                       literature_count (PubMed) – disk-cached, offline fixture fallback
    targets.py         score_targets (deterministic, weighted, explainable), record_finding
  agents/
    model.py           AnthropicModel (SDK tool use) and ScriptedModel (deterministic test double)
    loop.py            the tool-calling loop: validate -> run -> sanitize -> log -> feed back
    prompts.py         system prompts for planner / analyst / target / critic + stage prompts
    orchestrator.py    runs the roles in order, persists state.json, builds the report
  report/builder.py    report.md + report.html from facts + run log
  eval/
    score.py           evaluate_run(): state recovery, target rank, drug recall, false leads, decoys
    stability.py       top-k Jaccard and nomination overlap across runs
scripts/               make_synthetic, run_pipeline, run_baseline, run_eval
config/                default.yaml, synthetic.yaml, caa_template.yaml
data/ground_truth/     synthetic.yaml, caa_template.yaml
data/fixtures/         mini_hallmark.gmt (test stand-in), knowledge_fixtures.json (offline DB stand-in)
tests/                 35 tests: tools, privacy, registry, loop error handling, end-to-end + eval
```

Each run directory contains `report.md/.html`, `run.jsonl`, `state.json`, `config_used.json`,
`tables/*.csv` (full results), `figures/*.png`, `cluster_assignments.csv` (local only).

`--resume-from runs/<earlier>` reuses an earlier run's analysis (states, DE tables, clustering)
and runs only `--stages target critic` by default - re-run a nominator or critic (or try a
different model for them) without repeating the 15-minute analysis.

## Scoring

`score_targets` ranks candidates by a weighted sum (weights in `targets.weights`):

| component | source | default weight |
|---|---|---|
| effect size (log2FC, significance-weighted) | `pseudobulk_de`, else `find_markers` | 0.25 |
| specificity (log2 target-cluster mean / max other cluster) | `gene_specificity` | 0.25 |
| druggability (DGIdb class + Open Targets tractability) | `druggability_lookup` | 0.25 |
| known compounds (approved inhibitors) | `druggability_lookup` | 0.15 |
| consistency across disease samples | `gene_specificity` | 0.10 |

Flags (`low_specificity`, `shared_across_conditions`, `inconsistent_across_samples`,
`not_druggable`, `no_known_compounds`, `missing_evidence:*`) travel with each gene into the
report and the critic's review, and the first three also penalise the score (x0.6, x0.6, x0.8),
so a druggable pan-condition gene cannot outrank a clean state-specific one. Change the weights in the config, not in code, so the
choice is recorded in `config_used.json`.

## Known limitations / next steps

* The knowledge fixtures were hand-written for tests and are **illustrative only**. In this
  development sandbox the live APIs were unreachable, so the live clients are implemented
  and cached but **unverified** – run `druggability_lookup` on a few genes with
  `SCTARGET_OFFLINE` unset and check the `source` column says `live` the first time and
  `cache` the second. API schemas (DGIdb v5 GraphQL, Open Targets v4) change occasionally.
* Differential abundance uses per-sample proportions + Mann-Whitney; with 3 vs 6 samples the
  p-values are coarse by construction. Swap in scCODA or Milo via a new tool if you want a
  model-based test – the registry makes that a ~40-line addition.
* The ScriptedModel is a test double. Its adaptive steps (subcluster the largest cell type,
  then re-cluster globally) are heuristics; the real model gets the same guidance in
  `prompts.py` and can do better or worse – that difference is part of what you measure.
* Not yet implemented: LINCS/CMap signature reversal, DepMap essentiality, a sandboxed
  free-form code tool. All three fit the same `@tool` pattern.
* Baseline runs call tools directly (not through the loop), so their `run.jsonl` has no
  tool-call events; `state.json` is complete.

## Citation / provenance

Target biology this project is benchmarked against: Stephanou et al., *Neuroepithelial
reprogramming and ERBB vulnerability in canine acanthomatous ameloblastoma*, Mol Ther Oncol
2026. Knowledge sources: DGIdb, Open Targets Platform, Ensembl, PubMed/NCBI E-utilities,
MSigDB. Analysis: scanpy, anndata, PyDESeq2, harmonypy, leidenalg.
