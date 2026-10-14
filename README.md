# sc-target-agent

Agentic target discovery from single-cell / single-nucleus RNA-seq.
Deterministic analysis tools, orchestrated by an auditable LLM tool loop, scored against
ground truth you already know.

Built for the CoPAIR fellowship project: *can an agent independently rediscover the
disease-specific epithelial state and the ERBB4 / neratinib vulnerability in canine
acanthomatous ameloblastoma (CAA) snRNA-seq, and does it make fewer false leads than a
naive script?*

```
   .h5ad  ──►  planner ──► analyst ──► target nominator ──► critic ──► report.md/.html
               (plan)      (states)    (druggable genes)    (flags)    + run.jsonl audit log
                  │            │              │                │
                  └────────────┴──────────────┴────────────────┘
                         16 deterministic tools (scanpy, PyDESeq2, DGIdb, Open Targets …)
                         only aggregated results ever reach the model
```

## Why this design

* **Tools, not free-form code.** The model chooses *which* analysis to run and with what
  parameters; the numbers come from tested Python functions. Every call is logged with its
  arguments and the exact payload the model saw (`run.jsonl`), so any sentence in the report
  can be traced to a tool result.
* **Replicate-aware by default.** Disease states are found by per-*sample* differential
  abundance; state genes by pseudobulk DESeq2 against *sibling clusters of the same cell
  type* (so you get state markers, not cell-type markers). A `gene_specificity` tool and a
  critic role exist specifically to catch pan-tumour and cell-type decoys.
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
# 1. environment
conda env create -f environment.yml && conda activate sctarget     # or: pip install -e ".[full,dev]"

# 2. synthetic CAA-like dataset with a planted answer (ERBB4/PEG3 state, decoys EGFR/KRAS/...)
python scripts/make_synthetic.py --out data/synthetic

# 3. run the agent pipeline with the deterministic stand-in model, and the naive baseline
export SCTARGET_OFFLINE=1
python scripts/run_pipeline.py --config config/synthetic.yaml --model scripted --run-dir runs/agent
python scripts/run_baseline.py --config config/synthetic.yaml --run-dir runs/baseline

# 4. score both against the ground truth
python scripts/run_eval.py --ground-truth data/ground_truth/synthetic.yaml runs/agent runs/baseline

# 5. tests
pytest -q
```

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
                       gene_specificity, enrichment, plot_genes
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
tests/                 30 tests: tools, privacy, registry, loop error handling, end-to-end + eval
```

Each run directory contains `report.md/.html`, `run.jsonl`, `state.json`, `config_used.json`,
`tables/*.csv` (full results), `figures/*.png`, `cluster_assignments.csv` (local only).

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
