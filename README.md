# sc-target-agent

Auditable, agent-driven drug-target discovery from single-cell or single-nucleus RNA-seq with two
conditions (disease vs reference, treated vs control, knockout vs wild-type).

A four-role agent (planner, analyst, nominator, critic) works only through ~20 deterministic
single-cell tools; every tool call and the exact payload the model saw are logged; no per-cell
data or real sample identifier ever leaves your machine. If a target is already known for the
disease, a frozen answer key lets you score the method blind before trusting it on an open
question. Built as a CoPAIR fellowship project at Cornell (Communities of Practice for AI-enabled
Research); the method history is in `docs/EVALUATION_PROTOCOL.md` and the validation write-up in
`docs/REPORT.md`.

```
   .h5ad  ──►  planner ──► analyst ──► nominator ──► critic ──► report.md / report.html
               (plan)      (states)    (targets)     (flags)    + run.jsonl (every tool call)
                  │            │            │            │
                  └────────────┴────────────┴────────────┘
                     deterministic tools: QC, doublets, Harmony + Leiden, derived lineages,
                     per-sample differential abundance, pseudobulk DESeq2 (three contrasts),
                     specificity, condition attribution, ambient check, DGIdb / Open Targets /
                     PubMed, weighted score  —  only aggregated tables reach the model
```

---

## Contents

1. [What you need](#1-what-you-need)
2. [Install and check](#2-install-and-check)
3. [Walkthrough: a public dataset from GEO to a ranked target list](#3-walkthrough-a-public-dataset-from-geo-to-a-ranked-target-list)
4. [Your own dataset](#4-your-own-dataset)
5. [Choosing the model backend](#5-choosing-the-model-backend)
6. [Reading the report](#6-reading-the-report)
7. [The rules that keep a run honest](#7-the-rules-that-keep-a-run-honest)
8. [Validation so far](#8-validation-so-far)
9. [What is in the box](#9-what-is-in-the-box)
10. [How a candidate is scored](#10-how-a-candidate-is-scored)
11. [Troubleshooting](#11-troubleshooting)
12. [Limitations](#12-limitations)
13. [Data handling](#13-data-handling)
14. [Using it from Claude Code](#14-using-it-from-claude-code)
15. [Provenance](#15-provenance)

---

## 1. What you need

- **A Linux or macOS machine with conda** (Miniconda or Mambaforge). No GPU. 32 GB RAM is
  comfortable for 150k cells; 16 GB works to about 60k.
- **Data**: raw UMI counts for one dataset, as an `.h5ad`, per-sample 10x matrices, or a GEO
  "one matrix + metadata" deposit; a column that identifies the **library / biological replicate**
  (at least 3 per condition, or the per-sample tests are weak); a column that identifies the
  **condition**, with one value as the case and the rest as reference. A cell-type annotation is
  optional (lineages are derived from the data otherwise).
- **For the agent run only**: a Claude Code login (the default; no API key) or an Anthropic API
  key. The calibration run, the offline demo and the tests need neither.
- **Network** for gene-level lookups (DGIdb, Open Targets, NCBI E-utilities) during real runs;
  `SCTARGET_OFFLINE=1` keeps everything local with bundled illustrative fixtures.

## 2. Install and check

```bash
git clone https://github.com/AndreasStephanou/sc-target-agent
cd sc-target-agent
conda env create -f environment.yml
conda activate sctarget

python scripts/sctarget.py check-env     # 1 minute: exercises the sparse / large-data tool chain
python scripts/sctarget.py demo          # 5 minutes, offline: synthetic data with a planted answer
SCTARGET_OFFLINE=1 pytest -q             # 87 tests, ~7 minutes
```

`check-env` matters: a fresh install once passed every test and the demo yet could not cluster any
real dataset (a dependency drifted). It builds a 9k-cell sparse object and runs approximate
neighbours, Leiden, Scrublet, Harmony and DESeq2 through the installed packages. If it fails, the
message says which step and how to rebuild from the pinned `environment.yml`.

`demo` runs the whole pipeline with a deterministic scripted model and a naive baseline on a
synthetic dataset, then scores both against the planted answer. Expected: the scripted run finds the
planted state, ranks the planted target first with its drugs, nominates no decoys and PASSes; the
baseline finds the state but lets decoys through. Open `runs/demo_scripted/report.html`.

## 3. Walkthrough: a public dataset from GEO to a ranked target list

This is the sequence the skill follows, on a real public deposit. The example is GSE286057
(ovarian cancer omental metastases vs omentum, 2025, per-sample 10x files); the same steps apply to
any GEO series with per-sample `matrix.mtx` / `features.tsv` / `barcodes.tsv` files.

### 3.1 Look before you download

List the samples and check that the compartment you care about is actually in the data. Many
deposits are immune-sorted or lose their epithelium at dissociation; a target test aimed at a
compartment the data barely hold is wasted.

```bash
python scripts/fetch_geo_10x.py GSE286057 --out data/ovarian/raw --keep tumor --keep distant --keep proximal --dry-run
```

`--keep` is a case-insensitive **substring** match on the GEO sample title (mind that `tumor` also
matches `peritumoral`; check the printed sheet). Then download two libraries, one per condition,
build them, and look at the composition:

```bash
python scripts/fetch_geo_10x.py GSE286057 --out data/ovarian/probe --keep P14_tumor --keep P1_proximal
python scripts/build_from_10x_samples.py --sheet data/ovarian/probe/samples.csv --out data/ovarian/probe/probe.h5ad
python scripts/sctarget.py inspect data/ovarian/probe/probe.h5ad --name probe --force
#  -> read the "compartment composition" line: fraction of cells positive for immune / epithelial /
#     fibroblast / endothelial markers per condition. <5% epithelium everywhere = no epithelial target here.
```

### 3.2 Download, build, inspect

```bash
python scripts/fetch_geo_10x.py GSE286057 --out data/ovarian/raw --keep tumor --keep distant --keep proximal --keep distal --keep lesser_omentum
# edit data/ovarian/raw/samples.csv: set the `condition` column (two values), drop rows you do not want
python scripts/build_from_10x_samples.py --sheet data/ovarian/raw/samples.csv --out data/ovarian/ovarian.h5ad
python scripts/sctarget.py inspect data/ovarian/ovarian.h5ad --name ovarian
```

`inspect` guesses the sample, condition and cell-type columns and the species, reports whether `X`
holds integer counts, prints samples per condition and the compartment composition, and writes two
files to edit: `config/ovarian.yaml` and `data/ground_truth/ovarian.yaml` (a template). For a GEO
"one dense matrix + metadata table" deposit use `scripts/build_from_matrix.py --dry-run` first,
then build; for an `.h5ad` you already have, skip straight to `inspect`.

### 3.3 Verify the labels before any run

Public metadata is wrong more often than you expect (in one validation set tumour and benign were
swapped for 4 of 13 libraries). Check each library against a marker you trust, in both directions:

```bash
python - <<'EOF'
import anndata as ad, numpy as np, pandas as pd
a = ad.read_h5ad("data/ovarian/ovarian.h5ad")
def pos(g): return np.asarray((a.X[:, a.var_names.get_loc(g)] > 0).todense()).ravel()
m = {g: pos(g) for g in ["PTPRC", "EPCAM", "PAX8", "FOLR1", "MSLN"]}      # your markers here
rows = []
for s in sorted(a.obs["sample"].unique()):
    k = (a.obs["sample"] == s).to_numpy()
    rows.append({"sample": s, "condition": a.obs.loc[k, "condition"].iloc[0], "cells": int(k.sum()), **{g: round(float(v[k].mean()), 2) for g, v in m.items()}})
print(pd.DataFrame(rows).sort_values("condition").to_string(index=False))
EOF
```

A library that runs against its label is **not** silently dropped or relabelled. Decide now, in
writing, how it will be handled: keep it in the primary run, and if the result hinges on it, run a
second arm with `dataset.exclude_samples: [...]` and `dataset.exclude_samples_reason: "..."` in a
copy of the config. Both arms are reported. (The critic also runs a marker-free version of this
check, `check_reference_libraries`, during every run.)

### 3.4 Edit the config

`config/ovarian.yaml` after `inspect`, with the fields you will typically touch:

```yaml
dataset:
  path: data/ovarian/ovarian.h5ad
  species: homo_sapiens              # mouse / dog / ... are mapped to human symbols for drug lookups
  sample_key: sample                 # one value per library
  condition_key: condition
  celltype_key: null                 # null -> lineages are derived from the data; or an obs column
  disease: Tumor                     # the case condition
  reference: [Omentum]               # one or more reference conditions
  condition_descriptions:
    Tumor: "high-grade serous ovarian cancer metastases in the omentum (6 patients)"
    Omentum: "omentum without tumour: distant omentum from the same patients and 8 non-cancer donors"
  counts_layer: null                 # null = counts are in X; the loader keeps a copy before normalising
  forbidden_obs_columns: [title, patient]   # anything identifying: never shown to the model
  exclude_samples: []                # a disclosed second arm, with exclude_samples_reason
targets:
  objective: target                  # or `marker` for a perturbation / response study
  scoring_version: 3
```

### 3.5 Write and freeze the answer key (if you know one)

If a target is already known and functionally validated for this disease, write
`data/ground_truth/ovarian.yaml` **now**, commit it, tag it, and never edit it after the first run.
If it turns out to be wrong, say so in the write-up; do not fix it.

```yaml
expected_states:
  - name: HGSOC_tumour_epithelium
    marker_genes: [PAX8, WT1, EPCAM, KRT7, KRT8, KRT18, CLDN4, WFDC2, MUC16, FOLR1, MSLN, CDH1]
    min_marker_overlap: 5
expected_targets: [FOLR1]                      # mirvetuximab soravtansine
acceptable_targets: [MUC16, MSLN, TACSTD2, CLDN6, ERBB2, EPCAM]
expected_drugs: [mirvetuximab soravtansine]
decoy_targets: [MKI67, TOP2A, KDR, FLT1, PDGFRB, ACTA2, CSF1R]
per_cell_truth: null
```

```bash
git add config/ovarian.yaml data/ground_truth/ovarian.yaml && git commit -m "ovarian: key frozen before any run"
git tag ovarian-key-frozen
```

Two lessons from the validation datasets, written into the key before the run: the **expected target
must be transcriptionally distinctive** in the disease state (a protein-level target such as
HIF-2alpha/EPAS1 is invisible to an expression method and tests the key, not the method; put it
under `acceptable_targets`), and the **state markers should be genes you expect to be consistent
across patients**, not the textbook list. With no known answer, skip the key: the pipeline still
produces hypotheses with critic flags; only the evaluator does not apply.

### 3.6 Calibrate (free: no model, no tokens)

```bash
python scripts/sctarget.py calibrate --config config/ovarian.yaml --ground-truth data/ground_truth/ovarian.yaml
```

A fixed scripted procedure runs through the same tools (10 minutes for 30k cells, about an hour for
200k). Read `runs/ovarian_calibrate/report.md`, in this order: section 1, QC per lineage (a lineage
that lost more than half its cells means the thresholds do not fit the tissue; kidney tubule and
cardiomyocytes run 20–45% mitochondrial reads normally); section 3, which clusters are enriched
(tiny or immune-only ones mean the composition route is favouring consistent infiltrates over
patient-split malignant cells; look at `split_state_candidates`), doublet and dissociation-stress
suspects, and the recorded signatures; then the evaluation line. Iterate on the **config** only
(thresholds, caps, excluded libraries with a reason). Never change code or prompts to fit one
dataset.

### 3.7 Run the agent

```bash
python scripts/sctarget.py run --config config/ovarian.yaml                       # Claude Code login (default)
python scripts/sctarget.py run --config config/ovarian.yaml --model claude-sonnet-4-5   # API key in ANTHROPIC_API_KEY
```

Expect 1–2.5 hours and 1.5–2.5M input tokens for 150k cells; 20–90 minutes and 0.6–1.7M for
30k–200k cells on simpler designs. Run one dataset at a time (the druggability services
rate-limit). `run` refuses to start while the config still has scaffold placeholders and warns if no
key is committed. Run it in the background on a cluster:

```bash
nohup python scripts/sctarget.py run --config config/ovarian.yaml --run-dir runs/ovarian_agent_1 > runs/ovarian_agent_1.log 2>&1 &
tail -f runs/ovarian_agent_1.log
```

### 3.8 Score and tabulate

```bash
python scripts/run_eval.py --ground-truth data/ground_truth/ovarian.yaml runs/ovarian_agent_1
python scripts/sctarget.py report --ground-truth data/ground_truth/ovarian.yaml --out docs/results_ovarian.md runs/ovarian_*
```

The table has one row per run, with the arm label read from the run's own recorded configuration
(model, prompt and scoring versions, lineage source, excluded libraries), so every number in a
write-up traces to a run directory.

## 4. Your own dataset

The same four commands, from an `.h5ad` with raw counts:

```bash
python scripts/sctarget.py inspect   data/my.h5ad --name mydata      # columns, species, composition; writes config + key template
#   edit config/mydata.yaml; verify labels per library (3.3); write and commit the key if you know a target (3.5)
python scripts/sctarget.py calibrate --config config/mydata.yaml     # no tokens
python scripts/sctarget.py run       --config config/mydata.yaml     # the agent
python scripts/sctarget.py report    --ground-truth data/ground_truth/mydata.yaml runs/mydata_*
```

Export from Seurat or Scanpy with integer counts in `X` or a layer, one obs column per library and
one per condition; existing `leiden` / `X_umap` are reused if present and `already_processed: true`.
Non-human species: set `species` and the drug lookups go through human orthologs; unmapped genes
are marked `assumed_or_unverified` in the report.

## 5. Choosing the model backend

| backend | how | notes |
|---|---|---|
| `claude-code` (default) | `claude` CLI logged in to your claude.ai account (`claude auth status`) | No API key; billed to the subscription. Tool calls are emulated through a strict JSON protocol; the CLI's own tools are disabled. No seed or temperature control. Transient CLI errors (token refresh, rate limit) are retried. |
| API model id, e.g. `--model claude-sonnet-4-5` | `export ANTHROPIC_API_KEY=...` | Native tool use, `--seed`/temperature available. Never put the key in a file or a commit. |
| `--model scripted` | nothing | Deterministic stand-in used by `calibrate`, the demo and the tests. The floor the agent is compared against. |

Both model backends see identical tools, prompts and payloads; the backend is recorded in every run.

## 6. Reading the report

Open `runs/<run>/report.html` (self-contained) or `report.md`. Read **section 5, the critic's
flags**, first; then **section 4**, the evidence table and the per-library panels; then the
per-state tables.

- **Per-library panels**: one point per library (pseudobulk log-CPM), by lineage and recorded state,
  coloured by condition. A real target shows most disease libraries above the reference ones in its
  own lineage; a single high library is a patient effect; a shift in every lineage is ambient RNA.
- **Evidence table** per nominated gene: `State log2FC` (vs sibling clusters of the same lineage),
  `Lineage log2FC` (disease vs reference cells of the lineage) and `Axis` (which lineage contrast
  carried it), specificity, consistency across disease libraries, attribution, compounds, flags,
  critic notes.
- **`target_class`**: `state_specific` (defines this state among its lineage) or `disease_wide`
  (up in disease cells across the lineage). Both can be targets; they are different claims.
- **`attribution`**: `reference_driven` or `ambient_like` means the signal is not from these cells;
  `shared_with_reference` means a lineage or pan-condition gene; `reference_pool_in_cluster` says
  whether the in-cluster reference comparison was usable (`too_thin`, `rare_in_reference`) and
  `rare_reference_counterpart` marks a gene that a small normal population also expresses — an
  on-target consideration, not a disqualification. Attribution is the most conservative step; when
  it sets a plausible target aside, read the critic's reason, not just the class.
- **`proliferation_marker`**: a cell-cycle gene; it marks that cells divide, not what the disease is.
- **Direction**: druggability counts compounds and does not know whether inhibiting the target
  helps; the critic has caught agonists and HIF-stabilising "inhibitors" before. Ask it.
- **Library label check**: the critic's marker-free check of which libraries resemble the opposite
  condition, with correlation heatmaps per lineage. A flagged library is a prompt to look and, if it
  cannot be resolved, to run the disclosed exclusion arm.
- **Novel candidates** above a known target are hypotheses, not errors; two in the ccRCC validation
  (CD70, FOLH1) turned out to be clinical-stage targets the key had omitted.

Every number traces to `runs/<run>/run.jsonl`, which holds each tool call with the exact
`result_sent` payload the model saw. `tables/` holds the full result tables; `state.json` the facts.

## 7. The rules that keep a run honest

- The **key is frozen** before the first run and never edited after; a wrong key is reported, not fixed.
- Every change to prompts, scoring or analysis rules is a **labelled version** (`prompt_version`,
  `scoring_version`, recorded in `config_used.json`, the report header and the eval line) with a
  section in `docs/EVALUATION_PROTOCOL.md` naming the failure it fixed. Nothing is changed to make one
  dataset pass; config-level thresholds and disclosed exclusion arms are the only per-dataset knobs.
- Comparison arms: a naive baseline (`scripts/run_baseline.py`) and the scripted procedure, so the
  agent's contribution is measured, not assumed.
- `scripts/rescore.py --recompute-attribution --full-candidates` re-runs attribution and ranking on
  saved runs when a rule changes, so regressions are checked without new tokens.
- Candidates are hypotheses for IHC / functional follow-up. The report says so; so should you.

## 8. Validation so far

Public datasets, each run from the frozen key, with the number of runs and the final-arm outcome.
Full tables are regenerated from the run directories (`docs/results_*.md`); `docs/REPORT.md` is the
write-up.

| dataset | design | known target | outcome |
|---|---|---|---|
| ccRCC, GSE159115 | 7 tumour / 6 benign kidney libraries, 29k cells | CA9 (HIF programme; EPAS1 is protein-level and was wrongly expected — documented) | CA9 **ranked 1st** in every tumour state, with and without the authors' annotation; strict state recovery met; GEO tissue labels found swapped for 4 libraries and corrected from the authors' files |
| HCC, GSE149614 | 10 tumour / 8 non-tumour liver, 63k cells; reviewer-style from a fresh clone; rebuilt from GEO in round 0.4b | GPC3 (codrituzumab) | unscored with all libraries (a "normal" library expressed the tumour programme; padj 0.20); **1st in its state with codrituzumab** once that patient was excluded in a pre-registered arm (round 0.2: 4th); the 0.4b agent run scored it 2nd in a stromal state and the evaluator now says so (protocol 23i) |
| Ovarian omental metastases, GSE286057 (2025) | 6 tumour / 23 omentum, 124k cells after the 0.4b rebuild; reviewer-style, first model-driven runs | FOLR1 (mirvetuximab soravtansine) | found in the tumour state in every agent run (**14th-17th** within it, five runs over rounds 0.3-0.4b) and held back each time: two contaminated reference libraries in the primary run, then 59-151 benign FOLR1+ epithelial cells in normal omentum; the 0.4b label check shows that donor before any run; strict state recovery met four times; limitation stated |
| IPF lung, GSE135893 (round 0.4) | 12 IPF / 10 control donors, non-malignant disease vs healthy | SPP1+ profibrotic macrophage state; aberrant basaloid / ITGB6 | SPP1 macrophage state **recovered strictly, SPP1 1st in its state** under prompts v8; aberrant basaloid state not reached (no run subclustered the transitional cluster); no decoy nominated |
| Doxorubicin vs vehicle mouse heart, GSE285475 (round 0.4) | 4 vs 4 wild-type hearts, perturbation, `objective: marker` | p53 programme, Cdkn1a | cardiomyocyte stress / atrophy programme found in every run (Trim63, Nppb, Fbxo32, Nppa: acceptable genes); Cdkn1a not induced at two weeks (time-point key failure, recorded); critic asked for the pre-registered batch arm |
| Fibrolamellar carcinoma snRNA-seq (lab data, round 0.4b) | 9 primary / 9 normal liver libraries (7 metastatic set aside), 254k nuclei | PRKACA (DNAJB1-PRKACA fusion); CA12, ERBB2 acceptable | tumour state reached by the expression route (tumour nuclei co-clustered with normal hepatocytes); **PRKACA 2nd within its state behind CA12** in the floor and both agent runs, nominated, neratinib found; literature markers fail the state criterion (3/5) |

Across all runs that reached nomination, two decoys were ever nominated without a critic flag
(TOP2A before the proliferation flag; PDGFRB from a CAF programme before the microenvironment flag,
prompts v11). Nobody but the author has run the pipeline yet; that is the next test.

## 9. What is in the box

```
sctarget/
  config.py             config loading (config/default.yaml + per-dataset yaml + env)
  session.py            AnalysisSession: holds AnnData + results locally; tools operate on it
  privacy.py            sanitize(): truncation, per-cell rejection, sample pseudonymisation
  runlog.py             append-only JSONL audit log
  envcheck.py           the environment check behind `sctarget.py check-env`
  tools/
    registry.py         @tool decorator -> JSON schema, role permissions, argument validation
    analysis.py         load_dataset, run_qc (adaptive, preview), detect_doublets, integrate_and_cluster,
                        derive_lineages / label_lineages, subcluster, merge_clusters, cluster_composition,
                        differential_abundance (+ split-state candidates), disease_programs, find_markers
                        (stress + cycling flags), pseudobulk_de (sibling / within-cluster / compartment /
                        whole-lineage contrasts), gene_specificity, condition_attribution,
                        check_reference_libraries, enrichment, plot_genes
    knowledge.py        map_orthologs (Ensembl), druggability_lookup (DGIdb + Open Targets),
                        literature_count (PubMed); disk-cached, 429 backoff, offline fixture fallback
    targets.py          score_targets (deterministic, weighted, explainable; scoring v1-v3), record_finding
  agents/
    model.py            AnthropicModel (API tool use), ScriptedModel (deterministic floor)
    claude_code.py      ClaudeCodeModel (subscription login; JSON tool protocol; transient retry)
    loop.py             the tool loop: validate -> run -> sanitize -> log -> feed back
    prompts.py          role prompts (PROMPT_VERSION history at the top)
    orchestrator.py     runs the roles in order, persists state.json, builds the report
  report/               builder.py (markdown/HTML), figures.py (per-library panels, heatmaps)
  eval/                 score.py (state recovery, target rank, decoys, novel, drug recall), stability.py
scripts/
  sctarget.py           check-env | demo | inspect | calibrate | run | report
  fetch_geo_10x.py      GEO per-sample 10x files -> sample sheet
  build_from_10x_samples.py  sample sheet -> .h5ad
  build_from_matrix.py  GEO dense matrix + metadata -> .h5ad
  build_ccrcc.py        worked example: the ccRCC validation set (label correction from authors' files)
  new_dataset.py        behind `inspect`
  run_pipeline.py, run_baseline.py, run_eval.py, rescore.py, collect_results.py, make_synthetic.py
config/                 default.yaml, synthetic*.yaml, ccrcc.yaml, hcc*.yaml, ovarian*.yaml
data/ground_truth/      frozen keys: synthetic, ccrcc, hcc, ovarian
docs/                   REPORT.md, EVALUATION_PROTOCOL.md (22 sections), results_*.md
skills/sc-target-agent/ the Claude Code skill (also under .claude/skills/)
tests/                  87 tests
```

Each run directory holds `report.md`, `report.html`, `run.jsonl`, `state.json`, `config_used.json`,
`tables/*.csv`, `figures/*.png` and `cluster_assignments.csv` (local only). `--resume-from
runs/<earlier>` reuses an earlier run's analysis and runs only the nominator and critic.

## 10. How a candidate is scored

`score_targets` (scoring v3) ranks candidates by a weighted sum, weights in `targets.weights`:

| component | source | default weight |
|---|---|---|
| effect size: the larger of the state axis (vs sibling clusters) and the lineage axis (disease vs reference cells of the lineage; which contrast is used is recorded as `lineage_axis_source`) | `pseudobulk_de`, else `find_markers` | 0.25 |
| specificity (log2 target-cluster mean / max other non-disease cluster) | `gene_specificity` | 0.25 |
| druggability (DGIdb class + Open Targets tractability) | `druggability_lookup` | 0.25 |
| known compounds (approved inhibitors) | `druggability_lookup` | 0.15 |
| consistency across disease libraries | `gene_specificity` | 0.10 |

Flags travel with each gene into the report and the critic's review; several also multiply the
score: `low_specificity` and `shared_across_conditions` ×0.6, `inconsistent_across_samples` ×0.8,
`reference_driven` and `ambient_like` ×0.2, `shared_with_reference` and `weak` ×0.5,
`proliferation_marker` ×0.5. `disease_wide` and `rare_reference_counterpart` are informational.
Change weights in the config, not in code, so the choice is recorded in `config_used.json`.

## 11. Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `check-env` fails at neighbours or Scrublet with `'csr_array' object has no attribute 'getnnz'` | a newer `pynndescent` with an older `scanpy` | `conda env create -f environment.yml` (pinned), or `pip install 'pynndescent<0.7'` |
| calibration report says "Run aborted: tool X failed" | an environment or data-format problem, not a biological negative | read the error in the finding; `check-env`; check `counts_layer` and the obs columns |
| `run` says the config still has placeholders | `inspect` scaffold not edited | edit `config/<name>.yaml` (`SAMPLE_COLUMN`, `DISEASE`, `DESCRIBE`...) |
| `pseudobulk_de needs raw counts` | counts not integer, or not found | ensure `X` or `counts_layer` holds integer UMI counts (the loader stashes `X` automatically when it looks like counts) |
| the agent run dies at the first call with an OAuth / token-refresh message | two `claude` processes refreshing at once | it is retried automatically now; otherwise wait a minute and rerun |
| druggability shows `source: none` or empty Open Targets fields | HTTP 429 rate limit | retried with `Retry-After`; run one dataset at a time; rerun the nominator with `--resume-from` |
| `inspect` warns "fewer than 5% epithelial-marker-positive cells" | the deposit is immune-sorted or lost its epithelium | an epithelial target cannot be tested here; use the dataset for immune/stromal questions or pick another |
| `fetch_geo_10x.py` selects 0 samples, or too many | `--keep` is a substring match on titles | check the printed titles; `tumor` also matches `peritumoral` |
| a "normal" library expresses the disease programme | mislabelled or contaminated reference | keep the primary run; declare and run the `exclude_samples` arm with a written reason; report both |

## 12. Limitations

- Expression-based: targets regulated post-translationally are invisible.
- The composition route favours infiltrates consistent across patients over malignant cells that
  partition by patient; merging mitigates this but does not remove it.
- Condition attribution is conservative and depends on how clean the derived lineages are: a rare
  benign population expressing the target inside the tumour cluster can set it aside (ovarian), and
  coarse derived lineages decide whether the guards against that apply.
- The nominator over-nominates (precision 0.1–0.3) and rarely subclusters without being asked.
- The subscription backend has no seed or temperature control; stability is behavioural and should
  be estimated by repeat runs (`run_eval.py --stability`).
- Differential abundance uses per-sample proportions and Mann–Whitney; with few replicates the
  p-values are coarse by construction.
- Validated so far on tumour datasets only; the programme route for non-tumour disease has only
  been exercised on synthetic data. Drug knowledge is human-centric.

## 13. Data handling

To the model (Claude Code or API): role prompts, tool descriptions, and tool results — cluster-level
tables (fold changes, enrichment statistics, marker lists, proportions), gene symbols and drug names.
Sample identifiers are pseudonymised (`S01`, `S02`, ...; stable across arms of one dataset) before
anything is sent; columns listed under `forbidden_obs_columns` are never shown. The privacy layer
rejects any payload the size of a per-cell array. To DGIdb, Open Targets and NCBI E-utilities: gene
symbols. Nothing else; there are no uploads. Check for yourself: every `tool_call` event in
`run.jsonl` has a `result_sent` field that is exactly what the model saw. Whether a given dataset may
be sent to an external AI service is your and your institution's decision; `SCTARGET_OFFLINE=1`
keeps even gene-symbol lookups local. Never commit datasets, run directories or API keys.

## 14. Using it from Claude Code

Inside a clone, Claude Code picks up `.claude/skills/sc-target-agent/SKILL.md` automatically; in any
other project, copy `skills/sc-target-agent/` into its `.claude/skills/`. The skill walks a session
through the steps above with the lessons from the validation datasets built in (what to check in
QC, how to read a thin or rare reference pool, why an expression method cannot nominate a
protein-level target, how to pre-register an exclusion arm). The same skill is submitted to the
Cornell CoPAIR skills incubator (`skills/genomics/sc-target-agent`).

## 15. Provenance

Knowledge sources: DGIdb, Open Targets Platform, Ensembl, PubMed / NCBI E-utilities, MSigDB
Hallmark. Analysis: scanpy, anndata, PyDESeq2, harmonypy, leidenalg, scrublet. Validation datasets:
GSE159115 (Zhang et al., ccRCC), GSE149614 (Lu et al., HCC), GSE286057 (ovarian cancer omental
metastases). Developed by Andreas Stephanou (Cornell University) as a CoPAIR fellowship project.
MIT licence.
