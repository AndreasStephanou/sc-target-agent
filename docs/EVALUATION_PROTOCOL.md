# Evaluation protocol

The project's claim is not "the agent found ERBB4" – you already know ERBB4. The claim is
*"given only the data, the agent recovers the known answer, ranks it above plausible decoys,
and does so consistently, while a straightforward script does not."* The protocol below is
what makes that claim defensible.

## 1. Freeze the ground truth first

Fill `data/ground_truth/caa.yaml` from the published analysis **before** the first agent run
on the real data and commit it:

* `expected_states`: each CAA-specific epithelial state with 6-10 marker genes from the paper
  and `min_marker_overlap` (4 is a reasonable bar).
* `expected_targets`: ERBB4. `acceptable_targets`: PEG3 (diagnostic marker, not a drug target).
* `expected_drugs`: neratinib (plus other approved pan-ERBB inhibitors you would accept).
* `decoy_targets`: genes that are druggable but obviously not the answer – pan-tumour
  proliferation genes (MKI67, TOP2A), EGFR if it is also high in COSCC, strong cell-type
  markers of stromal/immune cells. Decoys sharpen the false-lead metric.

Do not edit the file after the first run. If you must (e.g. a typo), note it in the commit.

## 2. Blind run conditions

* The agent sees only the `.h5ad`, the config, and the tool outputs. Do not put the paper,
  gene names, or hints in the config or prompts. The system prompts in `agents/prompts.py`
  are generic; keep them so (the only disease-specific string is the condition label).
* Use a fixed model version string (e.g. `claude-sonnet-4-5`) and record it. `run.jsonl`
  logs it per call.
* Set `SCTARGET_OFFLINE` unset so druggability comes from the live databases, and keep the
  `cache/` directory so repeated runs are reproducible and cheap.

## 3. Runs

| arm | command | n |
|---|---|---|
| agent | `run_pipeline.py --config config/caa.yaml --model <model>` | 5 |
| baseline | `run_baseline.py --config config/caa.yaml` | 1 (deterministic) |
| agent, second model (optional) | same with a different model ID | 3-5 |

Then:

```bash
python scripts/run_eval.py --ground-truth data/ground_truth/caa.yaml --stability --json eval_caa.json runs/agent_*
python scripts/run_eval.py --ground-truth data/ground_truth/caa.yaml runs/baseline_*
```

## 4. Metrics to report

| metric | what it answers |
|---|---|
| state recovery, lenient (found / marker overlap / cell-level F1) | did the analysis reach the right cell state? (counts the cluster's DE and marker tables) |
| state recovery, strict (found_strict) | did the agent record the state's genes itself? the lenient-strict gap measures signal the tools produced but the agent dropped |
| rank of ERBB4; recall@5 | is the right target near the top? |
| decoys above target / novel candidates above target | known-wrong genes above the answer (the error metric) vs genes absent from the answer key (new hypotheses - report, do not penalise) |
| nominated precision; decoys nominated **unflagged** | did the critic catch the decoys that were nominated? |
| drug recall | did it reach the actionable compound? |
| mean Jaccard of top-10 across runs; rank-1 consistency | is the answer stable or a coin flip? |
| tool calls, tool errors, tokens | cost and robustness |

Report agent (mean ± range over runs) next to baseline. The interesting cells are usually
*false leads* and *decoys nominated unflagged*: that is where the sibling-cluster contrast,
`gene_specificity` and the critic earn their keep.

## 5. Generalisation (second dataset)

A result on one dataset you co-authored is suggestive, not conclusive. Add one public
dataset with a known, druggable, state-specific target (any disease-vs-normal design with
≥3 samples per group and a published target works) and repeat 1-4. If the agent fails
there, the failure mode is the finding.

## 6. Ablations worth running (cheap, informative)

Each is a config or one-line prompt change; report the metric deltas:

* critic off (`--stages planner analyst target`) → how many decoys get through?
* `pseudobulk_de` disallowed for the analyst role → does cell-level DE inflate false leads?
* `gene_specificity` disallowed → does EGFR-like pan-tumour gene rise in rank?
* scoring weights: druggability 0 → does the ranking become "marker list" rather than targets?

## 7. What to write up

Architecture diagram (README), the frozen ground truth, the metrics table, 2-3 example
critic flags verbatim from `report.md`, the audit trail excerpt (a chain from a report
sentence to the `run.jsonl` event that produced the number), and the ablations. State the
limitations plainly: n=3 per group, one disease, hand-tuned weights, and that the agent
nominates hypotheses – IHC and the kinase screen remain the validation.

## 8. Prompt / interface versions (post-freeze arms)

The frozen method (`method-freeze-v1`) is the analysis: tools, thresholds, scoring. The role
*interface* - what each role is told and what `record_finding` accepts - is versioned separately
as `PROMPT_VERSION` (in `agents/prompts.py`, written to `config_used.json`, the report header and
the eval output). Changes to it after the freeze are allowed only as a new, labelled arm, run with
its own set of runs and reported next to the previous version, never as a replacement.

| version | change | motivation (which runs) |
|---|---|---|
| 1 | freeze | - |
| 2 | nominator prompt carries `find_markers` top genes beside the pseudobulk list and asks for the union as candidate set; `record_finding(kind='disease_state')` rejects fewer than 5 genes | `caa_cc_1`, `caa_cc_2`: the analyst recorded both states with no `genes`, and the nominator (which is told "pseudobulk first, markers second") was never shown the marker tables, so a known target present only in the marker table was unreachable - a hand-off defect of the interface, not of the analysis |

