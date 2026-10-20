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

## 9. Scoring versions (post-freeze arms)

`targets.scoring_version` (default 1 = frozen) is recorded in every run as `scoring_version`, in
the report header and in the eval output, like `prompt_version`. A new version is a new arm.

| version | change | motivation |
|---|---|---|
| 1 | freeze: effect size = cluster vs sibling clusters (composition states) | - |
| 2 | second effect axis for composition states: disease vs reference pseudobulk within the *shared* clusters of the state's compartment (clusters of its dominant cell type minus the disease-enriched ones); the larger axis is used and recorded as `target_class` state_specific / disease_wide (disease_wide needs lineage log2FC >= 1 and > the state log2FC). Correction before first use: a first draft included the disease-enriched clusters in the compartment, which labelled nearly every candidate disease_wide because the state itself made up most of the disease compartment; no run was reported under that draft | v1 scores composition states for state-specificity only, so a gene up across all disease cells of the lineage and absent from the reference lineage (`disease_wide` in attribution) is detected and labelled but under-ranked; program states already had this axis - v2 removes the asymmetry |

`scripts/rescore.py` re-scores an existing run's candidate sets under another version without
re-running the agent (deterministic given the stored facts), writing `<run>_s<version>/`. It can
only re-rank genes the nominator already chose; a live run under the new version is still needed
to see candidates the nominator would pick with that scoring in place.

## 10. De novo clustering arm

All runs above reuse the dataset's existing clustering (`reuse_clustering`), i.e. the partition
from the original analysis. That partition encodes the original analysts' choices, so a state
found inside it is not fully independent evidence. The de novo arm strips the object to counts,
sample, condition and the broad lineage annotation (`data/caa_denovo.h5ad`, `config/caa_denovo.yaml`)
and lets the pipeline recompute HVGs, PCA, Harmony and Leiden. The broad annotation is kept (used
only to define sibling sets and compartments) and disclosed as inherited.

Tooling notes from the first attempt (`runs/caa_denovo_1`, not a valid arm):
* harmonypy >= 1 returns cells x PCs; scanpy's `harmony_integrate` wrapper transposes it and the
  result is rejected. `integrate_and_cluster` now calls harmonypy directly and normalises the
  orientation.
* The tool used to fall back to uncorrected PCA silently when Harmony failed. It now returns an
  error and tells the model it may call again with `use_harmony=false`; the choice is recorded.
  In the failed run the clustering split by patient and the critic correctly flagged every
  CAA-restricted cluster as single-patient and tied it to `harmony_used=false` - a useful
  robustness record, at 2.4M input tokens.
* `condition_attribution` now judges `disease_wide` on the lineage's shared clusters (other
  enriched states excluded, reported separately), the same definition `score_targets` v2 uses;
  the critic had flagged the label mismatch in `runs/caa_cc_v2s2_1`.

## 11. Attribution: thin reference-pool guard (post-freeze)

`condition_attribution` compares a gene's expression in the cluster's disease cells with its
expression in the reference cells *inside the same cluster* (the `reference_driven` and part of
the `shared_with_reference` tests). That comparison is only meaningful when the reference cells
are a non-trivial share of the cluster. A disease-specific state always attracts a few reference
cells that resemble it (doublets, misassigned nuclei, a rare real population), and a pool of a few
hundred cells among tens of thousands can carry a mean dominated by a handful of high expressers.
Since the fix, the in-cluster comparison decides the class only when reference cells are
>= `analysis.attribution.min_reference_frac_in_cluster` (default 2%) of the cluster; below that
the pool is reported in `reference_pool_in_cluster` but the sibling-cluster comparison carries
the judgement. Motivation: `runs/caa_denovo_1`, where ~150 reference cells (0.5%) inside a 31k-cell
CAA cluster were ERBB4-high and the gene was classed shared_with_reference and dropped before
scoring, while the same gene in the same population under the original clustering (5% reference
cells, ERBB4-low) was state-specific. The change names no gene; it mirrors the pseudobulk replicate
floor. Which samples such thin pools come from is worth a critic check (doublets vs a real rare
population).

`scripts/rescore.py --recompute-attribution --full-candidates` re-runs attribution under the
current rules and scores every gene the nominator examined (its `gene_specificity` calls), not
only the ones it chose to score: "what would this run have nominated with current rules and no
premature drop", at zero token cost. A live run remains the confirmation.

## 12. Annotation-free structure (post-freeze)

The structural logic (sibling sets, compartments, the ambient-control lineage) used to depend on
the dataset's cell-type annotation (`celltype_key`). It now runs on a *lineage key* that is the
annotation when one is configured and present, and otherwise a data-derived coarse clustering
(`derive_lineages`: Leiden at `analysis.lineages.resolution`, default 0.2, on the same embedding
as the fine clustering; each fine cluster belongs to the coarse group holding most of its cells).
Tools derive it lazily the first time they need it, so a dataset with no annotation runs
unchanged; `analysis.lineages.ignore_annotation: true` forces the derived hierarchy even when an
annotation exists (useful for a fully independent de novo arm: the CAA de novo runs still inherit
`Broad_anno` from the original analysis).

Names are separate from structure. `label_lineages` scores groups against a small marker table
(species-agnostic symbols, matched through the human-symbol column when present) and returns a
label with a confidence; nothing downstream reads labels. The rule is deliberate: a model's
annotation can vary between runs, the hierarchy cannot, and both ERBB4 losses in the agent runs
traced to how cell groups were drawn and compared.

Known costs: the coarse resolution is a knob (a singleton group gives its cluster no siblings and
is reported); a disease lineage that drifts far from its normal counterpart can form its own
coarse group, which weakens the state-vs-normal sibling contrast (the reference-aware compartment
axis of scoring v2 is the mitigation). Synthetic check: the e2e pipeline passes with
`celltype_key` removed (`tests/test_lineages.py`).

## 13. Any species, any two conditions

What is dataset-specific lives in the config and the answer key; the tools and prompts are generic.
`scripts/new_dataset.py <h5ad> --name X` inspects an object, guesses sample / condition / cell-type
columns, counts location and species from the gene symbols, and writes `config/X.yaml` plus a
ground-truth template with every guess marked for checking.

* Conditions: `disease` / `reference` may be written `case` / `control`; any two-condition design
  (treated vs vehicle, knockout vs wild-type, infected vs mock) maps onto them.
* Objective: `targets.objective: target` (default; druggability scored) or `marker` (a perturbation
  or response study: druggability reported but not scored).
* Species: human data skips the ortholog step (identity mapping, no network); other species map
  to human through `dataset.human_symbol_var` when provided, else Ensembl. Drug knowledge is
  human-centric; for non-model species the same-name assumption is reported and should be flagged.
* Cell-type annotation is optional (section 12).
* Before any model run on a new dataset: freeze and commit its key; run the scripted model once
  (free) as calibration and check how many clusters pass each gate (section 1-2); then one model run.
* The evaluator is for validating the method on datasets with a known answer; the pipeline runs
  without a key and produces hypotheses with critic flags.

## 14. Lessons from the second dataset's calibration (ccRCC, GSE159115)

Three generic defects surfaced on the first scripted run of a dissociated human tumour dataset,
none of which the snRNA-seq CAA data could have shown. All fixes are rule-level and name no gene.

* **Fixed QC thresholds do not transfer across tissues.** `max_mt_pct: 10` removed 61% of benign
  epithelium and ~99% of the distal nephron (thick ascending limb, DCT, CNT, intercalated cells,
  median mito fraction 25-45%, normal for those cells), i.e. the tumour's normal counterpart.
  `run_qc` is now adaptive per sample (median + n MADs, mito threshold floored at 20% and capped at
  50%; the lower gene bound stays fixed because an adaptive one cuts low-complexity minority
  lineages) and reports the fraction of each annotated type removed, warning when a lineage loses
  more than half of its cells.
* **A state split across clusters fails the composition test piecewise.** Malignant cells partition
  by patient: clusters 0/2/3 were each ~97% tumour and present in every tumour library, yet each
  missed the per-sample enrichment test (p 0.05-0.10) because every library's cells fall mostly in
  one piece. `differential_abundance` now reports `split_state_candidates` and `merge_clusters`
  tests their union (the mirror of `subcluster`). The scripted model merges automatically; the
  analyst prompt (v3) says when to.
* **Dissociation-stress clusters masquerade as disease states.** A 370-cell heat-shock cluster
  (DNAJB1, HSPA1A, HSPA6, BAG3, ATF3, JUN) topped the enrichment ranking. `find_markers` flags
  clusters whose top markers are dominated by the immediate-early / heat-shock programme
  (`stress_artifact_suspect`); the scripted model skips them; the prompt names the trap.
* **Ranking candidate states by fold-ratio alone favours small extreme clusters** (600 and 370
  cells) over the large consistent ones (2,000 cells). The scripted model now orders by
  log2 ratio x log10(cells); the prompt says the same in words.

Also from this dataset: the GEO `tissue` characteristics were swapped for four of thirteen
libraries (confirmed by CA9/NDUFA4L2); the builder takes tissue from the authors' annotation
files and keeps the GEO value for the record. The composition route, as a design, favours
infiltrates that are consistent across patients (every tumour has TAMs; benign kidney has few)
over malignant cells that partition by patient: the agent run characterised macrophage and T-cell
clusters, and the critic correctly refused both as lineage-wide. That bias is now stated.

## 15. QC as a recorded choice; doublets

`run_qc(apply=false)` previews what the thresholds would remove - overall, per sample (with each
sample's mito and gene-count distribution) and per annotated type - without changing the data.
The analyst is told to preview, adjust the caps if a lineage would lose most of its cells, then
apply. The data-driven rule stays the default, so the result is reproducible; the deviation is
what gets logged and justified.

`detect_doublets` runs Scrublet per sample and FLAGS doublets (doublet_score, predicted_doublet
on every cell) rather than removing them; `cluster_composition` reports each cluster's
doublet_frac and marks `doublet_suspect` clusters (>= 20% and >= 3x the dataset rate). A tumour-
immune or tumour-endothelial doublet cluster looks exactly like a disease-enriched, state-specific
population, so this is a guard for the composition route, not only hygiene. Removal is an explicit
option (remove=true). The scripted model flags and skips suspect clusters; it does not remove.

Refinements after the second ccRCC calibration run (`runs/ccrcc_dry2`):
* The stress flag rejected a 2,300-cell tumour cluster (98% tumour, CA9 in 84% of cells) because
  30% of its sibling-contrast markers were heat-shock genes: the dissociation programme was
  superimposed on real tumour cells and differed by condition because tumour and benign tissue
  were handled differently. `find_markers` now removes the stress genes from the marker list and
  calls a cluster an artefact only when fewer than 10 strong markers remain without them.
* Scrublet's automatic threshold found no doublets in any library while three clusters were
  visibly doublets (CA9+/PECAM1+, CA9+/PTPRC+, PECAM1+/PTPRC+). `detect_doublets` falls back to
  the score quantile at the expected rate when the automatic threshold finds nothing, and says so.
* The merged cluster the merge step produced was never characterised (two-state cap, one slot
  spent on a cluster dropped late as stress). Merged clusters are characterised first; cap 3.

