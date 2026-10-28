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

## 16. The disease-wide axis in the hand-off (prompt v4)

In `runs/ccrcc_dry3` the merged tumour cluster matched the key 5/5 by the lenient metric, yet its
recorded signature and candidate pool held none of the tumour-defining genes, because both came
from the sibling contrast and the siblings of a tumour cluster were mostly other tumour clusters:
tumour-vs-other-tumour is flat for CA9 (log2FC -0.06) and NDUFA4L2 (+0.04). The compartment
contrast (disease vs reference cells within the lineage's shared clusters, already the second
effect axis of scoring v2) carried them at +12 and +15. The analyst now runs that contrast for
every composition state and records the union of both top lists; the scripted model does the
same; the nominator's brief carries the compartment top genes. This is the same kind of hand-off
rule as prompt v2, and the same gap that lost ERBB4 in the CAA cluster-16 de novo run.

A key lesson: in the same tables EPAS1, the belzutifan target, is up only 2.5-fold on the
compartment axis (padj 0.09, baseMean 77) while its transcriptional targets are up thousands-
fold. HIF-2alpha is activated post-translationally (VHL loss), so an expression-based pipeline
cannot be expected to nominate it; the key's choice of EPAS1 as the expected target was a
drafting error about what "target" means for this method. The key stays frozen; the result is
reported with this explanation, and CA9 (acceptable in the key) is the expected expression-level
hit. For future keys: an expected target must be transcriptionally distinctive in the disease
state, or the test is of the key rather than the method.

## 17. Finishing guard

In `runs/ccrcc_agent_2` the analyst ran the full procedure (QC preview and apply, doublets,
clustering, sibling and compartment contrasts for three clusters, the programme route) and then
returned an empty reply - no tool calls, no text - which the backend took as "done". No state
was recorded; the nominator still ranked CA9 first from the tables, but by the evaluator's rules
an unrecorded state was never found. Two generic fixes: the Claude Code wrapper treats an empty
reply as a malformed one (retry with a notice), and the orchestrator gives an analyst that ends
with contrasts computed but no `disease_state` and no "nothing qualifies" summary one short
finishing turn (6 turns) that says exactly that. The guard is backend-independent and is logged
as a stage event.


## 18. Reviewer run on a third dataset (HCC, GSE149614): environment and false-negative guards

Before publishing, the author re-cloned the repository from GitHub into a fresh directory and a
fresh conda environment and followed only the README and the skill, on a dataset the method had
never seen: hepatocellular carcinoma, 10 primary-tumour and 8 non-tumour liver libraries (PVTT
and lymph-node libraries excluded at build), human, authors' six-class annotation in which
`Hepatocyte` covers normal and malignant cells alike. The key (`data/ground_truth/hcc.yaml`,
tag `hcc-key-frozen`) was written from the literature before any run: expected target GPC3
(absent from adult hepatocytes, re-expressed in most HCC, codrituzumab / CAR-T target, and
transcriptionally distinctive, the lesson of section 16); acceptable AFP, MDK; decoys spanning
proliferation, hepatocyte lineage, endothelium, cholangiocyte and EGFR.

Label check before the run (the skill's step 2): GPC3 in 25-98% of hepatocyte-class cells in
tumour libraries, 0-5% in non-tumour, for 17 of 18 libraries. HCC04N (labelled non-tumour) has
GPC3 77%, AKR1B10 60%, SPINK1 91%, MDK 57% while HCC04T is a GPC3-negative tumour. Unlike the
ccRCC swap (section 13) there is no second authors' file to settle it, so the data stay as
deposited and the anomaly is recorded here; an arm without patient 4 is the disclosed follow-up
if the result hinges on it.

The first calibration run produced "no disease-specific state found" with 7 tool errors. Three
generic defects, none visible on the two earlier datasets:

1. **Install not reproducible for real data.** `pynndescent 0.7.0` (released after the ccRCC
   runs) returns scipy `csr_array`; `scanpy 1.11.5` calls `.getnnz()` on it. The failure is only on
   the approximate-neighbours path (> ~8k cells) and inside Scrublet (which doubles each sample
   with simulated doublets), so the dense 2k-cell synthetic demo and all 66 tests passed on an
   install that could not cluster any real dataset. Fix: `pynndescent<0.7` and `umap-learn<0.6`
   pinned in `environment.yml` and `pyproject.toml`; a new `sctarget.envcheck` module
   (`python scripts/sctarget.py check-env`, also run first by `demo`) builds a sparse 9k-cell
   dataset and runs normalise, HVG, PCA, neighbours, Leiden, Scrublet, Harmony and DESeq2 through
   the installed packages; a `slow` test does the same.
2. **A crashed tool read as a negative result.** The scripted analyst treats an errored tool like
   an empty one and continued through subclustering, re-clustering and the programme route,
   narrating "no cluster is enriched". The calibration report read as a clean biological negative.
   Fix: `ScriptedModel` aborts the role when a critical tool (`load_dataset`, `run_qc`,
   `integrate_and_cluster`, `differential_abundance`, `merge_clusters`, `derive_lineages`) returns
   an error, records a "Run aborted: tool X failed" finding with the error text, and ends with an
   explicit abort message. Tested. (The agent backends already read the error and act on it.)
3. **No generic builder for the common GEO deposit.** The repository had only the ccRCC-specific
   builder; a reviewer with a "one dense count matrix + metadata table" deposit has to write code
   before step 1 of the skill. Fix: `scripts/build_from_matrix.py` (orientation detection, join on
   cell id, site filter, integer check, `sample`/`condition`/`author_celltype` columns), with a
   test, and referenced from the skill.

None of these changes the analysis; prompt and scoring versions are unchanged. The HCC
calibration is re-run on the identical object and key after the environment is rebuilt from the
pinned file, and its result is reported in `docs/results_hcc.md` whatever it is.

### 18b. Raw counts in X were destroyed by normalisation

The second HCC calibration (pinned environment) clustered, merged a patient-split tumour pair and
found five enriched clusters, three of them tumour hepatocytes present in all ten patients, but
every pseudobulk, compartment-contrast and scoring call failed with "pseudobulk_de needs raw
counts". Cause: the object held counts in `X` with `counts_layer: null` (what the scaffold writes
when `X` looks like counts), and `_ensure_lognorm` kept a raw copy only when a layer was
configured; the first clustering call normalised `X` in place. Both earlier datasets had a counts
layer, so the path was never exercised. The analyst then recorded the tumour clusters from
cell-level markers only: proliferation genes for the cycling sub-state, housekeeping and plasma-
protein genes for the others, 0/5 overlap with the key.

Fix: `load_dataset` stashes `X` into `layers['counts']` (and sets `dataset.counts_layer`) when `X`
holds integer counts and no layer is configured, picks an integer layer if one exists, and warns
loudly if no raw counts exist anywhere; the choice is in the load summary and the run facts.
Tested. The analysis itself is unchanged; the third calibration runs on the identical object.

## 19. Scoring v3 and prompts v5: the whole-lineage contrast and the proliferation flag

The third HCC calibration (`runs/hcc_calibrate_3`, pinned environment, counts preserved) ran the
full chain with no errors: adaptive QC kept 97% of cells, three tumour-hepatocyte clusters were
enriched in all ten patients (15: cycling; 18; 14), a patient-split pair (13+22) was merged, a
doublet cluster was flagged, tumour-endothelium and T/NK programmes were recorded. GPC3 was never
scored and TOP2A was nominated without a flag - the first unflagged decoy in three datasets.

**Why GPC3 was missed.** The sibling contrast is flat (-0.13: the siblings are other tumour
clusters, as in ccRCC). The compartment contrast of scoring v2 compares disease with reference
cells *within the lineage's shared, non-enriched clusters*; its result had GPC3 at +2.4 log2 but
padj 0.36, because in a tumour 90-100% of the malignant cells sit in the enriched clusters and the
disease side of that contrast was the atypical 9% (about 1,400 of 16,000 hepatocyte-class tumour
cells) that co-clustered with normal hepatocytes, against about 450 reference cells. The
attribution tool saw the thin reference pool (`too_thin, 11 cells` in cluster 15) but had no
other lineage-level contrast to fall back on. The tumour-vs-normal-hepatocyte comparison a
hepatologist would run first was never computed. In the same run GPC3, MDK and NDUFA4L2 appear
"up in disease" inside T-cell, myeloid and B clusters: tumour ambient RNA, the right genes in the
wrong cells.

**Why TOP2A got through.** Cluster 15 is a cycling tumour sub-state; its markers against the
other tumour clusters are cell-cycle genes, and neither scoring nor the critic had a notion of a
proliferation programme.

**Changes (scoring v3, prompts v5), all generic, no gene from any key named in code:**

1. `pseudobulk_de(contrast='disease_vs_reference_within_lineage')`: all disease cells of the
   lineage vs all its reference cells, enriched clusters included. Both lineage-level contrasts
   now report `lineage_coverage` (what fraction of the lineage's disease cells they used).
2. Scoring v3 takes the lineage axis from the shared-cluster contrast when it is representative
   (>= 25% of the lineage's disease cells and >= 3 reference replicates) and from the whole-lineage
   contrast otherwise; `lineage_axis_source` is recorded per gene and the decision with its
   reasons in the score_targets result. `target_class` semantics are unchanged.
3. A bundled cell-cycle list (Tirosh et al. 2016 S and G2/M genes plus recurrent cycling markers,
   135 genes) gives: `find_markers` -> `cell_cycle_fraction_top20` and `cycling_substate` (>= 50%),
   recorded in facts (`cycling_clusters`), markers NOT stripped (they are that sub-state's
   biology); `score_targets` -> `proliferation_marker` flag with penalty 0.5; the scripted critic
   and the critic prompt flag such candidates as "cell-cycle gene, not a disease-specific target".
4. The scripted analyst runs the whole-lineage contrast under scoring v3 and records the union of
   the sibling, whole-lineage and shared-cluster top lists; the nominator's pool and brief carry
   all three. Prompts v5 say the same to the agent.

**Regression protocol.** Before any HCC result is reported under v3, the scripted calibration is
re-run on ccRCC and CAA with `--scoring-version 3` and the ranks of CA9 and ERBB4 compared with the
v2 runs. All three results go in `docs/results_*.md`, whatever they show. The `config/ccrcc.yaml`
and the earlier runs keep `scoring_version: 2` as the record of that arm.
