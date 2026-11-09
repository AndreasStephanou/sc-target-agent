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

### 19b. HCC under scoring v3, and the pre-registered arm without patient 4

Regression (scripted calibration, `--scoring-version 3`): ccRCC CA9 rank 1, no unflagged decoy
(unchanged); CAA ERBB4 rank 4 (was 5), drugs 0.75, PASS. HCC (`runs/hcc_calibrate_4`): TOP2A no
longer nominated (proliferation flag), 20 nominations with zero decoys, SQLE / FDPS / IDI1
(cholesterol synthesis), NDUFA4L2, CYP17A1 and the tumour-endothelium set (FOLH1, PGF, ESM1,
THY1) on top; GPC3 still unscored.

The whole-lineage contrast did run as intended (14,366 tumour vs 1,471 normal hepatocyte-class
cells, 10 vs 7 replicates; the axis decision recorded "shared contrast held 10% of lineage
disease cells -> whole-lineage used"). In it GPC3 is up 6-fold (log2 2.63, baseMean 764) at padj
0.20: not significant. AKR1B10 (0.95) and SPINK1 (0.46) are flat. The reason is the one the
label check found before the first run: HCC04N, labelled non-tumour, expresses GPC3 (77%),
AKR1B10 (60%) and SPINK1 (91%) and contributes 22% of the reference hepatocyte-class cells,
while HCC04T is a GPC3-negative tumour. A replicate-aware test correctly refuses to call a gene
disease-specific when a "normal" replicate expresses it highly. This is a data-label result, not
a method defect, and it is exactly what the pre-run check predicted.

Pre-registered follow-up (section 18): one arm excluding both patient-4 libraries. To keep this a
disclosed config choice rather than a rebuilt dataset, `dataset.exclude_samples` and
`dataset.exclude_samples_reason` were added to the loader; the exclusion is recorded in the load
summary, the run facts and the report header. The primary (all-libraries) result stands and is
reported alongside.

### 19c. The arm without patient 4, and a marker-free label check

`runs/hcc_nop4_calibrate` (scoring v3, `dataset.exclude_samples: [HCC04N, HCC04T]`, identical key
and method): GPC3 is first in the tumour-vs-normal hepatocyte contrast (log2 6.45, padj 8e-12;
with the pair in: 2.63, padj 0.20), ranked 4th overall behind CA2, RAMP1 and PGF (none a decoy),
nominated by the scripted model, codrituzumab found (drug recall 1.0), zero decoys nominated,
TOP2A flagged. The state criterion still reads "missed": the recorded signature (cholesterol
synthesis SQLE / FDPS / IDI1, NDUFA4L2, CYP17A1 - consistent across patients) shares only GPC3
with the key's literature marker list (AFP, MDK, SPINK1, AKR1B10 ... are patient-heterogeneous
here and never significant). As with EPAS1, the key stays frozen and the lesson is recorded: a
key's state markers should be genes expected to be consistent across patients, not the textbook
list.

The exclusion was identified before the first run with a marker check that used the expected
target itself (GPC3), as CA9 was used for ccRCC. That is a weakness for a blind test and
unusable when no marker is known. A marker-free check was therefore added:
`check_reference_libraries` builds one pseudobulk per library within each lineage (log-CPM over
the 1,000 most variable genes), first asks whether the lineage shows any condition structure
(mean within-label correlation exceeds between-label by >= 0.05; otherwise the labels are not
testable there and it says so), and then flags any library whose nearest libraries mostly carry
the other condition and whose mean correlation with the other condition exceeds that with its own
by >= 0.05. No gene names. The scripted critic runs it once per run and records a major flag per
library; the critic prompt asks the agent to do the same. Tests: no flags on the clean synthetic
(no lineage-level structure, reported as such); a disease library relabelled as reference, in a
synthetic with an injected condition-wide programme, is flagged and the critic records it. The
HCC result of this check is reported in `docs/results_hcc.md`.

### 19d. Report figures (patch 37, report only)

Feature-plot UMAPs and dot plots were replaced by replicate-level figures, since the pipeline's
claims are about libraries: per-library pseudobulk expression panels (one point per library, by
lineage and recorded state, coloured by condition, boxes per condition; `plot_genes` now draws
these), a per-state bar of library contributions, the library-correlation heatmaps from
`check_reference_libraries`, and an evidence table for the nominated candidates at the top of the
ranking section (state and lineage effect sizes with the axis used, specificity, attribution,
compound counts with the inhibitor caveat, flags, critic notes). No analysis change; no version
bump. Reports regenerated from an existing run directory without the data object omit the figures.

## 20. results-v2

Tag `results-v2` marks the state after the HCC reviewer run: scoring v3, prompts v5, the
environment check, the counts-in-X guard, the scripted abort, `exclude_samples`,
`check_reference_libraries`, the GEO matrix builder, replicate-level report figures, and
deterministic pseudonyms (assigned over the full library set before any exclusion, so S-numbers
agree across arms of one dataset). `docs/results_hcc.md` holds all seven HCC calibration runs,
including the two that failed for infrastructure reasons; `docs/REPORT.md` section 5 is the
write-up. The HCC all-libraries result is a miss and is the primary result; the pre-registered
arm without patient 4 ranks GPC3 4th with codrituzumab. No agent (model-driven) run has been made
on HCC yet; that, independent testing by another fellow, and a non-tumour dataset are the next
tests.

## 21. Fourth dataset (ovarian cancer omental metastases, GSE286057) and the attribution representativeness guard

A second reviewer-style round: fresh clone, fresh environment, a dataset published in 2025 chosen
from GEO by search, and the first model-driven (Claude Code) runs since scoring v3. Two deposits were
screened and dropped on composition before any analysis (a 2025 gastric atlas, GSE270680, 47
libraries: 70% PTPRC+ and 1-2% EPCAM+ cells in both conditions, so no epithelial target can be
tested there; `inspect` now prints this and warns). The chosen deposit: HGSOC omental metastases,
6 tumour libraries (5 with substantial tumour content) vs 23 omentum libraries (6 distant omentum
from the same patients, 17 from 8 non-cancer donors); peritumoral libraries left out as a declared
second arm. Key frozen before any run: expected FOLR1 (mirvetuximab soravtansine), acceptable
MUC16, MSLN, TACSTD2, CLDN6, ERBB2, EPCAM, with the written prediction that MSLN would be
shared_with_reference (mesothelium) and that P14_distant (EPCAM 10%) and P12_distant (3%) were the
declared exclusion if the result hinged on them.

Runs: scripted calibration (state found 5/5 lenient; FOLR1 significant in the sibling pseudobulk,
log2 4.8, padj 5e-8, rank 39 of the significant genes, outside the fixed top-20 pool: unscored);
agent, all libraries (2.0M tokens, 77 calls: state found, strict 3/5, FOLR1 scored and ranked 16th
but set aside as shared_with_reference; the critic traced the cause to two reference libraries via
`check_reference_libraries` and asked for the exclusion arm by name - the same two libraries the
pre-run marker check had identified; no decoy nominated; agonist antibodies flagged for
direction); scripted no-leak arm (FOLR1 enters the pool, rank 96); agent no-leak arm (1.45M tokens:
state found and **strict 5/5 for the first time in an agent run**, 8 nominations, precision 0.25,
no decoys, FOLR1 **blocking** as reference_driven).

The blocking flag was traced on the raw object: 59 FOLR1+ reference cells in 94,235 (0.06%), 26 of
them from one healthy donor's lesser omentum, KRT19/CLDN3/WFDC2/EPCAM-positive and CALB2/UPK3B-
negative - a benign Müllerian-type epithelial inclusion, not mesothelium and not tumour spill-over,
expressing folate receptor alpha as normal Müllerian epithelium does. 151 such reference cells
were 7.5% of the 2,016-cell tumour cluster (above the 2% thin-pool guard) and expressed FOLR1 more
highly than the tumour cells, so the in-cluster rule read "the signal comes from reference cells".
The rule never asked how representative those cells were of the reference condition (0.16%); the
per-sample abundance test had already said the state is absent from reference libraries.

**Change (attribution rule; prompts v6; no scoring-version change):** `condition_attribution`
adds a representativeness guard: the in-cluster reference pool decides the class only if the
reference condition holds the cluster's LINEAGE at all - reference cells of that lineage must be
>= `analysis.attribution.min_reference_frac_of_condition` (default 0.5%) of all reference cells.
(A first, cell-relative version of the guard - in-cluster reference cells as a share of all
reference cells - was tried and withdrawn the same day: the rescore of the ovarian runs flipped
Treg-state genes from shared_with_reference to state_specific, because 240 reference Tregs are
0.25% of reference cells yet a representative counterpart; the reference holds T cells in bulk.
The lineage-relative criterion keeps the Treg verdicts and still frees FOLR1, whose lineage the
omentum reference barely holds.) Below that the status reads `rare_in_reference`, the class is decided by the other evidence,
and `rare_reference_counterpart` marks genes the rare population expresses; score_targets carries
it as an informational flag (no penalty), the scripted critic and the critic prompt report it as an
on-target consideration to identify. Tested (synthetic reconstruction of the case; the old
behaviour returns when the guard is set to 0). The HCC case is unaffected (normal hepatocytes in
the tumour cluster were ~3% of reference cells). Regression by `scripts/rescore.py
--recompute-attribution --full-candidates` on all saved runs.

Also from this round: `fetch_geo_10x.py` and `build_from_10x_samples.py` (GEO per-sample 10x
deposits; substring matching, with a warning that 'tumor' matches 'peritumoral'), the composition
line in `inspect`, 429 backoff with Retry-After in the knowledge clients (19 druggability failures
left Open Targets clinical fields empty), transient-error retry in the Claude Code backend (the
first agent run died on an OAuth token-refresh race before its first tool call), and a sentence in
the critic prompt that excluding a library is a human, config-level decision to request (the
critic worked this out unaided; it should not have to).


### 21b. Third agent run, and the decision to freeze

`runs/ovarian_agent_noleak_2` (prompts v6, lineage-relative guard in place; 2.3M tokens, 85 calls):
strict state recovery 5/5, four states recorded (tumour epithelium, Treg, pDC - then doubted by the
critic, CAF), five nominations, no decoys, FOLR1 again dropped as reference_driven. The guard did
not apply: in this clustering the derived epithelial lineage L1 held reference cells in bulk
(several reference libraries had 20-46% of their L1 cells in the tumour cluster), i.e. mesothelium
was absorbed into "epithelial" at resolution 0.2, so the reference "held the lineage" and the 151
benign cells were judged representative. A third rule (reference_driven requires that disease
cells in the cluster mostly do not express the gene) was drafted and not applied: three touches to
one rule in two days, each defeated by a different property of coarse derived lineages, is the
signal to stop and state the limitation. Rescore regression of the lineage-relative guard on the
earlier datasets: ccRCC CA9 1 -> 1; CAA de novo ERBB4 1 -> 3 within state (rank 7 overall;
state_specific -> shared_with_reference under recomputed attribution with full candidates - noted,
not acted on); no decoy unflagged anywhere.

## 22. results-v3

Tag `results-v3`: four datasets, the attribution representativeness guard (lineage-relative),
prompts v6, the GEO tools, composition check, 429 backoff and transient-retry. `docs/results_ovarian.md`
holds all seven ovarian runs including the two rescores; `docs/REPORT.md` section 6 is the write-up.
Method frozen here pending independent use by another CoPAIR fellow.

## 23. Validation round 0.4: the expression route must run on the base clustering

A new round began with the 0.4 skill (design classification, pre-run design and label checks,
no-key protocol), aimed at the two design types never run: non-malignant disease (IPF, GSE135893)
and perturbation (doxorubicin vs vehicle, GSE285475). Its first free step, the scripted calibration
on the program-mode synthetic dataset (`make_synthetic --mode program`, 400 cells/sample, seed 0),
ended with "No disease-specific state found" on data where every CAA keratinocyte runs a ~2.5x
program on 8 genes.

Cause (run.jsonl): with composition flat, the scripted analyst subclustered, then re-clustered
globally at resolution 2.0, and then called `disease_programs` on that clustering
(`sct_leiden_r2`, 15 clusters of 95-466 cells). The keratinocyte lineage was split into several
pieces; each piece showed 3 strong genes (SNAP25, GAP43, NEFL), below the absolute (30) and
relative (5, 3x any other cluster) criteria. The unit test (`tests/test_program_route.py`, 200
cells/sample, seed 2) passed because its route never reached the re-clustering step. The code
comment in the other branch already stated the rule ("on the base clustering, so subclusters do
not fragment the test"); it was applied in one branch and not the other.

**Change (tool default + prompts v7; no scoring change):** `integrate_and_cluster` records
`base_cluster_key` (the first clustering established in the session; preserved across subcluster,
merge and re-clustering). `disease_programs` defaults to it when `cluster_key` is not given and
reports `base_clustering_used`. The scripted route passes no key in either branch. The analyst
prompt says to call the tool without a key and why. Regression test added
(`test_disease_programs_defaults_to_base_clustering_after_recluster`). Generic: a program shift is
lineage-wide by definition; the right partition for testing it is the coarsest one, whatever the
dataset.

**Also:** `sctarget.py calibrate` now defaults to scoring version 3 like `run` (the calibration
above was silently scored under v1 because the config did not name a version), so the scripted
floor and the agent are scored alike.

Why this matters beyond the synthetic case: tumour datasets create new states and take the
composition route; a perturbation or a non-malignant disease shifts programs inside existing
lineages and takes exactly the path that was broken. The 0.4 round would have reported two false
negatives.

### 23b. Same-lineage clusters must not compete in the relative programme criterion

With the base clustering restored, the same calibration still failed: keratinocyte clusters 1 and
0 showed 9 and 6 strong genes (the whole planted set), but the relative criterion ("at least
small_program_min strong genes AND at least 3x any other cluster") compared each against the other
keratinocyte cluster, so 9 < 3 x 6 and neither was flagged. A programme shared by two clusters of
one lineage is evidence for a lineage-wide shift, not against it.

**Change (tool only):** `disease_programs` compares a cluster's strong-gene count against clusters
of OTHER lineages only (dominant lineage from the lineage key); same-lineage clusters are
excluded from the "any other cluster" maximum. The absolute criterion is unchanged. Regression test
extended (the flagged clusters must all be keratinocyte). Generic: whatever the tissue, the
lineage that runs the programme may be split over several clusters at the base resolution.

### 23c. `inspect` composition line was human-only; epithelial warning made conditional

On the mouse doxorubicin object `inspect` printed an empty compartment-composition line: the
canonical markers (PTPRC, EPCAM, COL1A1, PECAM1 ...) were matched exactly against Title-case mouse
symbols. Matching is now case-insensitive (`new_dataset.py`), and the "fewer than 5% epithelial
cells" warning says that it matters only when the question concerns epithelium; heart, brain and
blood have none by design. On the DOX object the line now reads endothelial 0.64 vehicle vs 0.46
treated and fibroblast 0.41 vs 0.48, a composition shift consistent with the deposit's own finding
of endothelial toxicity, visible before any run.

Also in this round (procedural, recorded for the skill): covariates that exist only in the deposit's
metadata (here the processing batch, two libraries aligned with STARsolo and six with Cell Ranger)
must be added to the sample sheet before building so they reach `obs`; `build_from_10x_samples.py`
already carries every sheet column. With `batch` in `obs` the design check shows it balanced across
arms (Cramer's V 0.00); the batch-only arm is pre-registered in `config/dox_batch2.yaml`.

### 23d. Programme-basis states: attribution, consistency and nomination must follow the basis

With 23/23b in place the synthetic programme state was found (7 of 8 planted genes) but ERBB4, the
planted target with the largest within-cluster effect (log2FC 1.65) and full druggability, ranked
7th behind six undruggable programme genes, and SNAP25 was nominated twice. Three causes:

1. `condition_attribution` classed ERBB4 `shared_with_reference` because reference cells of the
   cluster express it at >= 0.5x the disease cells. In a programme state the cluster exists in both
   conditions by definition, so that rule penalises every programme gene with a non-zero baseline;
   SNAP25 escaped only because its baseline is lower. **Change:** for a programme-basis cluster the
   class is read from the replicate-aware within-cluster contrast that found the programme
   (`disease_wide` if padj < threshold and log2FC >= min_log2fc; `shared_with_reference` = not
   induced); `reference_driven`, `weak` and `ambient_like` are unchanged. Composition states are
   unaffected.
2. `gene_specificity`'s `disease_sample_consistency` is "detected in >= 10% of the cluster's cells
   in that disease sample", a composition-state notion; and `score_targets` flagged `cons < 0.67`,
   which with three disease samples flags exactly 2 of 3 (0.667). **Change:** programme states use
   per-sample induction (disease sample's within-cluster mean above the median of the reference
   samples' means); the cut-off is "fewer than two thirds" exactly.
3. The scripted nominator recorded the same gene once per programme cluster of the same lineage.
   **Change:** a gene is nominated once across states; the evaluator counts unique nominations.

All three are generic: any dataset whose signal is a programme shift in an existing lineage
(perturbation, most non-malignant disease) goes through these rules. Tumour (composition) paths are
untouched; regression on the saved tumour runs is by `scripts/rescore.py --recompute-attribution`.

### 23e. The scripted floor characterises the top enriched cluster of every lineage

IPF floor (`runs/ipf_calibrate_1`): six clusters were disease-enriched across three lineages
(basal, MUC5B+ and ciliated epithelium; an IPF-only AT2 cluster; two macrophage clusters). The
scripted analyst characterised the top three by effect x size, all epithelial, so the macrophage
states (the SPP1+ expansion that two independent IPF studies report) never received a marker
table and the evaluator could not even match them. **Change (scripted route only):** the top
enriched cluster of each lineage is characterised before a second cluster of the same lineage, cap
raised from 3 to 5 composition states; the nominator and critic cover up to 6 recorded states
instead of 4. Tumour datasets, where the enriched clusters are one lineage, are unaffected in
ordering. The agent prompt already leaves the order to the model; nothing changes there.

Floor limitation stated rather than patched: the scripted route never subclusters an enriched or
programme cluster, so a rare state inside one (the 374 KRT5-/KRT17+ aberrant basaloid cells inside
the 5,559-cell transitional cluster) is invisible to it. The agent has `subcluster`.

### 23f. The candidate pool keeps the largest significant effects, not only the most significant statistics

DOX floor (`runs/dox_calibrate_1`, 4 vs 4 wild-type hearts, `objective: marker`): composition flat
as predicted; the expression route found one programme, in a cardiomyocyte cluster (377 genes up,
91 strong), whose within-cluster contrast carries the canonical cardiomyocyte doxorubicin response
(Trim63 padj 7e-6, Nppb, Fbxo32, Nppa - all acceptable genes of the frozen key). None reached the
candidate pool: `top_up` stored the 50 most significant genes by test statistic, the scripted route
took 40 of them, and the atrophy and natriuretic genes ranked 137th-650th by statistic because
abundant structural genes with small fold changes dominate the statistic in a 3-vs-3 pseudobulk.
**Change:** `top_up` = the 40 most significant genes plus up to 20 of the largest significant
effects (log2FC >= 1) not already listed, each row labelled `pool = stat | effect`; the scripted
candidate set takes up to 60; the nominator's brief shows the 30 most significant plus every
effect-pool gene. Generic; the tumour runs are unaffected except for a slightly larger pool.

Key-drafting lesson (recorded, key unchanged): the expected target Cdkn1a (p21) is not induced in
this dataset at two weeks (down in the cardiomyocyte programme, log2 -0.8; modest and
non-significant elsewhere), as the pre-run label check had hinted. The p53 read-out is a
time-point-dependent response; the lineage-level atrophy / natriuretic response is what the data
hold. This joins EPAS1 (protein-level target) and the HCC literature markers as the third
key-drafting failure mode: a read-out that is real but not at the sampled time.

**23f addendum.** The DOX floor re-run under 23f was unchanged, because `disease_programs` rebuilt
the programme's gene list from the 50 strongest genes by statistic and then replaced the stored
pool with the genes of its "clean programme", discarding the effect-pool genes again. Now the
programme's list follows the same pool rule (40 by statistic + 20 largest strong effects) and the
ambient step removes suspect genes from the stored pool instead of restricting it to the
programme head. Re-run as `dox_calibrate_3`.

### 23g. Target rank is judged within the state that carries it

ccRCC regression run under prompts v8 (`runs/ccrcc_agent_v04_1`): tumour state found (lenient 9/5,
strict 8/5), CA9 ranked 1st in both tumour states (scores 0.91 and 0.92, ahead of EGLN3), no decoy
nominated, EGLN3's direction flagged again, CD70 and FOLH1 again above the key - yet the evaluator
printed "CA9: 6". The merged ranking of all states put the endothelium and stroma programmes' genes
(FOLH1, VWF, ACE at 1.0) above the tumour state's top gene, because v8 now records the top enriched
cluster of every lineage. Every key is worded "ranked first in its state"; the merged rank
penalises exactly the lineage coverage the method should have.

**Change (evaluation definition; no analysis change):** `target_rank` and `acceptable_rank` are
the best rank within a recorded state's own ranking, with the state named (`target_state`);
the merged rank is kept as `target_rank_global` and reported, not scored. Genes "above the first
true target" are counted in the state where it ranks best. For runs with one recorded state the
numbers are unchanged; for earlier multi-state runs the within-state rank can only be equal or
better than the published global one. All results tables of this round are regenerated under the
new definition and say so in their header.

### 23h. The floor's effect pool is taken over the whole stored candidate list (round 0.4b)

Motivating failure: ovarian scripted floor on the rebuilt object (`runs/ovarian_calibrate_v04_1`).
The tumour cluster was found (5/5) and FOLR1 was significant in its sibling contrast (log2FC 4.8,
padj 5e-8), the 39th most significant gene and the 9th largest effect. `_candidate_pool` (23f)
stores the 40 most significant rows as `pool = 'stat'` and then up to 20 of the largest remaining
effects as `pool = 'effect'`; the scripted target role read the first 20 `stat` rows plus every
`effect` row. A gene at significance rank 21-40 is therefore in neither set, whatever its effect,
while a weaker effect outside the head is kept. FOLR1 was again unscored by the floor, as in round
0.3 (then: outside a fixed top-20).

Change (scripted floor only; agent backends read the full stored list): the floor's pool is the 20
most significant rows plus the 20 largest significant effects (log2FC >= 1) among **all other
stored rows**, whether the tool labelled them `stat` or `effect`. Tested
(`tests/test_candidate_pool_boundary.py`). The stored table and its `pool` column are unchanged, so
the agent's view and the report are unaffected. Verified by re-running the ovarian floor
(`runs/ovarian_calibrate_v04_2`); the result is recorded in `validation_0.4/NOTES.md` whatever it is.

### 23i. The evaluator names the lineage of the state that carries a target (round 0.4b)

Motivating case: HCC agent run on the rebuilt object, all libraries (`runs/hcc_agent_v04_1`). The
within-state rule (23g) reported GPC3 2nd and recall@5 = 1.0. The state that carried it was
cluster 12, a CAF-like stromal state (83.7% tumour cells, fibroblast-dominant); GPC3 reached its
ranking through the disease-wide contrast and was classed `disease_wide`. In the two tumour
hepatocyte states the key means (clusters 14 and 15) GPC3 was unscored, for the reason recorded
in round 0.2 (lineage contrast padj 0.20 with HCC04N in the reference). The number was right by
the rule and wrong as a reading: the expected target was not found where the key expects it.

Change (evaluator only; no method change, scoring v3 unchanged): the eval line prints the dominant
lineage of the carrying state next to its cluster (`state: {'GPC3': '12 (Fibroblast)'}`,
`target_state_lineage` in the JSON). An expected state in a key may carry an optional `lineage:`
field; when present, a rank in a state of another lineage is reported under
`target_rank_off_lineage` and not credited to `target_rank` or recall@k. Existing frozen keys have
no such field and keep their numbers; the HCC result above is read with the lineage shown and
recorded as "scored in a stromal state, unscored in the tumour hepatocyte states". New keys
should name the lineage. Tested (`tests/test_eval_lineage.py`).

### 23j. Prompts v9: a reference-driven candidate with per-sample support is reported as conditional, not dropped (round 0.4b)

Motivating case: ovarian no-leak agent run (`runs/ovarian_noleak_agent_v04_1`, prompts v8). State
found strict 5/5, FOLR1 17th within it, the critic named the 151 in-cluster reference cells as an
"epithelial/mesothelial-like population ... e.g. mesothelial or epithelial inclusions" and still
recorded a blocking flag, "not a safe tumour-state marker or target in this dataset". That is the
fourth agent run on this dataset to find FOLR1 in the right state and set it aside; the gene has
an approved agent (mirvetuximab soravtansine). Protocol 21b froze the attribution *rule* after three
touches and stated the limitation; the 0.4 skill then added a *reading* rule (a blocking
`reference_driven` verdict on a candidate that the per-sample and compartment evidence support is a
prompt to identify the reference cells; the candidate is reported as conditional on them, not
dropped). The repository's prompts still said "drop `reference_driven` genes outright" (target)
and "a blocking flag" (critic), so the method contradicted the skill that drives it.

Change (prompts v9; the attribution tool, its rule and scoring v3 are unchanged): the target prompt
keeps a `reference_driven` gene ranked with its flag when the state is absent from reference
libraries by per-sample abundance and the gene is specific and consistent across disease samples,
and reports it as conditional on the identity of the in-cluster reference cells; it drops the gene
when those cells are a representative counterpart of the lineage. The critic prompt asks for the
cells to be named (lineage, libraries, share of the reference condition) before any flag, records a
major flag that states the condition and what would settle it (reference-only subcluster, histology,
exclusion arm) in the conditional case, and keeps the blocking flag for the representative case.
Verified by re-running the no-leak arm under v9 (`runs/ovarian_noleak_agent_v04_2`); the outcome is
recorded in `validation_0.4/NOTES.md` whatever it is. The rule is general: it applies to any design
in which the reference holds a small population expressing the candidate (benign inclusions, rare
progenitors, a tissue-resident counterpart of a disease state).

### 23k. `subcluster(condition=...)`: a diagnostic split of one condition's cells inside a cluster (round 0.4b)

Motivating case: the v9 ovarian no-leak run (23j). The critic wrote that the `reference_driven`
verdict on FOLR1 "depends on reference cell identity", that the 151 cells "may be normal
fallopian-tube-like or mesothelial epithelium, or tumour-contaminated", and that "a reference-only
subcluster of cluster 24 would settle this" - and had no tool to do it: `subcluster` split all
cells of a cluster and switched the session's cluster key, which is not a move a critic should
make. The same question arises in every design where the reference holds a small population
expressing the candidate.

Change (tool; no method or scoring change): `subcluster` takes an optional `condition`. With it,
only the cells of that condition inside the chosen clusters are re-clustered (labels
`<cluster>.<condition>.<k>` in a new key `<key>_sub_<condition>`), the session cluster key and
`cluster_assignments.csv` are left unchanged, and the result tells the caller to run
`find_markers` on the new key. Available to the analyst and the critic (the tool's existing role
set). Tested (`tests/test_subcluster_condition.py`). The critic prompt (v10) names it as the
settling move for a representative-but-doubtful `reference_driven` verdict. Not re-run on the
ovarian arm in this round: the fix is recorded, its first use will be reported when it happens.

### 23l. Prompts v11: nominated genes from a non-disease-cell lineage carry a microenvironment flag (round 0.4b)

Motivating failure: FLC agent run 1 (`runs/flc_agent_1`, prompts v8). The tumour state was found
by the expression route and PRKACA ranked 2nd in it behind CA12 with ERBB2 2nd (neratinib found),
but PDGFRB, nominated from the CAF-like fibroblast programme (cluster 8), carried no critic flag -
the only decoy nominated unflagged in the round. The FLC key lists PDGFRB and COL1A1 as decoys
because the fibrous lamellae are the tumour's scaffold, not the tumour cell. The same pattern ran
unlabelled through earlier rounds (FOLH1 and ESM1 from tumour endothelium in HCC and ccRCC, CD70
and TREM2 from TAMs), where every key lists KDR, FLT1, PDGFRB, COL1A1 or CSF1R as decoys for that
reason. The critic prompt asked about proliferation and attribution, not about compartment.

Change (critic prompt only, v11; scoring unchanged): every nominated gene whose state is not of the
disease cells' lineage (stroma, vessels, immune in a tumour design) carries a flag naming the
lineage and saying it is a microenvironment candidate, not a target on the disease cells; the
nomination may stand as a microenvironment hypothesis. The scripted critic is not changed: it
has no notion of the disease cells' lineage (the config does not name it), and inventing one
from the data would be a second rule to validate. Verified by a second FLC agent run
(`runs/flc_agent_2`, also the dataset's first stability pair); the outcome is recorded in
`validation_0.4/NOTES.md` whatever it is.

### 23m. Prompts v12: a stage briefing never drops an item (round 0.4c, CAA)

Motivating failure: the CAA agent run (`runs/caa_agent_v04_1`, prompts v11). The critic's summary
ends "the candidate list I received was truncated after ERBB4, so any further candidates were not
reviewed". Its briefing was one JSON dump of the recorded states, the nominations, the top-8 ranked
rows per state and the abundance table, cut at a flat 9,000 characters; the claims ran 31,715.
Re-measured after the fact on every agent run of rounds 0.4 and 0.4b (`state.json` +
`config_used.json` through `stage_prompt`): the critic's claims ran 23-52k characters against the
cap in all of them, so in every run the critic was briefed on the first six or seven nominations
and had to find the rest, if at all, through its own tool calls. The nominator's briefing had the
same cap per section (6,000) and its state list exceeded it in the ovarian runs (seven states,
7.2-7.7k). Which nominations each run's critic never received is listed in
`validation_0.4/NOTES.md` (round 0.4c); FLC run 1's PDGFRB, the decoy that motivated 23l, is
among them, so part of what 23l attributed to the critic's reasoning was a briefing it did not have.

Change (prompt construction only, v12; tools, scoring and the scripted critic unchanged - the
scripted critic reads the facts directly): a briefing may shorten what it says about an item, never
drop one. An index line names every nomination and every recorded state (gene, cluster) and is
never shortened; free text per item (`detail`) is clipped progressively (900 -> 450 -> ... -> 120
characters) until the list fits its budget; gene lists are shortened evenly across keys, never
dropped by key; the critic is told that every indexed gene must appear in its review. Budgets are
set so that no earlier run's briefing is truncated (nominator 17-27k characters, critic 17-35k;
5-9k tokens, against 1.6-3.5M input tokens per run). Tests: `tests/test_prompt_briefing.py`.

The metric "decoys nominated unflagged" in rounds 0.3-0.4b measured the critic with an incomplete
briefing; its values stand as recorded (they are what the method did) but are not a measure of
the critic's rule-following for nominations past the cap. Verified prospectively by the CAA agent
run 2 (`runs/caa_agent_v04_2`, v12) and the CAA exclusion-arm agent run; outcomes in the notes.

### 23n. `map_orthologs`: an empty precomputed entry is "unknown", not "no ortholog" (round 0.4c)

Motivating failure: in the CAA agent run the ortholog tool returned `no_ortholog` for GABRB1 and
GPC5 because `var[human_symbol]` was an empty string for them, while the druggability lookup
used the same symbols as human and returned 92 interactions for GABRB1; the critic flagged the
inconsistency as "druggability unverified" (a correct reading of what it was shown). A biomaRt
table leaves a gene out when it has no human ortholog and when the query simply did not cover it;
the empty cell cannot tell the two apart. Change: an empty entry falls through to Ensembl like a
gene absent from the table (live: `one2one` or a true `no_ortholog`; offline: `assumed_same_symbol`,
flagged), and the result carries `var_column_empty: true` and names both sources. Test updated in
`tests/test_knowledge_parsing.py`.

Verified in CAA agent run 2 (v12): GPC5 now resolves through Ensembl (`ortholog_one2one`,
ENSG00000179399). GABRB1 still came back `no_ortholog`, this time from Ensembl itself: the dog
symbol resolves to the current-assembly gene ENSCAFG00845008573, for which Compara lists no
human homologue; the dataset's own identifier (ENSCAFG00000001889, CanFam3.1) is no longer served.
That is an annotation gap in Compara, not biology (GABRB1 is a one-to-one orthologue by every
other resource), and it will recur for any non-model species whose assembly moved. Second part of
the change: when Ensembl reports no human homologue, the client asks whether the query symbol is
itself a human gene (`lookup/symbol/homo_sapiens/<symbol>`); if it is, the mapping returns the
symbol as `same_symbol_human_gene` with `symbol_unverified: true`, so the druggability lookup
proceeds on a named assumption that the critic can see, and `assumed_or_unverified` lists the
gene. A true absence (no human gene of that name) stays `no_ortholog`.

### 23o. `run_eval` rewrites the report's evaluation section (round 0.4c)

The report is built at the end of the run, before `run_eval` scores it, so every scored run's
`report.md` carried the "no frozen answer key was scored" placeholder in section 6 while the
result sat in `eval.json` and the log. `run_eval` now rewrites section 6 of `report.md` and
`report.html` from `eval.json` (`sctarget.report.builder.refresh_eval_section`). Cosmetic, but
the report is the artefact a reader opens. Tests: `tests/test_report_eval_refresh.py`.

### 23p. A condition-exclusive lineage makes the within-lineage abundance test void, not failed (round 0.4c, de novo arm)

Motivating case: the CAA de novo floor (`runs/caa_denovo_calibrate_v04_1`; Harmony on the library,
lineages derived, nothing of the paper's processing reused). Harmony did not absorb the CAA state:
the 59k CAA keratinocytes are distinct enough from OSCC and gingival keratinocytes to form their
own coarse group, derived lineage L4 (99% disease cells, 645 reference cells), with the shared
epithelium in L0 and L1. The scripted critic re-ran `differential_abundance(within_celltype=True)`
as it must and flagged the two largest CAA clusters (16: 31k cells; 13) as "enrichment does not
hold within its parent cell type" (major), because within L4 the per-sample proportion of one CAA
sub-state against the others shows no shift. The flag is a misreading: the lineage is itself the
disease-associated population, so its composition test is the all-cells one; the within-lineage
test compares disease sub-states with each other and is void. Two consequences follow for the
same reason and were visible in the run: the sibling contrast inside L4 cancels every pan-CAA
gene (ERBB4 unranked in cluster 17, 2nd in 16, 1st in cluster 3 - the small CAA cluster inside the
shared lineage L0, where the siblings are OSCC and gingival keratinocytes); and attribution inside
L4 rests on a 645-cell reference pool. The same pattern will occur in any tumour-vs-other-tissue
design run with derived lineages when the tumour cells are transcriptionally far from their tissue
of origin (the HCC "pan-tumour gene cancels among tumour siblings" case, protocol 19/21, seen from
the lineage side).

Change: (i) scripted critic: before emitting "does not hold", it sums the disease fraction over
the clusters of the state's lineage (`cluster_composition`); at >= 0.9 it records one minor flag
per lineage, "Lineage L is condition-exclusive (x% disease cells): the within-lineage test is
void for its clusters", whose detail names the two consequences (read the disease-wide axis for
pan-state genes; thin reference pool for attribution). Below 0.9 the major flag stands unchanged.
(ii) Critic prompt (v13): the same check and the same wording rule. Tests:
`tests/test_critic_exclusive_lineage.py` (3). The de novo agent run 1 ran under v12; a second de
novo agent run under v13 serves as verification and as the arm's stability pair.

### 23q. Report figures draw a state under the label it has at report time, so a merged state is drawn as one population (round 0.4c; report only)

Motivating case: CAA de novo agent run 1 (`runs/caa_denovo_agent_v04_1`). The analyst recorded
states 16, 17 and 13; the critic then ran `merge_clusters(16, 17, 13)`, which wrote the label
`13+16+17` into a new active cluster column. Every cell-level report figure looked each recorded
state up by its recorded label in the active column, found no cell, and drew "state 16: 0 cells"
(state-share panels empty, no shaded state group in the per-library expression panels). FLC agent
run 1 (state 8, merged into 8+19 after the record) shows the same. Runs that recorded the state
under the merged label itself (ccRCC, IPF, FLC run 2) were not affected. Independently, the
clustering UMAP is drawn when the clusters are computed, so a merge or a subcluster step made
afterwards never appeared on any figure: the report could not show the merged clusters together.
The analysis tools were never affected (`_resolve_key` finds a recorded label in the base
clustering column); this is report code only, and no number in any report changes.

Change (report code only; `scripts/rebuild_report.py` applies it to an existing run):
(i) `figures.resolve_states` maps every recorded state to cells under the active key: the label
itself if present; else the merged label that lists it as a piece, with all recorded pieces of one
merge collapsed into one drawn state (title "state 13+16+17 ... recorded as states 16, 17, 13,
merged by merge_clusters"); else the label in another clustering column (the base clustering after
a subcluster step). `_groups` (expression panels) and `state_library_bars` use it; a candidate
panel for a merged piece says "drawn as merged state ...". (ii) New report-time figure `state_map`:
the recorded states on the embedding under the final cluster key, grey elsewhere, one colour and a
centroid label per resolved state. It is the report's one cell-level figure and a location map
only; the replicate-level figures carry the claims. (iii) `integrate_and_cluster` writes
`umap_coordinates.csv` (per cell, run directory only, never sent to the model) and `_resume_session`
reloads it, so a rebuilt report can draw the map when the dataset file has no embedding. Runs
made before this change have no coordinates file; their rebuilt reports get (i) and skip the map
unless the dataset carries `X_umap`. (iv) Found while rebuilding: a reloaded session re-derived
the data-derived lineages, which renumbered them (the rebuilt CAA de novo panels said L14, L8, L9
where the run and its report text say L1, L5, L6). `derive_lineages` assigns every fine cluster
wholly to one lineage and records the map in `facts.lineages.fine_to_lineage`, so the column is a
function of the cluster labels; `_resume_session` now restores it from that record
(`lineage_column_from_facts`; a merged or subcluster label through its first piece / parent)
instead of re-deriving. (v) Also found while rebuilding: the reloaded dataset file brings back the
cells the run's QC removed; they have no cluster assignment and formed a spurious "unassigned"
lineage group in the rebuilt FLC panels. `_resume_session` now drops cells without an assignment
and logs the count, so a rebuilt or resumed session sees exactly the run's cells. Tests:
`tests/test_report_merged_states.py` (6: pieces resolve to one drawn state with the merged cell
count; groups, shares, map and expression panel draw it; unmerged states unchanged; a state
recorded on the base clustering after a subcluster step; lineage column restored from the record;
a session resumed from a scripted run has the run's cells, labels, lineages and embedding).
Reports of the two affected runs rebuilt in place (section 6 refreshed by `run_eval`); no number in
either changed.

## 24. results-v4

Tag `results-v4`: round 0.4 (protocol 23-23g, prompts v7-v8, scoring unchanged at v3): two new
design types run blind from frozen keys (IPF GSE135893, non-malignant; doxorubicin heart GSE285475
wild-type arms, perturbation in `marker` mode), the expression route repaired in four places,
the candidate pool rule, lineage coverage in the floor and the prompt, within-state target rank in
the evaluator, the 0.4 report (executive summary, publication figures, PDF/SVG), `rebuild_report.py`,
`build_from_mtx.py`, and the 0.4.1 skill package in `skills/` and `.claude/skills/`. Results:
`docs/results_ipf_v04.md`, `docs/results_dox_v04.md`, `docs/results_ccrcc_v04.md`; write-up
`docs/REPORT.md` section 7; running notes `validation_0.4/NOTES.md`. 89 tests. Nothing pushed;
this is a working copy awaiting the author's review of each method change.
