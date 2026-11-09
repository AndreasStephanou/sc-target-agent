# Can an agent rediscover a known drug target from single-cell data, blind?

*sc-target-agent: validation on four datasets. Draft for the CoPAIR fellowship report. Numbers in
the tables are regenerated from the run directories by `scripts/collect_results.py`
(`docs/results_caa.md`, `docs/results_ccrcc.md`, `docs/results_hcc.md`, `docs/results_ovarian.md`); the method history is `docs/EVALUATION_PROTOCOL.md`.*

## 1. Question

Single-cell atlases of disease routinely end with a list of "potential targets" chosen by the
analysts. We asked whether an AI agent, given only the count matrix and the condition labels,
would recover a target that was independently validated, rank it above plausible decoys, do so
consistently, and do it better than a fixed script using the same tools. The known answers were
ERBB4 in canine acanthomatous ameloblastoma (CAA; our own snRNA-seq, neratinib sensitivity shown
in vitro), the HIF programme of clear cell renal cell carcinoma (ccRCC; public scRNA-seq,
CA9 and the belzutifan target HIF-2alpha), GPC3 in hepatocellular carcinoma (HCC; public
scRNA-seq, codrituzumab / CAR-T target), and FOLR1 in ovarian cancer omental metastases (public
scRNA-seq, 2025; mirvetuximab soravtansine). The third and fourth were run as a reviewer would run
them: fresh clone, fresh environment, README and skill only; the fourth also carried the first
model-driven runs since scoring v3.

## 2. Design

Four roles (planner, analyst, nominator, critic), each a fresh model conversation, act only
through ~20 deterministic tools (per-sample differential abundance, sample-blocked pseudobulk
DESeq2 against sibling clusters, a disease-vs-reference compartment contrast, specificity,
condition attribution with ambient and thin-pool guards, doublet and dissociation-stress flags,
DGIdb / Open Targets / PubMed lookups, a weighted deterministic score). State passes between
roles only as structured findings. Every tool call and the exact payload the model saw are
logged; no per-cell data or real sample identifiers reach the model.

Evaluation uses an answer key frozen and committed before the first model run: expected state
markers, expected and acceptable targets, expected drugs, decoys. Metrics: state recovery
(lenient: the cluster's tables carry the signature; strict: the agent itself recorded it),
target rank, decoys above the target and decoys nominated unflagged, drug recall, nominated
precision, stability across runs, cost. Comparison arms: a naive baseline script and the tool
chain under a fixed scripted procedure (also our free calibration tool). The analysis method was
frozen at a tag before any model run; every later change is a labelled version (prompts v1-v4,
scoring v1-v2, attribution rules) recorded in each run and motivated by a named failure.

## 3. CAA (canine, snRNA-seq, 205k nuclei, 6 CAA / 10 OSCC / 3 gingiva libraries)

| arm | state (lenient / strict) | ERBB4 rank | decoys | drugs | n |
|---|---|---|---|---|---|
| baseline script | wrong cluster / no | none | 2 unflagged | 0 | 1 |
| scripted tools | yes / no | 5 | none | 0.75 | 1 |
| agent, prompts v1 | yes / 0 genes recorded | none | none | 0 | 2 |
| agent, prompts v2 | yes / no | 2, 2, 2 | FGFR2, flagged 3/3 | 0.75 | 3 |
| agent, v2 + scoring v2 | yes / no | 1 (x4) | FGFR2, flagged | 0.75 | 3 rescored + 1 |
| de novo clustering, inherited annotation | yes / no | none, then 4 | none | 0, then 0.75 | 2 |
| de novo clustering, derived lineages | yes / no | 1 | none | 0.75 | 1 |

**The hand-off decided the outcome.** In every agent run the tools reached the right cluster.
Under v1 the analyst recorded the state without genes and the nominator never saw the marker
table where ERBB4 sat; two interface rules (a recorded state must carry genes; the nominator sees
both tables) took the result from 0/2 to 3/3 with the analysis unchanged.

**FGFR2.** Listed as a decoy (expressed in OSCC), it ranked first in every v2 run and the critic
flagged it each time. The later lineage axis showed it state-specific in this data, and FGFR2 is
recurrently mutated in human ameloblastoma; the key stays frozen and this is reported as an open
hypothesis.

**ERBB4 is disease-wide.** Modestly up against sibling CAA clusters, ~8-fold up against reference
epithelium within shared clusters: elevated across CAA epithelium, not only in the neuroepithelial
state. Under scoring v2 it ranked first in all runs.

**De novo.** Recomputing the embedding (after fixing a silent Harmony failure and making such
failures loud) found the state; ERBB4 was first lost to an attribution rule that let 150 reference
cells (0.5% of a 31k-cell cluster) overrule the rest, then recovered under a thin-pool guard on
the identical partition. With the lineage annotation also removed, the agent ranked ERBB4 first
through a smaller CAA population inside the shared epithelium.

## 4. ccRCC (human, dissociated scRNA-seq, GSE159115, 7 tumour / 6 benign libraries, 29k cells)

| arm | state (lenient / strict) | EPAS1 | CA9 (acceptable) | decoys | n |
|---|---|---|---|---|---|
| scripted, calibration passes 1-3 | no -> yes / no | none | - -> 1 | none | 3 |
| scripted, calibration pass 4 | yes / no | none | 1 in every tumour state | none | 1 |
| agent, authors' annotation | yes / **yes** (8/5) | none | 1 | none | 1 |
| agent, derived lineages | yes / **yes** (7/5) | none | 1 | none | 1 |

**Public metadata was wrong.** GEO tissue labels were swapped for 4 of 13 libraries; a CA9 check
per library confirmed the authors' annotation files and the builder uses them.

**Four generic defects surfaced by calibration**, none visible on nuclear data: fixed QC removed
61% of benign epithelium and ~99% of the distal nephron (now adaptive per sample with a
lineage-loss report and a preview mode); malignant cells split by patient and each piece failed
the per-sample test (now `split_state_candidates` and `merge_clusters`); a heat-shock
dissociation cluster topped the ranking (now a stress flag that strips stress genes and rejects
only if little remains); doublet detection failed silently on a missing dependency while three
clusters were plainly doublets (now flagged, with a fallback). Then the sibling contrast was flat
for every tumour-defining gene (CA9 log2FC -0.06) because the siblings were other tumour
sub-states; the compartment contrast carried them (+12) and now feeds the signature and the
candidate pool.

**The strict metric passed for the first time.** The agent recorded the tumour state with the
HIF programme as its signature (NDUFA4L2, CA9, ANGPTL4, EGLN3, HILPDA...), merged the split
tumour on its own, ranked CA9 first and EGLN3 second, nominated no decoys, with or without the
authors' annotation.

**EPAS1 was missed, and could not have been found.** Its transcript is up 2.5-fold (padj 0.09,
baseMean 77) while its targets are up thousands-fold; HIF-2alpha is activated post-
translationally. An expression-based method nominates transcriptionally distinctive genes. This
was a drafting error in the key, reported rather than corrected: an expected target must be
transcriptionally distinctive, or the test is of the key.

**The critic.** It refused a tumour/T-cell doublet cluster, named a patient-driven split,
rejected a programme in mixed cells, and twice flagged that EGLN3's approved drugs are PHD
inhibitors that stabilise HIF - the wrong therapeutic direction, a judgement the scoring cannot
make. Among "novel" candidates above the key, CD70 and FOLH1 are clinical-stage ccRCC targets the
key had omitted.

## 5. HCC (human, dissociated scRNA-seq, GSE149614, 10 tumour / 8 non-tumour libraries, 63k cells)

Run as a pre-publication review: the repository re-cloned from GitHub into a fresh directory and
conda environment, following only the README and the skill, on a tissue the method had never seen.
The key was written from the literature before any run (expected GPC3, transcriptionally
distinctive; acceptable AFP, MDK; decoys spanning proliferation, hepatocyte lineage, endothelium,
cholangiocyte, EGFR). The pre-run label check (GPC3 per library, the skill's step 2) agreed with
the deposited labels for 17 of 18 libraries; HCC04N, labelled non-tumour, carried GPC3 in 77% of
its hepatocyte-class cells (AKR1B10 60%, SPINK1 91%) while HCC04T was a GPC3-negative tumour. With
no second authors' file to settle it, the data were run as deposited and an arm without patient 4
was declared in advance as the follow-up if the result hinged on that library.

| arm | state (lenient / strict) | GPC3 | decoys | drugs | n |
|---|---|---|---|---|---|
| calibration, fresh install | run aborted (clustering crashed) | - | - | - | 1 |
| calibration, pinned env | no state genes (raw counts destroyed) | - | - | - | 1 |
| calibration, scoring v2 | no / no | none (unscored) | TOP2A **unflagged** | 0 | 1 |
| calibration, scoring v3, all libraries | no / no | none (unscored; log2 2.6, padj 0.20) | none | 0 | 2 |
| calibration, scoring v3, without patient 4 | no (1/5) / no | **4** (log2 6.5, padj 8e-12) | none | 1.0 (codrituzumab) | 2 |

**Three infrastructure defects first, none visible on the earlier datasets.** A fresh install
passed all tests and the offline demo yet could not cluster any real dataset (`pynndescent 0.7`
returns a scipy sparse-array type that `scanpy 1.11` cannot handle, only on the >8k-cell path and
inside Scrublet); the scripted analyst treated the crashed tool as an empty result and the report
read as a clean biological negative; and with counts in `X` and no counts layer configured, the
first normalisation destroyed the raw counts so every pseudobulk call refused. Fixes: pinned
environment plus `sctarget.py check-env` (which exercises exactly those paths), a hard abort with
"tool X failed" in the scripted model, and a loader that stashes raw counts before anything
touches `X`. Sections 18-18b of the protocol.

**Then a method limitation.** The scoring v2 compartment contrast compares disease with reference
cells within the lineage's shared, non-enriched clusters. In a tumour, 90-100% of malignant cells
sit in the enriched clusters, so that contrast compared the atypical 10% that co-cluster with
normal hepatocytes and GPC3 (+2.4 log2) did not reach significance; the sibling contrast was flat
(-0.13) as in ccRCC. A cycling tumour sub-state also put TOP2A through unflagged, the first
unflagged decoy in three datasets. Scoring v3 adds the plain whole-lineage contrast (all disease
cells of the lineage vs all its reference cells), chooses the representative lineage axis by a
stated rule and records the choice per gene, and flags cell-cycle genes (section 19). Regression:
ccRCC CA9 stays first; CAA ERBB4 moves from 5 to 4.

**The label problem decided the outcome.** Under v3 with all libraries, the whole-lineage contrast
had GPC3 at log2 2.63, baseMean 764, padj 0.20: a replicate-aware test refusing to call a gene
disease-specific while a "normal" library expresses it in 77% of its cells and contributes 22% of
the reference hepatocytes. With the patient-4 pair excluded by a one-line, recorded config change
(`dataset.exclude_samples`), GPC3 was first in the tumour-vs-normal hepatocyte contrast (log2 6.45,
padj 8e-12), ranked 4th overall behind CA2, RAMP1 and PGF (none a decoy; RAMP1 is also high in
normal fibroblasts and the critic flagged it), nominated, with codrituzumab found and no decoy
nominated. The state criterion still reads "missed": the recorded signature is cholesterol
synthesis (SQLE, FDPS, IDI1), NDUFA4L2 and CYP17A1, consistent across patients, while the key's
literature markers (AFP, MDK, SPINK1, AKR1B10) are patient-heterogeneous here and never
significant. As with EPAS1, the key stays frozen and the lesson is recorded: state markers in a
key should be genes expected to be consistent across patients.

**A marker-free label check.** The exclusion had been identified with the expected target itself,
which is a weakness for a blind test and useless without a known marker. `check_reference_libraries`
now asks, within each lineage and with no gene names, whether any library resembles the opposite
condition; it is gated on the lineage showing condition structure at all. On the raw object it
flagged HCC04N in the T/NK lineage (borderline in hepatocytes, 0.03 under the 0.05 margin); after
QC its hits shift, mostly immune-lineage crossovers in both arms, which says the immune compartment
separates poorly by condition here and that the check is sensitive to filtering. One tumour
library's hepatocytes correlate with the normals (a well-differentiated tumour), which dilutes
every tumour-vs-normal contrast. Reported as partial, independent support for the exclusion.

**Other candidates.** FOLH1 (PSMA; expressed in HCC neovasculature), PGF, ESM1, THY1 and CA2 come
from a tumour-endothelium programme the analyst recorded; SQLE / FDPS / IDI1 / SCD are the
cholesterol and lipid synthesis programme of HCC; MAGEA3 is a cancer-testis antigen. All are
hypotheses; the critic's direction caveat applies to each.

## 6. Ovarian cancer omental metastases (human, scRNA-seq, GSE286057, 6 tumour / 23 omentum libraries, 145k cells)

Chosen by searching GEO for 2024-2026 deposits with tumour and normal tissue and an approved drug
with a transcript-level target; a 2025 gastric atlas was downloaded first and dropped on
composition (70% immune, 1-2% epithelium in both conditions: nothing epithelial to test), which
`inspect` now reports. Design: 6 HGSOC omental implants (5 with substantial tumour content) vs 23
omentum libraries (6 distant omentum from the same patients, 17 from 8 non-cancer donors);
peritumoral libraries left out as a declared second arm. Omentum has no epithelium, so this is the
pure composition case: every epithelial gene is tumour-specific, and the key said in advance that
the test was whether FOLR1 is found and ranked among the other real ovarian-cancer surface
targets (MUC16, MSLN, TACSTD2, CLDN6, ERBB2, EPCAM, all acceptable), not whether it beats them;
that MSLN should come out shared with the mesothelium; and that two distant-omentum libraries with
tumour cells (EPCAM 10% and 3%) were the pre-registered exclusion.

| arm | state (lenient / strict) | FOLR1 | nominated / decoys | tokens |
|---|---|---|---|---|
| scripted calibration, all libraries | yes (5/5) / no | unscored (sibling log2 4.8, padj 5e-8, 39th of the significant genes, outside the fixed top-20 pool) | 20 / none | 0 |
| **agent**, all libraries | yes / no (3/5) | **16**; set aside as shared_with_reference | 11 / none | 2.0M |
| scripted, two leaking libraries excluded | yes / no | 96 | 20 / none | 0 |
| **agent**, two leaking libraries excluded | yes / **yes (5/5)** | blocked (reference_driven) | 8 / none | 1.45M |
| **agent**, same arm, prompts v6 | yes / **yes (5/5)** | blocked (reference_driven) | 5 / none | 2.3M |

**What the agent added over the script.** The script found the tumour state and had FOLR1
significant in its tables but never looked past rank 20. The agent recorded the state with its
real identity (EPCAM/KRT8/KRT18/KRT19, 84% tumour, 6/6 patients), reached FOLR1 and scored it, set
it aside as `shared_with_reference`, and the critic then traced the cause: `check_reference_libraries`
flagged two libraries labelled Omentum as tumour-like in up to six lineages and asked for the
`exclude_samples` arm by name - the same two libraries the pre-run marker check had identified.
It also called the macrophage state a composition shift (blocking), the TNFRSF4/9/18 antibodies
agonists (direction), and EPCAM pan-epithelial. In the exclusion arms it met the strict state
criterion both times, found a CAF state (INHBA, LOXL2, FAP, MMP11) and a Treg state, and nominated
no decoy in any of the three runs.

**Why FOLR1 ended blocked, and why we stopped there.** With the two leaking libraries gone, the
tumour cluster still held 151 reference-labelled cells, which expressed FOLR1 more highly than the
tumour cells, and the in-cluster rule read that as "the signal comes from reference cells". On the
raw object those cells are 59 FOLR1+ reference cells in 94,235 (0.06%), 26 from one healthy
donor's lesser omentum, KRT19/CLDN3/WFDC2/EPCAM-positive and CALB2/UPK3B-negative: a benign
Müllerian-type epithelial inclusion, which expresses folate receptor alpha as normal Müllerian
epithelium does. The pipeline was right about the cells and wrong about the verdict. Two generic
guards were tried (protocol section 21): a representativeness test on the in-cluster reference
pool, first cell-relative (withdrawn the same day: it also freed Treg genes from a legitimate
reference comparison), then lineage-relative (correct on rescore, but defeated in the third agent
run when the coarse derived lineage absorbed mesothelium and so "held" the epithelial lineage).
A third rule was drafted and not applied. The honest reading is that condition attribution is the
pipeline's most conservative step, that its verdicts depend on how clean the derived lineages are,
and that when it sets a known target aside the critic says exactly why - which is what a reviewer
needs. FOLR1 is reported as found and ranked 16th in the primary run, set aside with a stated
reason, and blocked in the exclusion arms by a rare benign population; the method stays frozen.

**Other findings of this round, all fixed:** the first agent run died before its first tool call
on an OAuth token-refresh race (transient CLI errors now retried); 19 druggability lookups failed
on HTTP 429 and left the Open Targets clinical fields empty (Retry-After backoff); the critic
reasoned its way to "I have no tool to exclude samples" (it is told now that this is a human,
config-level decision to request); and two GEO deposits in a row needed a generic per-sample 10x
fetcher and builder, now in `scripts/`.

## 7. What holds across the four datasets

The tools do most of the work; the model's contribution is visible at the interfaces: whether
findings carry between roles, whether a merged cluster is examined, whether a thin reference
pool is read as a caveat or a verdict. Nearly every failure was a hand-off or rule defect, and
every fix was a generic rule tied to a named run. The compartment contrast (disease vs reference
within the lineage's shared clusters) was the robust one; the sibling contrast is fragile when
disease cells drift far from their normal counterpart (both de novo arms). The critic caught a
metadata label swap, a dissociation artefact, doublet clusters, a label inconsistency between two
tools, thin reference pools and a drug-direction error; the marker-free library check now adds
doubtful labels. One decoy was nominated unflagged in 15 scripted and agent runs (TOP2A, HCC,
before the proliferation flag); none since.

The two reviewer rounds were worth more than the two development datasets: most of their
findings were infrastructure defects every new user would have hit, and the method change the
first forced (the whole-lineage contrast) is the comparison a biologist would have run first. A
reference set with a mislabelled, contaminated or rare-counterpart library is the common case, not
the exception; in all three public datasets the known target's fate was decided by what sat in the
reference. The pipeline now says so when it can see it, the exclusion arm is a disclosed config
choice, and the remaining limitation is stated rather than patched: attribution verdicts depend on
how clean the derived lineages are.

Across the 20-odd runs that reached nomination: the known target was ranked first in CAA and
ccRCC, fourth in HCC once a mislabelled library was excluded, and sixteenth in the primary ovarian
run; the agent met the strict state criterion in ccRCC (twice) and ovarian (twice); one decoy was
ever nominated unflagged (TOP2A, before the proliferation flag).

## 8. Limitations

Expression-based: protein-level targets are invisible. The composition route favours
infiltrates consistent across patients over malignant cells that partition by patient (mitigated
by merging, not removed). The nominator over-nominates (precision 0.1-0.3) and never subclusters
unprompted. The subscription backend has no seed or temperature control; stability is
behavioural. Four datasets, all tumours; the programme route for non-tumour disease has only
been exercised on synthetic data. Condition attribution is conservative and lineage-dependent: a
rare benign population expressing the target inside the tumour cluster can set it aside (ovarian),
and the derived lineages at coarse resolution decide whether the guards against that apply. The HCC result depends on excluding one patient whose "normal"
library expressed the tumour programme; the all-libraries result is the primary one and is a miss.
The marker-free label check is sensitive to cell filtering and cannot judge lineages without
condition structure. Drug knowledge is human-centric; non-human data rest on
ortholog mapping. All candidates are hypotheses for IHC / functional validation.

## 9. For the CoPAIR community

The repository runs on any two-condition dataset from a four-line config (`scripts/sctarget.py
check-env / inspect / calibrate / run / report`, plus `build_from_matrix.py`, `fetch_geo_10x.py` and `build_from_10x_samples.py` for GEO deposits), needs no API key (Claude Code login), derives lineages when
no annotation exists, and ships a skill (`skills/sc-target-agent/SKILL.md`) that walks a Claude
Code session through a new dataset with the lessons above built in. The protocol document is the
method history; `tests/` must pass after any change.

## 7. Round 0.4: a non-malignant disease and a perturbation (IPF, GSE135893; doxorubicin heart, GSE285475)

The 0.4 skill (design classification, pre-run design and label checks, a no-key protocol, labelled
limitations) was written on the assumption that the method generalises beyond tumours. This round
tested that on the two design types the protocol had named as never run, with the same rules as
before: keys frozen and tagged before any run, arms pre-registered, every method change a labelled
version motivated by a named failure, tumour datasets re-run as regression. Protocol sections
23-23g; results tables `docs/results_ipf_v04.md`, `docs/results_dox_v04.md`, `docs/results_ccrcc_v04.md`.

**What broke first.** The free step that opens every analysis, the scripted calibration on the
synthetic programme dataset, ended with "no disease-specific state found" on data where every case
keratinocyte runs a planted programme. Four generic defects sat on the expression route - the path
every perturbation and most non-malignant diseases take - and none had been visible on the tumour
datasets, which take the composition route: the programme test ran on a resolution-2.0
re-clustering that split each lineage into pieces too small to pass (23); two clusters of one
lineage that both carried the programme cancelled each other under "3x any other cluster" (23b);
attribution classed every programme gene with a non-zero baseline `shared_with_reference`,
consistency was detection-based and the two-thirds cut-off flagged exactly 2 of 3 samples, and the
nominator nominated one gene per cluster (23d); and the candidate pool kept the most significant
genes by statistic only, so the largest inductions never reached scoring (23f, three touches to
reach one rule end-to-end). After them the synthetic calibration passes (ERBB4 first, three unique
nominations) and the tumour regression is clean.

**IPF (12 donors vs 10, libraries pooled to donor; key: three states, ITGB6 secondary).** The
composition route dominates, not the programme route the design table had predicted: basal-cell
expansion, MUC5B+ secretory cells, ciliated cells (a composition shift the critic rejected), an
IPF-only AT2 state and two macrophage clusters are all enriched across donors, and the expression
route finds the senescence / aberrant-basaloid programme (CDKN2A, MMP7, KRT17, SOX9) spread over
the transitional and AT2 compartments. Neither the floor nor the v7 agent characterised the
macrophage clusters (three epithelial states filled the cap); prompts v8 - the top enriched cluster
of every lineage, the expression route in every analysis - turned the SPP1+ profibrotic macrophage
state into a strict recovery (9/5 and 6/5 markers; SPP1 ranked 1st within its state, 3rd merged).
The 374 KRT5-/KRT17+ aberrant basaloid cells sit inside a 5,559-cell transitional cluster that no
run subclustered, so that state and ITGB6 (significant only in the AT2 programme, log2 0.9) were
not reached; the v8 rule that asks for exactly that subclustering was not followed once, and is
left as the stated limitation rather than iterated on one dataset. Nominations (MMP1, MMP7, TNC,
ITGA2, SPP1, MMP9, CTSK ...) are published IPF biology the key's acceptable list did not contain;
no decoy was nominated in any run. The pre-run compartment-restricted label check and the critic's
marker-free check flagged the same donor (control-like macrophages in the largest IPF explant);
the pre-registered arm was not needed because the macrophage state was recovered with it in.

**Doxorubicin vs vehicle (4 vs 4 wild-type hearts, `objective: marker`; key: p53 programme,
Cdkn1a).** Composition is flat, as predicted; the expression route finds the cardiomyocyte
response (Trim63, Nppb, Fbxo32, Nppa, Ankrd1, Myh7 - acceptable genes of the key) in floors and
agent alike, and the agent's critic rejected a second "state" carried by 3 of 4 libraries and asked
for the pre-registered batch arm, which supports the atrophy response (Trim63 nominated at 3 vs
3). The expected read-out, p21, is not induced at two weeks of weekly dosing (down in the
cardiomyocyte programme, modest and non-significant elsewhere), as the pre-run label check had
hinted. This is the third key-drafting failure mode after EPAS1 (protein-level) and the HCC
literature markers (not consistent in the data): a real response gene sampled at the wrong time.
The key is unchanged and the result is reported as "programme found, expected target absent at
this time point, acceptable genes at ranks 9-12 of their state". Nominations were weak and the
critic said so.

**Regression and cost.** ccRCC under v8 and all 23x changes: tumour state found (lenient 9/5,
strict 8/5), CA9 first in both tumour states, EGLN3 second and flagged for direction, CD70 and
FOLH1 above the key, no decoy. The merged target rank fell to 6th because v8 records the TAM,
endothelium and stroma states too; the evaluator now judges the rank within the state that carries
the target (23g) and keeps the merged rank as information. v8 roughly doubles the analyst's tool
calls: 1.2-4.1M input tokens per run in this round against 0.6-2.3M before. Stability pair on
ccRCC under v8: both runs find the tumour state strictly (9/5 and 8/5) with CA9 first in its state
and no decoy; top-10 Jaccard 0.33, nominated-set Jaccard 0.42, rank-1 gene CA9 in both; the
variation is in the secondary states (TAM, endothelium, stroma) and their candidates, as before.

**What the round says about generality.** The code was design-agnostic; the method was not, in
four places, all on the expression route, all invisible on tumours, all generic once seen. The
skill's design table was right that perturbations are the pure programme case and wrong that
non-malignant disease is "mostly programme": IPF creates disease-only populations and the
composition route found them. Keys for non-tumour designs fail in a new way (time), and the
pre-run checks the skill added (compartment-restricted label check, design check with donor
nesting and covariates) found the same libraries the critic later flagged. Nobody but the author
has run any of this.

## 7b. Round 0.4b: the three earlier tumour datasets under the 0.4 skill (HCC, ovarian, fibrolamellar carcinoma)

Same rules as round 0.4. The HCC and ovarian objects had been deleted and were rebuilt from GEO
with the repository's builders (HCC identical; ovarian 124k cells for the same 29 libraries, with
region and donor columns in `obs`); the FLC object (the lab's Parse snRNA-seq, 354k nuclei, 9
primary / 9 normal / 7 metastatic libraries) and its key came from the 2026-11-02 test clone,
frozen there before any run. Keys unchanged; the pre-registered arms of earlier rounds carried
over. Protocol 23h-23l; results tables `docs/results_{hcc,ovarian,flc}_v04b.md`; running notes
`validation_0.4/NOTES.md`. Seven floors, seven agent runs, about 23M input tokens.

**Pre-run checks did their job.** `design_check` and `label_check` reproduced every library anomaly
the earlier rounds had found after the fact: HCC04N expressing the tumour programme, the two
tumour-bearing "distant omentum" libraries, FLC83_T and FLC57_T with low tumour content. A second
rule in the label check, added when the whole-library median-gap rule called every ovarian marker
uninformative (tumour cells are a minority of tumour libraries), also listed the healthy donor's
PAX8+ library that holds the benign epithelial inclusion behind FOLR1's `reference_driven` verdict
(section 6). Two other helper corrections: a covariate nested within the condition (tumour vs
omentum regions) is no longer called "perfectly confounded", and pairing is no longer reported as
pseudo-replication.

**HCC.** Floors as in round 0.2: with all libraries GPC3 is unscored (lineage contrast padj 0.20
with HCC04N in the reference; each tumour hepatocyte cluster's siblings are the other tumour
clusters, so the pan-tumour gene cancels); without patient 4, **GPC3 1st within its state**
(round 0.2: 4th), codrituzumab found, nothing above it. The agent run with all libraries recorded
six states and ranked GPC3 2nd - in a CAF-like stromal state, through the disease-wide contrast;
in the tumour hepatocyte states it was unscored as before. The within-state rule (23g) reported
that as recall 1.0, which was right by rule and wrong as a reading; the evaluator now prints the
lineage of the carrying state and honours an optional `lineage:` on a key's expected state (23i).
Frozen keys have no such field, so their numbers stand and this one is read as "scored in a
stromal state, unscored in the tumour states".

**Ovarian.** The floor found the tumour state (5/5) and left FOLR1 unscored although it was the 9th
largest effect of the sibling contrast: 39th by statistic, inside the stored 40-gene head but
outside the 20 the floor read, and not labelled an effect gene. 23h takes the effect set over the
whole stored list; the re-run scores FOLR1 24th. Agent, all libraries: state 6/5, **FOLR1 14th**
(0.3: 16th), the critic calling it "the best-supported druggable gene" and asking for the
exclusion arm by library, which is the reading the skill asks for. No-leak arm: floor 35th
(`reference_driven`); agent strict 5/5, FOLR1 17th, blocking flag, in two runs - the second under
prompts v9, which say a reference-driven candidate with per-sample support is reported as
conditional on the identity of the reference cells rather than dropped (23j; the skill had said so
since 0.4.0 and the prompt had said "drop outright"). The v9 critic named the cells ("fallopian-
tube-like or mesothelial epithelium, or tumour-contaminated"), named the settling experiment ("a
reference-only subcluster of cluster 24"), and still took the representative branch, because the
derived epithelial lineage absorbs mesothelium at this resolution (21b). The tool for the
settling move now exists (`subcluster(condition=...)`, 23k, prompts v10); it has not yet been
used in a run. Five agent runs on this dataset now agree: FOLR1 is found in the right state and
held back by 59-151 benign reference epithelial cells, visible before the run.

**Fibrolamellar carcinoma (first run of this dataset under any version with a result).** The
tumour state came through the **expression route**: with one condition per library, Harmony on
the library co-clustered tumour nuclei with normal hepatocytes, and the state is a programme
inside a shared hepatocyte cluster (2,009 genes up within the cluster: CA12, KRT7, CD24, IRF4).
The 0.4 rule that the route runs in every analysis is what reached it; the v6 floor on 2026-11-02
had missed it. Floor and both agent runs: **PRKACA 2nd within its state behind CA12** (both key
answers), nothing above them, ERBB2 2nd among acceptable genes with neratinib found, PRKACA
nominated by the agent in both runs; the two agent runs share the same top three in the tumour
state (nominated-set Jaccard 0.62, top-10 Jaccard 0.54). The key's literature markers again fail
the state criterion (3/5) because they split between the composition state and the programme
state - the same key-drafting failure as HCC. The first agent run nominated PDGFRB from the
CAF-like matrix programme without a flag, the only decoy nominated unflagged in the round; prompts
v11 (23l) ask the critic to flag every nomination from a non-disease-cell lineage as a
microenvironment candidate, and the second run flagged FAP, COL1A2 and VWF that way with no decoy
unflagged.

**What the round says.** The 0.4 skill found a state the previous method missed (FLC), moved every
floor without touching a key (GPC3 4th to 1st, FOLR1 unscored to scored, PRKACA 6th to 2nd), and
caught the library problems before tokens were spent. It did not change the two hard verdicts
(FOLR1 with benign reference epithelium; GPC3 with HCC04N), and it costs about twice as much per
agent run. Three prompt versions in one round, each verified once, is thin verification. The
next step is not another method change but an independent run by someone other than the author.

## 7c. Round 0.4c: the development dataset under the 0.4 skill (canine acanthomatous ameloblastoma)

The dataset the method was built on (lab snRNA-seq, 205k nuclei, 6 CAA / 10 OSCC / 3 gingiva
libraries; design row 3, tumour vs another tumour plus healthy tissue), last run under prompts v2
and scoring v2 (results-v3), now run as the skill says: both helper checks first, key tagged
(`caa-key-frozen-v04`, unchanged since its 0.3 freeze), floor, agent run, a pre-registered
exclusion arm, a second agent run, then (at the author's request) the de novo arm with Harmony and
derived lineages. Protocol 23m-23p; results `docs/results_caa_v04c.md`; notes
`validation_0.4/NOTES.md` (round 0.4c). Three floors, five agent runs, about 9.8M input tokens.

**Pre-run checks.** Design: three condition levels by design, 10-vs-3 reference imbalance, two dogs
paired (tumour and matched gingiva; after the 0.4b fix no longer reported as pseudo-replication),
BRAF status associated with condition (biology, not batch). Labels: PEG3 and RBFOX1 separate CAA
from both references in every library; nothing flagged; decision recorded before the run: all 19
libraries kept. The floor's marker-free check then flagged one OSCC library whose epithelial and
mesenchymal pseudobulks correlate with CAA; 14.5% of its cells sit in the two CAA-exclusive
clusters against 0.3-1.2% for the other OSCC libraries. A marker check cannot see a library that is
"a bit CAA-like" when the marker is 10% there and 60-87% in CAA; the lineage-wise correlation can.
Exclusion arm pre-registered (`config/caa_nos3.yaml`), the all-libraries run primary.

**Results.** Floor: both CAA-exclusive keratinocyte states found by composition (clusters 6 and
12, 95-96% disease cells, 6/6 CAA libraries); key state matched 9/5; **ERBB4 5th within its state**
behind GABRA1, GABRA4, GABRB1 (acceptable) and GABRB2; nominated; afatinib, dacomitinib and
neratinib found; no decoy ranked or nominated (0.3 floor: ERBB4 unranked, 24 novel genes above the
first acceptable gene). Agent run 1 (v11): strict state match 6/5, **ERBB4 4th**, nominated, drugs
found; FGFR2 nominated and flagged major by the critic as present in 86% of the cluster's
reference cells - the paper's own negative control, correctly set aside. Agent run 2 (v12): strict
5/5, **ERBB4 4th**, nominated, drugs found, no decoy nominated; top-10 Jaccard with run 1 0.67.
Exclusion-arm floor: ERBB4 4th, otherwise unchanged; the excluded library was not what held ERBB4
behind the GABA-A subunits. Exclusion-arm agent run (v12): strict 5/5, **ERBB4 3rd** with only
GABRA4 and GABRB1 above it, nominated, drugs found; FGFR2 nominated and flagged again ("pan-
epithelial, shared with reference"); the critic reviewed all ten nominations and noted that the
marker-free check now flags only endothelial libraries, so the arm had removed what it was meant
to remove. The GABA-A receptor subunits above ERBB4 in every run are the same
neuroepithelial programme the paper describes (GABRB1 is in its top-20 signature); they are druggable
only by anaesthetics and benzodiazepine-site ligands, which the critic said in both runs. The
critic in run 2 used `subcluster(condition=OSCC)` unprompted to characterise the 1,237 OSCC cells
inside cluster 6 and made the cluster-6 candidates conditional on them: the settling move built for
the omentum case (23k), applied to a reference that is another tumour.

**De novo arm** (Harmony on the library, Leiden recomputed, lineages derived, nothing of the paper's
processing reused). Harmony did not absorb the CAA state: the CAA keratinocytes (59k nuclei, 99%
disease cells, split into three Leiden clusters) formed their own derived lineage, with a smaller
CAA cluster inside the shared epithelial lineage and a cycling CAA cluster besides. Floor: key state
matched 7/5 to the shared-epithelium CAA cluster, **ERBB4 1st within its state** (global 4th),
nominated, drugs found, no decoy ranked or nominated. Agent run 1 (v12): **ERBB4 1st within its
state** in two states, global 2nd behind IGF2, nominated twice, drugs found; FGFR1 nominated from
one CAA piece and flagged major ("not disease-elevated, shared with reference"); the critic merged
the three CAA pieces into one population itself. Agent run 2 (v13): **ERBB4 1st within its state**,
global 3rd behind two stromal IGF-axis genes from a fibroblast cluster that the critic flagged as a
microenvironment candidate; nominated, drugs found, no decoy nominated; the critic wrote that the
CAA lineage "is condition-exclusive, so the within-lineage test is void, not failed" (23p). ERBB4
ranks 1st here and 3rd-5th with the paper's clustering because the GABA-A subunits spread across
the three CAA pieces and the shared-epithelium CAA cluster's siblings are OSCC and gingival
keratinocytes, not other CAA cells. Stability between the two de novo agent runs is low (top-10
Jaccard 0.18, nominated 0.19): with derived lineages there are five enriched clusters plus program
clusters and the two runs recorded different subsets, so the lists differ by state choice rather
than by disagreement about one state; the scored quantity (ERBB4 1st within its state) is the same
in the floor and both runs. Strict state recovery was not met on this arm (2/5, 2/5, 4/5): the key's
literature markers are spread across the CAA pieces.

**What the round found in the method.** The agent's own notes pointed at three defects, the
floor at a fourth. (i) The
critic wrote that the candidate list it received "was truncated after ERBB4". Every stage briefing
was one JSON dump cut at a flat character cap; re-measured on every agent run of rounds 0.4 and
0.4b, the critic's claims exceeded the cap in all of them, so in every run the critic had been
briefed on the first three to eight nominations and had to find the rest through its own tool
calls. PDGFRB in FLC run 1, the decoy that motivated the microenvironment flag, was among the
nominations that critic was never shown; GPC3 in the HCC agent run and PRKACA in both FLC runs
likewise. The "decoys nominated unflagged" numbers of rounds 0.3-0.4b stand as recorded but did not
measure the critic's rule-following on nominations past the cap. Prompts v12 (23m): an index line
names every nomination and state and is never shortened; free text is clipped per item until the
list fits; gene lists shrink evenly and never lose a key. Verified in run 2 ("all 7 nominated genes
have a recorded flag"). (ii) The ortholog tool returned "no ortholog" for GABRB1 and GPC5 because
the dataset's precomputed table was empty for them, while the druggability lookup used the same
symbols as human; an empty entry now falls through to Ensembl (GPC5 resolves), and when Ensembl
Compara itself lists no human homologue for a symbol that is a human gene - GABRB1 on the dog's
moved assembly - the symbol is kept on a flagged assumption rather than dropped (23n). (iii) The
report's evaluation section now reflects the scoring (23o). (iv) On the de novo arm the scripted
critic flagged two CAA clusters as "enrichment does not hold within its parent cell type" because
their derived lineage was itself 99% CAA, so the within-lineage test compared CAA sub-states with
each other. When a lineage is condition-exclusive the lineage is the disease population and the
test is void, not failed; the scripted critic now checks lineage composition first and the critic
prompt carries the rule (23p, prompts v13; verified in de novo run 2).

**What the round says.** On the development dataset the 0.4 skill reproduces the 0.3 answer
(ERBB4 found in the right state, nominated, with its drugs) with a more honest ranking (within the
state, behind the programme's own channel subunits), records both CAA-specific states as
distinct, flags the paper's negative control when it is nominated, and finds a reference library
anomaly the marker check could not. Without the paper's clustering and with batch correction the
answer is the same and ERBB4 ranks 1st within its state, at the cost of a less stable nomination
list across runs. The briefing defect is the most consequential finding of the
three 0.4 rounds: it was present in every agent run since the method was frozen, it was found by
reading what the critic said about its own inputs, and nothing in the evaluator could have shown
it. That is an argument for keeping the agent's notes in the report and reading them.
