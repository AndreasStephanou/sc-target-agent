# Can an agent rediscover a known drug target from single-cell data, blind?

*sc-target-agent: validation on three datasets. Draft for the CoPAIR fellowship report. Numbers in
the tables are regenerated from the run directories by `scripts/collect_results.py`
(`docs/results_caa.md`, `docs/results_ccrcc.md`, `docs/results_hcc.md`); the method history is `docs/EVALUATION_PROTOCOL.md`.*

## 1. Question

Single-cell atlases of disease routinely end with a list of "potential targets" chosen by the
analysts. We asked whether an AI agent, given only the count matrix and the condition labels,
would recover a target that was independently validated, rank it above plausible decoys, do so
consistently, and do it better than a fixed script using the same tools. The known answers were
ERBB4 in canine acanthomatous ameloblastoma (CAA; our own snRNA-seq, neratinib sensitivity shown
in vitro), the HIF programme of clear cell renal cell carcinoma (ccRCC; public scRNA-seq,
CA9 and the belzutifan target HIF-2alpha), and GPC3 in hepatocellular carcinoma (HCC; public
scRNA-seq, codrituzumab / CAR-T target), the third run as a reviewer would run it: fresh clone,
fresh environment, README and skill only.

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

## 6. What holds across the three datasets

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

The reviewer run on HCC was worth more than either earlier dataset: three of its five findings
were infrastructure defects that every new user would have hit, and the method change it forced
(the whole-lineage contrast) is the comparison a biologist would have run first. A dataset with a
mislabelled or contaminated reference library is probably the common case, not the exception;
the pipeline now says so when it can see it, and the exclusion arm is a disclosed config choice.

## 7. Limitations

Expression-based: protein-level targets are invisible. The composition route favours
infiltrates consistent across patients over malignant cells that partition by patient (mitigated
by merging, not removed). The nominator over-nominates (precision 0.1-0.3) and never subclusters
unprompted. The subscription backend has no seed or temperature control; stability is
behavioural. Three datasets, all tumours; the programme route for non-tumour disease has only
been exercised on synthetic data. The HCC result depends on excluding one patient whose "normal"
library expressed the tumour programme; the all-libraries result is the primary one and is a miss.
The marker-free label check is sensitive to cell filtering and cannot judge lineages without
condition structure. Drug knowledge is human-centric; non-human data rest on
ortholog mapping. All candidates are hypotheses for IHC / functional validation.

## 8. For the CoPAIR community

The repository runs on any two-condition dataset from a four-line config (`scripts/sctarget.py
check-env / inspect / calibrate / run / report`, plus `build_from_matrix.py` for GEO matrix deposits), needs no API key (Claude Code login), derives lineages when
no annotation exists, and ships a skill (`skills/sc-target-agent/SKILL.md`) that walks a Claude
Code session through a new dataset with the lessons above built in. The protocol document is the
method history; `tests/` must pass after any change.
