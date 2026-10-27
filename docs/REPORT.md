# Can an agent rediscover a known drug target from single-cell data, blind?

*sc-target-agent: validation on two datasets. Draft for the CoPAIR fellowship report. Numbers in
the tables are regenerated from the run directories by `scripts/collect_results.py`
(`docs/results_caa.md`, `docs/results_ccrcc.md`); the method history is `docs/EVALUATION_PROTOCOL.md`.*

## 1. Question

Single-cell atlases of disease routinely end with a list of "potential targets" chosen by the
analysts. We asked whether an AI agent, given only the count matrix and the condition labels,
would recover a target that was independently validated, rank it above plausible decoys, do so
consistently, and do it better than a fixed script using the same tools. The known answers were
ERBB4 in canine acanthomatous ameloblastoma (CAA; our own snRNA-seq, neratinib sensitivity shown
in vitro) and the HIF programme of clear cell renal cell carcinoma (ccRCC; public scRNA-seq,
CA9 and the belzutifan target HIF-2alpha).

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

## 5. What holds across both datasets

The tools do most of the work; the model's contribution is visible at the interfaces: whether
findings carry between roles, whether a merged cluster is examined, whether a thin reference
pool is read as a caveat or a verdict. Nearly every failure was a hand-off or rule defect, and
every fix was a generic rule tied to a named run. The compartment contrast (disease vs reference
within the lineage's shared clusters) was the robust one; the sibling contrast is fragile when
disease cells drift far from their normal counterpart (both de novo arms). The critic caught a
metadata label swap, a dissociation artefact, doublet clusters, a label inconsistency between two
tools, thin reference pools and a drug-direction error. No decoy was ever nominated unflagged.

## 6. Limitations

Expression-based: protein-level targets are invisible. The composition route favours
infiltrates consistent across patients over malignant cells that partition by patient (mitigated
by merging, not removed). The nominator over-nominates (precision 0.1-0.3) and never subclusters
unprompted. The subscription backend has no seed or temperature control; stability is
behavioural. Two datasets, both tumours; the programme route for non-tumour disease has only
been exercised on synthetic data. Drug knowledge is human-centric; non-human data rest on
ortholog mapping. All candidates are hypotheses for IHC / functional validation.

## 7. For the CoPAIR community

The repository runs on any two-condition dataset from a four-line config (`scripts/sctarget.py
inspect / calibrate / run / report`), needs no API key (Claude Code login), derives lineages when
no annotation exists, and ships a skill (`skills/sc-target-agent/SKILL.md`) that walks a Claude
Code session through a new dataset with the lessons above built in. The protocol document is the
method history; `tests/` must pass after any change.
