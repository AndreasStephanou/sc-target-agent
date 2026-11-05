# Validation results (summary)

Regenerated tables with every run's exact numbers are `docs/results_caa.md` and
`docs/results_ccrcc.md` in the repository, produced by `scripts/collect_results.py` from the run
directories; the arm labels come from each run's own recorded configuration.

## Dataset 1 — canine acanthomatous ameloblastoma (snRNA-seq, unpublished, De Vlaminck lab)

205k nuclei, 6 CAA / 10 OSCC / 3 gingiva libraries; OSCC and gingiva used as reference. Known
answer: a neuroepithelial CAA state and ERBB4, with neratinib sensitivity shown in vitro. No sample
identifiers or per-cell data are reproduced here.

| arm | state (lenient / strict) | ERBB4 rank | decoys | drug recall | n |
|---|---|---|---|---|---|
| baseline script | lenient only, wrong cluster / no | none | 2 ranked above, none nominated (no critic) | 0 | 1 |
| scripted tools (calibration) | yes / no | 5 | none | 0.75 | 1 |
| agent, prompts v1 | yes / 0 genes recorded | none | none | 0 | 2 |
| agent, prompts v2 | yes / no | 2, 2, 2 | FGFR2, flagged 3/3 | 0.75 | 3 |
| agent, v2 + scoring v2 | yes / no | 1 (×4) | FGFR2, flagged | 0.75 | 3 rescored + 1 |
| de novo clustering, inherited annotation | yes / no | none → 4 after thin-pool guard | none | 0 → 0.75 | 2 |
| de novo clustering, derived lineages | yes / no | 1 | none | 0.75 | 1 |

Reading: the tools reached the right cluster in every agent run; the differences between arms are
hand-off and rule defects, each fixed by a generic versioned rule. FGFR2 (a decoy in the key) ranked
first in the v2 runs, was flagged by the critic every time, is state-specific in this data and is
recurrently mutated in human ameloblastoma; the key stays frozen and this is an open hypothesis.

## Dataset 2 — human clear cell renal cell carcinoma (scRNA-seq, public, GSE159115)

29k cells, 7 tumour / 6 benign libraries from 7 patients. Known answer: the HIF programme; key
expected EPAS1 (a drafting error: protein-level target) with CA9 and VEGFA acceptable.

| arm | state (lenient / strict) | EPAS1 | CA9 (acceptable) | decoys | n |
|---|---|---|---|---|---|
| scripted, calibration passes 1–3 | no → yes / no | none | 3 → 1 | KRT8/18 ranked, none nominated | 3 |
| scripted, calibration pass 4 | yes / no | none | 1 in every tumour state | none | 1 |
| agent, authors' annotation | yes / **yes** (8/5) | none | 1 | none | 1 |
| agent, derived lineages | yes / **yes** (7/5) | none | 1 | none | 1 |

Reading: four generic defects surfaced by calibration on dissociated tissue (QC, patient split,
stress cluster, silent doublet failure) and one in the contrast design (sibling contrast flat for
tumour genes); all fixed as versioned rules. The strict state criterion passed for the first time.
EPAS1 was never found and could not have been. The critic flagged a tumour/T-cell doublet cluster,
a patient-driven split, a programme in mixed cells and EGLN3's wrong therapeutic direction. CD70 and
FOLH1 appeared above the key as "novel" and are clinical-stage ccRCC targets the key omitted.

## Dataset 3 — human hepatocellular carcinoma (scRNA-seq, public, GSE149614), run as a reviewer

Fresh clone, fresh conda environment, README and skill only. 10 primary-tumour / 8 non-tumour
libraries (PVTT and lymph-node libraries excluded at build), 63k cells, authors' six-class
annotation in which `Hepatocyte` covers normal and malignant cells. Key frozen before any run:
expected GPC3 (transcriptionally distinctive; codrituzumab / CAR-T target), acceptable AFP, MDK.
Pre-run label check: one library labelled non-tumour (HCC04N) expressed the tumour programme
(GPC3 77%, AKR1B10 60%, SPINK1 91% of hepatocyte-class cells); an arm without that patient was
declared in advance.

| arm | state (lenient / strict) | GPC3 | decoys | drugs | n |
|---|---|---|---|---|---|
| calibration, fresh install | run aborted (clustering crashed: library version conflict) | - | - | - | 1 |
| calibration, pinned env | no state genes (raw counts destroyed by normalisation) | - | - | - | 1 |
| calibration, scoring v2 | no / no | unscored | TOP2A **unflagged** | 0 | 1 |
| calibration, scoring v3, all libraries | no / no | unscored (log2 2.6, padj 0.20) | none | 0 | 2 |
| calibration, scoring v3, without patient 4 | no (1/5) / no | **4** (log2 6.5, padj 8e-12) | none | 1.0 (codrituzumab) | 2 |

Reading: three install-time defects every new user would have hit (unpinned environment; a
crashed tool reported as a biological negative; raw counts in `X` destroyed before pseudobulk),
one method limitation (the shared-cluster lineage contrast compared only the 10% of tumour cells
that co-cluster with normals; scoring v3 adds the whole-lineage contrast and a proliferation
flag; regression on ccRCC and CAA clean), and one data-label result (a replicate-aware test
correctly refuses GPC3 while a "normal" library expresses it in 77% of cells; with that patient
excluded by a recorded config line, GPC3 is first in tumour-vs-normal hepatocytes). The key's
literature state markers were patient-heterogeneous and never significant, so the state criterion
reads "missed" even in the successful arm: a key-drafting lesson recorded alongside EPAS1. A
marker-free library check (`check_reference_libraries`) was added and gives partial independent
support for the exclusion (flags HCC04N in the T/NK lineage on the raw object; borderline in
hepatocytes). Protocol sections 18-19.

## Dataset 4 — ovarian cancer omental metastases (scRNA-seq, public, GSE286057, 2025), run as a reviewer

Chosen by a GEO search for 2024-2026 deposits; a 2025 gastric atlas was downloaded first and
dropped on composition (70% immune, 1-2% epithelium in both conditions). 6 HGSOC omental implants
vs 23 omentum libraries (6 distant omentum from the same patients, 17 from 8 non-cancer donors);
peritumoral libraries left out as a declared arm. Key frozen before any run: expected FOLR1
(mirvetuximab soravtansine); acceptable MUC16, MSLN, TACSTD2, CLDN6, ERBB2, EPCAM, with the
written note that, with no epithelium in normal omentum, every epithelial gene would be
state-specific and the test was whether FOLR1 is found and ranked among them. The first
model-driven runs since scoring v3.

| arm | state (lenient / strict) | FOLR1 | nominated / decoys | tokens |
|---|---|---|---|---|
| scripted calibration, all libraries | yes / no | unscored (significant, 39th; outside the fixed top-20 pool) | 20 / none | 0 |
| agent, all libraries | yes / no (3/5) | **16**, set aside as shared_with_reference | 11 / none | 2.0M |
| scripted, two leaking libraries excluded | yes / no | 96 | 20 / none | 0 |
| agent, two leaking libraries excluded | yes / **yes (5/5)** | blocked (reference_driven) | 8 / none | 1.45M |
| agent, same arm, prompts v6 | yes / **yes (5/5)** | blocked (reference_driven) | 5 / none | 2.3M |

Reading: the agent reached FOLR1 where the script never looked, set it aside, and the critic
traced the reason to two reference libraries carrying tumour cells (the same two the pre-run label
check had flagged) and asked for the exclusion arm by name. In the exclusion arms the block came
from 59 FOLR1+ reference cells in 94,235 (0.06%), 26 of them from one healthy donor: a benign
Müllerian-type epithelial inclusion that expresses folate receptor alpha, as normal Müllerian
epithelium does. Two generic guards were tried and the second was defeated by coarse derived
lineages absorbing mesothelium; the method was frozen with the limitation stated rather than
patched again. Also found and fixed: a token-refresh race that killed the first agent run, HTTP 429
failures that emptied the Open Targets fields, and the lack of a generic GEO per-sample fetcher.
Protocol sections 21-22.

## Across the four

One decoy was nominated without a critic flag in ~22 runs (TOP2A, HCC, a cycling tumour sub-state, before the proliferation flag); none since. The known target's fate in all three public datasets was decided by what sat in the reference set (a swapped label, a contaminated library, a rare benign counterpart); the pipeline now says so when it can see it. The agent met the strict state criterion four times (ccRCC x2, ovarian x2); the script never did. The compartment contrast was the robust axis;
the sibling contrast is fragile when disease cells drift far from their normal counterpart.
Precision of nomination 0.1–0.3. Cost 0.6–1.7M input tokens per agent run on the Claude Code
backend.

## Not yet done

Independent testing by another CoPAIR fellow; the with-skill vs without-skill comparison; an
agent run on HCC; a non-tumour disease dataset (IPF, GSE135893); ablations (critic off,
specificity off, druggability weight 0).


## Dataset 5 — idiopathic pulmonary fibrosis (human scRNA-seq, public, GSE135893), validation round 0.4

First non-malignant design (row 4). 12 IPF donors vs 10 controls (19 + 12 libraries pooled to
donor; other ILD diagnoses dropped at build), 89k cells, authors' four compartments as lineages.
Key frozen before any run: three states (aberrant basaloid epithelium, SPP1+ profibrotic
macrophages, CTHRC1+ fibroblasts), expected target ITGB6 (secondary metric), 20 acceptable
clinical-stage targets, 19 decoys; pre-registered arms for site (Vanderbilt-only) and for the
largest IPF donor whose immune compartment was control-like in the pre-run check.

| arm | states found (lenient / strict) | ITGB6 | nominated / decoys | tokens |
|---|---|---|---|---|
| scripted floor, before 23e | 0/3 / 0/3 (macrophage clusters enriched, never characterised) | unranked | 18 / none | 0 |
| agent, prompts v7 | 0/3 / 0/3 (three epithelial states; macrophages not recorded; expression route not run) | unranked | 8 / none | 1.70M |
| agent, prompts v8 | **1/3 / 1/3** (SPP1 macrophages 9/5 and 6/5; SPP1 ranked 3rd) | unranked | 14 / none | 4.10M |

Reading: the composition route dominates (basal expansion, MUC5B+, ciliated, an IPF-only AT2
state, two macrophage clusters), as the corrected design table says; the expression route finds
the senescence / aberrant-basaloid programme spread over the transitional and AT2 compartments.
Prompts v8 (lineage coverage; expression route in every analysis) turned the macrophage state
from "enriched but never characterised" into a strict recovery. The 374 KRT5-/KRT17+ aberrant
basaloid cells sit inside a 5,559-cell transitional cluster; neither the floor nor the model
subclustered it, so that state and ITGB6 (significant only in the AT2 programme at log2 0.9) were
not reached. Nominations (MMP1, MMP7, TNC, ITGA2, SPP1, MMP9, CTSK ...) are published IPF biology
the key's acceptable list did not contain. Cost: v8 roughly doubles the analyst's tool calls.

## Dataset 6 — doxorubicin vs vehicle, mouse heart (scRNA-seq, public, GSE285475 wild-type arms), round 0.4

First perturbation design (row 5), `objective: marker`. 4 vs 4 animals, 76k cells, lineages
derived; batch (two alignment pipelines) balanced across arms, pre-registered batch-2 arm. Key
frozen before any run: expected read-out Cdkn1a (p53 programme), acceptable p53 / atrophy /
natriuretic / Nrf2 genes, decoys incl. the transgene of the excluded arms, cycle, IEG, sex and
lineage genes.

| arm | programme state | Cdkn1a | acceptable genes ranked | nominated / decoys | tokens |
|---|---|---|---|---|---|
| scripted floor, before 23f | cardiomyocyte programme found (Trim63, Nppb, Fbxo32, Nppa significant) but outside the statistic-ranked pool | not induced | none | 3 / none | 0 |
| scripted floor, 23f complete | same programme; Nppb 12th, Myh7 21st, Ankrd1 28th | not induced | Nppb | 3 / none | 0 |
| agent, prompts v8 | programme recorded (Mb, Cryab, Cish, Lmod2, Rorc, Irs2, Ankrd1, Nppb, Txnip); a second cardiomyocyte "state" rejected by the critic; batch arm requested by the critic | not induced | Nppb 12th, Ankrd1 9th (state 20); Myh7 9th (state 18) | 5 / none | 1.17M |

Reading: the pipeline found the lineage-level cardiomyocyte stress / atrophy / natriuretic
response and said how many libraries carry it; the key's expected read-out (p21) is not induced
at two weeks, which the pre-run label check had hinted. Recorded as the third key-drafting failure
mode (time-point-dependent read-out), key unchanged. Three method fixes came out of this dataset
(protocol 23, 23b, 23d, 23f): before them the expression route reported a false negative on the
synthetic programme data and the floor's pool could not hold the response genes.


## Regression and stability after round 0.4 (ccRCC, GSE159115, prompts v8, scoring v3)

Two agent runs under every change of the round: tumour state found strictly in both (9/5 and 8/5
markers), CA9 first within its state in both, EGLN3 second and flagged for direction, CD70 and
FOLH1 above the key as in earlier rounds, no decoy nominated. Top-10 Jaccard 0.33 and nominated-set
Jaccard 0.42 between the runs; the variation is in the secondary states v8 now records (tumour-
associated macrophages, endothelium and stroma programmes) and their candidates. Cost 2.2-2.9M
input tokens per run (0.6-1.7M before v8). The evaluator now judges a target's rank within the state
that carries it (protocol 23g); the merged rank is reported as `target_rank_global`.

## What has and has not been tested (added in skill 0.4.0)

| dimension | tested | not tested |
|---|---|---|
| design type | tumour vs tissue of origin (ccRCC, HCC); tumour vs tissue lacking the lineage (ovarian omentum); tumour vs another tumour + healthy (CAA); non-malignant disease vs healthy (IPF, round 0.4); perturbation treated vs vehicle (doxorubicin heart, round 0.4) | knockout vs wild-type; infection; paired-only; blood / bone marrow; brain |
| objective | `target`; `marker` (doxorubicin heart, one agent run and four floors) | - |
| species | human (4), dog (1), mouse (1) | others |
| modality | scRNA-seq (3), snRNA-seq (1) | spatial, multiome, CITE-seq |
| replicates per condition | 6-23 libraries | 2-3 libraries (the minimum the skill allows) |
| annotation | inherited (3) and derived (4) | - |
| backends | Claude Code login, API (sonnet-4-5), scripted | other model families |
| independent use | none (round 0.4 was run by the author with the skill driving the session) | the with-skill vs without-skill comparison; a run by another CoPAIR fellow |

The known target's fate in every public dataset was decided by what sat in the reference set. The
pre-run checks added in skill 0.4.0 (`scripts/design_check.py`, `scripts/label_check.py`) exist to
surface that before tokens are spent; they have been run on the synthetic and ccRCC objects but have
not yet changed the course of a real run.
