# Validation round 0.4 — running notes

Working copy of the repository (clone; remote removed; nothing pushed). Skill 0.4.0 installed in
`skills/` and `.claude/skills/`. All logs in this directory. Every method change is a protocol section
(`docs/EVALUATION_PROTOCOL.md` 23, 23b) and a labelled commit; every key was committed and tagged
before any run on its dataset.

## Step 1 — free checks (done)

| check | result |
|---|---|
| `check-env` | pass (`check_env.log`) |
| ccRCC pre-run checks | counts located; 7 vs 6 libraries; `library`/`tissue` skipped as aliases; with `--donor-key patient`: 5 patients with two libraries, 4 patients in both conditions (matched pairs); CA9 separates all 13 libraries; NDUFA4L2 lists two lower-content tumour libraries; EPCAM lists libraries in both directions (EPCAM is lost in ccRCC vs tubule: a bad label marker, which is the point of letting the researcher choose it); PTPRC declined as uninformative |
| CAA pre-run checks | 3 condition levels (by design); 10 vs 3 reference imbalance; 2 dogs in both conditions (tumour + matched gingiva); BRAF_status strongly associated with condition (biology: CAA is BRAF-driven); breed with 15 levels for 19 samples flagged as "too many levels to estimate" after the fix; PEG3 (the paper's diagnostic marker, not the target) separates CAA from both references in every library; ERBB4 (the target, run deliberately) lists two CAA libraries at 31% and one OSCC at 36% |
| synthetic `marker` mode calibration | **FAILED** (no state, no nominations) → protocol 23 / 23b, see below |

## Issue 1 (fixed, generic): expression route ran on a fragmented clustering — protocol 23

Scripted route: flat composition → subcluster → re-cluster at resolution 2.0 → `disease_programs` on
the 15 fragments (3 strong genes each). Fix: `base_cluster_key` tracked in facts; `disease_programs`
defaults to it; scripted route passes no key; analyst prompt v7 says why. Regression test added.

## Issue 2 (fixed, generic): same-lineage clusters competed in the relative criterion — protocol 23b

On the base clustering both keratinocyte clusters carried the programme (9 and 6 strong genes) and
"≥3× any other cluster" compared them against each other. Fix: compare against other lineages only.

## Issue 3 (fixed): `calibrate` defaulted to scoring v1 while `run` defaults to v3

The first calibration was silently scored under v1. `calibrate` now defaults to 3.

## Issue 4 (fixed in the skill's helper): covariates with ~one level per sample inflated Cramér's V

`design_check.py` now reports "too many levels for n samples" instead of warning.

## Issue 5 (skill text): whole-tissue label markers are uninformative in a non-malignant disease

IPF: MMP7, SPP1, KRT17 medians differed by <20 points across whole libraries, because the disease
shifts composition rather than replacing the tissue. Within the epithelial compartment MMP7 separates
IPF from control (36% vs 10%); within the immune compartment SPP1 does (29% vs 2%). Rule added to the
skill: in design rows 4-5 run the label check within the marker's compartment (`--celltype-key`).

## Issue 6 (to fix, generic): covariates known only from the deposit never reach `obs`

DOX: the GEO batch (two libraries processed with STARsolo, six with Cell Ranger) was not in the object,
so `design_check` could not see it. The sample-sheet builder carries every sheet column into `obs`, so
the fix is procedural: add such columns to `samples.csv` before building (done; object rebuilt).

## Issue 7 (to fix, generic): `inspect`'s compartment composition matched human symbols only

Empty composition line for the mouse DOX object. Fix: case-insensitive marker matching in `new_dataset.py`.

## Datasets built

| dataset | object | design | key | pre-registered arms |
|---|---|---|---|---|
| IPF, GSE135893 | `data/ipf/ipf.h5ad`, 89,326 cells, 22 donors (libraries pooled) | row 4; IPF 12 vs Control 10; other ILDs dropped at build | `ipf-key-frozen`: states = aberrant basaloid, SPP1 macrophages, CTHRC1 fibroblasts; target ITGB6; 20 acceptable; 19 decoys; primary metric = state recovery | `ipf_vanderbilt.yaml` (site), `ipf_no_vuild59.yaml` (control-like immune compartment in the largest IPF donor) |
| DOX, GSE285475 WT | `data/dox/dox.h5ad`, 75,809 cells, 4 vs 4 animals | row 5; `objective: marker` | `dox-key-frozen`: state = p53 response programme; target Cdkn1a; acceptable p53/atrophy/Nrf2 genes; decoys incl. Tfeb, cycle, IEG, sex, lineage; pre-run risk recorded (Cdkn1a+ fraction lower in DOX libraries) | `dox_batch2.yaml` (batch) |

## Runs

| run | arm | outcome |
|---|---|---|
| `synthetic_program_calibrate_v04` | before fixes | FAIL: no state |
| `synthetic_program_calibrate_v04_fix23` | base clustering only | FAIL: 9 and 6 strong genes in two keratinocyte clusters, neither 3x the other |
| `synthetic_program_calibrate_v04_fix23b` | both fixes | (pending) |
| `ipf_calibrate_1` | primary, scripted | (running) |

## Issue 8 (fixed, generic): programme-basis states were mis-scored — protocol 23d

Attribution classed induced programme genes `shared_with_reference` (reference cells of the same
cluster express them at baseline); consistency was detection-based and the `< 0.67` cut-off flagged
2 of 3 samples; the scripted nominator nominated one gene per cluster of the same lineage. Fixed:
basis-aware attribution from the within-cluster contrast, per-sample induction as consistency,
exact two-thirds, one nomination per gene. Synthetic `marker`-mode calibration now PASSES (ERBB4
rank 1, three unique nominations, consistency 1.0). Tumour paths untouched.

## IPF scripted floor (`runs/ipf_calibrate_1`, 11 min, 96 tool calls, 0 errors) — read 2026-11-05

Eval: all three key states "missed" by marker overlap; ITGB6 unranked; 18 nominations, precision 0;
no decoy nominated. Offline comparison with the authors' annotation (evaluation only; the column is
forbidden to the model):

| cluster | cells | IPF % | what it is (authors) | floor's treatment |
|---|---|---|---|---|
| 1 | 1,669 | 95 | basal cells (KRT5/KRT17/TP63): bronchiolisation | composition state; MMP1, TNC, ITGA2, PLAU, PTPRZ1 nominated |
| 2 | 5,559 | 77 | SCGB3A2+ / transitional AT2, and 341 of the 374 KRT5-/KRT17+ aberrant basaloid cells | programme state (CDKN2A, MMP7, KRT17, SOX9, CDKN1A, SOX4, PTGS2 all significant within the cluster); never subclustered, so the aberrant basaloid state stays hidden inside it |
| 4 | 7,127 | 42 | AT2 | programme state (MMP7, CDKN2A, ITGB6 log2 0.88 padj 0.011, SPP1 3.7) |
| 13 | 300 | 100 | AT2, IPF-only | composition state; CACNA2D2, ALOX15B, CSF3R, LRRK2 nominated |
| 0 | 2,533 | 85 | MUC5B+ secretory | enriched, NOT characterised (cap) |
| 15, 18 | 969, 4,305 | 87, 83 | macrophages (the SPP1+ expansion) | enriched, NOT characterised (cap of 3 composition states, all taken by epithelium) |
| 3+12+14 | 14,941 | 85 | ciliated | merged composition state; critic: enrichment does not hold within the lineage (log2 0.30); its five nominations all flagged shared_with_reference |
| 11 | 2,659 | 87 | PLIN2+ / myofibroblasts | log2 2.3 but MWU p 0.20 (donor heterogeneity): not enriched |

Library check (critic): S17 (IPF) immune compartment control-like with 6,276 cells — the same donor
the pre-run compartment check flagged (VUILD59); pre-registered arm exists. Also S20 (IPF) in
endothelial and immune, S15 and S03/S11 (controls) in epithelial/mesenchymal with small cell counts.

Reading: the composition route dominates in IPF, as the corrected design table says (basal, MUC5B+,
ciliated, an IPF-only AT2 state, two macrophage clusters); the expression route finds the
senescence / aberrant-basaloid programme spread over the transitional and AT2 compartments. The
floor's two generic weaknesses: (a) it characterises the top 3 enriched clusters overall, so in a
multi-lineage disease the immune states never get a marker table (eval cannot even match them);
(b) it never subclusters an enriched or programme cluster, so a rare state inside a programme
cluster is invisible. (a) is fixed as protocol 23e (lineage coverage); (b) is left to the agent, which
has `subcluster`, and is a stated floor limitation. ITGB6 reaches significance only in the AT2
programme at log2 0.88 (below the 1.0 strong-gene threshold): the key's expected target is demanding
at compartment resolution and will probably need the aberrant basaloid subcluster to rank.

## IPF agent run 1 (`runs/ipf_agent_1`, prompts v7, scoring v3, Claude Code backend) — 14.5 min, 1.70M input tokens, 72 tool calls, 2 recovered model errors

Eval: three key states missed (basaloid best match cluster 1 at 2/5: KRT17, MMP7); ITGB6 unranked;
8 unique nominations (MMP1, MMP7, TNC, ITGA2 from the basal cluster; SCTR, AGTR2 from the IPF-only
AT2 cluster; MUC5B, GPR87 from the secretory cluster), precision 0 against the key, no decoy
nominated. Critic: cluster 0 is a composition shift (within-lineage log2 0.54), S15 (control) looks
disease-like in epithelium, cluster 13 not distinct from normal AT2, MMP7 pan-epithelial.

Reading against the authors' annotation (evaluation only): the agent found the same five enriched
clusters as the floor (epithelial 1, 0, 13; immune 15, 18), recorded only the three epithelial
ones and never characterised the macrophage clusters — the same lineage-coverage gap fixed for the
scripted floor in 23e, now observed in the model. It never ran `disease_programs` (the prompt made it
optional when composition had changed), and never subclustered cluster 2, so the KRT5-/KRT17+
aberrant basaloid cells (341 of 374 inside cluster 2) stayed hidden; what it called "aberrant
basal-like" (cluster 1) is the basal-cell expansion of bronchiolisation, real IPF biology with the
wrong name. Its nominations are defensible biology (MMP1/MMP7/TNC in IPF basal/epithelium are
published) and the key's acceptable list did not contain them; they are novel candidates, not
errors. -> prompts v8: lineage coverage, run the expression route once in every analysis, and
subcluster a programme cluster of the case lineage when a rare state is plausible.

## DOX scripted floor (`runs/dox_calibrate_1`) — see protocol 23f

Composition flat; one programme (cardiomyocyte cluster 23: Trim63, Nppb, Fbxo32, Nppa significant)
found but its genes fell outside the statistic-ranked pool; Cdkn1a not induced at this time point
(key-drafting lesson, key unchanged). Rerun as `dox_calibrate_2` with the 23f pool.

## IPF agent run 2 (`runs/ipf_agent_2`, prompts v8, scoring v3) — 22 min, 4.10M input tokens, 119 tool calls, 1 recovered model error

Eval: **SPP1 profibrotic macrophage state FOUND (lenient 9/5, strict 6/5; cluster 18; SPP1 ranked 3rd,
CHI3L1 28th)**; aberrant basaloid still missed (best match cluster 1 at 2/5); ITGB6 unranked; 14
unique nominations (MMP1, MMP7, TNC, ITGA2 | SPP1, MMP9, CTSK, PLA2G7 | FN1 | LRRK2 | MUC5B | TRPV4,
PLAU, GPR87), precision 0.07 (SPP1), no decoy nominated. v8 changed what it was meant to change:
five composition states across two lineages (1, 13, 0 epithelial; 18, 15 immune) and the expression
route run once (7 programme clusters; cluster 2 recorded as a programme). It did NOT subcluster
cluster 2 despite the new rule, so the rare KRT5-/KRT17+ state remains inside it. Critic: cluster 0
blocked as a composition shift; cluster 15 and 13 identity genes shared with reference; immune
libraries S17/S20 resemble control (the pre-registered VUILD59 arm exists; the macrophage state was
recovered with them in, so the arm is not triggered). Cost note: v8 roughly doubles the analyst's
tool calls; 4.1M tokens per run on 89k cells.

Decision: no further prompt change on this dataset (the subclustering rule exists in v8 and was not
followed once; n = 1 is not a pattern); the gap is recorded as the IPF result's limitation.

## DOX floor re-runs 2-4 (protocol 23f and addenda)

`dox_calibrate_2` (pool rule in `_run_pseudobulk` only): unchanged — `disease_programs` rebuilt the
programme list from its own statistic head and the ambient step truncated the stored pool to it.
`dox_calibrate_3` (programme list follows the pool rule; ambient step removes suspects only): the
pool now holds Nppb, Ankrd1, Myh7, but the scripted nominator took the first 20 pool genes, all
statistic-ranked. `dox_calibrate_4` (nominator takes the statistic head plus every effect gene):
Nppb ranked 12th (score 0.22; `weak`: detected in < 10% of the cluster's disease cells), Myh7
21st, Ankrd1 28th; Ciart, Ces1d, Ndufa6, Mb lead. The floor now reports what the data hold: a
cardiomyocyte metabolic / structural programme with natriuretic and atrophy genes present but
sparse at two weeks; Cdkn1a not induced. Three touches to reach one pool rule end-to-end is the
cost of a pool that was built in three places; all three now share `_candidate_pool`.

## DOX agent run 1 (`runs/dox_agent_1`, prompts v8, scoring v3, objective marker) — 23 min, 1.17M input tokens, 63 tool calls, 0 errors

First perturbation design through the model. Composition flat (as predicted); expression route run;
the cardiomyocyte DOX programme recorded (cluster 20, L8: Mb, Cryab, Cish, Lmod2, Rorc, Irs2, Ankrd1,
Nppb, Txnip; 3 DOX vs 3 vehicle libraries because one DOX library has no cells in the cluster), plus a
second cardiomyocyte cluster recorded as a composition state and then rejected by the critic
(enrichment driven by 3 of 4 DOX libraries; within-lineage log2 0.91, p 0.2). Nppb ranked 12th and
Ankrd1 9th within the programme state, Myh7 9th in the other; Cdkn1a not induced (key-drafting
lesson, unchanged). Nominations Irs2, Rorc, Akt1, Ptgds, Dgat2, all flagged by the critic as sparse,
disease_wide or not real druggability; no decoy. The critic found S08 (a DOX library) vehicle-like
in lineage L4 and recommended an exclusion arm: the pre-registered batch-2 arm
(`config/dox_batch2.yaml`) answers it and is run as a free floor (`dox_batch2_calibrate_1`).

Reading: on a perturbation the agent now does what the design table says it should (programme
route, lineage-level response, replicate caveats); what it cannot do is rescue a key whose
expected read-out is not in the data at the sampled time. The honest DOX result is "state:
cardiomyocyte stress/atrophy programme found; expected target absent at 2 weeks; acceptable
natriuretic / atrophy genes present at ranks 9-12 of their state".

## DOX pre-registered batch-2 arm, scripted floor (`runs/dox_batch2_calibrate_1`, 3 vs 3)

Without the two batch-1 libraries the expression route flags eight programme clusters (the
batch-1 libraries inflated the reference dispersion); the cardiomyocyte programme (cluster 17:
Myl3, Myh6, Tnnt2, Hspb7, Pln, Csrp3 ...) is recorded and **Trim63 is nominated** (acceptable gene;
precision 0.14), with Bdh1, Macrod1, Lonrf1, Tspan4, Tsc22d3, Selenbp1. The arm the critic asked for
supports the atrophy response and shows the batch cost. Cdkn1a still absent.

## ccRCC regression run 1 (`runs/ccrcc_agent_v04_1`, prompts v8, scoring v3) — 11.5 min, 2.19M tokens, 87 calls, 0 errors

Tumour state found (lenient 9/5, strict 8/5); **CA9 ranked 1st in both tumour states** (0.91,
0.92; EGLN3 2nd and flagged for direction again); no decoy nominated; CD70 and FOLH1 above the key
as before; five states recorded (tumour x2, TAM, tumour endothelium programme, stroma programme).
The evaluator's merged rank put CA9 6th behind endothelium/stroma genes — a metric artefact of
lineage coverage, fixed as protocol 23g (within-state rank; merged rank kept as
`target_rank_global`). Cost 2.19M vs 0.6-1.7M in earlier rounds: v8's coverage and the expression
route add calls.

## ccRCC regression run 2 (`runs/ccrcc_agent_v04_2`) and stability — 11.3 min, 2.88M tokens, 89 calls, 0 errors

State found (lenient 8/5, strict 8/5; merged cluster 0+1+2+3); CA9 1st within its state, VEGFA 9th;
16 nominations, no decoy. Stability over the two v8 runs: top-10 Jaccard 0.33, nominated-set
Jaccard 0.42, rank-1 gene CA9 in both runs' tumour state (FOLH1 in both runs' endothelium
programme). The tumour result is unchanged by the round's seven method changes; the variation is in
secondary states and their candidates.

## Round closed: tag results-v4 (working copy only)

# Round 0.4b — the three earlier tumour datasets under skill 0.4 / prompts v8 (HCC, ovarian, FLC)

Started 2026-11-05 15:00 after `results-v4`. The HCC and ovarian objects had been deleted from the
machine, so both were rebuilt from GEO with the repository's own builders; the FLC object (the
lab's own Parse snRNA-seq, 354,510 nuclei x 41,352 genes, 17.7 GB) was copied from the test clone.
Keys unchanged from their frozen versions (HCC, ovarian: identical to the original clone; FLC:
verbatim copy of the test clone's `flc-key-frozen` key, committed here as `flc-key-frozen-v04`
before any run). Pre-registered arms carried over unchanged: `hcc_nop4` (drop patient 4),
`ovarian_noleak` (drop GSM8717127, GSM8717133), `flc_clean` (drop FLC83_T, FLC57_T, FLC26_ICBQ).

## Datasets rebuilt (validation_0.4/{hcc,ovarian,flc}/)

- **HCC** GSE149614, `build_from_matrix.py`: 63,101 cells x 25,712 genes, 10 tumour / 8 non-tumour
  libraries, counts integer in X. Identical dimensions to the round-0.2 object.
- **Ovarian** GSE286057, `fetch_geo_10x.py --drop peritumoral` + `build_from_10x_samples.py`:
  123,889 cells x 36,591 genes, 6 tumour / 23 omentum libraries, 29 GSMs as before. The round-0.3
  write-up quotes 145k cells for the same 29 libraries; the current deposit yields 124k with the
  same builder, so the earlier figure is either a pre-filter barcode count or an older file set.
  The sheet now carries `region` and `donor_group` so the design check can see them.
- **FLC** (lab data): design row 1, tumour vs tissue of origin; `tumor_status` has three levels,
  the 7 metastatic libraries are excluded by config as the declared second design.

## Pre-run checks

- **HCC design check:** paired (8 donors in both conditions), stage and virus not associated with
  the condition, no warnings. **Label check** (hepatocyte class; AKR1B10, SPINK1, MDK, GPC3,
  CYP3A4): HCC04N again expresses the tumour programme (GPC3 77%, SPINK1 91%, AKR1B10 60%) while
  HCC04T is GPC3-negative (2%): the round-0.2 anomaly reproduced on the rebuilt object. Kept as
  deposited; `hcc_nop4` is the arm.
- **Ovarian design check:** 23 vs 6 imbalance; five control donors contribute 2-3 omentum regions
  (pseudo-replication in the reference, stated, library kept as the sample unit so the frozen key's
  predictions apply); `region` nested within the condition (a sub-labelling: tumour vs omentum
  regions, not a batch). **Label check** (EPCAM, PAX8, KRT7, WFDC2, MSLN, whole library): every
  marker "uninformative" by the median-gap rule because tumour cells are a minority of tumour
  libraries; the new second rule (issue 9) lists GSM8717133 (P14_distant, EPCAM 10%), GSM8717127
  (P12_distant, EPCAM 2.5%) and GSM8717110 (P4 lesser omentum, PAX8 5%) - the first two are the
  pre-registered `ovarian_noleak` exclusions, the third is the healthy donor whose benign
  Müllerian-type inclusion blocked FOLR1 in round 0.3 (protocol 21), now visible before any run.
- **FLC design check:** 9 primary / 9 normal / 7 metastatic libraries, counts integer; no donor
  column in `obs` (pairing is in the library names: 57, 82, 83, 84, 89, 96 and 26, 34 have a T and
  an N library). **Label check** (CA12, SLC16A14, KRT7, CYP3A4, ALB, PRKACA): reproduces the test
  clone's pre-run notes - FLC83_T and FLC57_T express the tumour markers at normal levels with
  CYP3A4 at 52-70% (low tumour content), FLC26_ICBQ (Normal) has SLC16A14 in 40% of nuclei, FCF82_N
  has little hepatocyte signal; FLC24_MOYV and FLC84_T are CA12-low but SLC16A14-high (tumour
  heterogeneity, kept). Decisions as written in the test clone on 2026-11-02: all libraries in the
  primary run; `flc_clean` is the arm.

## Issue 9 (fixed in the skill's helper, generic): the label check was blind to a minority compartment

`label_check.py` compared each library with the two condition medians only for markers whose
medians differed by 20 points. When the case cells are a minority of every library (tumour nests in
omentum) no marker reaches the gap and the script reported nothing, although the per-library table
showed a reference library at 10% EPCAM against a 0.2% median. Second rule added, independent of
the gap: a library carrying the other condition's marker at >= 3x its own condition's median and
>= 2 points above it is listed (`--outlier-fold`, `--outlier-gap`). Verified on ovarian (finds the
two pre-registered exclusions and the benign inclusion), HCC (HCC04N; also lists small tumour
libraries with normal-hepatocyte CYP3A4, which is the expected admixture) and FLC (FLC26_ICBQ).

## Issue 10 (fixed in the skill's helper, generic): "perfectly confounded" for a sub-labelling of the condition

`design_check.py` called `region` (tumour / proximal / distal / lesser omentum / distant) perfectly
confounded with the condition and said any result may be a region effect. Every level sits inside
one condition because the column refines the condition; that is harmless when the levels are
biological sub-groups and confounding only when they are processing batches. The script now
reports such a column as NESTED WITHIN CONDITION with that distinction, and keeps PERFECTLY
CONFOUNDED for a column with as many levels as conditions.

## HCC scripted floors (`runs/hcc_calibrate_v04_1`, 10 min, 127 calls; `runs/hcc_nop4_calibrate_v04_1`, 7.5 min, 125 calls)

All libraries: six states recorded (three disease-only hepatocyte clusters 14, 15, 18 each present
in 10/10 patients, a myeloid state, an endothelial and a T/NK programme); key state 0/5 as in every
earlier HCC floor (the key's literature markers are patient-heterogeneous, the recorded limitation
of that key); **GPC3 unscored**: log2FC 2.6 in the whole-lineage contrast but padj 0.20 with
HCC04N in the reference, and ~0 in each hepatocyte cluster's sibling contrast because the siblings
are the other tumour hepatocyte clusters (the floor merged 13+22 only). Two proliferation decoys
(TOP2A, MKI67) ranked in the cycling hepatocyte cluster, nominated by neither role. Same picture as
round 0.2.

Without patient 4 (pre-registered arm): **GPC3 ranked 1st within its state** (cluster 13; global
8th), recall@5 1.0, codrituzumab found, nothing above it, no decoy nominated; key state again
missed by markers (1/5: GPC3 itself). Round 0.2's floor had GPC3 4th in this arm; the within-state
rank (23g) and the effect pool (23f) account for the difference. The HCC04N anomaly remains the
single determinant of whether the expected target is scored at all.

## HCC agent run 1, all libraries (`runs/hcc_agent_v04_1`, prompts v8, scoring v3) — 15.5 min, 3.61M input tokens, 111 calls, 0 tool errors -> protocol 23i

Six states recorded across hepatocyte (cycling tumour cluster 15; NDUFA4L2/SQLE/MAGEA cluster 14),
myeloid (C1Q+ lipid-associated macrophages), endothelial programme (THY1/ESM1), fibroblast
(CAF-like cluster 12) and T/NK (Treg-like, sub-threshold). The agent merged 13+22 (normal
hepatocytes) and ran `check_reference_libraries` twice; the critic's eleven flags are substantive
(cycling genes are not targets; cluster 15 genes are disease_wide, not state-specific; cluster 12's
enrichment is a lineage-composition effect, log2 -0.05 within the fibroblast lineage; T/NK tumour
libraries S01/S12 resemble normal; hepatocyte S08 flagged). Nominations: 15 genes, CYP17A1 and SQLE
first, GPC3 13th; codrituzumab found; no decoy nominated.

**The eval line said GPC3 2nd within its state and recall@5 = 1.0. The state was cluster 12, the
CAF-like stromal state** (GPC3 `disease_wide`, reached through the disease-wide contrast), not a
tumour hepatocyte state: in clusters 14 and 15 GPC3 was unscored as in every all-library HCC run
(lineage contrast padj 0.20 with HCC04N). The number was right by rule 23g and misleading as a
reading. Protocol 23i: the evaluator now prints the carrying state's dominant lineage
(`'12 (Fibroblast)'`) and honours an optional `lineage:` on a key's expected state, under which an
off-lineage rank is reported and not credited. Frozen keys have no such field, so their numbers
stand; this run is read as "GPC3 scored in a stromal state, unscored in the tumour hepatocyte
states". Compared with round 0.2 (GPC3 unscored everywhere with all libraries), the change is that
v8's lineage coverage characterised a stromal state in which the pan-tumour gene is detectable.

## Ovarian scripted floor (`runs/ovarian_calibrate_v04_1`, 16.5 min, 131 calls) -> protocol 23h

Tumour epithelium found (cluster 6, lenient 5/5; strict 1/5), EPCAM 2nd, CLDN6 23rd, ERBB2 57th;
FOLR1 **unscored** although significant in the sibling contrast (log2FC 4.8, padj 5e-8; 39th by
statistic, 9th by effect among 1,834 significant genes). Cause: the floor's pool read the 20 most
significant rows plus the rows the tool had labelled `effect`; a gene at statistic rank 21-40 was
in neither set. Fixed as protocol 23h (effect set taken over the whole stored list; test added);
floor re-run as `ovarian_calibrate_v04_2`.

**Re-run (`runs/ovarian_calibrate_v04_2`, 13 min):** state found 5/5 again; **FOLR1 now scored,
24th of 111 within the tumour state** (score 0.47, class state_specific, flag
`shared_with_reference`; global 58th), 22 genes above it all epithelial tumour genes (MMP7,
CDKN2A, WFDC2, RAB25 ...), none a decoy; EPCAM 2nd, CLDN6 25th, ERBB2 57th. The floor reaches
the gene for the first time; the `shared_with_reference` flag is the round-0.3 attribution
limitation (benign FOLR1+ epithelium in the omentum reference, protocol 21b) and is what the
pre-registered `ovarian_noleak` arm tests.

## Ovarian agent run 1, all libraries (`runs/ovarian_agent_v04_1`, prompts v8, scoring v3) — 18 min, 3.91M input tokens, 112 calls, 0 errors

Tumour epithelium found (lenient 6/5 including FOLR1; strict 4/5, 20 genes recorded); seven states
processed (five composition, one per lineage, plus fibroblast and T/NK programmes). **FOLR1 14th
within the tumour state** (round 0.3: 16th), classed state_specific with the
`shared_with_reference` flag, not nominated; EPCAM 2nd, ERBB2 15th; 15 nominations led by MMP7,
EPCAM, SLC34A2; no decoy. The critic's reading is the one the skill asks for: it called FOLR1 "the
best-supported druggable gene" (mirvetuximab soravtansine, live DGIdb), attributed the
shared-with-reference class to the tumour-like reference libraries S23/S27 that
`check_reference_libraries` flagged in up to six lineages (the same GSM8717127/GSM8717133 that the
pre-run label check and the `ovarian_noleak` arm name), said FOLR1 "needs the exclusion arm before
it is called tumour-specific", and recorded that the arm had not been run - rather than dropping
the gene. Pre-registered arm run as `ovarian_noleak_calibrate_v04_1` and `ovarian_noleak_agent_v04_1`.

## Ovarian no-leak arm, scripted floor (`runs/ovarian_noleak_calibrate_v04_1`, 14 min)

Without GSM8717127/GSM8717133: tumour epithelium found (cluster 24, lenient 6/5 incl. FOLR1;
strict 2/5); **FOLR1 35th of 83 within its state** (score 0.17, state_specific, flag
`reference_driven`, sample consistency 1.0; round 0.3's no-leak floor: 96th), EPCAM 4th, TACSTD2
18th; the top of the state is disease_wide epithelial genes (DEFB1, ASS1, MMP7). No decoy. The
`reference_driven` flag is the benign FOLR1+ Müllerian-type inclusion of a healthy donor
(protocol 21, 21b) - the library the new label-check rule lists as GSM8717110 (PAX8 5%) - which the
arm by construction does not remove. The limitation stands exactly as frozen in 0.3.

## Ovarian no-leak arm, agent (`runs/ovarian_noleak_agent_v04_1`, prompts v8) — 19 min, 3.47M input tokens, 105 calls, 0 errors

Tumour epithelium found **strict 5/5** (as in both 0.3 no-leak agent runs); **FOLR1 17th within its
state** (0.3: unranked after being set aside), EPCAM 2nd, TACSTD2 12th; 13 nominations led by MMP7
and EPCAM, no decoy. The critic's blocking flag is the 0.3 verdict again: FOLR1 `reference_driven`
(151 reference cells in the cluster, mean 5.5 vs 1.4 in disease cells; "an epithelial/
mesothelial-like population present in normal omentum ... e.g. mesothelial or epithelial
inclusions"). The critic identified the reference cells as a population and read the gene as not
tumour-exclusive, which is the 0.4 attribution reading rule halfway: it named the cells but still
withheld the candidate rather than reporting it as conditional on them. The 0.4 skill text says
"reported as conditional, not dropped"; the critic prompt (v8) does not yet say this in those
words. Noted as a prompt-text gap, not acted on in this round (one dataset, the known case).

Across the three ovarian runs of this round the picture is unchanged from 0.3 and stated more
precisely: FOLR1 is found in the right state every time (14th / 24th / 35th / 17th), is called
druggable, and is held back by 59-151 benign FOLR1+ reference epithelial cells that the new
label-check rule now shows before the run (GSM8717110).

## Ovarian no-leak arm, agent under prompts v9 (`runs/ovarian_noleak_agent_v04_2`) — 17 min, 2.83M input tokens, 93 calls, 0 errors -> protocol 23j outcome

Strict 5/5 again; FOLR1 17th within the state again (global 97th), EPCAM 3rd, TACSTD2 11th; six
nominations (CTLA4, IL2RA, MMP7, LOXL2, IL3RA, TNFSF13B; precision 0), no decoy. The critic applied
v9 as written and reached the **representative** branch: "normal omentum carries cluster 24 at up
to 0.46 of epithelium per library, mean 0.126, so the reference holds this epithelial population
in bulk and the lineage counterpart is representative" - blocking kept for FOLR1 and CD24. In the
same review it wrote that the verdict "depends on reference cell identity", that the 151 cells
"may be normal fallopian-tube-like or mesothelial epithelium, or tumour-contaminated", and that "a
reference-only subcluster of cluster 24 would settle this". So v9 produced the reasoning the skill
asks for (cells named, settling experiment stated) but not the conditional nomination, because the
representativeness test is answered by the derived lineage L1, which absorbs mesothelium at this
resolution (protocol 21b, the exact property that defeated the lineage-relative guard). **v9 stays**
(it is the correct rule and changed nothing else), the FOLR1 outcome is unchanged across five
agent runs on this dataset, and the limitation is restated: with derived lineages, "the reference
holds the lineage" cannot distinguish mesothelium from benign epithelium. The one step the critic
named and did not take - a reference-only subcluster of the tumour cluster - is the generic
settling move for this design and belongs in the skill as the next thing to try when the flag is
representative-but-doubtful.

## FLC scripted floor (`runs/flc_calibrate_1`, 46 min, 121 calls, 0 errors)

Primary vs normal liver, 9 vs 9 libraries, metastatic set aside (254k nuclei through QC and
Scrublet; the floor's slowest run). Tumour state recorded (cluster 6, derived lineage L1), key
markers 3/5 (CA12, KRT7, PRKACA; five needed, so "missed" by the key's own threshold);
**PRKACA 2nd within its state behind CA12** (both key answers; global 3rd), recall@5 1.0,
nothing above them, KRT7 7th, no decoy ranked or nominated; drugs: none of the trialled agents
(no PRKACA inhibitor exists; ENMD-2076/neratinib/everolimus not reached by the floor). The v6
floor in the test clone (2026-11-02) had PRKACA 6th with four novel genes above it and markers 1/5.

Worth recording: the tumour state came through the **expression route**, not composition. Cluster 6
(hepatocyte-like lineage L1) is present in both conditions with no abundance shift; the within-
cluster tumour-vs-normal contrast has 2,009 genes up at log2FC >= 1 (CA12, KRT7, CD24, IRF4 ...).
With one condition per library and Harmony correcting on the library, part of the tumour-vs-normal
difference in the hepatocyte lineage was integrated away and FLC tumour nuclei co-clustered with
normal hepatocytes. The 0.4 rule that the expression route runs in every analysis (23, 23e) is what
reached the state; the v6 floor, which ran it only as a fallback on a fragmented clustering,
missed it. A design in which every library is one condition cannot separate library from condition
in the integration; this is the general reason the expression route must always run, and it should
be said in the design-types reference.

## FLC agent run 1 (`runs/flc_agent_1`, prompts v8, scoring v3) — 37 min, 2.48M input tokens, 95 calls, 0 errors -> protocol 23l

Five states: tumour-lineage epithelial cluster 0 (L0), disease-enriched hepatocyte cluster 10,
tumour endothelium (cluster 20), and two programme states - **cluster 6, the hepatocyte/biliary-
like programme in a shared L1 cluster, which is where the key's answers sit**, and a CAF-like
matrix programme (cluster 8). **PRKACA 2nd within its state behind CA12; ERBB2 2nd among
acceptable (neratinib found); PRKACA nominated** (rank 8 of 11 nominations, "kinase, clinically
actionable, H-89/fasudil"); recall@5 1.0, nothing above the key answers, KRT7 7th; precision 0.27
(three key genes of eleven). The v6 run in the test clone had PRKACA unranked, 89 genes above the
first acceptable gene, and PDGFRB unflagged. Key state again "missed" at 3/5 markers (AKAP12,
CA12, TESC in cluster 10; the key needs 5 of 14 literature markers, which in this dataset split
across the composition state and the programme state).

**One decoy nominated unflagged: PDGFRB**, from the CAF-like fibroblast programme, recorded as "the
most druggable (14 approved inhibitors, but non-selective)". The critic reviewed attribution and
cycling and never asked which compartment the gene is in. Prompts v11 (23l): every nominated gene
from a non-disease-cell lineage carries a microenvironment flag naming the lineage. FLC agent run 2
under v11 serves as verification and as the dataset's stability pair.

## FLC agent run 2 (`runs/flc_agent_2`, prompts v11) — 43 min, 3.53M input tokens, 111 calls, 0 errors -> 23l verified; stability pair

**PRKACA 2nd within its state behind CA12 again** (global 3rd), nominated again; ERBB2 2nd,
KRT7 6th, neratinib found; 10 nominations, precision 0.30; **no decoy unflagged**: the critic
applied v11 and recorded "Microenvironment flag: clusters 20 and 8+19 are not tumour-cell states"
on VWF (endothelium) and on FAP and COL1A2 (CAF matrix programme), which is what the rule asks
for; PDGFRB was not nominated this time. The critic also tested cluster 0's enrichment within its
lineage (log2 0.69, p 0.18, does not hold) and flagged S22 as Normal-like in L0. Stability over
the two agent runs: nominated-set Jaccard 0.62, top-10 Jaccard 0.54; the tumour programme state's
top five are CA12, PRKACA, SLC7A5 in the same order in both runs (then HKDC1/IRF4 vs IRF4/MUC13).
Key state still 3/5 by markers in both runs (the key's literature markers split between the
composition and the programme state; a limitation of that key, recorded on 2026-11-02).

## Round 0.4b closed: tag results-v5 (working copy only)

Runs this round: 7 scripted floors, 7 agent runs (HCC 1, ovarian 3, FLC 2, plus the two FLC/ovarian
prompt verifications counted above), ~23M input tokens. Method changes: protocol 23h (floor pool),
23i (evaluator lineage), 23j (prompts v9), 23k (`subcluster(condition=...)`, prompts v10), 23l
(prompts v11); skill 0.4.2 helper fixes (issues 9-11) and text. Tests: 93 passed, 1 skipped.
Outcomes: HCC GPC3 1st in its state without patient 4, unscored in the tumour states with all
libraries (HCC04N), scored in a stromal state by the agent; ovarian FOLR1 found in the right state
in every run (14th-35th), held back by benign reference epithelium in every run, now visible
before the run; FLC PRKACA 2nd behind CA12 in all three runs, nominated by the agent, neratinib
found, no decoy unflagged under v11. Originals untouched; nothing pushed.

## Issue 11 (fixed in the skill's helper, generic): pairing was reported as pseudo-replication

The donor warning counted libraries per donor, so a paired design (one tumour and one normal
library per patient) was reported as "8 donors contribute more than one library ... over-counts
replication". Pseudo-replication is several libraries of one donor in the same condition; the
count is now per donor-condition pair. HCC: no warning (correctly: pairing only); ovarian: 5
control donors with 2-3 regions each (correctly flagged).

# Round 0.4c: canine acanthomatous ameloblastoma (CAA) under the 0.4 skill — started 2026-11-08

The development dataset (lab snRNA-seq, 204,967 nuclei, 19 libraries: 6 CAA, 10 OSCC, 3 healthy
gingiva; two dogs contribute a tumour and a matched gingiva library). Key `data/ground_truth/caa.yaml`
unchanged since commit d1e2d45 (frozen before the first 0.3 run); tagged `caa-key-frozen-v04` at
HEAD before any 0.4 run. Config `config/caa.yaml` as in 0.3 (reused leiden_1.0 clustering,
Broad_anno lineages, refined_anno forbidden). Last run on this dataset: prompts v2 / scoring v2
(results-v3); 0.4 has never been run here.

**Pre-run checks (skill 0.4.2 helpers, re-run because the 0.4a logs predate issues 9-11):**
`design_check`: three condition levels by design (CAA vs OSCC + HGIN); imbalance 10 vs 3 in the
reference; Apollo and Micah paired (tumour + gingiva), no longer reported as pseudo-replication
after issue 11; BRAF_status strongly associated with condition (biology: CAA is BRAF V595E-driven,
not a batch); breed has 15 levels for 19 samples, not estimable. `label_check` with PEG3, BNC2,
RBFOX1 (top of the paper's CAA-specific signature; not the target), ERBB4 (the target, run
deliberately) and PTPRC: PEG3 and RBFOX1 separate CAA from both references in every library
(PEG3 59-87% vs 3-18%); the second rule lists nothing. ERBB4 alone lists sample_7 and sample_9
(CAA, 31%) and sample_20 (OSCC, 36%); with PEG3/RBFOX1 clean in those libraries this is target
heterogeneity, not a label problem. **Decision (before the run): all 19 libraries kept; no
exclusion arm.** The reference cannot be cleaned with the expected target.

## CAA scripted floor (`runs/caa_calibrate_v04_1`, 15 min, 73 calls, 0 errors) — PASS

Both CAA-exclusive epithelial clusters found by composition: cluster 6 (log2 4.4, 96% disease
cells) and cluster 12 (log2 3.9, 95%), each present in all six CAA libraries; the paper's two
CAA-specific keratinocyte states. Key state matched to cluster 12 at 9/5 markers (ERBB4, FAM13A,
GABRB1, MAGI2, NAALADL2, NEO1, NFIA, PEG3, SESN3). **ERBB4 5th within cluster 12** (global 7th),
behind GABRA1, GABRA4, GABRB1 (acceptable), GABRB2; three novel genes above it, no decoy ranked
or nominated; ERBB4 nominated (11 nominations, precision 0.18); drugs afatinib, dacomitinib,
neratinib (recall 0.75). Four programme states besides (clusters 4, 11 recorded; 0, 1, 10, 14, 8,
15 listed). The 0.3 scripted baseline on this dataset (results-v1) had ERBB4 unranked and 24
novel genes above the first acceptable gene; the best 0.3 agent runs (prompts v2) had ERBB4 1st-2nd.

The GABA-A receptor subunits at the top are the same biology as the key's GABRB1 (neuroepithelial
programme of CAA; Fig. 3A of the paper) and are druggable only by anaesthetics and benzodiazepine-
site ligands; they are honest "novel" entries, not decoys, and the within-state table shows ERBB4
as the only kinase among the top five.

**Marker-free library check (floor critic, major):** S13 = sample_3 (OSCC, 5,284 cells) correlates
with CAA libraries in its Epithelial and Mesenchymal pseudobulks. Confirmed from the object:
14.5% of its cells sit in the CAA-exclusive clusters 6+12, against 0.3-1.2% for every other OSCC
library (sample_20: 7.0%) and 4-6% for healthy gingiva. The pre-run PEG3 check did not list it
(10.7%, highest OSCC but far from CAA). Also flagged: S11 = sample_24 (HGIN) in Endothelial Cells,
79 cells, ignored. **Pre-registered second arm `config/caa_nos3.yaml`** (sample_3 excluded, reason
in the config), to be run as floor and agent run after the primary agent run; the all-libraries
run stays primary.

Report section 6 says "no key was scored" because the report is written before `run_eval`
runs; the eval lives in `eval.json` and the log. Same in every 0.4 run; cosmetic, noted.

## CAA agent run 1 (`runs/caa_agent_v04_1`, prompts v11, scoring v3) — 15.7 min, 1.62M input tokens, 71 calls, 0 errors — PASS

Three states: clusters 6 and 12 by composition (the two CAA-specific keratinocyte states, 96% and
95% disease cells, present in 6/6 CAA libraries, confirmed within-lineage by the critic at log2 3.96
and 3.81) and one CAA-wide epithelial programme over the shared clusters (4, 11, 0, 1, 14). Key state
matched to cluster 12 at 13/5 markers, **strict match 6/5** (the agent itself recorded BNC2, ERBB4,
GABRB1, MAGI2, NAALADL2, PEG3 ... among 19 genes). **ERBB4 4th within cluster 12** (global 6th)
behind GABRA4, GABRB1 (acceptable), GABRB2; nominated with "approved pan-ERBB inhibitors";
afatinib / dacomitinib / neratinib found (recall 0.75); 10 unique genes nominated, precision 0.20.
**One decoy nominated, FGFR2, flagged major by the critic** ("not specific to the CAA state;
present in 86% of reference cells in cluster 6") - the paper's own negative control (infigratinib
minimal response), correctly set aside. No decoy unflagged. The critic also flagged S13 (sample_3)
and asked for the exclusion arm, flagged the 10-vs-3 OSCC/HGIN reference composition ("CAA vs
normal gingiva cannot be separated from odontogenic vs squamous lineage"), and called the
GABA-A compound counts inflated by anaesthetics and benzodiazepine-site ligands. The best 0.3 agent
runs (prompts v2) had ERBB4 1st-2nd globally with 7-9 nominations; the 0.3 runs did not score
within-state rank, did not record both keratinocyte states as distinct, and ran no library check.

Two remarks in the agent's own notes pointed at tool defects, both generic, both fixed this round:

## Issue 12 (fixed, generic; prompts v12, protocol 23m): stage briefings dropped items past a flat character cap

The critic wrote "the candidate list I received was truncated after ERBB4, so any further candidates
were not reviewed". `stage_prompt` built the critic's briefing as one JSON dump cut at 9,000
characters (the nominator's sections at 6,000 each). Re-measured on every agent run of rounds 0.4
and 0.4b from `state.json` + `config_used.json`: the critic's claims ran 23-52k characters in all
of them. Nominations the pre-v12 critic was never briefed on (title not within the first 9,000
characters; the critic could still reach them through its own tool calls, and sometimes did):

| run | nominations | briefed | never briefed |
|---|---|---|---|
| caa_agent_v04_1 | 12 | 7 | GABRB2, GABRA1, EPHA3, CA10, GPC5 |
| ccrcc_agent_v04_1 | 12 | 5 | APOC1, GPNMB, SLCO2B1, FOLH1, ANGPT2, CD276, PPARG |
| ccrcc_agent_v04_2 | 16 | 6 | TREM2, CD14, KIT, NTRK1, HPGDS, FOLH1, ANGPT2, CD276, COL6A3, FN1 |
| dox_agent_1 | 5 | 5 | - |
| flc_agent_1 | 13 | 5 | CA5A, VWF, ITGA6, CA12, **PRKACA**, FAP, **PDGFRB**, INHBA |
| flc_agent_2 | 13 | 4 | CA5A, VWF, F2RL3, FAP, COL1A2, CA12, SLC7A5, PRKACA |
| hcc_agent_v04_1 | 17 | 4 | CYP17A1, SCD, SQLE, MMP9, CCL18, SPP1, FOLH1, CA2, PGF, FN1, **GPC3**, LOXL2, PKM |
| ipf_agent_1 | 9 | 8 | GPR87 |
| ipf_agent_2 | 16 | 4 | SPP1, MMP9, CTSK, PLA2G7, FN1, LRRK2, MMP7, MUC5B, TRPV4, TNC, PLAU, GPR87 |
| ovarian_agent_v04_1 | 15 | 5 | TNFRSF4, TNFRSF9, SPP1, CTSL, TDO2, LOXL2, SERPINE1, EGFL6, EIF2AK2, OAS1 |
| ovarian_noleak_agent_v04_1 | 13 | 3 | IL2RA, TIGIT, CCR8, TDO2, INHBA, FAP, LOXL2, TNFSF13B, CTSL, STAT1 |
| ovarian_noleak_agent_v04_2 | 6 | 2 | MMP7, LOXL2, IL3RA, TNFSF13B |

Consequences for earlier conclusions: PDGFRB in FLC run 1 (the decoy that motivated 23l / v11) was
never in the critic's briefing, so part of what 23l attributed to the critic's reasoning was a
briefing it did not have; in FLC run 2 the critic flagged VWF, FAP and COL1A2 that it was also not
briefed on, so it reached them through `ranked_targets`. The metric "decoys nominated unflagged"
for rounds 0.3-0.4b stands as recorded but does not measure the critic's rule-following on
nominations past the cap. The floor's scripted critic is unaffected (it reads the facts directly).

Fix (prompt construction only): an index line names every nomination and every state and is never
shortened; per-item free text is clipped progressively until the list fits; gene lists are
shortened evenly, never dropped by key; the critic is told every indexed gene must appear in its
review. Measured on the earlier runs under v12: nominator 12-27k characters, critic 17-35k, no
truncation, every nominated gene present. Tests `tests/test_prompt_briefing.py` (5). Verified
prospectively by CAA agent run 2 (v12) and the exclusion-arm agent run.

## Issue 13 (fixed, generic; protocol 23n): `map_orthologs` called an empty precomputed entry "no ortholog"

GABRB1 and GPC5 came back `no_ortholog` because `var[human_symbol]` was an empty string for them,
while `druggability_lookup` used the same symbols as human (92 interactions for GABRB1); the critic
flagged "druggability unverified" - a correct reading of what it was shown. An empty table entry
now falls through to Ensembl like a gene absent from the table, with `var_column_empty: true` in
the result. Test updated.

## Issue 14 (fixed, cosmetic; protocol 23o): report section 6 said "no key was scored" on every scored run

The report is written before `run_eval`; `run_eval` now rewrites section 6 of report.md/html from
eval.json. Tests `tests/test_report_eval_refresh.py` (2). Both finished CAA reports refreshed.

Full suite after the three fixes: 100 passed, 1 skipped (`validation_0.4/full_tests_round_0.4c.log`).

## CAA exclusion arm, scripted floor (`runs/caa_nos3_calibrate_v04_1`, 13.8 min, 78 calls, 0 errors) — PASS

sample_3 removed (18 libraries). Same two composition states; key state cluster 12 at 9/5 markers;
**ERBB4 4th within cluster 12** (5th with all libraries; global 8th) behind GABRA1, GABRA4, GABRB1;
GABRB2 no longer above it; nothing else moves (PEG3 18th, drugs 0.75, 14 nominations, no decoy).
The marker-free check now flags S12 and S16 (OSCC) and S11 (HGIN) in the Endothelial lineage only;
endothelium carries little condition structure in this design (the critic's own r values are near
zero), so these are the check's noise floor, not label problems. The exclusion changes the result
by one rank: sample_3 was not what held ERBB4 at 5th; the GABA-A subunits are.

## CAA agent run 2 (`runs/caa_agent_v04_2`, prompts v12, scoring v3) — 38 min, 1.57M input tokens, 79 calls, 0 errors — PASS; 23m verified; stability pair

Same three states; key state cluster 12 at 13/5 markers, **strict 5/5**; **ERBB4 4th within
cluster 12 again** (global 6th) behind GABRA4, GABRB1, GABRB2; nominated; afatinib / dacomitinib /
neratinib found. 7 nominations (precision 0.29; run 1 had 10 at 0.20): KCNH8, CACNA2D3, ERBB4,
CA10, GABRB1, EPHA3, GABRA1. No decoy nominated (FGFR2 was not nominated this time).

**23m verified:** the critic opened with "all 7 nominated genes have a recorded flag" and every
nominated gene has a flag line in section 4 (run 1's critic had written that its list was cut
after ERBB4). **23k used unprompted by the critic:** it ran `subcluster(clusters=[6],
condition=OSCC)` and characterised the two largest OSCC-only subclusters of cluster 6 (ECM/EMT-like
programme), concluding that the 1,237 OSCC cells inside cluster 6 carry the state signature and
that cluster 6 candidates (KCNH8, CACNA2D3) are conditional on what those cells are - the ovarian
settling move, applied here to a reference that is another tumour. It again asked for the S13
exclusion arm (run as `caa_nos3_*`), re-ran within-celltype abundance (cluster 6 log2 3.96, 12 3.81),
and flagged ERBB4 as CAA-epithelium-wide rather than cluster-12-specific with non-selective drugs
(the paper's reading too: the ERBB4 dependency is a property of the CAA keratinocyte programme).

Stability over the two primary agent runs: top-10 Jaccard 0.67, nominated-set Jaccard 0.42,
rank-1 gene GABRA1 in both; ERBB4 4th within cluster 12 in both.

**Issue 13 follow-up:** GPC5 now resolves (`ortholog_one2one` via Ensembl, `var_column_empty`).
GABRB1 still `no_ortholog`, now from Ensembl itself: the dog symbol resolves to the current-assembly
gene ENSCAFG00845008573 with no human homologue in Compara; the dataset's CanFam3.1 identifier is no
longer served (HTTP 500). An annotation gap, not biology; it will recur for any non-model species
whose reference assembly moved. Second generic fix (23n, part 2): when Ensembl reports no homologue,
the client checks whether the symbol is a human gene; if so the mapping returns it as
`same_symbol_human_gene`, `symbol_unverified: true`, listed under `assumed_or_unverified`.

## De novo arm requested (2026-11-08): no reused clustering, Harmony on the library, derived lineages

At the author's request the pre-registered independent arm (protocol section on `caa_denovo`:
object stripped to counts, sample, condition and Broad_anno; `config/caa_denovo.yaml`, Harmony on
`sample`, `lineages.ignore_annotation: true`, nothing reused from the paper's processing) is run
under the 0.4 skill: floor `caa_denovo_calibrate_v04_1`, then agent run(s). Last run under prompts
v2 / scoring v2 (`caa_denovo_2`: ERBB4 4th global; `caa_denovo_3`, derived lineages: ERBB4 1st).
Same key, same arms otherwise. With one condition per library in 17 of 19 libraries, Harmony on the
library can absorb part of the tumour-vs-reference difference (the FLC lesson, 0.4b); the
expression route runs in every analysis and is what would recover a state integrated away.

## CAA exclusion arm, agent run (`runs/caa_nos3_agent_v04_1`, prompts v12) — 15 min, 1.71M input tokens, 72 calls, 0 errors — PASS

sample_3 excluded. Key state cluster 12 at 12/5, **strict 5/5**; **ERBB4 3rd within cluster 12**
(global 4th), only GABRA4 and GABRB1 (acceptable) above it; nominated; drugs 0.75. Ten nominations
(precision 0.20) including **FGFR2, flagged major again** ("pan-epithelial / shared with reference,
not state-defining"; 10 approved pan-FGFR inhibitors, narrow window). The critic: "every nominated
gene was reviewed and recorded" (23m holds in a second run); within-celltype log2 5.50 / 4.07 for
clusters 6 / 12; the marker-free check now flags only endothelial libraries (S11, S12, S16), "so a
second arm would only test endothelial contrasts" - the arm removed what it was meant to remove.
GABRB1 again `no_ortholog`: this run started before the 23n part-2 fallback (human-symbol check)
was in place; the homology answer is cached, the human lookup is new and will apply from the next run.

Across the three CAA agent runs ERBB4 is 4th, 4th, 3rd within its state and nominated each time;
FGFR2 nominated in two of three and flagged both times; no decoy unflagged.

## CAA de novo floor (`runs/caa_denovo_calibrate_v04_1`, 42 min, 117 calls, 0 errors) — PASS

Recomputed HVGs / PCA / **Harmony on the library** / Leiden r1.0 (20 clusters), lineages derived
(coarse Leiden r0.2, 8 groups), nothing of the paper's processing reused. Five enriched clusters,
all epithelial by the (undisclosed) annotation: **16 (31k cells, 99.5% CAA), 17 (25.6k, 99.3%),
13 (3.1k, 98%)**, 3 (3.0k, 85% CAA) and 12 (2.6k, 74%). **Harmony did not absorb the CAA state**
(unlike FLC): CAA keratinocytes are different enough from OSCC and gingival keratinocytes that
clusters 13, 16 and 17 form their own derived lineage L4 (59.5k cells; 69 HGIN and 576 OSCC cells),
while the shared epithelium is L0 (OSCC 48k / HGIN 3.2k / CAA 7.2k) and L1.

Key state matched to **cluster 3** (7/5: BNC2, CASK, ERBB4, FHIT, NAALADL2, NFIA, PEG3), a small CAA
cluster inside the shared epithelial lineage L0, where **ERBB4 is 1st within its state** (global
4th), nothing above it, PEG3 10th there; nominated; drugs 0.75; no decoy ranked or nominated; 19
nominations (precision 0.05; the floor over-nominates with seven states). In cluster 16 ERBB4 is
2nd behind GABRA4; in cluster 17 it is unranked: the siblings of 16 and 17 are each other and
cluster 13 (all CAA), so a pan-CAA-epithelium gene cancels in the sibling contrast (the HCC
pattern, protocol 19/21), and the within-lineage contrast for L4 has 645 reference cells. The 0.3
de novo runs (prompts v2): ERBB4 4th global (inherited lineages), 1st (derived lineages).

**Reading error in the scripted critic, generic:** it flagged "Cluster 16 enrichment does not hold
within its parent cell type" and the same for 13, because within L4 (99% CAA) the per-sample
proportions of 16 vs the rest of L4 show no shift. When a derived lineage is itself
condition-exclusive, the lineage *is* the disease state and the within-lineage test is void, not
failed. Recorded as issue 15 below. S13 (sample_3) flagged again (L1, L6) by the marker-free check
with the paper's processing removed; 63 of its 5,284 cells sit in clusters 3/16/17.

## Issue 15 (fixed, generic; protocol 23p, prompts v13): "enrichment does not hold within its parent cell type" on a condition-exclusive lineage

The scripted critic's within-lineage check assumed the lineage holds both conditions. With derived
lineages the CAA keratinocytes (clusters 13, 16, 17; 59k cells, 99% CAA) were their own lineage L4,
so the within-L4 test compared CAA sub-states with each other and "failed" for 16 and 13. Fix: the
critic sums the disease fraction over the lineage's clusters first; at >= 90% it records one minor
flag per lineage saying the test is void, the lineage is the disease population, pan-state genes
cancel in the sibling contrast (read the disease-wide axis) and attribution rests on a thin reference
pool. The LLM critic prompt carries the same rule (v13). Tests `tests/test_critic_exclusive_lineage.py`
(3). The de novo agent run 1 is under v12; de novo agent run 2 (v13) verifies and gives the arm its
stability pair.

## CAA de novo agent run 1 (`runs/caa_denovo_agent_v04_1`, prompts v12) — 29 min, 2.16M input tokens, 92 calls, 1 error — PASS

Harmony on the library, derived lineages, nothing reused. Five composition states recorded: the
three L4 CAA clusters (16, 17, 13 - the critic later merged them: "one L4 tumour population split
by Leiden", merged log2 6.50, 6/6 samples), cluster 3 (CAA cluster inside the shared epithelium L0)
and cluster 12 (L2, cycling). Key state matched to cluster 3 (7/5; strict 2/5 - the agent recorded
the state by its own markers, which overlap the key's literature list in two genes). **ERBB4 1st
within its state** (cluster 16; also 1st in cluster 3; global 2nd behind IGF2), nominated twice
(clusters 16 and 3), afatinib / dacomitinib / neratinib found; PEG3 2nd, CACNB2 and NEO1 5th among
acceptable. 7 unique nominations (precision 0.14); **FGFR1 (decoy) nominated from cluster 17 and
flagged major** ("not disease-elevated, shared with reference", within-cluster log2FC -0.34); no
decoy unflagged. The critic reviewed all eight nominations, recommended the S13 exclusion arm
(sample_3 again, in L1 and L6, with the paper's processing removed), and read ERBB4 in cluster 3
as the best-supported candidate because that cluster has an adequate in-cluster reference pool
while the L4 clusters have a "too thin" one (reference cells 52% ERBB4-positive there). It wrote
"16 and 13 are not enriched within cell type (log2 0.07 and -0.82)" - the reading that v13 / 23p
now corrects to "the lineage is condition-exclusive, the test is void" - but drew the right
conclusion anyway through the merge. The 0.3 de novo runs (prompts v2): ERBB4 4th (inherited
lineages) and 1st (derived) globally.

With Harmony and no inherited partition the answer is the same as with the paper's clustering
(ERBB4 found in the CAA keratinocyte state, nominated, with its drugs), ranked 1st instead of
3rd-5th because the GABA-A subunits spread across the three L4 pieces and the sibling contrast
of cluster 3 is against OSCC/gingival keratinocytes, not against other CAA cells.

## CAA de novo agent run 2 (`runs/caa_denovo_agent_v04_2`, prompts v13) — 55 min, 2.78M input tokens, 104 calls, 0 errors — PASS; 23p verified; stability pair

Four states recorded: cluster 17 (L4, the largest CAA piece, high confidence), 3 (L0), 12 (L2,
cycling) and 15 (L6 fibroblasts, by the expression route, low confidence); clusters 16 and 13 were
left as "not characterised separately" (the analyst recorded only the top enriched cluster of each
lineage). Key state matched to cluster 3 (6/5); strict missed (best 4/5 in cluster 12: the key's
literature markers are spread across the CAA pieces). **ERBB4 1st within its state** (cluster 3;
global 3rd behind IGF1 and EPHA3, both from the fibroblast cluster 15 and both flagged major by the
critic as "stromal target, not a tumour-cell target; reference almost all OSCC"), nominated,
afatinib / dacomitinib / neratinib found; PEG3 5th, CACNB2 7th among acceptable. 12 nominations
(precision 0.08); no decoy ranked above ERBB4 or nominated.

**23p verified.** The critic's own flag: "Lineage composition: L4 (clusters 16/17/13) is
condition-exclusive, so the within-lineage test is void, not failed ... L4 is itself the CAA tumour
population, so the all-cells test (log2 6.56, p = 0.00045) is the relevant one and it holds. Sibling
contrasts inside L4 compare CAA sub-states with each other and cancel pan-state genes." It named a
second degenerate case on its own: cluster 12 is the only cluster of L2, so its within-lineage
proportions are 1.0 in every sample and only the all-cells test counts (the scripted floor did not
flag cluster 12 either). It recommended the S13 exclusion arm again (sample_3, flagged in L1 and L6)
for cluster 15; not run on this arm: the primary arm's exclusion arm already removed that library,
and cluster 15 is a stromal program the critic set aside anyway. GABRB1 was not nominated in either
de novo run, so the 23n part-2 fallback was not exercised in a live run; its test is
`test_compara_gap_falls_back_to_human_symbol`.

**Stability pair (de novo):** top-10 Jaccard 0.18, nominated 0.19 (ERBB4, PNOC and FGF14 shared);
rank-1 genes IGF2 (run 1) and IGF1 (run 2). Much lower than the primary arm (0.67 / 0.42): with
derived lineages there are five composition-enriched clusters plus program clusters, and the two
runs recorded different subsets (run 1 all three L4 pieces, 3 and 12; run 2 only 17 of L4, plus 3,
12 and the fibroblast cluster 15), so the nomination lists differ by state choice, not by
disagreement about the same state. The scored quantity is stable: ERBB4 1st within its state in both
runs and in the floor, nominated with its drugs each time, no decoy unflagged. In both runs the
global top is a stromal IGF-axis gene that the critic flags as a microenvironment candidate.

## Round 0.4c closed: tag results-v6 (working copy only)

Three floors, five agent runs (9.8M input tokens), four generic fixes (issues 12-15; protocol
23m-23p; prompts v13), tests 104 passed / 1 skipped (`validation_0.4/full_tests_round_0.4c_close.log`),
results table `docs/results_caa_v04c.md` (+ .json) regenerated from the run directories. Not pushed;
the original clone, the incubator and GitHub untouched.

## Issue 16 (fixed, report-side, generic; protocol 23q): figures could not draw a state whose cluster was merged after it was recorded

Reported by the author after the round closed ("whenever clusters are merged the subsequent plots
aren't able to plot the merged clusters together"). Confirmed in `caa_denovo_agent_v04_1`: the
state-share figure shows "state 16: 0 cells", "state 17: 0 cells", "state 13: 0 cells" and the
per-library expression panels have no shaded state group for them, because the critic's
`merge_clusters(16, 17, 13)` relabelled those cells `13+16+17` in the active column and the report
looked the recorded labels up there. Checked every run that merged (12): the same defect in
`flc_agent_1` (state 8 recorded, then merged into 8+19); the other ten recorded the state under the
merged label or merged clusters that were not recorded states. The UMAP is drawn at clustering
time, so no figure ever showed a merge. Analysis tools were not affected (they resolve a label
against the base clustering), and no number changes.

Fix in the report code only: `resolve_states` (recorded label → the merged label that contains it,
pieces collapsed into one drawn state; or the base-clustering column after a subcluster step), used
by the state-share panels and the expression-panel groups; a report-time `state_map` figure of the
recorded states on the embedding under the final key; UMAP coordinates saved per run so rebuilt
reports can draw the map (older runs: only if the dataset carries `X_umap`). The first rebuild
showed a second rebuild defect: the reloaded session re-derived the data-derived lineages and
renumbered them (panels said L14, L8, L9 where the run says L1, L5, L6); the resume now restores
the lineage column from the run's own fine-cluster-to-lineage record instead. A third: the reloaded
FLC file brought back the QC-removed cells, which had no cluster assignment and formed an
"unassigned" lineage group; the resume now drops cells without an assignment. Tests
`tests/test_report_merged_states.py` (6). The two affected reports were rebuilt in place with
`scripts/rebuild_report.py` and their evaluation sections refreshed; see the rebuild logs
`validation_0.4/caa/caa_denovo_agent_v04_1_rebuild.log` and `validation_0.4/flc_agent_1_rebuild.log`.
Full suite after the fix: 110 passed, 1 skipped (`validation_0.4/full_tests_round_0.4c_23q.log`).
Committed as 9dd9dbc; the unpushed `results-v6` tag was moved onto it so the tag the skill manifest
names (protocol 23-23q) carries the report fix.

## CAA agent results read against the paper's fine annotation (`refined_anno`, never shown to the agent)

Asked by the author whether the CAA agent results make sense. Cross-tabulation of the agents'
clusters against the paper's labels (a check, not an input: the column is in
`forbidden_obs_columns`):

- **Clusters 6 and 12 are one paper state.** 97% of cluster 6 (29,449 of 30,237) and 84% of
  cluster 12 are "CAA-specific Basal KCs"; all six dogs contribute to both (S17 is 78% in cluster
  6, the others 18-34% in each). The agents' "EMT-like" (6) versus "PEG3 / GABA-high" (12) split is
  substructure within the paper's basal state at resolution 1.0, not the paper's two states. The
  exclusion-arm critic's reading ("one CAA epithelial spectrum") was right.
- **The paper's junctional state was missed in the primary arm.** The 8,358 CAA-specific junctional
  KCs sit in cluster 13 (1,558; the cluster is 2,300 cells, 70% CAA, 6/6 dogs, log2 1.59, MWU
  p 0.07, so it fails the enrichment rule on the p-value), cluster 12 (3,328), cluster 0 (1,946)
  and cluster 11 (648). No primary-arm run recorded it; the expression route flagged clusters 0, 11
  and 4 as CAA-altered programmes, which is where the junctional and the CAA basal cells inside
  shared clusters show up. The de novo arm's cluster 3 (3,049 cells: 1,558 junctional, 999
  suprabasal) is the only recorded state that isolates part of it, and it is the state where
  ERBB4 ranks 1st with siblings that are OSCC and gingival keratinocytes.
- **Cluster 4 (programme state in all three primary runs)** is 85% OSCC proliferating basal KCs
  with 2,012 CAA proliferating basal KCs inside; GABRA1 / EPHA3 / CA10 / GPC5 there are the
  CAA-wide neuro programme measured in cycling CAA cells against cycling OSCC cells. The critics
  said "CAA-wide, not a cluster-4 state" each time. The score still puts GABRA1 at rank 1 globally
  (score 1.000) in all three runs because a programme candidate's effect size is the
  within-cluster contrast and nothing discounts a programme drawn from a 13% minority of a shared
  cluster; the evaluator uses within-state ranks so this does not touch the scores, but a reader
  of section 3 sees a sparse GABA-A subunit above ERBB4.
- **De novo clusters:** 13+16+17 = 54,060 CAA basal + 4,461 junctional KCs (the merge was right);
  cluster 12 = proliferating basal KCs of both conditions (the critic's "cycling sub-state" is
  right); cluster 15 = fibroblasts, 2,287 CAA vs 6,980 OSCC (the critic's "stromal, not a
  tumour-cell target" is right).
- **Candidates:** GABA-A subunits above ERBB4 are the paper's neuroepithelial programme and the
  critics dismissed them as targets for the right reason; ERBB4 survives in every run at medium
  confidence with the right caveats (normal epithelium expresses it; approved pan-ERBB drugs act
  through EGFR/ERBB2); CA10 correctly dismissed (catalytically inactive, no inhibitors). FGFR2 is
  genuinely state-specific in the data (72% of cluster 6 cells) while the paper's kinase screen
  showed no response to FGFR inhibition: the exclusion-arm critic calling it "best-supported
  actionable" is where the expression data and the functional result disagree, and the key's
  decoy label reflects the functional result. KCNH8 (v12 run) is a novel marker with no literature.

Consequence for the write-up: REPORT 7c, README and the skill's validation text said "both
CAA-specific states recorded as distinct"; corrected to "the basal population as two clusters;
the junctional state missed in the primary arm". The key is frozen and scores one shared
signature, so this changes no score; recorded as a key-design limitation (one marker list per
reported state next time). Not a method change.

## Context audit and the model behind the runs (protocol 23r)

Asked by the author whether the skill runs out of context. The skill package itself is small:
SKILL.md ~5.8K tokens (287 lines, description 721 characters), all references together ~14K, so
~20K with everything open. The agent roles are where context grows. From the `model_call` usage
records of the 16 agent runs: per-role peak input tokens 14K planner, 80-186K analyst, 110-217K
nominator (ovarian no-leak run 1: three consecutive turns over 200K), 40-88K critic. Fixed material
per role is 7-9K tokens (system prompt 0.4-2K, protocol 0.2K, tool schemas 5-7K) plus Claude Code's
own system prompt (~5K); growth is tool results re-sent on every turn by the stateless CLI backend
(condition_attribution ~3.2K tokens per call, gene_specificity and score_targets ~3K,
check_reference_libraries ~4K). Every role conversation in every run ended with `end_turn`; no turn
budget and no context error anywhere.

Found while doing this: every agent run of rounds 0.4a-0.4c was launched as bare `--model
claude-code`, so the CLI chose its own default model and the run log recorded only the backend name
(`"model": "claude-code"`). Which model answered cannot be established from the logs. On this
machine the CLI's non-interactive default resolves today to `claude-sonnet-5-5` (checked with a
one-word `claude -p` call); the 217K-token turns succeeded, so the window was above 200K. The API
backend default in every dataset config (`agent.model: claude-sonnet-4-5`, a 200K model) would have
failed in the nominator on the larger datasets; left as is for now because `AnthropicModel` also
sends `temperature`, which the 5.x models reject, so moving that default is a separate change.

Fix (generic; protocol 23r; commit e212986): `ClaudeCodeModel` reads `modelUsage` from the CLI
result envelope and writes the model id into every `model_call` usage record (`usage.model`),
into `facts.model_ids`, and into the report header; `--model claude-code:<id>` pins the model.
The evaluator then broke on the new non-numeric usage field (`int('claude-fable-5-1')`): token
totals for `run_eval` and the report now go through one helper that sums numeric fields only
(`usage_totals_from_events`; test `tests/test_eval_usage_totals.py`). Full suite 110 passed,
1 skipped (`validation_0.4/full_tests_model_id.log`).

## CAA agent run 3, pinned to Claude Fable 5.1 (`runs/caa_agent_v04_3`, prompts v13, scoring v3) — 27 min, 2.13M input tokens, 94 calls, 1 tool error — PASS

Asked by the author to run the pipeline on the latest model. `--model claude-code:claude-fable-5-1`
(the CLI accepted `claude-fable-5-1` and `claude-opus-5-5` under this login), all libraries,
prompts v13; every model call records `claude-fable-5-1`. 27.2 min, 2.13M input tokens (1.6M for
the two primary runs on the CLI default), 94 tool calls (71 / 79), peak context 166K in the
nominator.

Evaluation (`validation_0.4/caa_agent_v04_3.log`, `runs/caa_agent_v04_3/eval.json`): PASS. State
found (cluster 12, 13/5 lenient); strict missed by one gene in each CAA cluster (cluster 6 recorded
BCL2, CACNB2, RBFOX1, COL8A1; cluster 12 GABRB1, PEG3, NAALADL2, ERBB4; 4/5 each: the analyst
recorded the sibling-contrast tops rather than the lineage-wide list the earlier runs recorded).
ERBB4 4th within its state (global 5th, the best global rank of the three primary runs); above it
GABRA4 and GABRB2 (novel, as in every run) and GABRB1 (acceptable). Drugs 0.75 (afatinib,
dacomitinib, neratinib). FGFR2 nominated and flagged ("disease-wide, 86% of the reference cells in
the cluster express it"); no unflagged decoy. Nominated precision 0.2: 10 nominations, five of them
microenvironment candidates from the three programme states (SELE, EPHA3, IGF1, KIT, TNFRSF11B),
each flagged "not the disease cells" or "no inhibitor". Stability with the v12 run: top-10 Jaccard
0.54, nominated-set 0.31 (v11 vs v12 was 0.67 / 0.42).

What differed with this model: (i) the critic ran the reference-only subcluster on the 1,237 OSCC
cells inside cluster 6 unprompted (as the v12 critic did) and read the largest piece as an
EMT/invasive carcinoma programme (MDGA2, SULF1, FN1, ZEB1, PTHLH; CACNA2D3 in 97% of its cells),
making CACNA2D3 conditional; the one tool error is its attempt to do the same for the 42 HGIN
cells (the tool refused, the critic moved on). (ii) It named the exclusion arm (S13) as the next
step, which the round had already run. (iii) More programme-state nominations, all flagged as
microenvironment, which lowers precision without touching the within-state rank. (iv) Its ERBB4
note ("disease-wide across both CAA epithelial states, shared with cluster 6") says what the
fine-annotation check found afterwards. Same answer, a third more tokens, a longer critic section.

Not a method change; no prompt or rule touched. Tag results-v7 (working copy only; nothing pushed).

## Frozen environment, container image, CI (protocol 23s; 2026-11-13)

The user asked whether to move from conda to Docker. Answer given: keep conda as the primary path,
add a lock file, ship Docker as the secondary path for reviewers and CI. Reasons: Docker is
installed on this machine but the daemon is not reachable for the account (rootless podman fails
too, no subuid range), Singularity is; the real gap was that `environment.yml` floats every version
but two and the versions behind results-v7 were recorded nowhere; the agent run needs the Claude
Code CLI login on the host.

Done: `scripts/freeze_env.sh` (explicit conda lock with md5 + pip layer), `environment/*.lock`
(182 conda packages, 57 pip), `Dockerfile` installing from the lock files only, `.dockerignore`
(no data, no runs), `.github/workflows/ci.yml` (build, `pytest -q`, offline demo, publish to GHCR
on master and tags). README sections 1, 2, 9, 11 and the skill's install notes updated. The lock
files were re-installed into a scratch prefix and the fast tests run against it (result in this
section's follow-up line).

Follow-up: the lock files were installed into a fresh prefix (`conda create --file` then the pip
layer with `--no-deps`, then `pip install --no-deps -e .`); `pip check` clean; the same versions
(scanpy 1.11.5, anndata 0.12.19, pydeseq2 0.5.4, numpy 2.4.6, pandas 2.3.3, pynndescent 0.6.0);
`pytest -m "not slow"`: 110 passed, 1 skipped, 1 deselected in 6.5 min. Committed and pushed; the
CI build of the image is the first full-suite run inside the container.

CI run 1 (9f4484a) failed: the image built, `pytest` died at collection with `libgomp.so.1` missing.
Cause: `~/.condarc` lists `defaults` before conda-forge with `channel_priority: disabled`, so 8
locked packages come from `defaults` (`_openmp_mutex`, `libxcrypt`, `libdrm`, `libexpat`, `tk`,
`xcb-util-image`, `xorg-libxdamage`, `python-dateutil`); the defaults `_openmp_mutex` brings no conda
libgomp, so scikit-learn and numba use the system one, which
the scratch prefix on this machine found and the Debian base image lacks. Fix: `apt-get install
libgomp1` in the Dockerfile; the lock is left as the record of what ran. Written into protocol 23s.
