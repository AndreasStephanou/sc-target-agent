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

## FLC scripted floor (`runs/flc_calibrate_1`, 46 min, 121 calls, 0 errors)

Primary vs normal liver, 9 vs 9 libraries, metastatic set aside (254k nuclei through QC and
Scrublet; the floor's slowest run). Tumour state recorded (cluster 6, derived lineage L1), key
markers 3/5 (CA12, KRT7, PRKACA; five needed, so "missed" by the key's own threshold);
**PRKACA 2nd within its state behind CA12** (both key answers; global 3rd), recall@5 1.0,
nothing above them, KRT7 7th, no decoy ranked or nominated; drugs: none of the trialled agents
(no PRKACA inhibitor exists; ENMD-2076/neratinib/everolimus not reached by the floor). The v6
floor in the test clone (2026-11-02) had PRKACA 6th with four novel genes above it and markers 1/5.

## Issue 11 (fixed in the skill's helper, generic): pairing was reported as pseudo-replication

The donor warning counted libraries per donor, so a paired design (one tumour and one normal
library per patient) was reported as "8 donors contribute more than one library ... over-counts
replication". Pseudo-replication is several libraries of one donor in the same condition; the
count is now per donor-condition pair. HCC: no warning (correctly: pairing only); ovarian: 5
control donors with 2-3 regions each (correctly flagged).
