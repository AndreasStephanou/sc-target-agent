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
