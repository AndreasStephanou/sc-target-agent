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
