# Design types: what changes between a tumour, a non-malignant disease and a perturbation

The pipeline's code accepts any two conditions (`dataset.disease` vs `dataset.reference`; the
comments in `config/default.yaml` say treated vs vehicle, knockout vs wild-type and infected vs mock
map onto them) and two objectives (`targets.objective: target | marker`). What differs between
designs is not the tool chain but **what the reference contains**, **which flags mean what**, and
**what a fair key looks like**. Classify the study before anything runs, and write the row into the
config (`dataset.condition_descriptions` and a comment) so the model and the report carry it.

| # | design | reference contains | objective | expect | what the test is | validated? |
|---|---|---|---|---|---|---|
| 1 | tumour vs its tissue of origin (ccRCC vs benign kidney, HCC vs non-tumour liver) | the tumour's normal counterpart lineage | `target` | `state_specific` and `disease_wide` candidates in the malignant lineage; attribution informative | known target ranked first in its state | yes (ccRCC, HCC) |
| 2 | tumour vs a tissue that lacks the lineage (ovarian metastasis vs omentum) | no normal counterpart; pure composition | `target` | every lineage gene is "specific"; attribution verdicts fragile; the real question is ranking among lineage genes | known target *found and ranked among* the real lineage targets (list them as acceptable) | yes (ovarian), with the attribution limitation |
| 3 | tumour vs another tumour or lesion (CAA vs OSCC) | a different disease sharing the tumour background | `target` | pan-tumour genes (proliferation, EMT, stress) are the decoys; case-exclusive states | known target ranked first; pan-tumour decoys flagged | yes (CAA) |
| 4 | non-malignant disease vs healthy (fibrosis, IBD, neurodegeneration, infection in tissue) | the same lineages, shifted in abundance and programme | `target` or `marker` | both routes: disease-only states do occur (IPF aberrant basaloid cells, SPP1 macrophages) and take the composition route; programme shifts in the remaining lineages take the expression route; proliferation flag may hit a genuinely expanding population | known disease gene ranked first in its lineage, or state recovery as the primary metric when no transcript-level approved target exists (IPF); composition-only findings labelled as such | **in progress** (IPF, GSE135893, round 0.4) |
| 5 | perturbation (treated vs vehicle, KO vs WT, infected vs mock, in vitro or in vivo) | the same cells without the perturbation | `marker` | condition-specific response genes per lineage; druggability irrelevant to the ranking; the pure expression-route case (protocol 23/23b fixed this path) | the perturbed pathway's known read-out ranked first; the perturbed gene itself is a decoy | **in progress** (doxorubicin vs vehicle, GSE285475, round 0.4) |
| 6 | paired / matched only (same donors in both conditions, no between-donor replication) | - | - | within-cluster contrasts are not blocked by donor; dispersion inflated | - | **not supported**: state the limitation; the tools will run but the tests are not what the design needs |
| 7 | more than two levels, time series, continuous dose | - | - | - | split into pairwise runs, each pre-registered | not supported in one run |

"Validated" means a frozen key, a blind run and a written result in `references/validation-results.md`.
For rows 4 and 5 the code paths exist (`objective: marker` reweights the score to effect size 0.40,
specificity 0.35, consistency 0.25 and zero on druggability; the critic and attribution run the same)
but no real dataset has been through them. Say so to the researcher.

## What the reference represents, and why it decides the reading

The known target's fate in all three public tumour datasets was decided by the reference set: a
swapped label (ccRCC), a "normal" library expressing the tumour programme (HCC), and a rare benign
population expressing the target inside the tumour cluster (ovarian). Generalised:

- **Reference holds the lineage** (rows 1, 4, 5): in-cluster attribution is meaningful; a
  `shared_with_reference` verdict usually means a lineage gene. `disease_wide` is the natural class
  for a programme shift.
- **Reference lacks the lineage** (row 2): in-cluster attribution has nothing to compare against;
  `too_thin` / `rare_in_reference` statuses are expected; any reference cells that do sit in the
  cluster are by definition unusual (benign inclusions, spill-over) and can dominate the verdict.
  Read `reference_driven` here as "identify these cells" rather than "discard the gene".
- **Reference is another disease** (row 3): the shared disease background is the decoy source;
  a gene up in both conditions versus a healthy third group is not a case-specific target.

## Label checks in non-malignant and perturbation designs

In a tumour-vs-normal design a lineage marker of the tumour reads over the whole library (CA9 in
ccRCC: 19-58% of cells in tumour libraries, under 1% in benign). In a non-malignant disease the
diseased tissue still contains every normal lineage, so whole-library marker fractions mostly track
composition and the condition medians sit close together (IPF, GSE135893: MMP7 13 points apart, SPP1
12, KRT17 7). Run `label_check.py` **within the compartment the marker belongs to** (`--celltype-key`
with the authors' compartment or the pipeline's derived lineage): within epithelium MMP7 separated
IPF from control 36% vs 10%; within immune cells SPP1 did 29% vs 2%, and listed the largest IPF
donor as control-like for its macrophages (regional heterogeneity, pre-registered as an arm, kept in
the primary run). In a perturbation the same applies, and a response gene's fraction-positive can
even move against the literature direction at the library level (doxorubicin, 2 weeks: Cdkn1a+
fraction lower in treated hearts) because composition and depth differ; record such an observation
in the key header as a risk, do not edit the key.

## Choosing the sample column

The per-sample tests treat each value of `sample_key` as an independent replicate. If donors
contribute several libraries (regions, sorts, lanes), the libraries are not independent:
`scripts/design_check.py --donor-key` reports it. Options, in order of preference: pool libraries
to the donor (set `sample_key` to the donor column; the pipeline pseudobulks per sample anyway);
keep libraries but state the inflation and interpret p-values as descriptive; or, for a paired
design, state that donor effects are not modelled.

## Keys and decoys per design type

A fair key names: expected state markers that are **consistent across replicates in this dataset
family** (6-20 genes, `min_marker_overlap` about 5), an **expected target that is transcriptionally
distinctive**, acceptable targets (real but not the test), expected drugs, and decoys. Decoys are
druggable-looking genes that a naive ranking would put up and that the method should flag:

| design | decoys to list |
|---|---|
| tumour (rows 1-3) | proliferation (MKI67, TOP2A, UBE2C), lineage markers of the tumour's cell of origin, pan-epithelial (KRT8/18, EPCAM when not the target), endothelial receptors of real anti-angiogenic drugs (KDR, FLT1), stromal (PDGFRB, ACTA2), pan-tumour stress/EMT genes when the reference is another tumour |
| non-malignant disease (row 4) | generic injury and interferon-response genes (IFIT1-3, ISG15), heat-shock and dissociation genes (HSPA1A/B, FOS, JUN), the markers of the lineage that expands (e.g. collagen genes in fibrosis are the state, not a target unless the question says so), complement and acute-phase genes |
| perturbation (row 5) | the perturbed gene itself (KO), its direct sensor or receptor when the ligand is the treatment, the sex and strain genes if groups differ, immediate-early genes |

Write in the key **what the test is** (ranked first; found and ranked among; recovered in the
expected lineage) and which pre-registered exclusion arm exists, before the first run.

Three key-drafting failure modes have now been seen, each recorded in the write-up and never fixed
by editing the key: a target regulated at the protein level (EPAS1, ccRCC); literature state
markers that are not consistent across replicates in the dataset (HCC); and a response gene that is
real but **not at the sampled time point** (Cdkn1a two weeks into weekly doxorubicin: not induced,
while the cardiomyocyte atrophy and natriuretic genes Trim63, Fbxo32, Nppa, Nppb were). For a
perturbation, ask when the read-out peaks relative to when the cells were collected, and put
time-robust read-outs (the lineage-level consequence) in the expected list and the acute ones in
`acceptable`.

## The no-key protocol

Most real users have no known answer. Without a key the evaluator does not apply, and nothing in
the pipeline tells the researcher whether the run is trustworthy. Do these instead, and report them
in the Evaluation section:

1. **Positive controls that must behave.** Lineage markers of the recorded state's lineage should
   appear as `shared_with_reference` or `low_specificity`, not as top candidates. If the condition
   has any published associated gene, look for it in `tables/*pseudobulk_de*` and `*ranked_targets*`
   and report where it sits (not as a score, as a sanity check).
2. **Negative controls that must not pass.** Proliferation genes, dissociation-stress genes and
   doublet clusters should carry their flags. One unflagged obvious artefact in the top 5 is a
   reason to re-read the calibration, not to trust the list.
3. **Stability.** Three agent runs; the top-5 overlap per state (Jaccard) and the set of nominated
   genes. `scripts/run_eval.py --stability` computes this but currently requires `--ground-truth`;
   without a key, compare `tables/*ranked_targets*.csv` across runs by hand and report the overlap.
4. **Replicate panels.** For each top candidate, the per-library panel in the report: most case
   libraries above the reference ones in the candidate's own lineage, not one library or every lineage.
5. **Direction and plausibility** for each top candidate, by the researcher.

Say in the write-up that these are sanity checks, not a validation, and that the candidates are
hypotheses with named follow-up experiments.

## What to expect in non-tumour designs (untested; stated as expectation)

- The strongest signal is often a **composition shift** (more macrophages, fewer alveolar type 2
  cells) rather than a new state. The critic's `differential_abundance(within_celltype=True)` check
  asks exactly this; a state whose enrichment vanishes within its cell type is a composition
  finding and should be reported as one, with the lineage-level `disease_wide` genes as the
  programme, not a "state".
- **`proliferation_marker`** will flag cycle genes in a genuinely expanding population. The right
  reading is "this population expands; describe its identity genes from the lineage contrasts",
  not "artefact". The score still halves for the cycle gene itself, which is correct: a cell-cycle
  gene is not a disease-specific target.
- **Interferon and stress programmes** are shared across many diseases and dissociation protocols;
  list them as decoys and expect the critic's `stress_artifact_suspect` flag to catch the
  dissociation part only.
- **`marker` mode** ranks by effect, specificity and consistency; the druggability columns are
  still printed. For a KO, the knocked-out gene should appear with a large negative effect: a
  positive control for the sign convention, and a decoy for nomination.
