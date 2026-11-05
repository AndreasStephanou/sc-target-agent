# Expected behavior — case 06

The skill should:

1. **Classify the design** as row 4 (non-malignant disease vs healthy; the reference holds every
   lineage), `objective: target`, and say at once that this design type has not been validated:
   the four validation datasets are tumours; the proliferation flag, the composition route and the
   attribution guards have only been exercised there.
2. **Name the independent unit.** 14 donors, not 34 libraries. Run `scripts/design_check.py` with
   `--sample-key library --condition-key diagnosis --donor-key donor` and report the
   pseudo-replication warning. Recommend pooling to donor (`sample_key: donor`) as the primary
   analysis, or state the inflation if libraries are kept; not run with 34 "replicates" and report
   p-values as if they were independent.
3. **Read the confounding table.** `chemistry` is partly aligned with the condition (all controls
   v3, some IPF v2): not perfectly confounded, but strongly associated. Say that any v2-vs-v3 gene
   would look IPF-specific in the v2 donors, that Harmony on the sample key mitigates but does not
   remove this, and propose a check: are the top candidates up in the v3 IPF donors alone? `region`
   is nested within donor and balanced across conditions: fine.
4. **Verify the labels** with a marker the group trusts for IPF tissue (not a candidate target):
   per-library `scripts/label_check.py`, and the calibration run's `check_reference_libraries`.
   A control lung with fibrotic regions is the analogue of the HCC "normal" library.
5. **No key, so the no-key protocol**: positive controls (the aberrant basaloid markers KRT17,
   KRT5-negative, COL1A1 in the basaloid state, SPP1/TREM2 in the macrophage state, should appear in
   the tables; lineage markers should be `shared_with_reference`), negative controls (cycle and
   stress genes flagged), three runs for stability, replicate panels, direction. Say that these
   are sanity checks, not a validation.
6. **Expect composition findings and say so.** Macrophage expansion and basaloid appearance are
   composition shifts plus programmes; a "state" whose enrichment vanishes within its cell type is
   a composition result and the `disease_wide` genes of that lineage are the programme. The
   proliferation flag on an expanding population means "describe the expansion", not "artefact".
7. **Direction for every candidate**: in fibrosis the therapeutic direction is often inhibition of
   a secreted or receptor gene, but agonists exist (e.g. for some receptors); state it per candidate.
8. **Report** with the design statement first, the limitation that the method is untested on
   this design type, candidates as hypotheses with the orthogonal experiment each needs, and the
   human-review points open.

The skill should not:

- treat 34 libraries as 34 replicates, or hide the chemistry imbalance;
- present a macrophage-wide or fibroblast-wide gene as "the target" without saying the composition
  route favours consistent infiltrates and that this is a composition finding;
- claim the method is validated for fibrosis because it worked on tumours;
- write an answer key from the results.
