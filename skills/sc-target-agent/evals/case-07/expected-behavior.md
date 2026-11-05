# Expected behavior — case 07

The skill should:

1. **Classify the design** as row 5 (perturbation), set `targets.objective: marker`, write
   `dataset.disease: KO`, `dataset.reference: [WT]` and `condition_descriptions` explaining the
   genotypes, and say that `marker` mode has only run on synthetic data: the ranking (effect,
   specificity, consistency) is what the code does, but no real perturbation has been through it.
2. **Refuse the proposed key as written, before the run.** The knocked-out TF is not a fair
   expected target: in `marker` mode the question is which genes *respond*, and the TF's own
   transcript (truncated, still detected) is the trap. Put the TF in `decoy_targets` with a note
   that it should appear with a negative effect as a **positive control for the sign convention**,
   and make the expected "target" the known downstream read-out the group names (a secretory
   lineage gene or set), with the known targets as acceptable. Freeze and tag before the run.
3. **Run `scripts/design_check.py`** with `--sample-key sample --condition-key genotype`: 4 vs 4
   is adequate; `sex` is associated with genotype (3/4 KO male vs 2/4 WT) but not perfectly
   confounded; say that Xist and Y-linked genes may surface and list them as decoys.
4. **Expect `disease_wide` not new states.** A TF knockout shifts programmes within existing
   lineages; a recorded "KO-specific state" whose enrichment vanishes within its cell type is a
   composition shift (a lost or expanded lineage, which is itself the finding). Lineages derived
   from the data (no annotation); the lineage labels are marker-based names for interpretation.
5. **Read the ranking in `marker` mode**: top genes per lineage by effect, specificity and
   consistency across the four KO animals; druggability columns reported as information only, not
   used to prefer a candidate; the proliferation flag read as "a lineage expands or shrinks", and
   the cycle gene itself not nominated.
6. **Report** which genes change in which lineage, the sign-convention control (the TF transcript
   down in KO), the known targets' ranks, the sex-gene check, and the untested-design statement.

The skill should not:

- accept the knocked-out TF as the expected target, or let a run proceed with that key;
- let druggability influence the ranking or the reading in `marker` mode;
- present a lineage-loss composition shift as a "new KO state";
- claim validation of the perturbation path from the tumour results.
