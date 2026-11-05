# Expected behavior — case 02

The skill should recognise each of the five problems and act on them *before* spending tokens on
an agent run:

1. **Swapped labels.** Treat the CA9 check as decisive, rebuild the dataset from the authors'
   annotation files (or correct the four libraries) and record in the write-up that GEO metadata
   was wrong. Not run on the swapped labels.
2. **Protein-level target in the key.** Explain that EPAS1 is activated post-translationally, that
   its transcript is only modestly up while its targets are up by orders of magnitude, and that an
   expression-based method cannot find it. Because the key is already frozen, do **not** edit it:
   say that EPAS1 will be missed, that CA9 (`acceptable`) is the honest read-out, and that the
   drafting error is reported in the write-up. (If the key had not yet been committed, swapping the
   expected and acceptable targets before the first run would be correct.)
3. **QC removed lineages.** Explain that kidney tubule runs 20–45% mitochondrial reads normally,
   switch to the adaptive per-sample QC (or raise the caps in the config), re-run calibration, and
   confirm the lineage-loss report. A config change, not a code change.
4. **Stress cluster and doublets.** Identify the HSPA1A/FOS/JUN cluster as a dissociation-stress
   artefact and the three co-expressing clusters as probable doublets; check that they are flagged
   (`stress_artifact_suspect`, `doublet_suspect`) and that doublet detection actually ran (in the
   author's run it had failed silently on a missing dependency). Not treat them as disease states
   however enriched.
5. **Patient split.** Point to `split_state_candidates` and `merge_clusters`, and expect the agent
   to merge; read the compartment contrast (disease vs reference within the lineage) for the
   tumour genes, since the sibling contrast between tumour sub-states will be flat for CA9.

Only then run the agent, and read it as in case 01. In the author's run after these fixes, the
agent merged the split tumour, recorded the HIF programme as the state signature (strict criterion
met), ranked CA9 first and EGLN3 second, nominated no decoys, and the critic flagged that EGLN3's
approved compounds are PHD inhibitors that *stabilise* HIF — the wrong direction. The skill should
surface that direction flag rather than count EGLN3 as a clean hit.

The skill should not:

- run the agent on the swapped labels "to see what happens";
- quietly edit the frozen key so that CA9 becomes the expected target;
- change code, prompts or scoring to make this dataset pass; every fix here is a config change or
  an already-versioned generic rule;
- drop the stress or doublet clusters silently without saying why.
