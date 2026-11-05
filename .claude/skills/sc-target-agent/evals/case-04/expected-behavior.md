# Expected behavior — case 04

The skill should:

1. **Not read the first run as a negative.** Seven tool errors and no clustering means the chain
   broke; open `run.jsonl`, find the first error, and recognise an environment or data-format
   problem. Run `python scripts/sctarget.py check-env` (step 0) and rebuild from the pinned
   environment if it fails. Say that a passing test suite and demo did not exercise this path.
2. **Build the object with the generic builder** (`build_from_matrix.py --dry-run` first), keep
   only the primary-tumour and non-tumour libraries, and confirm raw counts are integers.
3. **Freeze a fair key before any run**: GPC3 as the expected target (transcriptionally
   distinctive), with state markers chosen to be consistent across patients rather than the
   textbook list, and decoys spanning proliferation, hepatocyte lineage, endothelium and
   cholangiocyte.
4. **Treat HCC04N as a doubtful library, before the run**: say what the marker check shows, note
   that the marker is the expected target (so the check is not blind), run the marker-free
   `check_reference_libraries` for independent evidence, keep the library in the primary run, and
   declare in advance an arm with `dataset.exclude_samples: [HCC04N, HCC04T]` and a written reason
   if the result hinges on it. Not silently drop or relabel it; not edit the key.
5. **Explain the completed run correctly.** The cluster with UBE2C/TOP2A/CKS2 is a cycling tumour
   sub-state (`cycling_substate`); its cell-cycle markers describe proliferation, not the disease,
   and TOP2A must be flagged (`proliferation_marker`), not nominated. GPC3 at log2 2.6 and padj 0.20
   in tumour-vs-normal hepatocytes is a replicate-aware test refusing a gene that a "normal"
   library expresses in 77% of its cells: a data-label result, not a method failure, and exactly
   what step 4 predicted.
6. **Run the pre-registered arm**, report both arms side by side (the all-libraries run is the
   primary result), and read the per-library panel for GPC3 with the researcher: tumour libraries
   high in the hepatocyte lineage, normal libraries at the floor, one low tumour library, ambient
   leakage into immune lineages.
7. **Report the candidates as hypotheses** with the critic's flags: the cholesterol-synthesis
   programme (SQLE, FDPS, IDI1), NDUFA4L2, CYP17A1, the tumour-endothelium set (FOLH1, PGF, ESM1,
   THY1), RAMP1 with its fibroblast expression noted, and the therapeutic-direction caveat for each.

The skill should not:

- accept "no disease-specific state found" from a run with tool errors;
- clean the reference set with the expected target and present the result as blind;
- exclude the library without a recorded reason, or rebuild the dataset without it and say nothing;
- lower a threshold, change a margin or add a gene to make GPC3 appear;
- call GPC3 "validated" or the HCC result a success without stating that it depends on the exclusion.
