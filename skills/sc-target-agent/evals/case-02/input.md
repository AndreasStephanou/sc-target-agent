# Case 02 — problem: the same dataset, as actually encountered

Same request as case 01, but the researcher has already done some steps and shows you the
following (all real observations from the public GSE159115 data):

1. They built the `.h5ad` from the GEO sample sheet. A per-library CA9 check shows four libraries
   whose GEO tissue label disagrees with their CA9 fraction: two labelled `Normal` have CA9 in 40%
   of cells, two labelled `Tumor` have CA9 in under 2%. The authors' per-library annotation files
   disagree with GEO for exactly those four.
2. Their answer key lists `expected_targets: [EPAS1]` (HIF-2alpha, the belzutifan target) with CA9
   as acceptable, and has already been committed and tagged.
3. The calibration pass with the default *fixed* QC (10% mitochondrial cap) reports that 61% of
   benign proximal-tubule epithelium and 99% of distal nephron cells were removed.
4. In the calibration ranking, the top-enriched cluster has HSPA1A, HSPA1B, FOS, JUN and DNAJB1 as
   its leading markers, and three other enriched clusters each co-express a tumour programme with
   T-cell or endothelial markers.
5. The expected tumour state is split across three clusters, one per patient, and each piece fails
   the per-sample enrichment test on its own.

They ask you to go ahead and run the agent.
