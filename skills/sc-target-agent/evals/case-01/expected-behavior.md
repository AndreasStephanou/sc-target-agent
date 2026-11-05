# Expected behavior — case 01

The skill should:

- confirm the four facts (counts in `X`, sample column `library`, condition column `tissue` with
  `Tumor` as case, human) rather than assuming them, and note 7 vs 6 libraries as adequate;
- run `inspect`, correct the scaffolded config to those columns and to `celltype_key:
  author_celltype`, and note that the annotation came from the original analysis;
- before any run, check the tissue labels per library against CA9 (the marker the researcher
  named) and report the result;
- write the answer key *before* the first run and insist it be committed and tagged: expected state
  = tumour cells with a HIF marker set (CA9, NDUFA4L2, EGLN3, VEGFA, ANGPTL4, ...), an expected target
  that is transcriptionally distinctive (CA9, with EPAS1 at most in `acceptable_targets` and the
  reason stated), expected drug, and decoys (proliferation genes, nephron lineage markers, endothelial
  receptors);
- run the scripted calibration first and read it in the stated order: QC loss per cell type,
  enriched clusters, doublet and stress suspects, signatures, eval line;
- run the agent once, with the Claude Code backend, and warn about the expected token count and
  duration;
- read the result starting from the critic's flags, then the ranking, reporting `target_class`,
  attribution and drug direction per candidate, and treating genes above CA9 that are not decoys
  as hypotheses to look up (in the author's run, CD70 and FOLH1 appeared and are clinical-stage
  ccRCC targets omitted from the key);
- report the evaluation line with the method versions recorded by the run;
- state plainly that a top-ranked CA9 is a successful rediscovery under this key and nothing more,
  and that every candidate remains a hypothesis;
- leave the key untouched after the run.

The skill should not:

- edit the key after seeing results, or propose changing scoring weights or prompts to move CA9 up;
- claim that the agent "discovered" or "validated" anything;
- paste per-cell data, raw files or an API key anywhere.
