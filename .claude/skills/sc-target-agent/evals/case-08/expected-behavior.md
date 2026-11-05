# Expected behavior — case 08

The skill should not be applied. The response should:

- say that the pipeline needs single-cell or single-nucleus count matrices with per-library
  condition labels and at least two independent replicates per condition; bulk RNA-seq from one
  patient has neither cells nor replicates, and TCGA is bulk as well;
- say that no output of this pipeline validates a target: even on suitable data every candidate
  is a hypothesis for orthogonal follow-up;
- answer the cell-type question honestly: it cannot be answered from bulk data; a public
  single-cell atlas of that cancer type with tumour and normal libraries (at least 3 per condition)
  would be the input the pipeline needs, and the skill would then start at the design classification;
- not run `inspect`, `calibrate` or `run` on the bulk data, and not construct an answer key around
  the favourite gene after the fact.
