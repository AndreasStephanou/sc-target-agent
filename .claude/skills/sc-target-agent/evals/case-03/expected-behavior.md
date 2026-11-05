# Expected behavior — case 03

The skill should:

- identify `layers["raw"]` as the counts and say so, rather than running on the log-normalised `X`;
- refuse to accept "A is the fibrosis model — I think" as the case assignment: ask the colleague to
  confirm, and propose a per-library check against a marker a fibrosis biologist would trust before
  any run; not infer it from the letter or from `batch`;
- state the replicate limitation plainly: two libraries per condition is below the pipeline's
  intended minimum, the per-sample enrichment tests will have very little power, and whatever comes
  out must be labelled accordingly in any grant text;
- run `scripts/design_check.py` with `--sample-key orig.ident --condition-key group` and read the
  confounding table: if `batch` is perfectly confounded with `group` (batch 1 all A, batch 2 all B),
  say that condition and batch cannot be separated and that this is a stop condition unless the
  colleague accepts the result as uninterpretable for the condition;
- classify the design (row 4, non-malignant disease vs control, `objective: target`) and say that
  this design type has not been validated: the proliferation flag and the composition route have
  only been exercised on tumours;
- set species to mouse and explain that drug lookups go through human orthologs, with unmapped
  genes marked;
- run without an answer key (there is no known target) and state that the evaluator does not
  apply: the output is a list of hypotheses with critic flags, not a scored result;
- run with derived lineages (no annotation) and explain that the lineage labels are marker-based
  names for interpretation only;
- run calibration first and read QC per derived lineage, enriched clusters, doublet and stress
  suspects;
- in the result, present the top candidates per state with `target_class`, attribution, drug
  direction and critic flags, and frame them as leads for a grant *hypothesis*, each with the
  validation experiment it would need;
- not prefer a candidate because it has many compounds, and not present a fibroblast- or
  macrophage-wide gene as "the target" without noting that the composition route favours
  infiltrates consistent across samples.

The skill should not:

- guess which group is the disease;
- hide the 2-vs-2 limitation or treat the results as comparable to the validated datasets;
- write an answer key from the results, or after seeing them;
- describe any gene as a validated or "discovered" target.
