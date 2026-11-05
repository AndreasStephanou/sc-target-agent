# Example output

*Illustrative structure and level of caution for the request in `example-request.md`. Numbers are
rounded from the author's development runs; a real response carries the exact values from the run
directory.*

## Before running

Four facts confirmed with you: counts are integers in `layers["counts"]`; the sample column is
`library` (19 values, 6 case / 13 reference once the second carcinoma is treated as reference);
the condition column is `tumour_type` with case = the ameloblastoma label; species is dog, so
drug lookups use human orthologs and the report marks unmapped genes.

Answer key written to `data/ground_truth/<name>.yaml` from your in-vitro result and the published
markers, committed and tagged `<name>-key-frozen` before any run. I did not open it again.

Calibration (scripted, no tokens): QC kept 90% of nuclei with no lineage losing more than 15%;
two epithelial clusters enriched in the case condition across all six libraries; no doublet or
stress suspects (nuclear data); the eval line found the expected state and ranked the expected
receptor fifth, with a druggable decoy above it. That means the fixed script already reaches the
state; the agent run measures what reasoning adds.

## Agent run (Claude Code backend, prompts v4, scoring v2, derived lineages, 0.9M input tokens)

**Data summary.** 184k nuclei after adaptive QC, 19 libraries, lineages derived from the data
(annotation ignored), seven lineages labelled by markers for interpretation only.

**Disease-associated states.** One state recorded: a case-enriched epithelial cluster present in
6/6 case libraries and <1% of reference libraries. Signature (sibling + compartment contrast)
includes neural-crest and neuroepithelial genes and the expected receptor. No patient-split
caveat; the critic noted that the state is small (about 2% of case epithelium) and that a second,
larger case-wide signal exists across the shared epithelial clusters.

**Candidate ranking (state 1, top 5).**

| rank | gene | target_class | attribution | drugs (direction) | critic |
|---|---|---|---|---|---|
| 1 | *expected receptor* | disease_wide | disease_specific | 2 approved inhibitors (inhibit) | none |
| 2 | FGFR2 | state_specific | disease_specific | several inhibitors | flagged: listed as expressed in the reference carcinoma |
| 3 | gene C | state_specific | shared_with_reference | 1 compound | flagged: lineage gene |
| 4 | gene D | disease_wide | disease_specific | none | — |
| 5 | gene E | state_specific | ambient_like | 1 compound | flagged: ambient-like |

**Critic flags (verbatim, abridged).** FGFR2 is expressed in the reference carcinoma and may be a
pan-oral-tumour gene rather than specific to this state; gene C is a lineage marker; gene E's
signal looks ambient; the state is small and the disease-wide axis carries most of the receptor's
evidence.

**Evaluation against the frozen key.** State found (lenient and strict); expected target rank 1;
decoys above target 0; decoys nominated unflagged 0; drug recall 0.75 (three of four expected
compounds); one "novel" candidate above a decoy (FGFR2, which is a decoy in the key and was
flagged). Versions: prompt_version 4, scoring_version 2.

## What this does and does not show

The agent recovered the receptor you expected and ranked it first, and the critic flagged every
decoy it nominated. That is a successful blind rediscovery under this key. It does not validate the
receptor further than your in-vitro result already does, and the FGFR2 result is an open question:
your key treats it as a decoy, but in this data it is state-specific and FGFR2 is recurrently
mutated in the human counterpart of this tumour. I have left the key unchanged and reported this as
a hypothesis.

## Recommended next checks

IHC for the receptor and FGFR2 on the tumour sections; repeat the agent run twice for a stability
estimate; compare with the annotation-inherited arm.

## Human-review points

- Treating the second carcinoma as reference rather than a third condition (your decision; it makes
  the test harder and is stated in the report).
- Whether the small state size (2% of case epithelium) is biologically meaningful or a subcluster
  the original analysis would have merged.
- Whether FGFR2 belongs in the decoy list; if the key was wrong, that goes in the write-up, not
  into the key.
- Therapeutic direction of each compound listed (all inhibitors here, consistent with your result).
