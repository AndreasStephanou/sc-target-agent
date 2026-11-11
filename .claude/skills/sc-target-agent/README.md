# sc-target-agent

Status: Experimental (version 0.4.0, draft)

## Overview

Single-cell atlases of disease usually end with a list of "potential targets" chosen by the
analysts. This skill runs an AI agent that, given only a count matrix and per-library condition
labels, finds condition-associated cell states, ranks candidates with evidence, and has a critic
flag the weak points, all through deterministic tools with a complete audit log. If a target is
already known for the condition, a frozen answer key lets a lab test the method blind before
trusting it on an open question.

Version 0.4.0 rewrites the procedure so that it is accurate about what has been validated and
usable beyond tumour-vs-normal designs. The researcher classifies the study design first
(`references/design-types.md`), chooses the objective (`target` for druggable targets, `marker`
for condition-specific response genes), checks replication, confounding and label integrity with
two small scripts before spending tokens, and reads every flag with a glossary that separates the
general reading from the tumour-specific one (`references/flags.md`). The pipeline code itself is
unchanged (repository tag `results-v3`).

The code is in the `sc-target-agent` repository (`references/repository.md`). This skill is the
procedure for using it well: what to check before a run, how to read a calibration pass, how to
read a set-aside candidate, and which decisions stay with the researcher.

## Intended users

Researchers with a two-condition single-cell or single-nucleus RNA-seq dataset: disease vs
reference tissue, tumour vs another lesion, treated vs control, knockout vs wild-type, infected vs
mock. Comfort with running Python scripts from a terminal is enough. Validation so far covers
four tumour datasets; non-malignant disease and perturbation designs are supported by the code but
untested, and the skill says so at every step where it matters.

## Inputs

- An `.h5ad` with integer counts (or per-sample 10x matrices, or a GEO matrix + metadata deposit).
- The sample column (one value per library) and the independent biological unit (donor, animal).
- The condition column, which value is the case, and what the reference represents.
- Species.
- The design type and objective (`references/design-types.md`).
- Optional: a cell-type annotation; a known target and marker set for an answer key.

## Outputs

A run directory with `report.md` / `report.html` (data summary, states, ranking with
`target_class`, attribution and druggability, critic flags), `state.json`, tables and figures, and
`run.jsonl`, the audit log of every tool call and the exact payload the model saw. With an answer
key, an evaluation line (state recovery lenient/strict, target rank, decoys above the target,
decoys nominated unflagged, drug recall). Without one, the no-key protocol's sanity checks and a
stability estimate. The skill's own output format adds a design statement and a per-candidate
therapeutic-direction line.

## Example use

A lab has snRNA-seq from six tumours and three normal tissues and an in-vitro hit on one receptor
tyrosine kinase. They classify the design (tumour vs another lesion plus healthy tissue), write and
freeze the key, run the design and label checks, run the free calibration pass, fix the QC
thresholds for their tissue, then run the agent and see whether it ranks the kinase at the top,
what else it ranks above it, and what the critic flags. A second example (`examples/`) shows a
knockout-vs-wild-type request in `marker` mode and what an honest response looks like for a design
the method has not been validated on.

## Human review

The researcher must classify the design and name the independent replicate; verify the condition
labels and pre-register how doubtful libraries are handled; decide whether a confounded or
under-replicated design should be run at all; judge whether the answer key is fair; decide whether
QC thresholds fit the tissue; decide whether a recorded state is biology, artefact or a
composition shift; judge each candidate's therapeutic direction and, for a candidate set aside by
attribution, whether the reference population identified changes the claim; and treat every
candidate as a hypothesis for orthogonal follow-up. Changes to code, prompts or scoring are method
changes and must be versioned and justified, never made to fit one dataset. Details in `SKILL.md`.

## Limitations

Expression-based (protein-level targets are invisible); replicate-aware tests cannot repair nested
libraries or a confounded covariate; druggability has no direction; attribution can set a real
target aside when a rare reference population expresses it; the composition route favours
infiltrates consistent across patients; nomination precision 0.1-0.3; no seed control on the
subscription backend; validated on four tumour datasets only, with `marker` mode run on synthetic
data alone; human-centric drug knowledge with ortholog mapping for other species.

## Dependencies

Runs code. The `sc-target-agent` repository with its conda environment (`environment.yml`:
scanpy, harmonypy, pydeseq2, scrublet, pyyaml, anthropic). The agent run needs a Claude Code
login or an `ANTHROPIC_API_KEY`; calibration and the demo run offline. Gene-level lookups go to
DGIdb, Open Targets and NCBI E-utilities. The two helper scripts in `scripts/` need only anndata,
h5py, numpy, pandas and scipy and send nothing anywhere.

## Data considerations

The pipeline's privacy guard blocks per-cell matrices and real sample identifiers from reaching
the model; cluster-level summary tables, marker lists and gene symbols are sent to Anthropic, and
gene symbols to the druggability services. Whether a given dataset may be used this way is the
researcher's and their institution's decision; do not run it on restricted human data without
checking. Evaluation material here is public (GSE159115, GSE149614, GSE286057), synthetic, or
fabricated for reasoning tests. Never commit the dataset, the run directories or an API key.

## Testing

`evals/` holds eight cases: a straightforward public dataset with a frozen key (ccRCC); the same
dataset as actually encountered (swapped metadata, tissue-inappropriate QC, a protein-level target
in the key); an ambiguous mouse request with two replicates per condition and no known answer; the
HCC reviewer run (a passing install that cannot cluster, a contaminated reference library, a
cycling sub-state); the ovarian reviewer run (dataset choice by composition, a known target set
aside by attribution, knowing when to stop); and two new fabricated cases for designs the method
has not been validated on: a non-malignant lung-fibrosis study with libraries nested in donors,
and a knockout-vs-wild-type study in `marker` mode; and one out-of-scope request (bulk RNA-seq,
one sample per arm) that the skill must decline. `evals/test-cases.json` lists them for the
CoPAIR checker. The helper scripts were run on the repository's synthetic dataset and on the
public ccRCC object (`scripts/README.md`). The revised workflow has not yet been exercised
end-to-end on a new dataset, by the author or by anyone else; that is the next test.
`references/validation-results.md` summarises the four validation datasets and lists what has
and has not been tested.

## Maintainers

Andreas Stephanou (@AndreasStephanou), De Vlaminck lab and Sethupathy lab, Cornell University.

## Changelog

See `CHANGELOG.md`.

## Acknowledgement

Skill structure follows the CoPAIR incubator conventions; other skill examples are at
https://github.com/anthropics/skills.
