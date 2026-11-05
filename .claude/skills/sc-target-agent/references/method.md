# Method summary

Condensed from `docs/REPORT.md` and `docs/EVALUATION_PROTOCOL.md` in the repository, which remain
the authoritative sources.

## Design

Four roles, each a fresh model conversation, act only through deterministic tools and pass state to
each other only as structured findings:

- **planner** — reads the dataset summary and sets the plan;
- **analyst** — QC (adaptive per sample, preview then apply), doublet detection (Scrublet per
  sample, with a fallback), clustering and Harmony integration, lineage derivation when there is no
  annotation (coarse Leiden, marker-based labels for interpretation only), per-sample differential
  abundance with split-state candidates and cluster merging, sibling pseudobulk DESeq2 and a
  disease-vs-reference compartment contrast within the lineage's shared clusters, dissociation-
  stress and doublet flags; records disease states with a signature (a state must carry genes);
- **nominator** — sees both the sibling and compartment tables plus marker lists, runs specificity,
  condition attribution (with ambient and thin-reference-pool guards), druggability (DGIdb, Open
  Targets, PubMed), and the deterministic score; nominates candidates;
- **critic** — re-examines every nominated state and candidate and writes flags (doublet, stress,
  patient split, label inconsistency, thin pool, wrong therapeutic direction, decoy-like).

Every tool call and the exact payload the model saw are logged to `run.jsonl`. No per-cell data or
real sample identifiers reach the model.

## Scoring

Scoring v2 combines the sibling axis (does the gene define this state among its lineage siblings)
and the compartment axis (is it up in disease cells across the lineage, against reference cells in
shared clusters) with specificity, attribution and druggability. `target_class` is
`state_specific` or `disease_wide` (the latter requires lineage-level log2FC ≥ 1). The objective
can be `target` (default) or `marker`.

## Evaluation

Against a key frozen and committed before the first run: state recovery (lenient — the cluster's
tables carry the signature; strict — the agent itself recorded it), target rank, acceptable-target
rank (reported, not scored), decoys above the target, novel genes above the target (reported
separately, not penalised), decoys nominated without a critic flag, drug recall, nominated
precision, stability across repeat runs, tokens. Comparison arms: a naive baseline script and the
same tools under a fixed scripted procedure (also the free calibration pass).

## Rules added after the freeze, and the failure that motivated each

| version / rule | failure it fixed |
|---|---|
| prompts v2: a recorded state must carry genes; the nominator sees both tables | v1 analyst recorded the state with 0 genes; ERBB4 never reached the nominator |
| scoring v2: compartment axis, `target_class` | ERBB4 was modestly up vs siblings but 8-fold up vs reference epithelium; sibling-only scoring ranked it 2nd behind a decoy |
| thin-reference-pool guard in attribution | 150 reference cells (0.5% of a 31k cluster) overruled the rest and ERBB4 was lost in the de novo arm |
| Harmony called directly, loud failure | harmonypy 2.x shape change made the scanpy wrapper fall back to PCA silently |
| derived lineages | annotation-free runs had no lineage axis |
| adaptive per-sample QC, lineage-loss report, preview mode | fixed 10% mito removed 61% of benign kidney epithelium and 99% of distal nephron |
| `split_state_candidates`, `merge_clusters` | malignant cells split by patient; each piece failed the per-sample test |
| stress flag (strip stress genes; reject only if <10 strong markers remain) | heat-shock dissociation cluster topped the ranking; a first version also rejected a real 98%-tumour cluster |
| doublet detection with fallback | Scrublet failed silently on a missing dependency while three clusters were plainly doublets |
| prompts v4: compartment contrast in signature, candidate pool and nominator brief | sibling contrast was flat for every tumour gene (CA9 log2FC −0.06) because the siblings were other tumour sub-states |
| analyst finishing guard; empty-reply retry in the Claude Code backend | the analyst ended a run with an empty reply and no recorded state |

Each has a numbered section in `docs/EVALUATION_PROTOCOL.md`.

## Objectives and design types (added in skill 0.4.0)

`targets.objective: target` (default) scores effect size 0.25, specificity 0.25, druggability 0.25,
known compounds 0.15, consistency across samples 0.10. `targets.objective: marker` reweights to
effect 0.40, specificity 0.35, consistency 0.25 and zero on the two drug terms; druggability is
still looked up and printed. Everything else (states, contrasts, attribution, critic) is identical.
The code takes any two conditions; the validation covers tumour designs only. `design-types.md`
maps each design to what the reference contains, which flags are informative and what a fair key
is; `flags.md` is the glossary.

Two steps of the method are the ones most likely to behave differently outside tumours and are
flagged as such in the skill: the proliferation flag (motivated by cycling tumour sub-states; a
genuinely expanding population in a non-malignant design will carry it) and condition attribution
(its verdict depends on what reference cells sit in the cluster, which a design without a normal
counterpart lineage cannot supply).
