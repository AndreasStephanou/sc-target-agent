# Flags, classes and statuses: what each one means and how to read it

Names are the ones written by the repository's tools (`score_targets`, `condition_attribution`,
`find_markers`, `differential_abundance`, `check_reference_libraries`) and by the critic. "Read as"
gives the generalised reading; "tumour reading" the one validated so far.

## Candidate classes (`target_class`, from `score_targets`)

| value | meaning | read as |
|---|---|---|
| `state_specific` | defines this state among its lineage siblings (sibling contrast) and is not shared with the reference | a marker of a case-specific cell state |
| `disease_wide` | up in case cells across the lineage versus the reference counterpart (lineage contrast, log2FC >= 1) | a programme shift in an existing lineage; the natural class for non-malignant disease and perturbations |
| `lineage_axis_source` | which lineage contrast carried the gene: the shared-cluster (compartment) one, or the whole-lineage one used when nearly all case cells sit in enriched clusters (tumours) | tells you whether the comparison was case cells vs their normal counterpart in the same clusters, or the whole lineage |

## Score flags (`flags` column of `*ranked_targets*.csv`)

| flag | effect on score | meaning |
|---|---|---|
| `proliferation_marker` | halves | a cell-cycle gene. Tumour reading: "the cells divide". Non-malignant reading: the population expands; describe its identity from the lineage contrasts, do not nominate the cycle gene |
| `low_specificity` | penalty | below `targets.min_specificity_log2fc` against other states |
| `inconsistent_across_samples` | penalty | not up in most case replicates; one or two libraries carry it |
| `disease_wide` | informational | see class above |
| `shared_with_reference` | penalty | expressed by reference cells of the lineage or by both conditions: a lineage or pan-condition gene |
| `reference_driven` | **blocking** | the in-cluster reference cells express it at least as highly as the case cells: the signal is theirs. Before accepting: read `reference_pool_in_cluster`, then identify the cells (lineage, libraries, count, markers) |
| `ambient_like` / `ambient_suspect` | blocking / major | the pattern matches ambient RNA from a dominant population (a shift in every lineage) |
| `not_druggable` / `no_known_compounds` | lowers druggability / compound terms | from DGIdb and Open Targets; `source=fixture` means the bundled offline fixture, unverified |
| `missing_evidence:<field>` | term treated as 0 | the lookup failed or returned nothing (HTTP 429 once emptied the Open Targets fields; re-run the lookup) |
| `rare_reference_counterpart` | none (informational) | a small reference population (< `min_reference_frac_of_condition`, default 0.5% of reference cells) also expresses it: an on-target consideration to identify, not a disqualification |

## Attribution statuses (`reference_pool_in_cluster`, from `condition_attribution`)

| status | meaning | consequence |
|---|---|---|
| (usable) | >= 2% of the cluster's cells are reference and the reference holds the lineage | the in-cluster comparison decided the class |
| `too_thin` | fewer than `min_reference_frac_in_cluster` (2%) reference cells in the cluster | in-cluster comparison not used; class decided by the other evidence |
| `rare_in_reference` | reference cells of this lineage are < 0.5% of all reference cells | the reference barely holds the lineage; in-cluster comparison not used |

Known failure (ovarian, row 2 of the design table): 151 benign epithelial cells (7.5% of the
cluster, above the 2% guard) expressed FOLR1 more highly than the tumour cells and the gene was
`reference_driven`; on the raw object they were 59 FOLR1+ cells out of 94k reference cells, a
benign Müllerian-type inclusion. The lineage-relative guard freed it on rescore but a coarser
derived lineage that absorbed mesothelium defeated the guard in the next run. **Rule for the
reader:** a blocking attribution verdict on a candidate that the per-sample abundance test and the
compartment contrast both support is a prompt to look at the reference cells, and the candidate is
reported as conditional on them, not dropped.

## State and cluster flags (from `find_markers`, `differential_abundance`, the critic)

| flag | meaning |
|---|---|
| `doublet_suspect` | co-expresses two lineages' markers or Scrublet-high; not a disease state however enriched |
| `stress_artifact_suspect` | heat-shock / immediate-early dominated (HSPA1A/B, FOS, JUN, DNAJB1); dissociation artefact. Stress genes are stripped and the cluster rejected only if < 10 strong markers remain |
| `cycling_substate` | a sub-state defined by cell-cycle genes; real biology, but its identity genes come from the lineage contrasts |
| `split_state_candidates` | clusters that are the same state split by replicate; `merge_clusters` is the remedy |
| composition shift (`differential_abundance(within_celltype=True)` loses the enrichment) | the "state" is a change in cell-type proportions, not a new state within a type; report as composition |
| library flagged by `check_reference_libraries` | within a lineage, that library correlates more with the opposite condition (no gene names used): swapped label, contamination, or low purity. A prompt to look; exclusion is a human config decision with a written reason |

## Critic severities

`blocking` (the candidate or state should not be reported as a lead without resolving it),
`major` (report with the caveat), `minor` (note). The critic also writes a summary stating which
candidates survive and which wet-lab experiment would settle each doubt; quote it.

## Therapeutic direction (not a flag; a human check)

Druggability counts compounds. The critic has caught EGLN3 (PHD inhibitors stabilise HIF: wrong
direction for a HIF-driven tumour) and agonist antibodies on TNFRSF4/9/18. For every top candidate
write: *compounds listed act by [inhibition | activation | ADC payload delivery]; the therapeutic
hypothesis for this design is [inhibit | activate | deliver]; consistent / inconsistent.*
