# Example output (perturbation design)

*Illustrative structure for the request in `example-request-perturbation.md`. No run of this kind
has been performed; the numbers below are placeholders in brackets, not results. A real response
carries the exact values from the run directory.*

## Design statement

Design type: row 5, perturbation (treated vs vehicle, in vitro), `targets.objective: marker`.
The reference is the same cells without the treatment. Independent replicates: 3 wells per
condition, one library each (`design_check.py`: no nested libraries, no other sample-level
covariates in the object; ask whether the wells came from one animal's marrow, in which case the
biological unit is the animal and the replicates are technical). Species mouse: drug lookups go
through human orthologs and are reported as information only in `marker` mode.

**This design type has not been validated.** The four validation datasets are tumours; `marker`
mode has run on synthetic data only. The ranking the code produces (effect size 0.40, specificity
0.35, consistency across wells 0.25) is what will be read; the proliferation flag and the
attribution guards were motivated by tumour cases.

Labels verified with the cytokine's immediate receptor-induced gene the group trusts (not one of
the expected read-outs): all three treated wells high, all vehicle wells at the floor.

Answer key frozen before the run: expected state = the responding macrophage sub-population with
the canonical read-out set [6-10 genes the group named]; expected target = the canonical read-out
gene; acceptable = the other read-outs; decoys = the cytokine's receptor itself (expressed in both
conditions), immediate-early genes (FOS, JUN, EGR1), cell-cycle genes, and Xist / Y-linked genes if
the wells differ by sex. The test is "ranked first in the responding sub-population".

## Calibration (scripted, no tokens)

[QC per derived lineage; enriched clusters; whether a treated-only cluster appears or the shift is
within existing sub-populations; `check_reference_libraries` on 3 vs 3; eval line.]

## Agent run

**Condition-associated states.** [One or two states; for each, whether
`differential_abundance(within_celltype=True)` keeps the enrichment (a real sub-state) or loses it
(a composition shift between sub-populations, reported as such).]

**Ranking (marker mode, top 5 per state).** [gene, effect, specificity, consistency 3/3,
`target_class`, attribution, critic flag.] Druggability columns are listed and not used.

**Sign-convention control.** [The receptor-induced marker is up in treated; any gene the group
expects down is down.]

**Critic flags.** [Verbatim.]

## Evaluation against the frozen key

[State found lenient/strict; expected read-out rank; decoys above it; decoys nominated unflagged;
versions.] Drug recall is not applicable in `marker` mode.

## What this does and does not show

A top-ranked canonical read-out would be a successful rediscovery under this key on a design the
method has not been validated on, and nothing more. Candidates beyond the known read-outs are
hypotheses for perturbation follow-up (knockdown of the candidate, blocking antibody to the
receptor).

## Human-review points

- The biological unit (wells vs animals) and therefore the replicate count.
- Whether the key's read-out set is fair (consistent across wells, transcript-level).
- Whether the responding "state" is a sub-population or a composition shift.
- Sex of the wells, if the marrow came from different animals.
