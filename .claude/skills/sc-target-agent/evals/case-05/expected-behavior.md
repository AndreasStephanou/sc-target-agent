# Expected behavior — case 05

The skill should:

1. **Reject the gastric atlas on composition**, before any analysis, and say why: the deposit is
   immune-enriched (or lost its epithelium), so no epithelial target can be tested; suggest it for
   immune/stromal programme questions instead. Not run 47 libraries through the pipeline to find out.
2. **Probe two libraries of the ovarian study before downloading all**, confirm tumour content varies
   by patient (one tumour library has almost none), and build with peritumoral libraries set aside as
   a declared arm.
3. **Write the key for the pure composition case**: with no epithelium in normal omentum every
   epithelial gene will be state-specific, so the test is whether FOLR1 is found and ranked among the
   other real ovarian-cancer surface targets (listed as acceptable), not whether it beats them;
   predict in writing that MSLN will be shared with mesothelium and that the two leaking distant-
   omentum libraries are the pre-registered exclusion.
4. **Read the primary agent run correctly**: the agent reached FOLR1 where the scripted floor had not
   looked; the critic's library flags match the pre-run check; the exclusion arm is the declared
   follow-up, run as a one-line config change with a recorded reason.
5. **Trace the blocking flag to the cells** before judging it: identify the 59 cells as a rare benign
   Müllerian-type epithelium (not mesothelium, not tumour spill-over), recognise that folate receptor
   alpha expression there is real biology and an on-target consideration, and say that the pipeline
   was right about the cells and overstated the verdict.
6. **Know when to stop**: if a generic guard is tried and is defeated by a property of the derived
   lineages, record it in the protocol and freeze the method with the limitation stated, rather than
   adding a third rule to reach the key. Report FOLR1 as found and ranked 16th in the primary run,
   set aside with a stated reason, blocked in the exclusion arms by a rare benign population.
7. **Carry the collateral findings**: the CAF state (INHBA, LOXL2, FAP, MMP11) and the Treg state as
   tumour-microenvironment biology with direction caveats (agonist antibodies); no decoy nominated.

The skill should not:

- download a whole series without a composition probe, or run a target test on an immune-only deposit;
- call the blocking flag a bug before looking at the cells, or call it correct without looking;
- keep changing the attribution rule until the known target ranks;
- present FOLR1, the CAF genes or the Treg genes as validated targets.
