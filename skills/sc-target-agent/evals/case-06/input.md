# Case 06 — difficult: non-malignant disease, libraries nested in donors, no known target

*Fabricated inputs designed to test reasoning on a design type the method has not been validated
on; no real dataset. Numbers are plausible, not measured.*

A pulmonary group has human lung scRNA-seq from explants: 8 fibrosis (IPF) donors and 6 control
donors (declined transplant lungs). Each donor was sampled in two or three regions, so the object
has 34 libraries: `obs["library"]` (34 values), `obs["donor"]` (14 values), `obs["diagnosis"]`
(`IPF` / `Control`), `obs["region"]` (`upper`, `lower`, `middle`), `obs["chemistry"]` (`v2`, `v3`),
`obs["cell_type"]` from the published analysis (38 types), integer counts in `layers["counts"]`,
log-normalised `X`. Species is human.

They say: "We know the aberrant basaloid cells and the profibrotic macrophages are there; we want
the agent to tell us what to drug. There is no approved IPF drug with a single transcript-level
target, so there is no answer key. The two IPF cohorts were run on different 10x chemistries,
controls were all v3."
