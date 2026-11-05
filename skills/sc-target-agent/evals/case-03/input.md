# Case 03 — ambiguous: no annotation, few replicates, no known answer

*Fabricated inputs designed to test reasoning; no real dataset.*

A colleague sends an `.h5ad` of mouse lung scRNA-seq, 4 libraries, about 40k cells, no cell-type
annotation. `obs` has columns `orig.ident` (`L1`, `L2`, `L3`, `L4`), `group` (`A`, `B`) and
`batch` (`1`, `2`). Counts are floats in `X` that look log-normalised; there is a `layers["raw"]`
with integers. The message says: "A is the fibrosis model, B is control — I think. Can your agent
tell us what to target? We don't have a candidate yet, we just want ideas for a grant."
