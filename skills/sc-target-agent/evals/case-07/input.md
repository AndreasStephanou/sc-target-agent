# Case 07 — ambiguous: knockout vs wild-type, `marker` mode, the knocked-out gene is the trap

*Fabricated inputs designed to test reasoning on a perturbation design; no real dataset.*

A mouse group has scRNA-seq of small-intestine epithelium from 4 knockout and 4 wild-type animals
(one library each), about 60k cells, integer counts in `X`, `obs["sample"]` (8 values),
`obs["genotype"]` (`KO` / `WT`), `obs["sex"]` (3 KO male / 1 female; 2 WT male / 2 female), no
cell-type annotation. The knocked-out gene is a transcription factor they name; the KO is a
whole-body deletion and the transcript is still detectable (a truncated message).

They ask: "Which genes change in which cell type when we knock it out? We don't care about drugs.
Can the agent run as a blind test: it should find the TF's known targets in the secretory lineage
without being told. And please run it with the TF as the expected target."
