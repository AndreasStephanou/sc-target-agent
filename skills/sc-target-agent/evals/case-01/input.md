# Case 01 — straightforward: public ccRCC with a known answer

I downloaded GSE159115 (Zhang et al. 2021, PNAS): scRNA-seq of clear cell renal cell carcinoma,
13 libraries from 7 patients, tumour and matched benign kidney, with the authors' cell-type
annotation files. I built a single `.h5ad` (about 29k cells) with the tissue label per library in
`obs["tissue"]` (`Tumor` / `Normal`), the library in `obs["library"]`, the authors' cell type in
`obs["author_celltype"]`, and integer counts in `X`. Species is human.

ccRCC is driven by VHL loss and HIF activation; CA9 is the classic marker and belzutifan targets
HIF-2alpha. I want to use this as a blind test of the pipeline: does it find the tumour-cell state
and rank CA9 or the HIF programme at the top, without being told? Then I want to run it on my own
data.

I have a Claude Code login on the cluster and the `sctarget` conda environment is built.
