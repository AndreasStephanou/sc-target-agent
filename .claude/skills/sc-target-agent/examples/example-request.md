# Example request

I have single-nucleus RNA-seq from a canine oral tumour study: 6 tumour libraries, 10 libraries of
a different oral carcinoma, and 3 normal gingiva libraries, about 200k nuclei after the original
QC, in one `.h5ad` with a cell-type annotation from the published analysis. We have an in-vitro
result showing the tumour cells are sensitive to a receptor tyrosine kinase inhibitor and we think
we know which receptor is responsible.

I want to know whether an agent, given only the counts and the condition labels, would find the
tumour-specific epithelial state and nominate that receptor, without being told what we expect,
and what else it would rank above it. Can you set this up as a blind test, run it on the cluster
with my Claude Code login, and help me read the result?
