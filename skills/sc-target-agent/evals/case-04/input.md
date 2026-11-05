# Case 04 — problem: HCC (GSE149614), as encountered in the reviewer run

A researcher has cloned the repository fresh, built the conda environment today, and downloaded
GSE149614 (hepatocellular carcinoma, 10 primary tumour and 8 non-tumour liver libraries, plus PVTT
and lymph-node libraries). GEO provides one dense count matrix and a metadata table (`Cell`,
`sample`, `site`, `patient`, `celltype`). The authors' `celltype` has a single `Hepatocyte` class.
They know GPC3 is the textbook HCC target (codrituzumab, CAR-T) and want a blind test, then their
own data.

They show you, in order:

1. `pytest -q` passes and the synthetic demo passes. The first calibration run on HCC ends after
   500 s with "No disease-specific state found" and 7 tool errors in the log.
2. A per-library check of GPC3 in hepatocyte-class cells: 25–98% in nine tumour libraries, 2% in
   one tumour library (HCC04T), 0–5% in seven non-tumour libraries, and 77% in HCC04N, labelled
   non-tumour. AKR1B10 and SPINK1 behave the same way in HCC04N.
3. A later calibration run that completes: three tumour-hepatocyte clusters enriched in all ten
   patients, one of them with UBE2C, TOP2A, CKS2, PTTG1 as its top markers; GPC3 is not in the
   ranking at all; TOP2A is nominated; in the tumour-vs-normal hepatocyte contrast GPC3 shows
   log2FC 2.6 at padj 0.20.
