# Case 05 — difficult: a 2025 deposit chosen by search (GSE286057), as encountered

A researcher wants to test the pipeline on a dataset published in the last two years with an
approved drug whose target is transcriptionally distinctive. A GEO search returns, among others, a
2025 gastric cancer atlas (47 tumour and matched-normal libraries, per-sample 10x files) and a 2025
ovarian cancer study (omental metastases from 7 patients with distant omentum, plus omentum from 8
non-cancer donors). They show you:

1. For the gastric atlas, after downloading all 47 libraries: PTPRC-positive cells 70% and
   EPCAM-positive cells 1-2% in both conditions.
2. For the ovarian study, a two-library probe: P9_tumor has 2% EPCAM+ cells; P13/P14/P15_tumor have
   10-47%; omentum libraries have ~0% EPCAM/PAX8/FOLR1 but 10-70% MSLN; two "distant omentum"
   libraries from cancer patients have 3% and 10% EPCAM+ cells.
3. After the primary agent run: the tumour state is found; FOLR1 ranks 16th and is set aside as
   shared_with_reference; the critic flags two reference libraries as tumour-like and asks for an
   exclusion arm.
4. After the exclusion arm: FOLR1 is blocking as reference_driven because 151 reference cells inside
   the 2,016-cell tumour cluster express it more highly than the tumour cells. On the raw object
   those are 59 FOLR1+ cells in 94,235 reference cells, 26 from one healthy donor, KRT19/CLDN3/
   WFDC2-positive and CALB2/UPK3B-negative.
