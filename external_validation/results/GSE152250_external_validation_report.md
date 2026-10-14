# GSE152250 external validation

Primary validation time point was prespecified as 24 h to match GSE97358. GSE152250 culture replicates are analyzed as independent/unpaired; no donor pairing is assumed. Cross-study effect sizes are not pooled.

Decision: **GO**

| Time | CTRL mean FPKM | TGFb1 mean FPKM | log2 ratio of means | one-sided exact p | two-sided exact p | all TGF > all CTRL |
|---:|---:|---:|---:|---:|---:|:---:|
| 6 h | 98.4793 | 101.491 | 0.0435 | 0.4000 | 0.8000 | False |
| 24 h | 108.605 | 171.000 | 0.6549 | 0.0500 | 0.1000 | True |
| 48 h | 76.4247 | 91.5926 | 0.2612 | 0.0500 | 0.1000 | True |

Boundary: this tests independent within-study directional replication of the CTHRC1 response. It does not convert three culture replicates into donor-level evidence and does not pool FPKM-scale magnitudes with the donor-paired GSE97358 estimand.
