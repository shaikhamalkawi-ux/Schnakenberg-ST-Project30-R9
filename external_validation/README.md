# GSE152250 independent external validation

## Purpose

This directory contains the reproducible external-validation layer for the manuscript's primary 24-h TGF-beta1 -> CTHRC1 cardiac-fibroblast result.

The validation was prespecified before inspecting the target result:

- primary validation time point: 24 h, because it matches the GSE97358 primary endpoint;
- secondary temporal context: 6 h and 48 h;
- primary criterion: mean CTHRC1 expression under TGF-beta1 exceeds control at 24 h;
- GSE152250 culture replicates are treated as independent/unpaired;
- no donor pairing is assumed;
- no cross-study numerical pooling with GSE97358 is performed.

## Official source inputs

The executable workflow downloads these files directly from NCBI GEO:

- `GSE152250_Fibroblasts_FPKM.txt.gz`
- `GSE152250_family.soft.gz`

Raw third-party GEO files are not redistributed in this repository. Their SHA-256 hashes are recorded in `results/GSE152250_input_sha256.txt`.

## Accepted result

At 24 h, mean CTHRC1 was 108.605 FPKM in control and 171.000 FPKM after TGF-beta1, a log2 ratio of means of +0.655. All three treated replicates exceeded all three controls. The exact directional 3-versus-3 permutation value is 0.05; the two-sided value is 0.10 because only 20 label allocations exist.

At 48 h the direction remained positive (76.425 vs 91.593 FPKM; log2 ratio +0.261; all treated above all controls). At 6 h the contrast was near-null (98.479 vs 101.491 FPKM; log2 ratio +0.043).

The accepted execution was GitHub Actions run 38585174863. The analysis artifact digest was `sha256:44b9dbe24e376a20a73100df9ef2ecc7d12f744492b05f2191fcc38c808a882e`.

## Claim boundary

This result is independent human primary cardiac-fibroblast directional/time-course replication. It does not convert three culture replicates into donor-level evidence, does not identify an onset time or myocardial hard delay, and does not justify pooling FPKM-scale magnitudes with the GSE97358 donor-paired estimand.
