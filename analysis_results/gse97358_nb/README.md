# GSE97358 donor-blocked count-model sensitivity — ED4

Completed 18 October 2026. This directory is a post-review statistical sensitivity on the same 168 samples / 84 donor pairs as the inherited primary analysis. It is not an independent cohort, new clinical evidence, or a fitted cardiac Schnakenberg mechanism. The main branch and prior manuscript files were not overwritten.

## Result
CTHRC1 log2 fold change = 0.883434548215682, 95% Wald interval [0.805045680498737, 0.961823415932626], two-sided Wald P = 4.06859974783245e-108. The inherited paired transformed mean remains 0.882 [0.801, 0.962], with 82/84 positive pairs. These are different estimands and different inference procedures, not interchangeable coefficients or evidence of improved biological significance.

All seven previously specified genes are reported in seven_gene_sensitivity.csv. q_BH_secondary6 applies only to POSTN, FN1, LOX, COL1A1, COL3A1 and ASPN; its blank primary entry is intentional. padj is a separate diagnostic BH adjustment over the tested gene universe, not the primary reporting family. ASPN remains unsupported at 0.05 (secondary q=0.432027361652092), and no gene or donor was removed. COL1A1 is a dispersion outlier, for which DESeq2 retained the gene-wise dispersion; it is not a dropped gene or a Cook's-distance exclusion.

## Execution
- Statistical specification: sensitivity/gse97358_nb/run_analysis.R, immutable original commit eb846acfca6155fb68029a39ca63a23a1675815b.
- Complete run commit: 39e7c60c9632f531eab25484061c80bca452b0b0.
- Successful GitHub Actions run: https://github.com/shaikhamalkawi-ux/Schnakenberg-ST-Project30-R9/actions/runs/39146108992 .
- R 4.3.3; DESeq2 1.42.0; optimized single-thread OpenBLAS.
- Full source matrix: 64,253 rows x 168 columns; 15,363 rows retained at counts >=10 in >=84 samples; all seven targets retained.
- Design ~ donor + treatment: rank 85, residual df 83. Batch is constant within each pair and absorbed by donor effects.
- Median-ratio factors estimated from the full matrix before filtering; parametric mean-dispersion trend; unshrunken treatment coefficients; no count replacement; independent filtering disabled.
- All 15,363 model fits converged; no missing P values; no target Cook's flags at the recorded F(85,83) 99th-percentile cutoff.
- The first attempt stopped at its 20-minute limit during dispersion estimation and did not yield completed NB results. The rerun changed runtime/BLAS only; the R script was byte-identical.

## Verification and limits
A separate Python NB GLM fit, conditional on DESeq2's estimated dispersions and independently reproduced size factors, agreed with all seven log2 coefficients within 1.833e-5 and SEs within 1.052e-9 (pre-set tolerances 1e-4). This is conditional numerical verification, not independent dispersion estimation or biological validation. Six-target BH values were independently checked. Original transformed means, median, positive counts and leave-one-donor-out bounds were reconstructed; historical bootstrap intervals were not rerun.

Original-source gzip SHA-256:
- Counts: f07c0230c8824b91cb73c28c0d544df989073a38172a2506167692cb6fb460d1
- Sample metadata: 0c1bb7564e7c943ab68a4ec9a486da473a44d92a89081b9a5cd001510afcf747

Original investigator-deposited HTSeq gene counts are available from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97358 . They are not redistributed in this repository; FASTQ alignment was not repeated. Exact donor-level outputs, complete model diagnostics, software records, verification scripts and the all-donor figure generator accompany the ED4 manuscript source package.

This analysis closes a specific statistical/reproducibility concern. It does not establish mathematical novelty, an in-vivo delay, physical diffusivity, or a clinical fibrosis biomarker. Exact ED4 author approval and journal submission remain separate steps.
