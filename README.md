# Food-Allergy-Biomarker-Statistical-Replication

This project investigates possible genetic differences associated with food allergy using publicly available DNA methylation and RNA-seq datasets. Rather than focusing on classification, the project examines whether candidate biomarkers show consistent statistical and biological evidence across independent, age-differentiated cohorts.

## Motivation

Food allergy results from a complex interaction between the immune system, genetics, and environmental factors. DNA methylation provides one way to investigate molecular changes that may be associated with allergic disease.

I was interested in whether methylation signals in resting, naive CD4+ T cells from one food-allergy cohort could also be observed in an independent cohort with a different age group, and whether those signals could be supported by changes in gene expression. Most prior computational biomarker work in this space (e.g. Alag 2019, Kilicarslan 2025) validates a candidate on a single cohort. This project examines two cohorts and two types of genetic data.

## Research Question

Which DNA methylation differences are associated with food allergy, and which candidate signals show evidence of consistency across independent cohorts and independent data types?

The analysis uses three GEO datasets:

* GSE114134: infant DNA methylation data (n=59)
* GSE189148: adolescent resting-state DNA methylation data (n=88 loaded, n=43 after QC/label restriction)
* GSE189149: adolescent RNA-seq data (n=48)

## Methods

* Quality control: detection p-value filtering, funnorm normalization, probe filtering, PCA-based covariate screening (chip, array position, sex, scan date, allergy status) with sex-linked probes removed before re-testing
* Sample identity verification: SNP-fingerprint correlation used to check all samples for mislabeling
* Differential methylation: limma-family testing plus bumphunter regional analysis, with permutation-based FWER correction matched across both cohorts (same permutation count used in each)
* Systematic genome-wide concordance scan: every bumphunter region in the infant cohort checked against every region in the adolescent cohort to see how many other loci show the same cross-cohort, direction-matched pattern, and where ISG15 ranks among them
* Cross-omics ranking: candidates with methylation evidence checked against independent RNA-seq differential expression 
* Specificity control: the identical pipeline re-run on a 14-gene panel of unrelated interferon-response genes (380 probes total) as a negative control to see if ISG15 is just part of general inflammation or a promising signal
* Artifact screening: the single strongest genome-wide hit by raw p-value checked against known imprinted-locus behavior before being treated as a candidate

## Results

### SNP check
SNP-fingerprint correlation confirmed one declared technical replicate, and flagged a few possible identity issues: GSM5695343 and GSM5695320/GSM5695319 showed SNP-correlation evidence inconsistent with their declared labels. Corrections are recorded in `data/GSE189148_sample_metadata_snp_annotated.csv,` but the underlying GEO sample labels were left unmodified. Everything is in the notebooks.

### Systematic concordance scan
A genome-wide scan (not a single-candidate test) found 4,392 direction-matched, cross-cohort concordant methylation regions. Within that full ranking:

* ISG15 (cg08469540, cg25610492) ranks 215th of 4,392 on methylation evidence alone, with p=0.0175 (GSE114134) and p=0.00589 (GSE189148). The region did not survive family-wise error correction (FWER=1 in both cohorts). However, this is consistent with EWAS power constraints at this sample size (Tsai & Bell, 2015 found that detecting even a 10% methylation difference at genome-wide significance requires roughly 112 subjects per group). Genome-wide significance was not a realistic bar for cohorts this size, which is why cross-cohort and cross-omics concordance act as the primary evidence for why ISG15 is a promising candidate for biomarker verification.
* RGS14 ranks 18th of 4,392 on methylation alone, and is the strongest methylation-only candidate in the scan, though not previously reported in a food-allergy context, so ISG15 is the central focus.
* Once cross-omics evidence (RNA-seq) is folded in, 1,584 candidates have both direction-matched methylation and RNA-seq support. On that combined ranking, RGS14 is 1st of 1,584 and ISG15 is 49th of 1,584 (top ~3%). 

### Specificity control
The same pipeline run on 14 unrelated interferon-response genes (MX1, IFI44L, PARP9, IFI27, IFITM1, IFIT1, IFIT3, STAT1, IRF7, ISG20, OAS2, OAS3, PSMB8, EPSTI1 — 380 probes total) found no comparable cross-cohort signal in 13 of the 14 genes. IFIT1 was the one exception, with one probe reaching nominal significance in both cohorts. However, IFIT1 did not show RNA-seq backing.

### Artifact screening
The infant cohort's single strongest genome-wide hit by raw p-value (chr5, the nc886/VTRNA2-1 locus) was checked and excluded because it was confirmed unrelated to allergy status by Fisher's exact test (p=0.585).

## Conclusions

The results suggest that the region surrounding ISG15 is a promising candidate for further investigation in food allergy, with the strongest evidence coming from cross-cohort and cross-omics concordance rather than genome-wide statistical significance, which was not achievable at this sample size. RGS14 is the systematic scan's top-ranked candidate on combined evidence and a genuinely novel candidate in this disease context, but lacks the independent replication ISG15 has and requires further validation before being treated as more than exploratory.

Taken together, this project's contribution is as much methodological as biological: a validation framework that ranks and specificity-checks candidates systematically, catches real sample identity issues and known artifact loci before they can be mistaken for findings, and reports negative and null results (13/14 specificity genes, FWER failure, RGS14's lack of external replication) with the same weight as positive ones.

## Repo Structure
```
├── notebooks/
│   ├── 1.DataPreparation.ipynb                    — initial data loading and setup
│   ├── 2.GSE189148_Loading_IDATs.ipynb            — raw IDAT loading, adolescent cohort
│   ├── 3.GSE189148_Quality_Control.ipynb          — QC, covariate screening, SNP identity verification
│   ├── 4.GSE114134_Analysis.ipynb                 — infant cohort processing and differential methylation
│   ├── 6.Followup_Analysis_(DMR_Check).ipynb      — ISG15/nc886 region-level check, artifact screening
│   ├── 7.Transcriptomic_ISG15_Check.ipynb         — RNA-seq validation of the ISG15 signal
│   ├── 8.ISG_Family_Check.ipynb                   — 14-gene specificity control panel
│   ├── 9.Methods_QC_and_Concordance_Check.ipynb   — systematic genome-wide concordance scan (RGS14 discovery)
│   └── 10.Transcriptomic_Candidate_Check.ipynb    — cross-omics ranking of all concordant candidates
├── scripts/
│   └── rerun_bootstrap.R                          — bootstrap nullMethod rerun for the adolescent cohort
├── data/                                          — sample metadata, SNP verification results, candidate lists
└── results/                                       — saved model objects (bumphunter, DMR results) and QC plots
```

## References
* Alag, A. (2019). *PLOS ONE* — 18-CpG food allergy methylation signature (single-cohort validation)
* Kilicarslan, S. et al. (2025). *Biomedicines* — ML classifiers on food allergy methylation data
* Tsai, P.C. & Bell, J.T. (2015). *International Journal of Epidemiology* — power and sample size estimation for EWAS

## Project Status
This repository is a focused continuation of my original food-allergy biomarker project. The classification component has been removed so that the project can focus on statistical biomarker analysis and biological interpretation. AI assistance was used to write the one script, debug code, and copy the necessary parts of the previous, more messy repo here.
