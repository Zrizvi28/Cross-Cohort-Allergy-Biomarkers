# Cross Cohort Food Allergy Biomarkers

Investigates possible genetic differences in resting, naive CD4 cells associated with food allergy using publicly available DNA methylation and RNA-seq datasets.

## Data:

* GSE114134: infant DNA methylation data (n=59)
* GSE189148: adolescent resting-state DNA methylation data (n=88 loaded, n=43 after QC/label restriction)
* GSE189149: adolescent RNA-seq data (n=48)

## Results
A genome-wide scan found 4,392 direction-matched, cross-cohort concordant methylation regions. Within that full ranking:

* ISG15 (cg08469540, cg25610492) ranks 215th of 4,392 on methylation evidence alone, with p=0.0175 (GSE114134) and p=0.00589 (GSE189148) 
* RGS14 ranks 18th of 4,392 on methylation alone, and is the strongest methylation-only candidate in the scan, though not previously reported in a food-allergy context, so ISG15 is the central focus
* Once cross-omics evidence (RNA-seq) is folded in, 1,584 candidates have both direction-matched methylation and RNA-seq support. On that combined ranking, RGS14 is 1st of 1,584 and ISG15 is 49th of 1,584 (top ~3%)
  
The same pipeline run on 14 unrelated interferon-response genes (MX1, IFI44L, PARP9, IFI27, IFITM1, IFIT1, IFIT3, STAT1, IRF7, ISG20, OAS2, OAS3, PSMB8, EPSTI1) found no comparable cross-cohort signal in 13 of the 14 genes. IFIT1 was the one exception, with one probe reaching nominal significance in both cohorts. However, IFIT1 did not show RNA-seq backing, so ISG15 is in the clear.

## References
* Alag, A. (2019)
* Kilicarslan, S. et al. (2025)
* Tsai, P.C. & Bell, J.T. (2015)

## Project Info
This repository is a focused continuation of my original food-allergy biomarker project. The classification component has been removed because it was pretty messy. AI assistance was used to write one script, debug code, and copy the necessary parts of the previous, more messy repo here. All thinking and the actual code was the author's work.
