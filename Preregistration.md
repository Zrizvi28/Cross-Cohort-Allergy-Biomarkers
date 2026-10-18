## Primary Analysis: cross-cohort directional concordance

**Hypothesis:** Among the infant cohort's most strongly allergy-associated probes, the direction of effect in the adolescent cohort agrees more often than chance.

**Inputs:**
- Infant M-values: data/mVals_114134_shared.rds (n = 59; 39 allergic, 20 control)
- Adolescent results: data/tt_189148.rds (n = 43; 28 allergic, 15 control)
- 751,107 probes shared by both cohorts

**Procedure:**
1. Fit limma (~ allergy status) on infant M-values; rank probes by p-value.
2. For the top k probes, concordance = fraction with the same logFC sign in
   adolescents
3. Permute infant allergy labels 1,000 times (seed 42), repeating step
   1–2 each time with adolescent results
4. One-sided p = (number of permutations with concordance ≥ observed + 1) /
   (1,000 + 1)

**Primary test:** k = 100.
**Decision rule:** p < 0.05 → shared signal detected ("yes")
p ≥ 0.05 → no detectable shared signal ("no")

**Sensitivity analyses:** k = 50, 200, 500. Reported here as support for or against the stability of the primary result, but do not affect primary result

**Secondary (descriptive):** reverse direction (adolescent top 100 checked in
infants), observed concordance only

**Disclosure:** The parameters above (k values,
permutation count, seed, one-sided test, decision rule) were fixed before
running this test. Infant samples are assumed to be exchangeable under the null (no technical
batch aligned with allergy status).

## Addition 1: Sex-based control analyses
Date: 2026-10-17 Primary result at time of writing: commit 9d199b7 (k = 100, concordance 0.67 vs null mean 0.517, one-sided p = 0.012)

**Timing and Disclosure**
Controls were decided shortly after primary result. Adolescent sex distribution (25 M, 18 F) observed prior to analysis. 60% male for allergic, 57% for healthy controls. Infant cohort sex distribution not yet known.

**Scope**
Results do not change primary analysis but outcomes are reported here to educate proper intepretation.

**1. Positive Control**
Question: Can concordance test at k=100 detect a known shared signal across cohorts?
Procedure: Identical to original, sex replacing allergy status.
Probe set: Autosomal probes only, after also removing probes flagged as cross-reactive or SNP-affected (Zhou et al. 2017 EPIC MASK_general list; Chen et al. 2013). This narrower probe set differs from the primary analysis. However, this control only tests the method's sensitivity, not the primary result directly. A one-sided p<0.05 is a pass.

**2. Negative Control**
Question: Is the primary concordance at k=100 driven by sex-related rather than allergy-related outcomes?
Procedure: Same as original analysis, but fit limma on sex status for the adolescent cohort and allergy for the infant. 
Test: Two-sided p = (number of permutations with |concordance − null mean| ≥ |observed − null mean| + 1) / (1,000 + 1).
Pass: Two-sided p>0.05

**Sex-adjusted analysis**
Next step if:
1. Infant sex is unbalanced by allergy status, defined as difference of >=10 percentage points in % male between allergic and healthy controls
2. Confound control is found to be significant
Procedure: Rerun with sex added as a covariate in both cohorts.

**References**
Chen Y-a et al. (2013). Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray. Epigenetics. https://doi.org/10.4161/epi.23470
Zhou W, Laird PW, Shen H (2017). Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. Nucleic Acids Research. https://doi.org/10.1093/nar/gkw967
Singmann P et al. (2015). Characterization of whole-genome autosomal differences of DNA methylation between men and women. Epigenetics & Chromatin. https://doi.org/10.1186/s13072-015-0035-3
Almstrup K et al. (2016). Pubertal development in healthy children is mirrored by DNA methylation patterns in peripheral blood. Scientific Reports. https://doi.org/10.1038/srep28657
Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010). Negative controls: a tool for detecting confounding and bias in observational studies. Epidemiology. https://doi.org/10.1097/EDE.0b013e3181d61eeb
Smyth GK et al. limma User's Guide (Bioconductor). https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf





