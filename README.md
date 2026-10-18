# Cross Cohort Food Allergy Signal

Tests whether DNA methylation differences associated with allergy in naive CD4 cells in infants transfer to a separate adolescent cohort using public cohorts and a preregistered analysis.

## Data:

* GSE114134: infant naive CD4 cell DNA methylation EPIC data (n=59, 39 allergic and 20 healthy)
* GSE189148: Adolescent resting naive CD4+ T-cell methylation (EPIC) (n=43 (28 allergic, 15 control) after QC and label restriction)
* GSE189149: Adolescent CD4+ T-cell RNA-seq (n=48, only for secondary analyses)

751,107 autosomal probes are shared by both methylation cohorts.

## Preregistration
All parameters and methodology were fixed ahead of each test in Pregistration.md

## Results
| Analysis | Concordance | Null mean | p | Outcome | Result file |
|:---|---:|---:|---:|:---|:---|
| **Primary (preregistered, k = 100)** | **0.67** | 0.517 | **0.012** | Shared signal detected | [concordance_permutation_n1000.csv](results/concordance_permutation_n1000.csv) |
| Positive control (sex vs sex) | 0.96 | 0.502 | < 0.001 | Pass | [control_positive_sex_sex.csv](results/control_positive_sex_sex.csv) |
| Negative control (allergy vs sex) | 0.30 | 0.486 | 0.007 (two-sided) | Significant | [control_negative_allergy_sex.csv](results/control_negative_allergy_sex.csv) |
| Sex-adjusted primary | 0.66 | 0.515 | 0.019 | Robust to sex | [control_sex_adjusted_primary.csv](results/control_sex_adjusted_primary.csv) |
| Regional (promoters, k = 100) | | | | *Not yet run* | |

All runs used 1000 permuations of infant labels with seed 42. The full methods and decision rules are in [Preregistration.md](Preregistration.md).

### Primary Concordance Test

Among the infant cohort's 100 most allergy-associated probes, does the direction of effect in adolescents agree more often than chance?

| k | Observed | Null mean | Null 95% range | One-sided p | Role |
|---:|---:|---:|:---:|---:|:---|
| **100** | **0.67** | 0.517 | 0.38–0.64 | **0.012** | **Primary** |
| 50 | 0.64 | 0.520 | 0.36–0.68 | 0.085 | Sensitivity |
| 200 | 0.59 | 0.514 | 0.41–0.62 | 0.080 | Sensitivity |
| 500 | 0.574 | 0.512 | 0.42–0.60 | 0.088 | Sensitivity |

**Result:** a shared directional signal was detected at the preregistered k = 100. The sensitivity values of k point the same way but do not reach p < 0.05.

**Probe level:** this is a set-level result. Among the top 100, no individual probe reaches FDR < 0.05 in adolescents (Benjamini–Hochberg across the 100). In infants, 10 probes reach FDR < 0.05 across all probes. The ISG15 probe cg08469540 ranks 2nd in infants and is direction-concordant. Full table: [top100_probe_table.csv](results/top100_probe_table.csv).

### Addition 1: sex-based controls

| Check | Result | Interpretation |
|:---|:---|:---|
| Infant sex × allergy | 46.2% male (allergic) vs 45.0% (control) | Balanced. Adolescents are also balanced (60% vs 57%). [infant_sex_by_allergy.csv](results/infant_sex_by_allergy.csv) |
| Positive control (infant sex vs adolescent sex) | 0.96, p < 0.001 | **Pass:** the test detects a known shared signal. |
| Negative control (infant allergy vs adolescent sex) | 0.30, two-sided p = 0.007 | **Significant:** triggers the sex-adjusted rerun. |
| Sex-adjusted primary (sex in both models, all probes) | 0.66, p = 0.019 | **Robust:** the primary result holds after sex adjustment. |

Both controls use autosomal probes with the Zhou et al. MASK_general list removed (741,453 probes). Sex is balanced by allergy status in both cohorts, and the primary result survives sex adjustment, so sample-composition confounding does not explain the negative control. Its below-chance concordance means that probes associated with infant allergy tend to show the opposite-signed sex difference in adolescents. One possible explanation is overlap between allergy-associated and sex-differential CpG sites (future direction).

### Addition 2: regional concordance
 
This analysis is preregistered: promoters (TSS1500 + TSS200, one region per gene) are the primary set, CpG islands with shores are the secondary set, and k = 100. It has not been run yet.

## Exploratory analyses (not preregistered)
 
The earlier work in `notebooks/` included bumphunter and DMR scans, a genome-wide ranking of direction-matched regions, RNA-seq cross-checks, and an interferon-gene panel. These analyses came before the preregistration and tried several methods, so they are secondary results.
 
* ISG15 (cg08469540, cg25610492) ranked 215th of 4,392 direction-matched regions, with nominal p in both cohorts. One ISG15 probe survived FDR in infant cohort.
* Of 14 other interferon-response genes, 13 showed no comparable cross-cohort signal (IFIT1 was the exception).
* RGS14 ranks 18th of 4,392 on methylation alone, and is the strongest methylation-only candidate in the scan (but no individual probes cleared significance in either cohort, so not the primary result)
* Once cross-omics evidence (RNA-seq) is folded in, 1,584 candidates have both direction-matched methylation and RNA-seq support

## Limitations
 
- The cohorts are small (59 and 43 samples), and no individual probe replicates at FDR < 0.05 in adolescents.
- Infant effect sizes for the top 100 are probably inflated, because the probes were selected for their strong infant signal (the winner's curse).
- The two cohorts differ in age, recruitment and processing, so concordance across them is a conservative test.
- The exploratory analyses came before the preregistration.
  
## References
 
- Martino D et al. (2018). Epigenetic dysregulation of naive CD4+ T-cell activation genes in childhood food allergy. *Nature Communications.* https://doi.org/10.1038/s41467-018-05608-4
- Zhou W, Laird PW, Shen H (2017). Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes. *Nucleic Acids Research.* https://doi.org/10.1093/nar/gkw967
- Lipsitch M, Tchetgen Tchetgen E, Cohen T (2010). Negative controls. *Epidemiology.* https://doi.org/10.1097/EDE.0b013e3181d61eeb
- Phipson B, Smyth GK (2010). Permutation p-values should never be zero. *Statistical Applications in Genetics and Molecular Biology.* https://doi.org/10.2202/1544-6115.1585
- Smyth GK et al. limma User's Guide. https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf
- Alag A (2019); Kilicarslan S et al. (2025); Tsai PC & Bell JT (2015)

**AI disclosure:** AI assistance (Claude, Anthropic) was used to edit `scripts/concordance_test.R`, `scripts/sex_controls.R` and `scripts/top100_table.R`, to debug code, and to help structure the preregistration and this README. The research questions, analysis parameters, preregistered decision rules and interpretation are the author's.
