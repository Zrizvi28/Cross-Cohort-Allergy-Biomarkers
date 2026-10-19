# Results

Every file here was produced by a script in `scripts/`. All permutation tests used 1,000 permutations of infant labels with seed 42. Analysis definitions are in [Preregistration.md](../Preregistration.md).

## Preregistered

| File | Analysis | Script |
|:---|:---|:---|
| `concordance_permutation_n1000.csv`, `concordance_log.txt` | **Primary:** directional concordance, k = 50/100/200/500 | `concordance_test.R` |
| `infant_sex_by_allergy.csv` | Addition 1: infant sex × allergy check | `sex_controls.R table` |
| `control_positive_sex_sex.csv`, `_null.rds`, `control_positive_log.txt` | Addition 1: positive control | `sex_controls.R positive` |
| `control_negative_allergy_sex.csv`, `_null.rds`, `control_negative_log.txt` | Addition 1: negative control | `sex_controls.R negative` |
| `control_sex_adjusted_primary.csv`, `_null.rds`, `control_adjusted_log.txt` | Addition 1: sex-adjusted primary | `sex_controls.R adjusted` |
| `regional_concordance.csv`, `regional_*_null.rds`, `regional_log.txt` | Addition 2: regional concordance (promoters, islands) | `regional_concordance.R` |

## Post hoc sensitivity (not preregistered)

| File | Analysis | Script |
|:---|:---|:---|
| `sensitivity_masked_primary.csv`, `_null.rds`, `sensitivity_masked_log.txt` | Primary rerun on the MASK_general probe set | `sensitivity_masked.R` |
| `top100_masked_flags.csv` | Which top-100 probes are masked (none) | `sensitivity_masked.R` |
| `infant_ancestry_by_allergy.csv` | Infant ancestry stratum × allergy | `sensitivity_masked.R` |
| `sensitivity_ancestry_primary.csv`, `_null.rds`, `sensitivity_ancestry_log.txt` | Primary rerun with infant ancestry as a covariate | `sensitivity_ancestry.R` |

## Descriptive

| File | Contents | Script |
|:---|:---|:---|
| `top100_probe_table.csv` | Top 100 infant probes: annotation, effects, p, FDR, β in both cohorts | `top100_table.R` |
| `regional_promoters_top100.csv`, `regional_islands_top100.csv` | Top 100 regions per set | `regional_concordance.R` |

## Exploratory (`exploratory/`)

These outputs come from analyses run **before** the preregistration (`notebooks/` and `scripts/exploratory/`): bumphunter and DMR scans, the interferon-gene panel, RNA-seq cross-checks, and gene-level probe statistics. They are kept for transparency and are hypothesis-generating only.
