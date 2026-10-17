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
