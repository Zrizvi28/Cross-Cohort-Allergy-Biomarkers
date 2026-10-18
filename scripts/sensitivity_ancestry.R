#!/usr/bin/env Rscript
# POST HOC sensitivity (not preregistered): primary concordance test with infant
# ancestry stratum added as a covariate. Adolescent side unchanged (tt_189148;
# no ancestry data available for adolescents).
# Usage (repo root):
#   N_PERM=5 Rscript scripts/sensitivity_ancestry.R   # smoke test: hides results
#   Rscript scripts/sensitivity_ancestry.R            # full run (1,000 permutations)

suppressPackageStartupMessages(library(limma))

data_dir <- here::here("data")
res_dir  <- here::here("results")
KS     <- c(50, 100, 200, 500)
N_PERM <- as.integer(Sys.getenv("N_PERM", "1000"))
SEED   <- 42
SMOKE  <- N_PERM < 1000

m_inf    <- readRDS(file.path(data_dir, "mVals_114134_shared.rds"))
meta_inf <- readRDS(file.path(data_dir, "meta_114134_baseline.rds"))
m_ado    <- readRDS(file.path(data_dir, "mVals_189148_shared.rds"))
meta_ado <- read.csv(file.path(data_dir, "GSE189148_resting_sample_metadata.csv"),
                     stringsAsFactors = FALSE)
tt_ado   <- readRDS(file.path(data_dir, "tt_189148.rds"))    # same adolescent fit as primary

stopifnot(setequal(colnames(m_inf), rownames(meta_inf)))
meta_inf    <- meta_inf[colnames(m_inf), ]
allergy_inf <- factor(meta_inf$`allergy status:ch1`, levels = c("control", "allergic"))
anc_col     <- grep("ancestry", colnames(meta_inf), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(anc_col))
ancestry    <- factor(meta_inf[[anc_col]])
stopifnot(!anyNA(allergy_inf), !anyNA(ancestry), nlevels(ancestry) >= 2)
message("Ancestry covariate: ", anc_col, " (", paste(levels(ancestry), collapse = ", "), ")")

# Put tt_189148 on the same allergic - control scale as the infant fit
ado_gsm     <- sub("_.*$", "", colnames(m_ado))
meta_ado    <- meta_ado[match(ado_gsm, meta_ado$gsm), ]
allergy_ado <- factor(ifelse(meta_ado$allergy_status == "control", "control", "allergic"),
                      levels = c("control", "allergic"))
chk  <- head(rownames(m_ado), 20000)
raw  <- rowMeans(m_ado[chk, allergy_ado == "allergic"]) - rowMeans(m_ado[chk, allergy_ado == "control"])
flip <- cor(tt_ado[chk, "logFC"], raw) < 0
ado_fc <- setNames(tt_ado[rownames(m_inf), "logFC"] * ifelse(flip, -1, 1), rownames(m_inf))

# Allergy is coefficient 2; ancestry is held fixed in every permutation
conc_at_k <- function(lab) {
  tt  <- topTable(eBayes(lmFit(m_inf, model.matrix(~ lab + ancestry))), coef = 2,
                  number = Inf, sort.by = "none")
  ord <- rownames(tt)[order(tt$P.Value)]
  vapply(KS, function(k) {
    top <- ord[seq_len(k)]
    mean(sign(tt[top, "logFC"]) == sign(ado_fc[top]))
  }, numeric(1))
}

observed <- conc_at_k(allergy_inf)
set.seed(SEED)
null <- t(vapply(seq_len(N_PERM), function(i) {
  if (i %% 50 == 0) message("permutation ", i, "  ", Sys.time())
  conc_at_k(sample(allergy_inf))
}, numeric(length(KS))))

if (SMOKE) {
  message("Smoke test OK (", N_PERM, " permutations). Result hidden.")
} else {
  p <- vapply(seq_along(KS), function(j) (sum(null[, j] >= observed[j]) + 1) / (N_PERM + 1),
              numeric(1))
  res <- data.frame(analysis = "posthoc_ancestry_adjusted_primary", k = KS,
                    n_probes = nrow(m_inf), observed = observed, null_mean = colMeans(null),
                    null_2.5 = apply(null, 2, quantile, 0.025),
                    null_97.5 = apply(null, 2, quantile, 0.975),
                    p_one_sided = p, n_perm = N_PERM, seed = SEED)
  write.csv(res, file.path(res_dir, "sensitivity_ancestry_primary.csv"), row.names = FALSE)
  saveRDS(null, file.path(res_dir, "sensitivity_ancestry_primary_null.rds"))
  print(res, row.names = FALSE)
}
