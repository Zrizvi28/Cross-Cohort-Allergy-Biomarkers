# concordance_test.R
# Run from the project root:  Rscript scripts/concordance_test.R
# Implements the primary analysis in Preregistration.md.
# AI disclosure: drafted with Claude (Anthropic); parameters, methodology, and decision criteria by Zamin Rizvi
library(limma)

K_VALUES  <- c(50, 100, 200, 500)
K_PRIMARY <- 100
N_PERM    <- 1000     
SEED      <- 42

data_dir <- here::here("data")
m_inf    <- readRDS(file.path(data_dir, "mVals_114134_shared.rds"))
tt_adol  <- readRDS(file.path(data_dir, "tt_189148.rds"))
status   <- readRDS(file.path(data_dir, "allergy_status.rds"))
batch    <- readRDS(file.path(data_dir, "batch.rds"))
y_inf    <- droplevels(status[batch == "GSE114134"])
stopifnot(length(y_inf) == ncol(m_inf), identical(rownames(m_inf), rownames(tt_adol)))
cat("Infants:", ncol(m_inf), "| probes:", nrow(m_inf), "\n"); print(table(y_inf))
rm(status, batch); invisible(gc())

concordance <- function(y) {
  fit <- eBayes(lmFit(m_inf, model.matrix(~ y)))
  t   <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  ord <- order(t$P.Value)
  sapply(K_VALUES, function(k) mean(sign(t$logFC[ord[1:k]]) == sign(tt_adol$logFC[ord[1:k]])))
}

observed <- concordance(y_inf)
cat("Observed concordance:", round(observed, 3), "\n")

set.seed(SEED)
null <- matrix(NA, N_PERM, length(K_VALUES))
for (i in seq_len(N_PERM)) {
  null[i, ] <- concordance(sample(y_inf))
  if (i %% 50 == 0 || i == N_PERM) { cat("permutation", i, format(Sys.time()), "\n"); invisible(gc()) }
}

p_one <- sapply(seq_along(K_VALUES), function(j) (sum(null[, j] >= observed[j]) + 1) / (N_PERM + 1))

results <- data.frame(k = K_VALUES, observed = observed, null_mean = colMeans(null),
                      null_2.5 = apply(null, 2, quantile, 0.025),
                      null_97.5 = apply(null, 2, quantile, 0.975),
                      p_one_sided = p_one, primary = K_VALUES == K_PRIMARY, n_perm = N_PERM)
print(results)
write.csv(results, here::here("results", paste0("concordance_permutation_n", N_PERM, ".csv")), row.names = FALSE)
saveRDS(null, file.path(data_dir, paste0("concordance_null_n", N_PERM, ".rds")))