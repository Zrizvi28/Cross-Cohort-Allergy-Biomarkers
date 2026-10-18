#!/usr/bin/env Rscript
# POST HOC sensitivity checks (not preregistered):
#   1. Infant ancestry strata by allergy status (descriptive)
#   2. How many top-100 / concordant probes are flagged by Zhou MASK_general (descriptive)
#   3. Primary concordance test rerun on the masked probe set (same k, permutations, seed)
# Usage (repo root):
#   N_PERM=5 Rscript scripts/sensitivity_masked.R   # smoke test: checks pipeline, hides step 3 result
#   Rscript scripts/sensitivity_masked.R            # full run (1,000 permutations)

suppressPackageStartupMessages(library(limma))

data_dir <- here::here("data")
res_dir  <- here::here("results")
KS     <- c(50, 100, 200, 500)
N_PERM <- as.integer(Sys.getenv("N_PERM", "1000"))
SEED   <- 42
SMOKE  <- N_PERM < 1000

m_inf    <- readRDS(file.path(data_dir, "mVals_114134_shared.rds"))
meta_inf <- readRDS(file.path(data_dir, "meta_114134_baseline.rds"))
tt_ado   <- readRDS(file.path(data_dir, "tt_189148.rds"))       # same adolescent fit as primary
m_ado    <- readRDS(file.path(data_dir, "mVals_189148_shared.rds"))
meta_ado <- read.csv(file.path(data_dir, "GSE189148_resting_sample_metadata.csv"),
                     stringsAsFactors = FALSE)

stopifnot(setequal(colnames(m_inf), rownames(meta_inf)))
meta_inf    <- meta_inf[colnames(m_inf), ]
allergy_inf <- factor(meta_inf$`allergy status:ch1`, levels = c("control", "allergic"))
ado_gsm     <- sub("_.*$", "", colnames(m_ado))
meta_ado    <- meta_ado[match(ado_gsm, meta_ado$gsm), ]
allergy_ado <- factor(ifelse(meta_ado$allergy_status == "control", "control", "allergic"),
                      levels = c("control", "allergic"))
stopifnot(!anyNA(allergy_inf), !anyNA(allergy_ado))

# ---------- 1. Infant ancestry by allergy ----------
anc_col <- grep("ancestry", colnames(meta_inf), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(anc_col)) {
  anc_tab <- table(allergy = allergy_inf, ancestry = meta_inf[[anc_col]])
  cat("\n== Infant", anc_col, "by allergy status ==\n"); print(anc_tab)
  print(round(100 * prop.table(anc_tab, 1), 1))
  cat("Fisher's exact p (descriptive):",
      signif(fisher.test(anc_tab, simulate.p.value = TRUE, B = 1e4)$p.value, 3), "\n")
  write.csv(as.data.frame(anc_tab), file.path(res_dir, "infant_ancestry_by_allergy.csv"),
            row.names = FALSE)
} else {
  cat("\nNo ancestry column found in infant metadata\n")
}

# ---------- Masked probe set (same as the controls) ----------
masked_probe_set <- function() {
  man <- read.delim(file.path(data_dir, "EPIC.hg19.manifest.tsv.gz"), stringsAsFactors = FALSE)
  mask_files <- file.path(data_dir, c("EPIC.hg19.mask.tsv.gz", "EPIC.hg38.mask.tsv.gz"))
  src <- read.delim(mask_files[file.exists(mask_files)][1], stringsAsFactors = FALSE)
  src_id  <- intersect(c("probeID", "Probe_ID"), colnames(src))[1]
  man_id  <- intersect(c("probeID", "Probe_ID"), colnames(man))[1]
  gen_col <- grep("general", colnames(src), ignore.case = TRUE, value = TRUE)[1]
  masked    <- src[[src_id]][as.logical(src[[gen_col]]) %in% TRUE]
  autosomal <- man[[man_id]][!is.na(man$CpG_chrm) & !(man$CpG_chrm %in% c("chrX", "chrY"))]
  list(keep = intersect(rownames(m_inf), setdiff(autosomal, masked)), masked = masked)
}
mp <- masked_probe_set()
message("Masked probe set: ", length(mp$keep), " of ", nrow(m_inf))

# ---------- 2. Masked probes among the primary top 100 ----------
top_tab <- read.csv(file.path(res_dir, "top100_probe_table.csv"), stringsAsFactors = FALSE)
top_tab$masked <- top_tab$probe %in% mp$masked
cat("\n== MASK_general among primary top 100 ==\n")
print(table(masked = top_tab$masked, concordant = top_tab$concordant))
write.csv(top_tab[, c("rank", "probe", "gene", "island", "concordant", "masked")],
          file.path(res_dir, "top100_masked_flags.csv"), row.names = FALSE)

# ---------- 3. Primary test on masked probes ----------
probes <- mp$keep
# Put tt_189148 on the same allergic - control scale as the infant fit
chk <- head(probes, 20000)
raw <- rowMeans(m_ado[chk, allergy_ado == "allergic"]) - rowMeans(m_ado[chk, allergy_ado == "control"])
flip <- cor(tt_ado[chk, "logFC"], raw) < 0
ado_fc <- setNames(tt_ado[probes, "logFC"] * ifelse(flip, -1, 1), probes)

m_use <- m_inf[probes, ]
conc_at_k <- function(lab) {
  tt  <- topTable(eBayes(lmFit(m_use, model.matrix(~ lab))), coef = 2, number = Inf,
                  sort.by = "none")
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
  message("Smoke test OK (", N_PERM, " permutations). Step 3 result hidden.")
} else {
  p <- vapply(seq_along(KS), function(j) (sum(null[, j] >= observed[j]) + 1) / (N_PERM + 1),
              numeric(1))
  res <- data.frame(analysis = "posthoc_masked_primary", k = KS, n_probes = length(probes),
                    observed = observed, null_mean = colMeans(null),
                    null_2.5 = apply(null, 2, quantile, 0.025),
                    null_97.5 = apply(null, 2, quantile, 0.975),
                    p_one_sided = p, n_perm = N_PERM, seed = SEED)
  write.csv(res, file.path(res_dir, "sensitivity_masked_primary.csv"), row.names = FALSE)
  saveRDS(null, file.path(res_dir, "sensitivity_masked_primary_null.rds"))
  cat("\n== Primary test on masked probes (POST HOC) ==\n"); print(res, row.names = FALSE)
}
