#!/usr/bin/env Rscript
# Addition 1 sex-based controls (Preregistration.md, commit ea7018a)
# Usage (from the repo root):
#   Rscript scripts/concordance_controls.R table      # Step 1: infant sex x allergy (run first)
#   Rscript scripts/concordance_controls.R positive   # positive control: infant sex vs adolescent sex
#   Rscript scripts/concordance_controls.R negative   # negative control: infant allergy vs adolescent sex
#   Rscript scripts/concordance_controls.R adjusted   # sex-adjusted rerun (only if triggered)
# Smoke test with few permutations first:  N_PERM=5 Rscript scripts/sex_controls.R positive

suppressPackageStartupMessages(library(limma))

args     <- commandArgs(trailingOnly = TRUE)
step     <- if (length(args) > 0) args[1] else "table"
data_dir <- here::here("data")
res_dir  <- here::here("results")
K        <- 100
N_PERM   <- as.integer(Sys.getenv("N_PERM", "1000"))
SEED     <- 42

# ---------- Inputs ----------
m_inf    <- readRDS(file.path(data_dir, "mVals_114134_shared.rds"))
m_ado    <- readRDS(file.path(data_dir, "mVals_189148_shared.rds"))
meta_inf <- readRDS(file.path(data_dir, "meta_114134_baseline.rds"))
meta_ado <- read.csv(file.path(data_dir, "GSE189148_resting_sample_metadata.csv"),
                     stringsAsFactors = FALSE)
scan_dates <- readRDS(file.path(data_dir, "scan_dates_189148.rds"))

# Line metadata up with matrix columns (stop loudly if they can't be matched)
stopifnot(setequal(colnames(m_inf), rownames(meta_inf)))
meta_inf <- meta_inf[colnames(m_inf), ]
ado_gsm  <- sub("_.*$", "", colnames(m_ado))
meta_ado <- meta_ado[match(ado_gsm, meta_ado$gsm), ]
stopifnot(identical(meta_ado$gsm, ado_gsm))

# Same factor coding in both cohorts so logFC signs are comparable
allergy_inf <- factor(meta_inf$`allergy status:ch1`, levels = c("control", "allergic"))
sex_inf     <- factor(ifelse(tolower(meta_inf$`Sex:ch1`) %in% c("male", "m"), "M", "F"),
                      levels = c("F", "M"))
allergy_ado <- factor(ifelse(meta_ado$allergy_status == "control", "control", "allergic"),
                      levels = c("control", "allergic"))
sex_ado     <- factor(meta_ado$Sex, levels = c("F", "M"))
scan_ado    <- factor(scan_dates[meta_ado$Basename])

stopifnot(!anyNA(allergy_inf), !anyNA(sex_inf), !anyNA(allergy_ado),
          !anyNA(sex_ado), !anyNA(scan_ado),
          length(allergy_inf) == 59, length(allergy_ado) == 43)

# ---------- Helpers ----------
fit_coef <- function(m, design) {
  fit <- eBayes(lmFit(m, design))
  topTable(fit, coef = 2, number = Inf, sort.by = "none")[, c("logFC", "P.Value")]
}

concordance <- function(inf_tt, ado_logfc, k = K) {
  top <- rownames(inf_tt)[order(inf_tt$P.Value)][seq_len(k)]
  mean(sign(inf_tt[top, "logFC"]) == sign(ado_logfc[top]))
}

run_test <- function(m_use, labels, covar, ado_logfc, sides, name) {
  design_for <- function(lab) {
    if (is.null(covar)) model.matrix(~ lab) else model.matrix(~ lab + covar)
  }
  observed <- concordance(fit_coef(m_use, design_for(labels)), ado_logfc)
  set.seed(SEED)
  null <- vapply(seq_len(N_PERM), function(i) {
    if (i %% 50 == 0) message(name, ": permutation ", i, "  ", Sys.time())
    concordance(fit_coef(m_use, design_for(sample(labels))), ado_logfc)
  }, numeric(1))
  null_mean <- mean(null)
  p <- if (sides == 1) {
    (sum(null >= observed) + 1) / (N_PERM + 1)
  } else {
    (sum(abs(null - null_mean) >= abs(observed - null_mean)) + 1) / (N_PERM + 1)
  }
  out <- data.frame(test = name, k = K, n_probes = nrow(m_use), observed = observed,
                    null_mean = null_mean,
                    null_2.5 = quantile(null, 0.025), null_97.5 = quantile(null, 0.975),
                    p = p, sides = sides, n_perm = N_PERM, seed = SEED, row.names = NULL)
  write.csv(out, file.path(res_dir, paste0("control_", name, ".csv")), row.names = FALSE)
  saveRDS(null, file.path(res_dir, paste0("control_", name, "_null.rds")))
  print(out)
  invisible(out)
}

# Autosomal probes minus Zhou et al. MASK_general (EPIC hg19)
# Chromosomes come from the manifest; the general mask comes from the manifest
# if it has one (legacy format), otherwise from the separate mask file.
masked_probe_set <- function() {
  man <- read.delim(file.path(data_dir, "EPIC.hg19.manifest.tsv.gz"),
                    stringsAsFactors = FALSE)
  man_id <- intersect(c("probeID", "Probe_ID"), colnames(man))[1]
  src <- man
  if (!any(grepl("general", colnames(man), ignore.case = TRUE))) {
    mask_files <- file.path(data_dir, c("EPIC.hg19.mask.tsv.gz", "EPIC.hg38.mask.tsv.gz"))
    mask_file  <- mask_files[file.exists(mask_files)][1]
    stopifnot(!is.na(mask_file))
    message("Using mask file: ", basename(mask_file))
    src <- read.delim(mask_file, stringsAsFactors = FALSE)
  }
  src_id  <- intersect(c("probeID", "Probe_ID"), colnames(src))[1]
  gen_col <- grep("general", colnames(src), ignore.case = TRUE, value = TRUE)[1]
  stopifnot(!is.na(man_id), !is.na(src_id), !is.na(gen_col),
            "CpG_chrm" %in% colnames(man))
  message("Using mask column: ", gen_col)

  masked    <- src[[src_id]][as.logical(src[[gen_col]]) %in% TRUE]
  autosomal <- man[[man_id]][!is.na(man$CpG_chrm) &
                             !(man$CpG_chrm %in% c("chrX", "chrY"))]
  probes <- intersect(rownames(m_inf), setdiff(autosomal, masked))
  message("Masked autosomal probe set: ", length(probes), " of ", nrow(m_inf))
  stopifnot(length(probes) > 0.5 * nrow(m_inf))  # catches probe-ID format mismatches
  probes
}

# ---------- Steps ----------
if (step == "table") {
  tab      <- table(allergy = allergy_inf, sex = sex_inf)
  pct_male <- 100 * prop.table(tab, 1)[, "M"]
  diff_pp  <- abs(pct_male["allergic"] - pct_male["control"])
  print(tab)
  print(round(pct_male, 1))
  cat("Difference in % male:", round(diff_pp, 1), "percentage points\n")
  cat("Fisher's exact p (descriptive only):", signif(fisher.test(tab)$p.value, 3), "\n")
  cat("Sex-adjusted rerun triggered (>= 10 pp):", diff_pp >= 10, "\n")
  write.csv(as.data.frame(tab), file.path(res_dir, "infant_sex_by_allergy.csv"),
            row.names = FALSE)

} else if (step == "positive") {
  probes  <- masked_probe_set()
  ado_sex <- fit_coef(m_ado, model.matrix(~ sex_ado + scan_ado))
  run_test(m_inf[probes, ], sex_inf, NULL,
           setNames(ado_sex[probes, "logFC"], probes), sides = 1, name = "positive_sex_sex")

} else if (step == "negative") {
  probes  <- masked_probe_set()
  ado_sex <- fit_coef(m_ado, model.matrix(~ sex_ado + scan_ado))
  run_test(m_inf[probes, ], allergy_inf, NULL,
           setNames(ado_sex[probes, "logFC"], probes), sides = 2,
           name = "negative_allergy_sex")

} else if (step == "adjusted") {
  ado_all <- fit_coef(m_ado, model.matrix(~ allergy_ado + sex_ado + scan_ado))
  run_test(m_inf, allergy_inf, sex_inf,
           setNames(ado_all$logFC, rownames(ado_all)), sides = 1,
           name = "sex_adjusted_primary")

} else {
  stop("Unknown step: ", step, " (use table, positive, negative or adjusted)")
}
