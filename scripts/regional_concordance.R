#!/usr/bin/env Rscript
# Addition 2: regional concordance (Preregistration.md, commit ed8e1e9)
# Usage (repo root):
#   N_PERM=5 Rscript scripts/regional_concordance.R   # smoke test: checks the pipeline, hides results
#   Rscript scripts/regional_concordance.R            # preregistered run (1,000 permutations)

suppressPackageStartupMessages({ library(limma); library(minfi); library(Matrix) })

data_dir <- here::here("data")
res_dir  <- here::here("results")
KS       <- c(50, 100, 200, 500)
N_PERM   <- as.integer(Sys.getenv("N_PERM", "1000"))
SEED     <- 42
SMOKE    <- N_PERM < 1000
MIN_PROBES <- 2

# ---------- Inputs (all shared probes, no masking) ----------
m_inf    <- readRDS(file.path(data_dir, "mVals_114134_shared.rds"))
m_ado    <- readRDS(file.path(data_dir, "mVals_189148_shared.rds"))
meta_inf <- readRDS(file.path(data_dir, "meta_114134_baseline.rds"))
meta_ado <- read.csv(file.path(data_dir, "GSE189148_resting_sample_metadata.csv"),
                     stringsAsFactors = FALSE)
scan_dates <- readRDS(file.path(data_dir, "scan_dates_189148.rds"))
stopifnot(identical(rownames(m_inf), rownames(m_ado)))

stopifnot(setequal(colnames(m_inf), rownames(meta_inf)))
meta_inf <- meta_inf[colnames(m_inf), ]
ado_gsm  <- sub("_.*$", "", colnames(m_ado))
meta_ado <- meta_ado[match(ado_gsm, meta_ado$gsm), ]
stopifnot(identical(meta_ado$gsm, ado_gsm))

allergy_inf <- factor(meta_inf$`allergy status:ch1`, levels = c("control", "allergic"))
allergy_ado <- factor(ifelse(meta_ado$allergy_status == "control", "control", "allergic"),
                      levels = c("control", "allergic"))
scan_ado    <- factor(scan_dates[meta_ado$Basename])
stopifnot(!anyNA(allergy_inf), !anyNA(allergy_ado), !anyNA(scan_ado))

# ---------- Annotation ----------
ann <- getAnnotation(readRDS(file.path(data_dir, "grSet_funnorm.rds")))
missing <- setdiff(rownames(m_inf), rownames(ann))
message("Shared probes without annotation: ", length(missing))
ann <- ann[intersect(rownames(m_inf), rownames(ann)), ]

# Promoters: TSS1500 or TSS200, one region per gene
nm  <- strsplit(ann$UCSC_RefGene_Name, ";")
gp  <- strsplit(ann$UCSC_RefGene_Group, ";")
ok  <- lengths(nm) == lengths(gp)
message("Probes with mismatched gene/group annotation (skipped): ", sum(!ok & lengths(nm) > 0))
probe_rep <- rep(rownames(ann)[ok], lengths(nm)[ok])
gene_u    <- unlist(nm[ok]); grp_u <- unlist(gp[ok])
keep      <- grp_u %in% c("TSS1500", "TSS200")
promoter_pairs <- unique(data.frame(probe = probe_rep[keep], region = gene_u[keep],
                                    stringsAsFactors = FALSE))

# Islands: Island + N_Shore + S_Shore, grouped by Islands_Name (shelves excluded)
isl <- ann$Relation_to_Island %in% c("Island", "N_Shore", "S_Shore") & ann$Islands_Name != ""
island_pairs <- data.frame(probe = rownames(ann)[isl], region = ann$Islands_Name[isl],
                           stringsAsFactors = FALSE)

# ---------- Region values: mean M of member probes, per sample ----------
region_values <- function(pairs, m) {
  n_per  <- table(pairs$region)
  pairs  <- pairs[pairs$region %in% names(n_per)[n_per >= MIN_PROBES], ]
  regs   <- sort(unique(pairs$region))
  w      <- sparseMatrix(i = match(pairs$region, regs), j = match(pairs$probe, rownames(m)),
                         x = 1 / as.numeric(n_per[pairs$region]),
                         dims = c(length(regs), nrow(m)))
  out <- as.matrix(w %*% m)
  dimnames(out) <- list(regs, colnames(m))
  attr(out, "n_probes") <- as.integer(n_per[regs])
  out
}

# ---------- Test ----------
fit_coef <- function(m, design) {
  fit <- eBayes(lmFit(m, design))
  topTable(fit, coef = 2, number = Inf, sort.by = "none")[, c("logFC", "P.Value")]
}
conc_at_k <- function(inf_tt, ado_logfc) {
  ord <- rownames(inf_tt)[order(inf_tt$P.Value)]
  vapply(KS, function(k) {
    top <- ord[seq_len(k)]
    mean(sign(inf_tt[top, "logFC"]) == sign(ado_logfc[top]))
  }, numeric(1))
}

run_set <- function(pairs, set_name) {
  r_inf <- region_values(pairs, m_inf)
  r_ado <- region_values(pairs, m_ado)
  stopifnot(identical(rownames(r_inf), rownames(r_ado)))
  message(set_name, ": ", nrow(r_inf), " regions with >= ", MIN_PROBES, " probes")

  ado_tt  <- fit_coef(r_ado, model.matrix(~ allergy_ado + scan_ado))    # fit once
  ado_fc  <- setNames(ado_tt$logFC, rownames(ado_tt))
  inf_tt  <- fit_coef(r_inf, model.matrix(~ allergy_inf))
  observed <- conc_at_k(inf_tt, ado_fc)

  set.seed(SEED)
  null <- t(vapply(seq_len(N_PERM), function(i) {
    if (i %% 100 == 0) message(set_name, ": permutation ", i, "  ", Sys.time())
    conc_at_k(fit_coef(r_inf, model.matrix(~ sample(allergy_inf))), ado_fc)
  }, numeric(length(KS))))

  if (SMOKE) {
    message(set_name, ": smoke test OK (", N_PERM, " permutations). Results hidden.")
    return(invisible(NULL))
  }

  p <- vapply(seq_along(KS), function(j) (sum(null[, j] >= observed[j]) + 1) / (N_PERM + 1),
              numeric(1))
  res <- data.frame(set = set_name, k = KS, n_regions = nrow(r_inf), observed = observed,
                    null_mean = colMeans(null),
                    null_2.5 = apply(null, 2, quantile, 0.025),
                    null_97.5 = apply(null, 2, quantile, 0.975),
                    p_one_sided = p,
                    primary = set_name == "promoters" & KS == 100,
                    n_perm = N_PERM, seed = SEED)
  saveRDS(null, file.path(res_dir, paste0("regional_", set_name, "_null.rds")))

  # Descriptive: top 100 regions
  top <- rownames(inf_tt)[order(inf_tt$P.Value)][1:100]
  desc <- data.frame(rank = 1:100, region = top,
                     n_probes = attr(r_inf, "n_probes")[match(top, rownames(r_inf))],
                     inf_logFC = inf_tt[top, "logFC"], inf_p = inf_tt[top, "P.Value"],
                     ado_logFC = ado_tt[top, "logFC"], ado_p = ado_tt[top, "P.Value"],
                     concordant = sign(inf_tt[top, "logFC"]) == sign(ado_tt[top, "logFC"]))
  write.csv(desc, file.path(res_dir, paste0("regional_", set_name, "_top100.csv")),
            row.names = FALSE)
  if (set_name == "promoters") {
    cat("ISG15 promoter rank among", nrow(inf_tt), "promoters:",
        match("ISG15", rownames(inf_tt)[order(inf_tt$P.Value)]), "\n")
  }
  res
}

prom <- run_set(promoter_pairs, "promoters")
isl  <- run_set(island_pairs, "islands")

if (!SMOKE) {
  out <- rbind(prom, isl)
  write.csv(out, file.path(res_dir, "regional_concordance.csv"), row.names = FALSE)
  print(out, row.names = FALSE)
}
