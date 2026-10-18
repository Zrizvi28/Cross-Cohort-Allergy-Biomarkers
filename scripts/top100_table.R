#!/usr/bin/env Rscript
# Descriptive master table for the primary analysis's top-100 infant probes.
# No new inferential test. Usage (repo root):  Rscript scripts/top100_table.R
# Output: results/top100_probe_table.csv

suppressPackageStartupMessages({ library(limma); library(minfi) })

data_dir <- here::here("data")
res_dir  <- here::here("results")
K <- 100

# ---------- Inputs ----------
m_inf    <- readRDS(file.path(data_dir, "mVals_114134_shared.rds"))
m_ado    <- readRDS(file.path(data_dir, "mVals_189148_shared.rds"))
meta_inf <- readRDS(file.path(data_dir, "meta_114134_baseline.rds"))
meta_ado <- read.csv(file.path(data_dir, "GSE189148_resting_sample_metadata.csv"),
                     stringsAsFactors = FALSE)
tt_ado   <- readRDS(file.path(data_dir, "tt_189148.rds"))   # same fit as the primary
scan_dates <- readRDS(file.path(data_dir, "scan_dates_189148.rds"))
stopifnot(all(c("logFC", "P.Value") %in% colnames(tt_ado)))

# Align metadata to matrix columns
stopifnot(setequal(colnames(m_inf), rownames(meta_inf)))
meta_inf <- meta_inf[colnames(m_inf), ]
ado_gsm  <- sub("_.*$", "", colnames(m_ado))
meta_ado <- meta_ado[match(ado_gsm, meta_ado$gsm), ]
stopifnot(identical(meta_ado$gsm, ado_gsm))

allergy_inf <- factor(meta_inf$`allergy status:ch1`, levels = c("control", "allergic"))
allergy_ado <- factor(ifelse(meta_ado$allergy_status == "control", "control", "allergic"),
                      levels = c("control", "allergic"))
sex_ado  <- factor(meta_ado$Sex, levels = c("F", "M"))
scan_ado <- factor(scan_dates[meta_ado$Basename])
stopifnot(!anyNA(allergy_inf), !anyNA(allergy_ado), !anyNA(sex_ado), !anyNA(scan_ado))

# ---------- Infant fit (all probes) and top 100 ----------
fit_inf <- eBayes(lmFit(m_inf, model.matrix(~ allergy_inf)))
tt_inf  <- topTable(fit_inf, coef = 2, number = Inf, sort.by = "none")  # allergic - control
tt_inf$FDR_all <- p.adjust(tt_inf$P.Value, method = "BH")               # FDR over all probes
top <- rownames(tt_inf)[order(tt_inf$P.Value)][seq_len(K)]

# ---------- Adolescent: make sign convention allergic - control ----------
raw_diff <- rowMeans(m_ado[top, allergy_ado == "allergic"]) -
            rowMeans(m_ado[top, allergy_ado == "control"])
flip <- cor(tt_ado[top, "logFC"], raw_diff) < 0
if (flip) message("tt_189148 is coded control - allergic; flipping its sign for the table")
ado_logFC <- tt_ado[top, "logFC"] * ifelse(flip, -1, 1)
ado_p     <- tt_ado[top, "P.Value"]

# ---------- Beta values from M-values ----------
to_beta  <- function(m) 2^m / (1 + 2^m)
grp_mean <- function(m, g, lvl) rowMeans(to_beta(m[top, g == lvl, drop = FALSE]))
inf_b_all <- grp_mean(m_inf, allergy_inf, "allergic"); inf_b_ctl <- grp_mean(m_inf, allergy_inf, "control")
ado_b_all <- grp_mean(m_ado, allergy_ado, "allergic"); ado_b_ctl <- grp_mean(m_ado, allergy_ado, "control")

# ---------- Post hoc: adolescent sex effect for the same probes ----------
fit_sex <- eBayes(lmFit(m_ado[top, ], model.matrix(~ sex_ado + scan_ado)))
sex_logFC <- topTable(fit_sex, coef = 2, number = Inf, sort.by = "none")[top, "logFC"]  # M - F

# ---------- Annotation ----------
ann <- getAnnotation(readRDS(file.path(data_dir, "grSet_funnorm.rds")))
ann <- ann[top, ]

# ---------- Assemble ----------
tab <- data.frame(
  rank = seq_len(K), probe = top,
  chr = ann$chr, pos = ann$pos,
  gene = ann$UCSC_RefGene_Name, island = ann$Relation_to_Island,
  inf_logFC = tt_inf[top, "logFC"], inf_p = tt_inf[top, "P.Value"],
  inf_FDR_all = tt_inf[top, "FDR_all"],
  inf_beta_allergic = inf_b_all, inf_beta_control = inf_b_ctl,
  inf_delta_beta = inf_b_all - inf_b_ctl,
  ado_logFC = ado_logFC, ado_p = ado_p,
  ado_FDR_top100 = p.adjust(ado_p, method = "BH"),             # FDR over these 100
  ado_beta_allergic = ado_b_all, ado_beta_control = ado_b_ctl,
  ado_delta_beta = ado_b_all - ado_b_ctl,
  concordant = sign(tt_inf[top, "logFC"]) == sign(ado_logFC),
  posthoc_ado_sex_logFC_MminusF = sex_logFC,
  posthoc_matches_sex_sign = sign(tt_inf[top, "logFC"]) == sign(sex_logFC),
  row.names = NULL, stringsAsFactors = FALSE
)

# Must reproduce the preregistered result
conc <- mean(tab$concordant)
cat("Concordance:", conc, "(primary reported 0.67)\n")
stopifnot(isTRUE(all.equal(conc, 0.67)))

cat("cg08469540 (ISG15) rank:", match("cg08469540", tab$probe), "\n")
cat("Adolescent FDR < 0.05 within top 100:", sum(tab$ado_FDR_top100 < 0.05), "\n")
cat("Infant FDR < 0.05 (all probes):", sum(tab$inf_FDR_all < 0.05), "\n")

write.csv(tab, file.path(res_dir, "top100_probe_table.csv"), row.names = FALSE)
cat("Saved results/top100_probe_table.csv\n")
