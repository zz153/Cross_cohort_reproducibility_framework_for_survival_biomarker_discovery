# ============================================================
# 07_meta_permutation_and_candidate_validation.R
# Permutation validation and candidate gene meta-analysis
# ============================================================
# Description:
#   1️⃣ Performs meta-analysis of per-gene Cox models across TCGA and CGGA.
#   2️⃣ Uses permutation of survival outcomes to estimate empirical discovery rate.
#   3️⃣ Summarizes results for a defined candidate gene set (e.g., top prognostic lncRNAs).
#
# Inputs:
#   - tcga_surv_input, cgga_surv_input : survival data with OS, Censor, and gene columns
#   - candidate gene IDs (ENSEMBL)
#
# Outputs:
#   - Figure7_meta_observed_all_genes.csv ........ All meta-analysis results
#   - Figure7_meta_permutation_counts_candidates.csv ..... Null-distribution counts
#   - Figure7_Permutation_DiscoveryHistogram_candidates.png
#   - CandidateGenes_Cox_Meta_Table.csv
#   - Figure7_CandidateGenes_CoxMeta_Table.png
#
# ============================================================

set.seed(42)
alpha <- 0.05
B <- 1000  # number of permutations
gene_start_col <- 11

# ============================================================
# 1️⃣ Identify common genes
# ============================================================
common_genes <- intersect(
  colnames(tcga_surv_input)[gene_start_col:ncol(tcga_surv_input)],
  colnames(cgga_surv_input)[gene_start_col:ncol(cgga_surv_input)]
)

has_var <- function(x) sd(x, na.rm = TRUE) > 0
common_genes <- common_genes[
  sapply(common_genes,
         function(g) has_var(tcga_surv_input[[g]]) && has_var(cgga_surv_input[[g]]))
]
cat("Common genes used:", length(common_genes), "\n")

# Clean inputs
tcga0 <- tcga_surv_input %>%
  dplyr::select(patient, OS, Censor, dplyr::all_of(common_genes)) %>%
  dplyr::filter(!is.na(OS), !is.na(Censor))
cgga0 <- cgga_surv_input %>%
  dplyr::select(patient, OS, Censor, dplyr::all_of(common_genes)) %>%
  dplyr::filter(!is.na(OS), !is.na(Censor))

# ============================================================
# 2️⃣ Core functions
# ============================================================

# Cox model per gene
cox_one_gene <- function(dat, gene) {
  f <- as.formula(paste0("Surv(OS, Censor) ~ ", gene))
  fit <- try(coxph(f, data = dat), silent = TRUE)
  if (inherits(fit, "try-error") || length(coef(fit)) == 0)
    return(tibble(logHR = NA, SE = NA, p = NA))
  ctab <- summary(fit)$coef[1, ]
  tibble(
    logHR = unname(ctab["coef"]),
    SE    = unname(ctab["se(coef)"]),
    p     = unname(ctab["Pr(>|z|)"])
  )
}

# Random-effects meta between TCGA and CGGA
meta_one_gene <- function(g, tcga, cgga) {
  tc <- cox_one_gene(tcga, g); cg <- cox_one_gene(cgga, g)
  if (any(is.na(tc)) || any(is.na(cg)))
    return(tibble(gene = g, logHR_meta = NA, SE_meta = NA, p_meta = NA, I2 = NA))
  m <- try(rma.uni(yi = c(tc$logHR, cg$logHR),
                   sei = c(tc$SE, cg$SE),
                   method = "DL"), silent = TRUE)
  if (inherits(m, "try-error"))
    return(tibble(gene = g, logHR_meta = NA, SE_meta = NA, p_meta = NA, I2 = NA))
  tibble(
    gene = g,
    logHR_meta = as.numeric(m$b),
    SE_meta = as.numeric(m$se),
    p_meta = as.numeric(m$pval),
    I2 = m$I2
  )
}

# Wrapper
run_meta_all_genes <- function(genes, tcga, cgga) {
  map_dfr(genes, ~meta_one_gene(.x, tcga, cgga)) %>%
    mutate(FDR = p.adjust(p_meta, method = "BH"))
}

# ============================================================
# 3️⃣ Observed meta-analysis
# ============================================================
meta_obs <- run_meta_all_genes(common_genes, tcga0, cgga0)
write_csv(meta_obs, "Figure7_meta_observed_all_genes.csv")

obs_sig <- meta_obs %>% filter(FDR < alpha)
cat("Observed significant meta-genes (FDR <", alpha, "):", nrow(obs_sig), "\n")

write_csv(obs_sig, "Figure7_meta_observed_significant.csv")

# ============================================================
# 4️⃣ Permutation baseline for candidate genes
# ============================================================
candidate_genes <- c("ENSG00000105997","ENSG00000128710",
                     "ENSG00000128713","ENSG00000253552",
                     "ENSG00000250133")

permute_surv <- function(df) {
  idx <- sample.int(nrow(df))
  df %>% mutate(OS = OS[idx], Censor = Censor[idx])
}

perm_counts <- numeric(B)

for (b in seq_len(B)) {
  tcga_b <- permute_surv(tcga0)
  cgga_b <- permute_surv(cgga0)
  meta_b <- run_meta_all_genes(candidate_genes, tcga_b, cgga_b)
  perm_counts[b] <- sum(meta_b$FDR < alpha, na.rm = TRUE)
  if (b %% 10 == 0) cat("Permutation", b, ":", perm_counts[b], "discoveries\n")
}

write_csv(data.frame(iter = 1:B, n_sig = perm_counts),
          "Figure7_meta_permutation_counts_candidates.csv")

# ============================================================
# 5️⃣ Empirical significance
# ============================================================
meta_obs_cand <- run_meta_all_genes(candidate_genes, tcga0, cgga0)
n_obs <- sum(meta_obs_cand$FDR < alpha, na.rm = TRUE)
emp_p <- (sum(perm_counts >= n_obs) + 1) / (B + 1)

cat("Observed candidate discoveries:", n_obs, "\n")
cat("Mean null discoveries:", mean(perm_counts), "±", sd(perm_counts), "\n")
cat("Empirical p-value:", signif(emp_p, 3), "\n")

# Histogram
png("Figure7_Permutation_DiscoveryHistogram_candidates.png", 1200, 900, res = 150)
hist(perm_counts, breaks = 20, col = "skyblue", border = "grey40",
     xlab = "# significant candidate genes (FDR < 0.05)",
     main = "Permutation baseline (candidate genes)")
abline(v = n_obs, lwd = 3, col = "red")
legend("topright",
       legend = c(paste("Observed =", n_obs),
                  paste("Mean null =", round(mean(perm_counts),1)),
                  paste("Empirical p =", signif(emp_p,3))),
       bty = "n")
dev.off()

# ============================================================
# 6️⃣ Detailed candidate table (TCGA, CGGA, pooled)
# ============================================================
cox_one <- function(df, gene) {
  f <- as.formula(paste0("Surv(OS, Censor) ~ `", gene, "`"))
  fit <- coxph(f, data = df)
  s <- summary(fit); ci <- exp(confint(fit))
  tibble(
    gene = gene,
    HR = exp(s$coef[1, "coef"]),
    HR_lo = ci[1, 1], HR_hi = ci[1, 2],
    p = s$coef[1, "Pr(>|z|)"],
    beta = s$coef[1, "coef"]
  )
}

meta_from_two <- function(tc, cg) {
  m <- rma.uni(yi = c(tc$beta, cg$beta), sei = c(tc$SE, cg$SE), method = "DL")
  tibble(
    gene = tc$gene,
    meta_HR = exp(m$b), meta_lo = exp(m$ci.lb), meta_hi = exp(m$ci.ub),
    meta_p = m$pval
  )
}

# Run for each candidate
tcga_cox <- map_dfr(candidate_genes, ~cox_one(tcga0, .x)) %>% rename_with(~paste0("TCGA_", .), -gene)
cgga_cox <- map_dfr(candidate_genes, ~cox_one(cgga0, .x)) %>% rename_with(~paste0("CGGA_", .), -gene)

names(tcga_cox)
names(cgga_cox)

# --- Standardize TCGA results ---
tcga_cox <- tcga_cox %>%
  dplyr::mutate(
    beta = TCGA_beta,
    SE   = (log(TCGA_HR_hi) - log(TCGA_HR_lo)) / (2 * 1.96),
    HR   = TCGA_HR,
    pval = TCGA_p
  ) %>%
  dplyr::select(gene, beta, SE, HR, pval)

# --- Standardize CGGA results ---
cgga_cox <- cgga_cox %>%
  dplyr::mutate(
    beta = CGGA_beta,
    SE   = (log(CGGA_HR_hi) - log(CGGA_HR_lo)) / (2 * 1.96),
    HR   = CGGA_HR,
    pval = CGGA_p
  ) %>%
  dplyr::select(gene, beta, SE, HR, pval)


meta_tbl <- map_dfr(candidate_genes, function(g) {
  tc <- tcga_cox %>% filter(gene == g)
  cg <- cgga_cox %>% filter(gene == g)
  meta_from_two(tc, cg)
})
names(tcga_cox)

tab <- tcga_cox %>%
  left_join(cgga_cox, by = "gene") %>%
  left_join(meta_tbl, by = "gene") %>%
  mutate(
    meta_FDR = p.adjust(meta_p, method = "BH"),
    dir_consist = ifelse(sign(beta.x) == sign(beta.y), "Yes", "No")
  ) %>%
  mutate(
    across(c(HR.x, HR.y, meta_HR), \(x) round(x, 2)),
    across(matches("(lo|hi)$"), \(x) round(x, 2)),
    across(matches("(_p|_FDR)$"), \(x) signif(x, 3))
  ) %>%
  transmute(
    Gene = gene,
    `TCGA HR (95% CI)` = sprintf("%0.2f", HR.x),
    `TCGA p` = pval.x,
    `CGGA HR (95% CI)` = sprintf("%0.2f", HR.y),
    `CGGA p` = pval.y,
    `Meta HR (95% CI)` = sprintf("%0.2f (%0.2f–%0.2f)", meta_HR, meta_lo, meta_hi),
    `Meta p` = meta_p,
    `Meta FDR` = meta_FDR,
    `Same direction?` = dir_consist
  )

write_csv(tab, "CandidateGenes_Cox_Meta_Table.csv")

