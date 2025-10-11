# ============================================================
# 06_meta_cox_analysis_and_forest_plots.R
# Cross-cohort Cox meta-analysis and visualization
# ============================================================
# Description:
#   Performs per-gene Cox proportional hazards analysis in TCGA and CGGA cohorts,
#   combines results via random-effects meta-analysis (DerSimonian–Laird),
#   and generates forest and funnel plots for key prognostic genes.
#
# Inputs:
#   - tcga_surv_input: TCGA survival dataframe (OS, Censor, expression matrix)
#   - cgga_surv_input: CGGA survival dataframe
#   - gene_info (optional): data.frame linking Ensembl IDs ↔ gene symbols
#
# Outputs:
#   - "Figure7_MetaCox_results.csv" ........ table of pooled HRs, CI, FDR, heterogeneity
#   - "Figure7_Forest_<gene>.png" .......... forest plot per gene (TCGA + CGGA + pooled)
#   - "Figure7_Funnel_<gene>.png" .......... funnel plot (bias visualization)
#
# Dependencies:
#   Packages: dplyr, purrr, tidyr, survival, metafor, ggplot2, stringr, readr
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr)
  library(survival); library(metafor)
  library(ggplot2); library(stringr); library(readr)
})

cat("\n📈 Running cross-cohort Cox meta-analysis...\n")

# ============================================================
# 1️⃣ Helper Functions
# ============================================================

# Safe formula builder
backtick <- function(x) paste0("`", x, "`")

# Run univariate Cox model for one gene in one cohort
fit_cox_one_gene <- function(df, gene, time_col = "OS", status_col = "Censor") {
  if (!gene %in% colnames(df)) return(NULL)
  f <- as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", backtick(gene)))
  fit <- coxph(f, data = df)
  s <- summary(fit)
  
  logHR <- unname(s$coef[1, "coef"])
  se    <- unname(s$coef[1, "se(coef)"])
  HR    <- exp(logHR)
  ci    <- exp(confint(fit))
  
  tibble(
    gene = gene,
    HR = HR,
    HR_lo = ci[1, 1],
    HR_hi = ci[1, 2],
    logHR = logHR,
    SE = se,
    z = unname(s$coef[1, "z"]),
    p = unname(s$coef[1, "Pr(>|z|)"]),
    n = nrow(df),
    events = sum(df[[status_col]] == 1, na.rm = TRUE)
  )
}

# ============================================================
# 2️⃣ Run per-cohort Cox models
# ============================================================

genes_of_interest <- colnames(tcga_surv_input)[11:ncol(tcga_surv_input)]

cat("🧪 Running Cox models in TCGA and CGGA...\n")
cox_tcga <- map_dfr(genes_of_interest, ~ fit_cox_one_gene(tcga_surv_input, .x)) %>%
  mutate(cohort = "TCGA")

cox_cgga <- map_dfr(genes_of_interest, ~ fit_cox_one_gene(cgga_surv_input, .x)) %>%
  mutate(cohort = "CGGA")

cox_long <- bind_rows(cox_tcga, cox_cgga)

# ============================================================
# 3️⃣ Random-effects meta-analysis per gene
# ============================================================

meta_one_gene <- function(g) {
  dat <- cox_long %>%
    filter(gene == g, is.finite(logHR), is.finite(SE), SE > 0)
  if (nrow(dat) < 2) return(NULL)
  
  m <- rma.uni(yi = logHR, sei = SE, method = "DL", data = dat)
  
  tibble(
    gene = g,
    k = m$k,
    pooled_logHR = as.numeric(m$b),
    pooled_SE = as.numeric(m$se),
    pooled_HR = exp(m$b),
    pooled_lo = exp(m$ci.lb),
    pooled_hi = exp(m$ci.ub),
    Q = m$QE, Q_p = m$QEp,
    I2 = m$I2, tau2 = m$tau2,
    meta_z = as.numeric(m$zval),
    meta_p = as.numeric(m$pval),
    method = m$method
  )
}

cat("📊 Performing random-effects meta-analysis...\n")
meta_tbl <- map_dfr(genes_of_interest, meta_one_gene) %>%
  filter(!is.na(meta_p)) %>%
  mutate(meta_FDR = p.adjust(meta_p, method = "BH")) %>%
  arrange(meta_FDR)

write_csv(meta_tbl, "Figure7_MetaCox_results.csv")

# ============================================================
# 4️⃣ Select genes for visualization
# ============================================================

top_genes_for_plots <- meta_tbl %>%
  slice_head(n = 8) %>%
  pull(gene)

cat("📎 Selected top genes for forest/funnel plots:\n",
    paste(top_genes_for_plots, collapse = ", "), "\n")

# ============================================================
# 5️⃣ Forest plot per gene
# ============================================================

plot_forest_one_gene <- function(g, gene_info_df = NULL) {
  dat <- cox_long %>% filter(gene == g)
  if (nrow(dat) < 2) return(NULL)
  
  m <- rma.uni(yi = dat$logHR, sei = dat$SE, method = "DL")
  
  # Map Ensembl → symbol if available
  title_gene <- g
  if (!is.null(gene_info_df)) {
    sym <- gene_info_df$external_gene_name[match(g, gene_info_df$ensembl_gene_id)]
    if (!is.na(sym)) title_gene <- paste0(sym, " (", g, ")")
  }
  
  slab <- paste0(dat$cohort, " (n=", dat$n, ", events=", dat$events, ")")
  
  png(paste0("Figure7_Forest_", g, ".png"), width = 1600, height = 1200, res = 180)
  par(mar = c(4, 5, 3, 2))
  forest(
    m,
    slab = slab,
    xlab = "Hazard Ratio (log scale)",
    atransf = exp,
    at = log(c(0.5, 0.75, 1, 1.5, 2, 3)),
    ilab = cbind(round(dat$HR, 2),
                 paste0("[", round(dat$HR_lo, 2), ", ", round(dat$HR_hi, 2), "]")),
    ilab.xpos = c(-3.2, -1.5),
    header = "Study                 HR (95% CI)",
    main = paste0("Meta-analysis: ", title_gene)
  )
  abline(v = log(1), lty = 2)
  dev.off()
}

invisible(lapply(top_genes_for_plots, plot_forest_one_gene, gene_info_df = gene_info))

# ============================================================
# 6️⃣ Funnel plot and leave-one-out diagnostics
# ============================================================

plot_funnel_one_gene <- function(g) {
  dat <- cox_long %>% filter(gene == g)
  if (nrow(dat) < 2) return(NULL)
  m <- rma.uni(yi = dat$logHR, sei = dat$SE, method = "DL")
  
  png(paste0("Figure7_Funnel_", g, ".png"), width = 1200, height = 900, res = 160)
  funnel(m, xlab = "log(HR)")
  dev.off()
}

# Leave-one-out influence diagnostics
loo_one_gene <- function(g) {
  dat <- cox_long %>% filter(gene == g)
  if (nrow(dat) < 2) return(NULL)
  m <- rma.uni(yi = dat$logHR, sei = dat$SE, method = "DL")
  as.data.frame(influence(m))
}

# Example diagnostic for top gene
plot_funnel_one_gene(top_genes_for_plots[1])

cat("\n✅ Meta-analysis and forest plots completed successfully.\n")
