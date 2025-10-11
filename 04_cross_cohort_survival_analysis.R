# ============================================================
# 04_cross_cohort_survival_analysis.R
# Univariate and Multivariate Cox + Kaplan–Meier Reproducibility Assessment
# ============================================================
# Clinical columns present in both inputs
clin_cols <- c("patient","Gender","Age","OS","Censor",
               "Radio_status","Chemo_status",
               "IDH_mutation_status","MGMTp_methylation_status","Subtype")

# Genes = everything after clinical columns (you already harmonized order)
genes <- setdiff(colnames(tcga_surv_input), clin_cols)

# Define covariate sets (edit as needed)
covar_sets <- list(
  none        = character(0),                                    # placeholder (will map to univariate)
  basic       = c("Age","Gender"),
  treatment   = c("Age","Gender","Radio_status","Chemo_status"),
  molecular   = c("Age","Gender","IDH_mutation_status","MGMTp_methylation_status","Subtype")
)

# ---- Helpers ----
safe_median <- function(x) stats::median(x, na.rm = TRUE)

logrank_p <- function(time, event, group) {
  # group must be factor/2-level
  if (length(unique(group[!is.na(group)])) < 2) return(NA_real_)
  sd <- survdiff(Surv(time, event) ~ group)
  pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
}

logrank_p_tertile <- function(time, event, expr) {
  qs <- quantile(expr, probs = c(1/3, 2/3), na.rm = TRUE)
  grp <- cut(expr,
             breaks = c(-Inf, qs[1], qs[2], Inf),
             labels = c("Low","Mid","High"))
  # Compare only High vs Low
  grp2 <- droplevels(grp[grp != "Mid"])
  time2 <- time[grp != "Mid"]
  event2 <- event[grp != "Mid"]
  if (length(unique(grp2)) < 2) return(NA_real_)
  sd <- survdiff(Surv(time2, event2) ~ grp2)
  pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
}


fit_cox <- function(df, time_col = "OS", event_col = "Censor", rhs_terms) {
  # rhs_terms: character vector of terms appearing RHS of formula
  fml <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                           paste(rhs_terms, collapse = " + ")))
  fit <- tryCatch(coxph(fml, data = df, ties = "efron"), error = function(e) NULL)
  if (is.null(fit) || !is.finite(fit$loglik[2])) return(NULL)
  fit
}

extract_hr_p <- function(fit, term) {
  s <- tryCatch(summary(fit), error = function(e) NULL)
  if (is.null(s)) return(c(HR = NA_real_, pval = NA_real_, beta = NA_real_))
  # Handle factor expansion (e.g., termSubType_B vs baseline)
  coefs <- as.data.frame(s$coef)
  coefs$term <- rownames(coefs)
  # Keep rows that start with the term (for factors, they appear as termLEVEL)
  sub <- coefs[grepl(paste0("^", term), coefs$term), , drop = FALSE]
  if (nrow(sub) == 0) {
    # numeric gene term should match exactly
    sub <- coefs[coefs$term == term, , drop = FALSE]
  }
  if (nrow(sub) == 0) return(c(HR = NA_real_, pval = NA_real_, beta = NA_real_))
  # For numeric genes there should be one row; if multiple (shouldn't), take first
  HR   <- sub[1, "exp(coef)"]
  pval <- sub[1, "Pr(>|z|)"]
  beta <- sub[1, "coef"]
  c(HR = unname(HR), pval = unname(pval), beta = unname(beta))
}

fit_cindex <- function(fit, dd) {
  tryCatch({
    lp <- predict(fit, type = "lp")   # linear predictor
    conc <- survConcordance(Surv(dd$OS, dd$Censor) ~ lp)
    conc$concordance
  }, error = function(e) NA_real_)
}

# ---- Inside analyze_one_cohort ----
analyze_one_cohort <- function(dd, cohort_name, genes, covar_sets) {
  dd <- dd %>%
    mutate(
      Gender = as.factor(Gender),
      Radio_status = as.factor(Radio_status),
      Chemo_status = as.factor(Chemo_status),
      IDH_mutation_status = as.factor(IDH_mutation_status),
      MGMTp_methylation_status = as.factor(MGMTp_methylation_status),
      Subtype = if ("Subtype" %in% colnames(dd)) as.factor(dd$Subtype) else NA
    )
  
  # --- KM (median split + tertile split) + Univariate Cox ---
  uni_tbl <- map_dfr(genes, function(g) {
    x <- dd[[g]]
    if (all(is.na(x))) return(NULL)
    
    # Median split
    med <- median(x, na.rm = TRUE)
    grp <- factor(ifelse(x >= med, "High", "Low"), levels = c("Low","High"))
    p_lr_med <- logrank_p(dd$OS, dd$Censor, grp)
    
    # Tertile split
    p_lr_tertile <- logrank_p_tertile(dd$OS, dd$Censor, x)
    
    # Univariate Cox (continuous)
    fit_u <- fit_cox(dd, rhs_terms = g)
    hrp   <- if (!is.null(fit_u)) extract_hr_p(fit_u, g) else c(HR=NA,pval=NA,beta=NA)
    
    tibble(
      cohort = cohort_name,
      gene   = g,
      model  = "cox_univariate",
      HR     = as.numeric(hrp["HR"]),
      pval   = as.numeric(hrp["pval"]),
      beta   = as.numeric(hrp["beta"]),
      logrank_p_median = p_lr_med,
      logrank_p_tertile = p_lr_tertile   # << new column
    )
  })
  
  # --- Multivariable Cox (with covariate sets) ---
  mv_tbl <- imap_dfr(covar_sets, function(covs, setname) {
    if (setname == "none") return(NULL)
    valid_covs <- covs[covs %in% colnames(dd)]
    if (length(valid_covs) == 0) return(NULL)
    map_dfr(genes, function(g) {
      rhs <- c(g, valid_covs)
      fit <- fit_cox(dd, rhs_terms = rhs)
      hrp <- if (!is.null(fit)) extract_hr_p(fit, g) else c(HR=NA,pval=NA,beta=NA)
      tibble(
        cohort = cohort_name,
        gene   = g,
        model  = paste0("cox_adj_", setname),
        HR     = as.numeric(hrp["HR"]),
        pval   = as.numeric(hrp["pval"]),
        beta   = as.numeric(hrp["beta"]),
        logrank_p_median = NA_real_,
        logrank_p_tertile = NA_real_   # not done for adjusted
      )
    })
  })
  
  bind_rows(uni_tbl, mv_tbl) %>%
    group_by(cohort, model) %>%
    mutate(FDR = p.adjust(pval, method = "fdr")) %>%
    ungroup() %>%
    mutate(direction = case_when(
      is.na(beta) ~ NA_character_,
      beta > 0    ~ "risk_up",
      beta < 0    ~ "risk_down",
      TRUE        ~ NA_character_
    ))
}
head(tcga_surv_input)
head(cgga_surv_input)

# ---- Run all analyses ----

# Analyze original cohorts
tcga_res_orig <- analyze_one_cohort(tcga_surv_input, "TCGA_orig", genes, covar_sets)
cgga_res_orig <- analyze_one_cohort(cgga_surv_input, "CGGA_orig", genes, covar_sets)

library(dplyr)

combined_res <- tcga_res_orig %>%
  rename_with(~paste0("TCGA_orig_", .), -c(gene, model)) %>%
  full_join(
    cgga_res_orig %>% rename_with(~paste0("CGGA_orig_", .), -c(gene, model)),
    by = c("gene", "model")
  )

alpha <- 0.05
colnames(combined_res)
# --- Cox significant in both cohorts (nominal p) ---

uni_res <- combined_res %>%
  filter(model == "cox_univariate")

cox_both_nominal <- uni_res %>%
  mutate(sig_TCGA = TCGA_orig_pval < alpha,
         sig_CGGA = CGGA_orig_pval < alpha,
         same_dir = TCGA_orig_direction == CGGA_orig_direction) %>%
  filter(sig_TCGA & sig_CGGA)

n_cox_both_nominal <- nrow(cox_both_nominal)

# --- KM (median split) significant in both ---
km_both_median <- uni_res %>%
  filter(TCGA_orig_logrank_p_median < alpha,
         CGGA_orig_logrank_p_median < alpha)

# --- KM (tertile split) significant in both ---
km_both_tertile <- uni_res %>%
  filter(TCGA_orig_logrank_p_tertile < alpha,
         CGGA_orig_logrank_p_tertile < alpha)

# --- Overlaps ---
overlap_km_cox_median <- intersect(cox_both_nominal$gene, km_both_median$gene)
overlap_km_cox_tertile <- intersect(cox_both_nominal$gene, km_both_tertile$gene)

# --- Summary ---
cat("Cox (nominal p < 0.05) significant in both:", n_cox_both_nominal, "\n")
cat("KM median significant in both:", nrow(km_both_median), "\n")
cat("KM tertile significant in both:", nrow(km_both_tertile), "\n")
cat("Overlap Cox+KM median:", length(overlap_km_cox_median), "\n")
cat("Overlap Cox+KM tertile:", length(overlap_km_cox_tertile), "\n")

# ---- Build robust gene summary ----
# Genes significant in Cox + KM (median or tertile) across both cohorts
robust_genes <- intersect(overlap_km_cox_median, overlap_km_cox_tertile)

robust_genes <- Reduce(intersect, list(
  cox_both_nominal$gene,
  km_both_median$gene,
  km_both_tertile$gene
))

# Filter to univariate Cox results for these genes
robust_summary <- uni_res %>%
  dplyr::filter(gene %in% robust_genes, model == "cox_univariate") %>%
  dplyr::select(
    gene,
    TCGA_orig_HR, TCGA_orig_pval, TCGA_orig_direction,
    CGGA_orig_HR, CGGA_orig_pval, CGGA_orig_direction,
    TCGA_orig_logrank_p_median, CGGA_orig_logrank_p_median,
    TCGA_orig_logrank_p_tertile, CGGA_orig_logrank_p_tertile
  ) %>%
  # Add concordance column
  mutate(
    concordance = ifelse(TCGA_orig_direction == CGGA_orig_direction,
                         "Same", "Opposite")
  ) %>%
  # Round for readability
  mutate(
    across(ends_with("HR"), ~ round(., 2)),
    across(ends_with("pval"), ~ signif(., 3)),
    across(starts_with("TCGA_orig_logrank_p"), ~ signif(., 3)),
    across(starts_with("CGGA_orig_logrank_p"), ~ signif(., 3))
  )

# View final table
print(robust_summary)

# If you want to export
write.csv(robust_summary, "robust_survival_genes_summary.csv", row.names = FALSE)

# Base R
write.csv(robust_summary, "robust_survival_genes_summary.csv", row.names = FALSE)

# Or with readr (nicer formatting, no row names column)
readr::write_csv(robust_summary, "robust_survival_genes_summary.csv")

# Example: replace with your actual gene vectors
tcga_triple <- intersect(
  intersect(
    uni_res$gene[uni_res$TCGA_orig_pval < 0.05],
    uni_res$gene[uni_res$TCGA_orig_logrank_p_median < 0.05]
  ),
  uni_res$gene[uni_res$TCGA_orig_logrank_p_tertile < 0.05]
)

cgga_triple <- intersect(
  intersect(
    uni_res$gene[uni_res$CGGA_orig_pval < 0.05],
    uni_res$gene[uni_res$CGGA_orig_logrank_p_median < 0.05]
  ),
  uni_res$gene[uni_res$CGGA_orig_logrank_p_tertile < 0.05]
)

summary_tbl <- tibble(
  Method = c("Cox regression", "KM (median)", "KM (tertile)", "Triple overlap"),
  TCGA_sig = c(
    sum(uni_res$TCGA_orig_pval < 0.05, na.rm = TRUE),
    sum(uni_res$TCGA_orig_logrank_p_median < 0.05, na.rm = TRUE),
    sum(uni_res$TCGA_orig_logrank_p_tertile < 0.05, na.rm = TRUE),
    length(tcga_triple)
  ),
  CGGA_sig = c(
    sum(uni_res$CGGA_orig_pval < 0.05, na.rm = TRUE),
    sum(uni_res$CGGA_orig_logrank_p_median < 0.05, na.rm = TRUE),
    sum(uni_res$CGGA_orig_logrank_p_tertile < 0.05, na.rm = TRUE),
    length(cgga_triple)
  ),
  Overlap = c(
    length(cox_both_nominal$gene),
    nrow(km_both_median),
    nrow(km_both_tertile),
    length(intersect(tcga_triple, cgga_triple))
  )
) %>%
  mutate(
    Perc_TCGA = round(100 * Overlap / TCGA_sig, 1),
    Perc_CGGA = round(100 * Overlap / CGGA_sig, 1)
  )

print(summary_tbl)

# Save
write.csv(summary_tbl, "survival_overlap_summary.csv", row.names = FALSE)

# Number of same-direction concordant genes across cohorts
n_concordant <- robust_summary %>%
  filter(concordance == "Same") %>%
  nrow()

# Add as extra row
summary_tbl <- summary_tbl %>%
  add_row(
    Method = "Concordant overlap (same direction)",
    TCGA_sig = NA,
    CGGA_sig = NA,
    Overlap = n_concordant,
    Perc_TCGA = NA,
    Perc_CGGA = NA
  )

print(summary_tbl)
write.csv(summary_tbl, "survival_overlap_summary_with_concordant.csv", row.names = FALSE)



##----------------------------------------------------------------------##
##                             Figure 3 plot                            ##
##----------------------------------------------------------------------##

library(ggplot2)

summary_counts <- data.frame(
  Category = c("Cox (both cohorts)", 
               "KM median (both cohorts)", 
               "KM tertile (both cohorts)", 
               "Cox + KM median overlap", 
               "Cox + KM tertile overlap"),
  Count = c(n_cox_both_nominal, 
            nrow(km_both_median), 
            nrow(km_both_tertile), 
            length(overlap_km_cox_median), 
            length(overlap_km_cox_tertile))
)

ggplot(summary_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat="identity") +
  geom_text(aes(label = Count), vjust = -0.5, size=4) +
  theme_bw() +
  labs(title="Summary of Survival Analysis Across TCGA and CGGA",
       y="Number of Genes", x="") +
  theme(axis.text.x = element_text(angle=25, hjust=1),
        legend.position="none")

# Three-way overlap
overlap_all_three <- Reduce(intersect, 
                            list(cox_both_nominal$gene, 
                                 km_both_median$gene, 
                                 km_both_tertile$gene))
length(overlap_all_three)

summary_counts <- data.frame(
  Category = c("Cox (both cohorts)", 
               "KM median (both cohorts)", 
               "KM tertile (both cohorts)", 
               "Cox ∩ KM median (both cohorts)", 
               "Cox ∩ KM tertile (both cohorts)",
               "Cox ∩ KM median ∩ KM tertile (both cohorts)"),
  Count = c(n_cox_both_nominal, 
            nrow(km_both_median), 
            nrow(km_both_tertile), 
            length(overlap_km_cox_median), 
            length(overlap_km_cox_tertile),
            length(overlap_all_three))
)


gene_sets <- list(
  Cox_TCGA = uni_res$gene[uni_res$TCGA_orig_pval < 0.05],
  Cox_CGGA = uni_res$gene[uni_res$CGGA_orig_pval < 0.05],
  KM_median_TCGA = uni_res$gene[uni_res$TCGA_orig_logrank_p_median < 0.05],
  KM_median_CGGA = uni_res$gene[uni_res$CGGA_orig_logrank_p_median < 0.05],
  KM_tertile_TCGA = uni_res$gene[uni_res$TCGA_orig_logrank_p_tertile < 0.05],
  KM_tertile_CGGA = uni_res$gene[uni_res$CGGA_orig_logrank_p_tertile < 0.05]
)

upset(fromList(gene_sets), 
      nsets = 6, 
      order.by = "freq", 
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size")

library(ggVennDiagram)

venn_tcga <- ggVennDiagram(
  list(
    Cox = uni_res$gene[uni_res$TCGA_orig_pval < 0.05],
    KM_median = uni_res$gene[uni_res$TCGA_orig_logrank_p_median < 0.05],
    KM_tertile = uni_res$gene[uni_res$TCGA_orig_logrank_p_tertile < 0.05]
  ),
  label_alpha = 0, 
  label = "count"
) + ggtitle("TCGA: Survival Significance Overlap")

print(venn_tcga)

venn_tcga <- ggVennDiagram(
  list(
    Cox = uni_res$gene[uni_res$TCGA_orig_pval < 0.05],
    KM_median = uni_res$gene[uni_res$TCGA_orig_logrank_p_median < 0.05],
    KM_tertile = uni_res$gene[uni_res$TCGA_orig_logrank_p_tertile < 0.05]
  ),
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "steelblue") +   # lighter blue
  ggtitle("TCGA: Survival Significance Overlap") +
  theme(legend.position = "none")

print(venn_tcga)

venn_cgga <- ggVennDiagram(
  list(
    Cox = uni_res$gene[uni_res$CGGA_orig_pval < 0.05],
    KM_median = uni_res$gene[uni_res$CGGA_orig_logrank_p_median < 0.05],
    KM_tertile = uni_res$gene[uni_res$CGGA_orig_logrank_p_tertile < 0.05]
  ),
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "darkred") +   # you can pick a contrasting color
  ggtitle("CGGA: Survival Significance Overlap") +
  theme(legend.position = "none")

print(venn_cgga)

venn_cgga <- ggVennDiagram(
  list(
    Cox = uni_res$gene[uni_res$CGGA_orig_pval < 0.05],
    KM_median = uni_res$gene[uni_res$CGGA_orig_logrank_p_median < 0.05],
    KM_tertile = uni_res$gene[uni_res$CGGA_orig_logrank_p_tertile < 0.05]
  ),
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "#FFEEEE", high = "#FF9999") +  # pastel range
  ggtitle("CGGA: Survival Significance Overlap") +
  theme(legend.position = "none")

print(venn_cgga)

# Venn diagram of TCGA vs CGGA triple-overlap sets
venn_meta <- ggVennDiagram(
  list(
    TCGA_triple = tcga_triple,
    CGGA_triple = cgga_triple
  ),
  label_alpha = 0,
  label = "count"
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Overlap of Triple-Significant Genes (TCGA vs CGGA)") +
  theme(legend.position = "none")

print(venn_meta)

# You can also directly inspect overlap genes:
overlap_triple <- intersect(tcga_triple, cgga_triple)
length(overlap_triple)
head(overlap_triple)

# Find overlapping significant genes
overlap_triple <- intersect(tcga_triple, cgga_triple)

# Count them
length(overlap_triple)

# Show the first few
head(overlap_triple)

library(VennDiagram)
library(grid)

# Define overlaps
tcga_triple <- intersect(
  intersect(
    uni_res$gene[uni_res$TCGA_orig_pval < 0.05],
    uni_res$gene[uni_res$TCGA_orig_logrank_p_median < 0.05]
  ),
  uni_res$gene[uni_res$TCGA_orig_logrank_p_tertile < 0.05]
)

cgga_triple <- intersect(
  intersect(
    uni_res$gene[uni_res$CGGA_orig_pval < 0.05],
    uni_res$gene[uni_res$CGGA_orig_logrank_p_median < 0.05]
  ),
  uni_res$gene[uni_res$CGGA_orig_logrank_p_tertile < 0.05]
)

# Proper 2-set Venn
venn.plot <- draw.pairwise.venn(
  area1 = length(tcga_triple),
  area2 = length(cgga_triple),
  cross.area = length(intersect(tcga_triple, cgga_triple)),
  category = c("TCGA triple", "CGGA triple"),
  fill = c("#4C9F70", "#E07B91"),   # pastel green / pink
  alpha = 0.6,
  cat.pos = c(180, 0),              # TCGA label left, CGGA label right
  cat.dist = c(0.1, 0.1),           # push labels outwards
  cat.cex = 1.4,                    # label font size
  cex = 1.4,                        # numbers font size
  fontface = "bold"
)

grid.draw(venn.plot)


venn.plot <- draw.pairwise.venn(
  area1 = length(tcga_triple),
  area2 = length(cgga_triple),
  cross.area = length(intersect(tcga_triple, cgga_triple)),
  category = c("TCGA triple", "CGGA triple"),
  fill = c("#4C9F70", "#E07B91"),
  alpha = 0.6,
  scaled = FALSE,   # <---- force equal-size circles
  cat.pos = c(180, 0),
  cat.dist = c(0.1, 0.1),
  cat.cex = 1.4,
  cex = 1.4,
  fontface = "bold"
)
grid.newpage()
grid.draw(venn.plot)

library(VennDiagram)

venn.plot <- draw.pairwise.venn(
  area1 = length(tcga_triple),
  area2 = length(cgga_triple),
  cross.area = length(intersect(tcga_triple, cgga_triple)),
  category = c("TCGA triple", "CGGA triple"),
  fill = c("#4682B4", "#E57373"),   # Blue for TCGA, Red for CGGA
  alpha = 0.6,
  scaled = FALSE,                   # force equal-size circles
  cat.pos = c(180, 0),              # TCGA left, CGGA right
  cat.dist = c(0.1, 0.1),
  cat.cex = 1.4,                    # label size
  cex = 1.4,                        # numbers size
  fontface = "bold"
)

grid.newpage()
grid.draw(venn.plot)

# Genes of interest (6 overlap)
overlap_genes <- c("ENSG00000105996","ENSG00000105997","ENSG00000128710",
                   "ENSG00000128713","ENSG00000170178","ENSG00000228630")


robust_genes <- union(overlap_km_cox_median, overlap_km_cox_tertile)
head(overlap_km_cox_median)
head(overlap_km_cox_tertile)
print(robust_genes)


# Genes that are in both overlap_genes and robust_genes
common_genes_list <- intersect(overlap_genes, robust_genes)

# Count and view
length(common_genes_list)
print(common_genes_list)

# Define your sets
tcga_triple <- intersect(
  intersect(
    uni_res$gene[uni_res$TCGA_orig_pval < 0.05],
    uni_res$gene[uni_res$TCGA_orig_logrank_p_median < 0.05]
  ),
  uni_res$gene[uni_res$TCGA_orig_logrank_p_tertile < 0.05]
)

cgga_triple <- intersect(
  intersect(
    uni_res$gene[uni_res$CGGA_orig_pval < 0.05],
    uni_res$gene[uni_res$CGGA_orig_logrank_p_median < 0.05]
  ),
  uni_res$gene[uni_res$CGGA_orig_logrank_p_tertile < 0.05]
)

# Make the Venn diagram
venn.plot <- draw.pairwise.venn(
  area1      = length(cgga_triple),                          # CGGA on left
  area2      = length(tcga_triple),                          # TCGA on right
  cross.area = length(intersect(tcga_triple, cgga_triple)),  # overlap
  category   = c("CGGA triple", "TCGA triple"),
  fill       = c("#E57373", "#4682B4"),   # Red for CGGA, Blue for TCGA
  alpha      = 0.6,
  scaled     = FALSE,                   # equal circle sizes
  cat.pos    = c(-20, 20),              # adjust labels (left, right)
  cat.dist   = c(0.05, 0.05),           
  cat.cex    = 1.4,
  cex        = 1.4,
  fontface   = "bold"
)

grid.newpage()
grid.draw(venn.plot)

colnames(robust_summary)
cox_both <- cox_both_nominal %>%
  mutate(concordance = ifelse(same_dir, "Same", "Opposite"))


same_dir <- cox_both %>% filter(concordance == "Same") %>% pull(gene) 
opp_dir <- cox_both %>% filter(concordance == "Opposite") %>% pull(gene)

venn.plot <- draw.pairwise.venn(
  area1      = length(same_dir) + length(opp_dir),  # TCGA risk_up (all significant)
  area2      = length(same_dir) + length(opp_dir),  # CGGA risk (all significant)
  cross.area = length(same_dir),                   # overlap = same direction
  category   = c("TCGA Cox significant", "CGGA Cox significant"),
  fill       = c("#4682B4", "#E57373"),  # Blue for TCGA, Red for CGGA
  alpha      = 0.6,
  scaled     = FALSE,                    # force equal-size circles
  cat.pos    = c(180, 0),
  cat.dist   = c(0.05, 0.05),
  cat.cex    = 1.4,
  cex        = 1.4,
  fontface   = "bold"
)

grid.newpage()
grid.draw(venn.plot)
head(robust_summary)

# Filter robust_summary for same-direction genes only
robust_summary_same <- robust_summary %>%
  dplyr::filter(concordance == "Same")

# Check results
print(robust_summary_same)

# Save table
write.csv(robust_summary_same, "robust_survival_genes_same_direction.csv", row.names = FALSE)
head(robust_summary_same)

##----------------------------------------------------------------------------##
##                                SANITY CHECK.                               ##
##----------------------------------------------------------------------------##
robust_summary %>% 
  dplyr::filter(gene %in% overlap_genes) %>%
  dplyr::select(gene, TCGA_orig_HR, TCGA_orig_pval, CGGA_orig_HR, CGGA_orig_pval, concordance)

robust_summary %>% 
  dplyr::filter(gene %in% overlap_genes) %>%
  dplyr::select(gene,
                TCGA_orig_logrank_p_median, CGGA_orig_logrank_p_median,
                TCGA_orig_logrank_p_tertile, CGGA_orig_logrank_p_tertile)
