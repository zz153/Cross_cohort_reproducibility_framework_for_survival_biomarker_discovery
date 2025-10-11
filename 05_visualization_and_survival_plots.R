# ============================================================
# 05_visualization_and_survival_plots.R
# Expression Visualization and Kaplan–Meier Validation
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tibble)
  library(tidyr); library(survival); library(survminer)
  library(gridExtra)
})

cat("\n📊 Generating expression and survival plots...\n")

# ============================================================
# Helper: Build expression dataframe for a gene
# ============================================================
make_expr_df <- function(gene_id) {
  tcga_tum <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, tumor_samples], group = "TCGA_tumor")
  tcga_norm <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, normal_samples], group = "TCGA_normal")
  gtex_norm <- data.frame(gene = gene_id, expr = gtex_cpm[gene_id, ], group = "GTEx_normal")
  cgga_tum  <- data.frame(gene = gene_id, expr = cgga_cpm[gene_id, ], group = "CGGA_tumor")
  
  bind_rows(tcga_tum, cgga_tum, tcga_norm, gtex_norm) %>%
    mutate(group = factor(group,
                          levels = c("TCGA_tumor","CGGA_tumor",
                                     "TCGA_normal","GTEx_normal")))
}

# Common theme and colors
plot_colors <- c(
  "TCGA_tumor"="#D55E00",
  "CGGA_tumor"="#CC79A7",
  "TCGA_normal"="#009E73",
  "GTEx_normal"="#56B4E9"
)
plot_theme <- theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

overlap_genes <- robust_summary$gene

# ============================================================
# 1️⃣ Violin plots
# ============================================================
dir.create("expression_violin", showWarnings = FALSE)

for (g in overlap_genes) {
  df <- make_expr_df(g)
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    labs(title = paste0("Expression of ", g),
         y = "Expression (log2 CPM)", x = "") +
    plot_theme + scale_fill_manual(values = plot_colors)
  
  ggsave(file.path("expression_violin", paste0("violin_", g, ".pdf")), p, width = 5, height = 5)
}

# ============================================================
# 2️⃣ Bar plots (mean ± SEM)
# ============================================================
dir.create("expression_barplots", showWarnings = FALSE)

for (g in overlap_genes) {
  df <- make_expr_df(g)
  plot_df <- df %>%
    group_by(group) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              se_expr = sd(expr, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  
  p <- ggplot(plot_df, aes(x = group, y = mean_expr, fill = group)) +
    geom_col(width = 0.6) +
    geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                  width = 0.2) +
    labs(title = paste0("Expression of ", g),
         y = "Mean expression (log2 CPM)", x = "") +
    plot_theme + scale_fill_manual(values = plot_colors)
  
  ggsave(file.path("expression_barplots", paste0("bar_", g, ".pdf")), p, width = 5, height = 5)
}

# ============================================================
# 3️⃣ Box plots
# ============================================================
dir.create("expression_boxplots", showWarnings = FALSE)

for (g in overlap_genes) {
  df <- make_expr_df(g)
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
    labs(title = paste0("Expression of ", g),
         y = "Expression (log2 CPM)", x = "") +
    plot_theme + scale_fill_manual(values = plot_colors)
  
  ggsave(file.path("expression_boxplots", paste0("box_", g, ".pdf")), p, width = 5, height = 5)
}

# ============================================================
# 4️⃣ Kaplan–Meier Survival (tertiles + high/low)
# ============================================================
dir.create("survival_tertile_plots", showWarnings = FALSE)

plot_survival_tertile <- function(df, gene_id, cohort_name) {
  if (!gene_id %in% colnames(df)) {
    warning(paste("Gene", gene_id, "not found in", cohort_name))
    return(NULL)
  }
  
  df <- df %>%
    mutate(OS = as.numeric(OS), Censor = as.numeric(Censor))
  
  qs <- quantile(df[[gene_id]], probs = c(1/3, 2/3), na.rm = TRUE)
  df$expr_group <- cut(df[[gene_id]], breaks = c(-Inf, qs[1], qs[2], Inf),
                       labels = c("Low","Mid","High"))
  
  # 3-group survival
  if (length(unique(na.omit(df$expr_group))) >= 2) {
    fit <- survfit(Surv(OS, Censor) ~ expr_group, data = df)
    p1 <- ggsurvplot(fit, data = df, pval = TRUE, risk.table = TRUE,
                     palette = c("#1b9e77", "#7570b3", "#d95f02"),
                     legend.title = "Expression",
                     title = paste0(cohort_name, " survival (", gene_id, ", tertiles)"),
                     ggtheme = theme_bw(base_size = 14))
    ggsave(file.path("survival_tertile_plots",
                     paste0("KM_", cohort_name, "_", gene_id, "_tertiles.pdf")),
           p1$plot, width = 6, height = 6)
  }
  
  # High vs Low
  df2 <- df %>% filter(expr_group %in% c("Low","High"))
  if (length(unique(df2$expr_group)) == 2) {
    fit2 <- survfit(Surv(OS, Censor) ~ expr_group, data = df2)
    p2 <- ggsurvplot(fit2, data = df2, pval = TRUE, risk.table = TRUE,
                     palette = c("#1b9e77", "#d95f02"),
                     legend.title = "Expression",
                     title = paste0(cohort_name, " survival (", gene_id, ", High vs Low)"),
                     ggtheme = theme_bw(base_size = 14))
    ggsave(file.path("survival_tertile_plots",
                     paste0("KM_", cohort_name, "_", gene_id, "_HighLow.pdf")),
           p2$plot, width = 6, height = 6)
  }
}

# ============================================================
# Run survival plots for all genes and cohorts
# ============================================================
overlap_genes <- robust_summary$gene
cohorts <- list(TCGA = tcga_surv_input, CGGA = cgga_surv_input)

for (gene in overlap_genes) {
  for (cohort_name in names(cohorts)) {
    message("Plotting ", gene, " in ", cohort_name)
    plot_survival_tertile(cohorts[[cohort_name]], gene, cohort_name)
  }
}

cat("\n✅ All expression and survival plots generated successfully.\n")

