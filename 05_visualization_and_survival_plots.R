# ============================================================
# 05_visualization_and_survival_plots.R
# Expression Visualization and Kaplan–Meier Validation
# ============================================================
##----------------------------------------------------------------##
##                       Violin plots                             ##
##----------------------------------------------------------------##
# Build long expression dataframe for one gene at a time
make_expr_df <- function(gene_id) {
  # TCGA tumor
  tcga_tum <- data.frame(
    gene = gene_id,
    expr = tcga_cpm[gene_id, tumor_samples],
    group = "TCGA_tumor"
  ) %>% rownames_to_column("sample")
  
  # TCGA normal
  tcga_norm <- data.frame(
    gene = gene_id,
    expr = tcga_cpm[gene_id, normal_samples],
    group = "TCGA_normal"
  ) %>% rownames_to_column("sample")
  
  # GTEx normal
  gtex_norm <- data.frame(
    gene = gene_id,
    expr = gtex_cpm[gene_id, ],
    group = "GTEx_normal"
  ) %>% rownames_to_column("sample")
  
  # CGGA tumor
  cgga_tum <- data.frame(
    gene = gene_id,
    expr = cgga_cpm[gene_id, ],
    group = "CGGA_tumor"
  ) %>% rownames_to_column("sample")
  
  # Combine all
  bind_rows(tcga_tum, cgga_tum, tcga_norm, gtex_norm) %>%
    mutate(group = factor(group, 
                          levels = c("TCGA_tumor","CGGA_tumor",
                                     "TCGA_normal","GTEx_normal")))
}

# Output directory
dir.create("expression_plots", showWarnings = FALSE)

# Loop through genes and plot
for (g in overlap_genes) {
  df <- make_expr_df(g)
  
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
    labs(
      title = paste0("Expression of ", g),
      y = "Expression (log2 CPM)", x = ""
    ) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    scale_fill_manual(values = c(
      "TCGA_tumor"="#D55E00",
      "CGGA_tumor"="#CC79A7",
      "TCGA_normal"="#009E73",
      "GTEx_normal"="#56B4E9"
    ))
  
  ggsave(filename = file.path("expression_plots", paste0("expr_", g, ".pdf")),
         plot = p, width = 5, height = 5)
}

##--------------------------------------------------------------------##
##                                    Bar plot.                       ##
##--------------------------------------------------------------------##
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)


# Build long expression dataframe for one gene at a time
make_expr_df <- function(gene_id) {
  tcga_tum <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, tumor_samples], group = "TCGA_tumor")
  tcga_norm <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, normal_samples], group = "TCGA_normal")
  gtex_norm <- data.frame(gene = gene_id, expr = gtex_cpm[gene_id, ], group = "GTEx_normal")
  cgga_tum <- data.frame(gene = gene_id, expr = cgga_cpm[gene_id, ], group = "CGGA_tumor")
  
  bind_rows(tcga_tum, cgga_tum, tcga_norm, gtex_norm) %>%
    mutate(group = factor(group, 
                          levels = c("TCGA_tumor","CGGA_tumor",
                                     "TCGA_normal","GTEx_normal")))
}

# Output dir
dir.create("expression_barplots", showWarnings = FALSE)

# Loop
for (g in overlap_genes) {
  df <- make_expr_df(g)
  
  # Run tests
  anova_res <- summary(aov(expr ~ group, data = df))
  kruskal_res <- kruskal.test(expr ~ group, data = df)
  
  cat("Gene:", g, "\n")
  print(anova_res)
  print(kruskal_res)
  cat("\n-------------------------\n")
  
  # Bar plot with mean + SEM
  plot_df <- df %>%
    group_by(group) %>%
    summarise(mean_expr = mean(expr, na.rm = TRUE),
              se_expr = sd(expr, na.rm = TRUE)/sqrt(n()),
              .groups = "drop")
  
  p <- ggplot(plot_df, aes(x = group, y = mean_expr, fill = group)) +
    geom_col(position = position_dodge(), width = 0.6) +
    geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr), 
                  width = 0.2, position = position_dodge(0.6)) +
    labs(title = paste0("Expression of ", g),
         y = "Mean expression (log2 CPM)", x = "") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    scale_fill_manual(values = c(
      "TCGA_tumor"="#D55E00",
      "CGGA_tumor"="#CC79A7",
      "TCGA_normal"="#009E73",
      "GTEx_normal"="#56B4E9"
    ))
  
  ggsave(filename = file.path("expression_barplots", paste0("bar_", g, ".pdf")),
         plot = p, width = 5, height = 5)
}

##--------------------------------------------------------------------##
##                                    Box plot.                       ##
##--------------------------------------------------------------------##

library(ggplot2)
library(dplyr)

# Build long expression dataframe
make_expr_df <- function(gene_id) {
  tcga_tum <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, tumor_samples], group = "TCGA_tumor")
  tcga_norm <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, normal_samples], group = "TCGA_normal")
  gtex_norm <- data.frame(gene = gene_id, expr = gtex_cpm[gene_id, ], group = "GTEx_normal")
  cgga_tum <- data.frame(gene = gene_id, expr = cgga_cpm[gene_id, ], group = "CGGA_tumor")
  
  bind_rows(tcga_tum, cgga_tum, tcga_norm, gtex_norm) %>%
    mutate(group = factor(group, 
                          levels = c("TCGA_tumor","CGGA_tumor",
                                     "TCGA_normal","GTEx_normal")))
}

# Output dir
dir.create("expression_boxplots", showWarnings = FALSE)

# Loop
for (g in overlap_genes) {
  df <- make_expr_df(g)
  
  # Run tests
  anova_res <- summary(aov(expr ~ group, data = df))
  kruskal_res <- kruskal.test(expr ~ group, data = df)
  
  cat("Gene:", g, "\n")
  print(anova_res)
  print(kruskal_res)
  cat("\n-------------------------\n")
  
  # Boxplot with points
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +   # clean boxplot
    labs(title = paste0("Expression of ", g),
         y = "Expression (log2 CPM)", x = "") +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none") +
    scale_fill_manual(values = c(
      "TCGA_tumor"="#D55E00",
      "CGGA_tumor"="#CC79A7",
      "TCGA_normal"="#009E73",
      "GTEx_normal"="#56B4E9"
    ))
  
  ggsave(filename = file.path("expression_boxplots", paste0("box_", g, ".pdf")),
         plot = p, width = 5, height = 5)
}

library(ggplot2)
library(dplyr)



# Build expression dataframe for one gene at a time
make_expr_df <- function(gene_id) {
  tcga_tum <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, tumor_samples], group = "TCGA_tumor")
  tcga_norm <- data.frame(gene = gene_id, expr = tcga_cpm[gene_id, normal_samples], group = "TCGA_normal")
  gtex_norm <- data.frame(gene = gene_id, expr = gtex_cpm[gene_id, ], group = "GTEx_normal")
  cgga_tum <- data.frame(gene = gene_id, expr = cgga_cpm[gene_id, ], group = "CGGA_tumor")
  
  bind_rows(tcga_tum, cgga_tum, tcga_norm, gtex_norm) %>%
    mutate(group = factor(group, 
                          levels = c("TCGA_tumor","CGGA_tumor",
                                     "TCGA_normal","GTEx_normal")))
}

# Output directory
dir.create("expression_plots", showWarnings = FALSE)

# Loop through genes and plot boxplots
for (g in overlap_genes) {
  df <- make_expr_df(g)
  
  p <- ggplot(df, aes(x = group, y = expr, fill = group)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5, width = 0.6) +
    labs(
      title = paste0("Expression of ", g),
      y = "Expression (log2 CPM)", x = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none"
    ) +
    scale_fill_manual(values = c(
      "TCGA_tumor"="#D55E00",
      "CGGA_tumor"="#CC79A7",
      "TCGA_normal"="#009E73",
      "GTEx_normal"="#56B4E9"
    ))
  
  ggsave(
    filename = file.path("survival_tertile_plots", paste0("KM_", cohort_name, "_", gene_id, "_tertiles.pdf")),
    plot = print(p1),
    width = 6, height = 6
  )
}


overlap_genes %in% colnames(tcga_surv_input)
table(tcga_surv_input$Censor)
tcga_surv_input$Censor <- ifelse(tcga_surv_input$Censor == 1, 0, 1)



library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(gridExtra)

plot_survival_tertile <- function(df, gene_id, cohort_name) {
  if (!gene_id %in% colnames(df)) {
    warning(paste("Gene", gene_id, "not found in", cohort_name))
    return(NULL)
  }
  
  if (!dir.exists("survival_tertile_plots")) dir.create("survival_tertile_plots")
  
  df <- df %>%
    mutate(OS = as.numeric(OS),
           Censor = as.numeric(Censor))
  
  # Tertiles
  qs <- quantile(df[[gene_id]], probs = c(1/3, 2/3), na.rm = TRUE)
  df$expr_group <- cut(df[[gene_id]],
                       breaks = c(-Inf, qs[1], qs[2], Inf),
                       labels = c("Low", "Mid", "High"))
  
  # 3-group survival
  if (length(unique(na.omit(df$expr_group))) >= 2) {
    fit <- survfit(Surv(OS, Censor) ~ expr_group, data = df)
    p1 <- ggsurvplot(
      fit, data = df,
      pval = TRUE, risk.table = TRUE,
      palette = c("#1b9e77", "#7570b3", "#d95f02"),
      legend.title = "Expression",
      title = paste0(cohort_name, " survival (", gene_id, ", tertiles)"),
      ggtheme = theme_bw(base_size = 14)
    )
    ggsave(file.path("survival_tertile_plots", paste0("KM_", cohort_name, "_", gene_id, "_tertiles.pdf")),
           plot = p1$plot, width = 6, height = 6)
  }
  
  # High vs Low survival
  df2 <- df %>% filter(!is.na(expr_group) & expr_group %in% c("Low","High"))
  if (length(unique(df2$expr_group)) == 2) {
    fit2 <- survfit(Surv(OS, Censor) ~ expr_group, data = df2)
    p2 <- ggsurvplot(
      fit2, data = df2,
      pval = TRUE, risk.table = TRUE,
      palette = c("#1b9e77", "#d95f02"),
      legend.title = "Expression",
      title = paste0(cohort_name, " survival (", gene_id, ", High vs Low)"),
      ggtheme = theme_bw(base_size = 14)
    )
    ggsave(file.path("survival_tertile_plots", paste0("KM_", cohort_name, "_", gene_id, "_HighLow.pdf")),
           plot = p2$plot, width = 6, height = 6)
  }
}

print(robust_summary)

# Genes of interest
overlap_genes <- robust_summary$gene


# Cohorts to run
cohorts <- list(
  TCGA = tcga_surv_input,
  CGGA = cgga_surv_input
)

# Loop over genes and cohorts
for (gene in overlap_genes) {
  for (cohort_name in names(cohorts)) {
    message("Plotting ", gene, " in ", cohort_name)
    plot_survival_tertile(cohorts[[cohort_name]], gene, cohort_name)
  }
}


plot_survival_tertile <- function(df, gene_id, cohort_name) {
  if (!gene_id %in% colnames(df)) {
    warning(paste("Gene", gene_id, "not found in", cohort_name))
    return(NULL)
  }
  
  if (!dir.exists("survival_tertile_plots")) dir.create("survival_tertile_plots")
  
  df <- df %>%
    mutate(
      OS = as.numeric(OS),
      Censor = as.numeric(Censor)
    )
  
  # Tertiles
  qs <- quantile(df[[gene_id]], probs = c(1/3, 2/3), na.rm = TRUE)
  df$expr_group <- cut(
    df[[gene_id]],
    breaks = c(-Inf, qs[1], qs[2], Inf),
    labels = c("Low", "Mid", "High")
  )
  
  results <- list()
  
  # 3-group log-rank
  if (length(unique(na.omit(df$expr_group))) >= 2) {
    fit <- survfit(Surv(OS, Censor) ~ expr_group, data = df)
    sd <- survdiff(Surv(OS, Censor) ~ expr_group, data = df)
    pval <- pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
    
    results$tertile_p <- pval
    
    p1 <- ggsurvplot(
      fit, data = df,
      pval = FALSE, risk.table = TRUE,  # suppress p-value
      palette = c("#1b9e77", "#7570b3", "#d95f02"),
      legend.title = "Expression",
      title = paste0(cohort_name, " survival (", gene_id, ", tertiles)"),
      ggtheme = theme_bw(base_size = 14)
    )
    ggsave(
      filename = file.path("survival_tertile_plots", paste0("KM_", cohort_name, "_", gene_id, "_tertiles.pdf")),
      plot = p1, width = 6, height = 6
    )
  }
  
  # High vs Low
  df2 <- df %>% filter(expr_group %in% c("Low","High"))
  if (length(unique(na.omit(df2$expr_group))) == 2) {
    fit2 <- survfit(Surv(OS, Censor) ~ expr_group, data = df2)
    sd2 <- survdiff(Surv(OS, Censor) ~ expr_group, data = df2)
    pval2 <- pchisq(sd2$chisq, df = length(sd2$n) - 1, lower.tail = FALSE)
    
    results$highlow_p <- pval2
    
    p2 <- ggsurvplot(
      fit2, data = df2,
      pval = FALSE, risk.table = TRUE,  # suppress p-value
      palette = c("#1b9e77", "#d95f02"),
      legend.title = "Expression",
      title = paste0(cohort_name, " survival (", gene_id, ", High vs Low)"),
      ggtheme = theme_bw(base_size = 14)
    )
    ggsave(
      filename = file.path("survival_tertile_plots", paste0("KM_", cohort_name, "_", gene_id, "_HighLow.pdf")),
      plot = p2, width = 6, height = 6
    )
  }
  
  return(results)
}

