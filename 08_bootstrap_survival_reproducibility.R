tcga_surv_input <- read_csv("TCGA_Survival_Input.csv")
cgga_surv_input <- read_csv("CGGA_Survival_Input.csv")

genes_tp <- c("ENSG00000105997","ENSG00000128710","ENSG00000128713",
              "ENSG00000253552","ENSG00000250133")      # robust meta-validated
genes_fp <- false_pos_genes$gene                        # internal false positives
candidate_genes <- c(genes_tp, genes_fp)

# ------------------------------------------------------------
# 3️⃣ Bootstrapped Cox regression
# ------------------------------------------------------------
bootstrap_cox <- function(df, gene, B = 1000) {
  n <- nrow(df)
  res <- replicate(B, {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    fit <- try(coxph(Surv(OS, Censor) ~ df[idx, gene], data = df), silent = TRUE)
    if (inherits(fit, "try-error") || length(coef(fit)) == 0)
      return(c(HR = NA, p = NA))
    s <- summary(fit)
    c(HR = unname(s$coef[1, "exp(coef)"]),
      p = unname(s$coef[1, "Pr(>|z|)"]))
  })
  res <- as.data.frame(t(res))
  list(
    meanHR   = mean(res$HR, na.rm = TRUE),
    medianHR = median(res$HR, na.rm = TRUE),
    CI       = quantile(res$HR, probs = c(0.025, 0.975), na.rm = TRUE),
    sig_rate = mean(res$p < 0.05, na.rm = TRUE)
  )
}

boot_summary <- function(df, cohort, genes) {
  map_dfr(genes, function(g) {
    r <- bootstrap_cox(df, g, B = 500)  # increase to 1000 for final analysis
    tibble(Gene = g, Cohort = cohort,
           MeanHR = r$meanHR, MedianHR = r$medianHR,
           CI_lo = r$CI[1], CI_hi = r$CI[2], SigRate = r$sig_rate)
  })
}

cat("Running bootstraps (500x per gene)...\n")
tcga_boot <- boot_summary(tcga_surv_input, "TCGA", candidate_genes)
cgga_boot <- boot_summary(cgga_surv_input, "CGGA", candidate_genes)
boot_all <- bind_rows(tcga_boot, cgga_boot)
write_csv(boot_all, "Bootstrap_Cox_Summary.csv")

# ------------------------------------------------------------
# 4️⃣ Visualize bootstrap reproducibility
# ------------------------------------------------------------
boot_all <- boot_all %>%
  mutate(Type = case_when(
    Gene %in% genes_tp ~ "True Positive",
    Gene %in% genes_fp ~ "False Positive",
    TRUE ~ "Other"
  ))

boot_all$Gene <- factor(boot_all$Gene, levels = c(genes_tp, genes_fp))

p <- ggplot(boot_all, aes(x = Gene, y = SigRate, fill = Cohort)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Type, scales = "free_y", ncol = 1, strip.position = "right") +
  scale_fill_manual(values = c("TCGA" = "#1f78b4", "CGGA" = "#e31a1c")) +
  coord_flip() +
  theme_bw(base_size = 13) +
  labs(
    title = "Bootstrap reproducibility of prognostic vs false-positive genes",
    x = "Gene", y = "% significant bootstraps"
  ) +
  theme(strip.text.y.right = element_text(angle = 0))

ggsave("Bootstrap_SignificanceRates_TPvsFP.png", p, width = 8, height = 7)

cat("✅ Bootstrapping complete. Results and plots saved.\n")

# ------------------------------------------------------------
# 5️⃣ Annotate Ensembl → HGNC
# ------------------------------------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"),
                 filters = "ensembl_gene_id",
                 values = genes_fp, mart = mart)
print(mapping)
write_csv(mapping, "FalsePositive_HGNC_mapping.csv")

cat("✅ Mapping complete. Pipeline finished successfully.\n")
