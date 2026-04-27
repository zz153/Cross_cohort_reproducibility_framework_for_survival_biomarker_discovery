# ============================================================
# 03_differential_expression.R
# Differential Gene Expression Analysis (Tumor vs. Normal)
# ============================================================

# ------------------------------
# Load inputs
# ------------------------------
tcga_counts <- readRDS("TCGA_GBM_harmonized_counts.rds")
gtex_counts <- readRDS("GTEx_Brain_harmonized_counts.rds")
gene_info   <- readRDS("gene_info.rds")
tcga_meta   <- colData(readRDS("TCGA_GBM_data.rds"))

tumor_samples  <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])
tcga_normals   <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])

# ------------------------------
# Differential Expression
# ------------------------------
tcga_tumor_expr  <- tcga_counts[, colnames(tcga_counts) %in% tumor_samples]
tcga_normal_expr <- tcga_counts[, colnames(tcga_counts) %in% tcga_normals]
gtex_expr        <- gtex_counts

all_expr <- cbind(tcga_tumor_expr, tcga_normal_expr, gtex_expr)

# ------------------------------
# Pool TCGA + GTEx normals
# ------------------------------
group <- factor(c(
  rep("Tumor", ncol(tcga_tumor_expr)),
  rep("Normal", ncol(tcga_normal_expr) + ncol(gtex_expr))
), levels = c("Normal", "Tumor"))

# You can still keep batch to adjust for TCGA vs GTEx effects
batch <- factor(c(
  rep("TCGA", ncol(tcga_tumor_expr)+ncol(tcga_normal_expr)),
  rep("GTEx", ncol(gtex_expr))
))

# ------------------------------
# edgeR/voom pipeline
# ------------------------------
dge <- DGEList(all_expr)
keep <- filterByExpr(dge, group)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ group + batch)   # Tumor vs Normal, adjusting for batch
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# ------------------------------
# Get Tumor vs All_Normals results
# ------------------------------
res_tcga_pooled <- topTable(fit, coef="groupTumor", number=Inf, sort.by="P") %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info[,c("ensembl_gene_id","external_gene_name")], by="ensembl_gene_id") %>%
  mutate(
    status = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

write.csv(res_tcga_pooled, "TCGA_DEG_all_Tumor_vs_AllNormals.csv", row.names=FALSE)

dim(res_tcga_pooled)
head(res_tcga_pooled)

# ------------------------------
# Save voom object
# ------------------------------
saveRDS(v, "TCGA_voom_Tumor_vs_AllNormals.rds")

# To reload later:
v_tcga_pooled <- readRDS("TCGA_voom_Tumor_vs_AllNormals.rds")

# The log2CPM matrix is stored here:
expr_tcga_pooled <- v_tcga_pooled$E
dim(expr_tcga_pooled)
head(expr_tcga_pooled[,1:5])

dim(expr_tcga_pooled)
##---------------------------------------------------------------------------##
##                    CLINICAL DATA EXTRACTION                               ##
##---------------------------------------------------------------------------##
# Clinical metadata
tcga_meta <- as.data.frame(colData(tcga_data))

# Restrict to tumor samples only
tumor_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])

# Subset to tumor clinical info
clin_tcga <- tcga_meta[tumor_samples, ]
# Pull out the columns, deceased and time to survival

head(clin_tcga)

# Make sure CGGA_ID is character, not factor
cgga_gbm_clin$CGGA_ID <- as.character(cgga_gbm_clin$CGGA_ID)

# Set row names
rownames(cgga_gbm_clin) <- cgga_gbm_clin$CGGA_ID

# Optionally drop the CGGA_ID column (since now it's in rownames)
cgga_gbm_clin$CGGA_ID <- NULL

# ============================================================
# Script: fig2_differential_expression_CGGA.R
# Purpose: Run DE (CGGA Tumor vs Normal: TCGA+GTEx)
# ============================================================

# ------------------------------
# Load harmonized counts
# ------------------------------
cgga_counts <- readRDS("CGGA_GBM_harmonized_counts.rds")
tcga_counts <- readRDS("TCGA_GBM_harmonized_counts.rds")
gtex_counts <- readRDS("GTEx_Brain_harmonized_counts.rds")
gene_info   <- readRDS("gene_info.rds")
tcga_data   <- readRDS("TCGA_GBM_data.rds")
tcga_meta   <- as.data.frame(colData(tcga_data))

# Identify TCGA normal samples
tcga_normal_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])

# ------------------------------
# Build combined matrix: CGGA vs Normals
# ------------------------------
cgga_expr  <- cgga_counts
tcga_norms <- tcga_counts[, colnames(tcga_counts) %in% tcga_normal_samples]
gtex_expr  <- gtex_counts

all_GBM_dge_expr_cgga <- cbind(cgga_expr, tcga_norms, gtex_expr)

group <- factor(c(
  rep("CGGA_Tumor", ncol(cgga_expr)),
  rep("TCGA_Normal", ncol(tcga_norms)),
  rep("GTEx_Normal", ncol(gtex_expr))
), levels = c("TCGA_Normal","GTEx_Normal","CGGA_Tumor"))

batch <- factor(c(
  rep("CGGA", ncol(cgga_expr)),
  rep("TCGA", ncol(tcga_norms)),
  rep("GTEx", ncol(gtex_expr))
))

# ------------------------------
# edgeR + limma pipeline
# ------------------------------
dge_CGGA <- DGEList(all_GBM_dge_expr_cgga)
keep <- filterByExpr(dge_CGGA, group)
dge_CGGA <- dge_CGGA[keep,, keep.lib.sizes=FALSE]
dge_CGGA <- calcNormFactors(dge_CGGA)

design <- model.matrix(~ group + batch)
v_cgga <- voom(dge_CGGA, design, plot=TRUE)
saveRDS(v_cgga, "CGGA_voom.rds")
expr_cgga <- v_cgga$E   # CGGA log2CPM matrix
head(v_cgga)
cgga_data<-readRDS("CGGA_voom.rds")
fit_cgga <- lmFit(v_cgga, design)
fit_cgga <- eBayes(fit_cgga)

# Extract DE results: CGGA_Tumor vs Normals
res_cgga <- topTable(fit_cgga, coef="groupCGGA_Tumor", number=Inf, sort.by="P") %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info[,c("ensembl_gene_id","external_gene_name")], 
            by="ensembl_gene_id") %>%
  mutate(
    status = case_when(
      adj.P.Val < 0.05 & logFC > 1  ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Save results
write.csv(res_cgga, "DEG_CGGA_Tumor_vs_Normals.csv", row.names=FALSE)

cat("✅ CGGA vs Normal DGE complete! Results saved.\n")

res_cgga_unique <- res_cgga %>%
  group_by(ensembl_gene_id) %>%
  arrange(adj.P.Val) %>%
  slice(1) %>%  # keep the most significant entry
  ungroup()

dim(res_cgga_unique)
head(res_cgga_unique)

deg_tcga <- res_tcga_pooled %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)   # significant DEGs
dim(deg_tcga)
deg_cgga <- res_cgga_unique %>% filter(adj.P.Val < 0.05 & abs(logFC) > 1)
head(deg_cgga)

# Save significant DEGs for CGGA
write.csv(deg_cgga, 
          "CGGA_DEG_sig_Tumor_vs_Normals_unique.csv", 
          row.names = FALSE)

cat("✅ Significant DEG tables saved for TCGA and CGGA!\n")

# Extract Ensembl IDs by direction
up_tcga   <- deg_tcga %>% filter(logFC > 1)  %>% pull(ensembl_gene_id)
down_tcga <- deg_tcga %>% filter(logFC < -1) %>% pull(ensembl_gene_id)

up_cgga   <- deg_cgga %>% filter(logFC > 1)  %>% pull(ensembl_gene_id)
down_cgga <- deg_cgga %>% filter(logFC < -1) %>% pull(ensembl_gene_id)

# Overlaps
common_up   <- intersect(up_tcga, up_cgga)
common_down <- intersect(down_tcga, down_cgga)

# Counts
cat("Upregulated genes: TCGA =", length(up_tcga), 
    "CGGA =", length(up_cgga), 
    "Overlap =", length(common_up), "\n")

cat("Downregulated genes: TCGA =", length(down_tcga), 
    "CGGA =", length(down_cgga), 
    "Overlap =", length(common_down), "\n")

# Save overlap lists
write.csv(data.frame(ensembl_gene_id = common_up), 
          "Overlap_Up_TCGA_CGGA.csv", row.names=FALSE)

write.csv(data.frame(ensembl_gene_id = common_down), 
          "Overlap_Down_TCGA_CGGA.csv", row.names=FALSE)

cat("✅ Overlap lists saved: Overlap_Up_TCGA_CGGA.csv & Overlap_Down_TCGA_CGGA.csv\n")

common_up_annot <- data.frame(ensembl_gene_id = common_up) %>%
  left_join(gene_info[,c("ensembl_gene_id","external_gene_name")], 
            by="ensembl_gene_id")
write.csv(common_up_annot, "Overlap_Up_TCGA_CGGA_annotated.csv", row.names=FALSE)

common_down_annot <- data.frame(ensembl_gene_id = common_down) %>%
  left_join(gene_info[,c("ensembl_gene_id","external_gene_name")], 
            by="ensembl_gene_id")
write.csv(common_down_annot, "Overlap_Down_TCGA_CGGA_annotated.csv", row.names=FALSE)
head(common_up_annot)
head(common_down_annot)
##--------------------------------------------------------
##                     SURVIVAL ANALYSIS                  
##--------------------------------------------------------
head(cgga_gbm_clin)
head(clin_tcga)
dim(clin_tcga)

# Build a simplified survival table from clin_tcga

surv_tcga <- clin_tcga %>%
  as_tibble() %>%
  mutate(
    OS = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    Censor = ifelse(vital_status == "Dead", 1, 0),
    Age = round(age_at_diagnosis / 365.25),   # convert to years
    Gender = ifelse(gender == "male", "Male", "Female"),
    IDH_mutation_status = ifelse(paper_IDH.status == "Mutant", "Mutant", "Wildtype"),
    MGMTp_methylation_status = case_when(
      paper_MGMT.promoter.status == "Methylated" ~ "methylated",
      paper_MGMT.promoter.status == "Unmethylated" ~ "un-methylated",
      TRUE ~ NA_character_
    ),
    Radio_status = ifelse(grepl("Radiation", treatments, ignore.case=TRUE), 1, 0),
    Chemo_status = ifelse(grepl("Temozolomide", treatments, ignore.case=TRUE), 1, 0)
  ) %>%
  dplyr::select(patient, Gender, Age, OS, Censor,
                Radio_status, Chemo_status,
                IDH_mutation_status, MGMTp_methylation_status, 
                paper_Original.Subtype) %>%
  rename(Subtype = paper_Original.Subtype)
dim(surv_tcga)

# Check the cleaned table
head(surv_tcga)
dim(surv_tcga)
colnames(cgga_gbm_clin)

surv_cgga <- cgga_gbm_clin %>%
  dplyr::rename(
    Censor = `Censor (alive=0; dead=1)`,
    Radio_status = `Radio_status (treated=1;un-treated=0)`,
    Chemo_status = `Chemo_status (TMZ treated=1;un-treated=0)`
  ) %>%
  dplyr::select(Gender, Age, OS, Censor, 
                Radio_status, Chemo_status,
                IDH_mutation_status, MGMTp_methylation_status,
                `1p19q_codeletion_status`) %>%
  tibble::rownames_to_column("patient")
dim(surv_cgga)

write.csv(as.data.frame(expr_tcga_pooled), "TCGA_voom_log2CPM.csv")
write.csv(as.data.frame(expr_cgga), "CGGA_voom_log2CPM.csv")

common_genes <- unique(c(common_up_annot$ensembl_gene_id,
                         common_down_annot$ensembl_gene_id))

expr_tcga_common <- expr_tcga_pooled[rownames(expr_tcga_pooled) %in% common_genes, ]
dim(expr_tcga_common)
head(expr_tcga_common)
# Keep only columns starting with "CGGA_"
expr_cgga_only <- expr_cgga[, grepl("^CGGA_", colnames(expr_cgga))]

expr_cgga_common <- expr_cgga_only[rownames(expr_cgga) %in% common_genes, ]

dim(expr_cgga_common)
dim(expr_tcga_common)
dim(surv_tcga)
head(surv_tcga)
head(surv_cgga)
head(surv_tcga)
head(expr_tcga_common)
# --- Prepare TCGA survival input ---

# Make sure TCGA expression colnames are truncated
colnames(expr_tcga_common) <- substr(colnames(expr_tcga_common), 1, 12)

expr_tcga_t <- as.data.frame(t(expr_tcga_common))
rownames(expr_tcga_t)

expr_tcga_t$patient <- rownames(expr_tcga_t)
head(surv_tcga)

# Replace - with . in surv_tcga patient IDs
surv_tcga <- surv_tcga %>%
  mutate(patient = gsub("-", ".", patient))

# Merge with clinical
tcga_surv_input <- surv_tcga %>%
  distinct(patient, .keep_all = TRUE) %>%   # remove duplicates if any
  inner_join(expr_tcga_t, by = "patient")

# --- Prepare CGGA survival input ---
expr_cgga_t <- as.data.frame(t(expr_cgga_common))
expr_cgga_t$patient <- rownames(expr_cgga_t)

cgga_surv_input <- surv_cgga %>%
  distinct(patient, .keep_all = TRUE) %>%
  inner_join(expr_cgga_t, by = "patient")

# --- Save as CSV ---
write.csv(tcga_surv_input, "TCGA_Survival_Input.csv", row.names = FALSE)
write.csv(cgga_surv_input, "CGGA_Survival_Input.csv", row.names = FALSE)

# Check dimensions
dim(tcga_surv_input)
dim(cgga_surv_input)
head(cgga_surv_input$patient)

head(tcga_surv_input)

# Gene columns only (drop metadata)
genes_tcga <- colnames(tcga_surv_input)[11:ncol(tcga_surv_input)]
genes_cgga <- colnames(cgga_surv_input)[11:ncol(cgga_surv_input)]

# Compare
identical(genes_tcga, genes_cgga)   # should return TRUE

## FIGURE 2A PANEL ##

library(ggrepel)

plot_volcano <- function(res, title){
  top_labels <- res %>% filter(status!="NS") %>% arrange(adj.P.Val) %>% slice(1:10)
  ggplot(res, aes(x=logFC, y=-log10(P.Value), color=status)) +
    geom_point(alpha=0.6, size=1.2) +
    scale_color_manual(values=c("Up"="firebrick","Down"="royalblue","NS"="grey80")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    ggrepel::geom_text_repel(data=top_labels, aes(label=external_gene_name),
                             size=3, max.overlaps=10) +
    labs(title=title, x="log2 Fold Change", y="-log10 P") +
    theme_bw(base_size=14) +
    theme(legend.position="bottom")
}


p_vol_tcga <- plot_volcano(res_tcga_pooled, "TCGA GBM Tumor vs All Normals")
p_vol_cgga <- plot_volcano(res_cgga_unique, "CGGA Tumor vs Normals")

ggsave(file.path("Fig2A_Volcano_TCGA.png"), p_vol_tcga, width=5, height=5)
ggsave(file.path("Fig2A_Volcano_CGGA.png"), p_vol_cgga, width=5, height=5)

# ------------------------------
# Voom log2CPM matrices
# ------------------------------
v_tcga <- readRDS("TCGA_voom_Tumor_vs_AllNormals.rds")
expr_tcga_v <- v_tcga$E

v_cgga <- readRDS("CGGA_voom.rds")
expr_cgga_v <- v_cgga$E

# Clean IDs
rownames(expr_tcga_v) <- sub("\\..*", "", rownames(expr_tcga_v))
rownames(expr_cgga_v) <- sub("\\..*", "", rownames(expr_cgga_v))

top_tcga <- res_tcga_pooled %>% arrange(adj.P.Val) %>% slice_head(n=30) %>% pull(ensembl_gene_id)
top_cgga <- res_cgga_unique %>% arrange(adj.P.Val) %>% slice_head(n=30) %>% pull(ensembl_gene_id)

# Genes that are in both matrices
common_union <- intersect(c(top_tcga, top_cgga),
                          intersect(rownames(expr_tcga_v), rownames(expr_cgga_v)))

expr_tcga_top_z <- t(scale(t(expr_tcga_v[rownames(expr_tcga_v) %in% top_tcga, ])))
expr_cgga_top_z <- t(scale(t(expr_cgga_v[rownames(expr_cgga_v) %in% top_cgga, ])))
expr_union_top_z <- t(scale(t(expr_tcga_v[rownames(expr_tcga_v) %in% common_union, ])))

map_symbols <- function(ids) {
  symbols <- gene_info$external_gene_name[match(ids, gene_info$ensembl_gene_id)]
  symbols[is.na(symbols)] <- ids
  return(symbols)
}
rownames(expr_tcga_top_z) <- map_symbols(rownames(expr_tcga_top_z))
rownames(expr_cgga_top_z) <- map_symbols(rownames(expr_cgga_top_z))
rownames(expr_union_top_z) <- map_symbols(rownames(expr_union_top_z))

ann_tcga <- data.frame(Group = ifelse(colnames(expr_tcga_top_z) %in% tumor_samples, "Tumor", "Normal"))
rownames(ann_tcga) <- colnames(expr_tcga_top_z)

ann_cgga <- data.frame(Group = ifelse(grepl("^CGGA", colnames(expr_cgga_top_z)), "Tumor", "Normal"))
rownames(ann_cgga) <- colnames(expr_cgga_top_z)

pheatmap(expr_tcga_top_z,
         annotation_col = ann_tcga,
         cluster_cols = TRUE, cluster_rows = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Figure 2C: TCGA Tumor vs Normal (Top 30 DEGs)")

pheatmap(expr_cgga_top_z,
         annotation_col = ann_cgga,
         cluster_cols = TRUE, cluster_rows = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Figure 2D: CGGA Tumor vs Normal (Top 30 DEGs)")

pheatmap(expr_union_top_z,
         annotation_col = ann_tcga,
         cluster_cols = TRUE, cluster_rows = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         main = "Figure 2E: Union Top DEGs (Common across TCGA & CGGA)",
         fontsize_row = 6)   # try 6 or 7 instead of default 10


##---------------------------------------------------------------##
##                     VENN DIAGRAM                              ##
##---------------------------------------------------------------##
map_symbols <- function(ids) {
  symbols <- gene_info$external_gene_name[match(ids, gene_info$ensembl_gene_id)]
  symbols[is.na(symbols)] <- ids
  return(symbols)
}

up_tcga_symbols   <- map_symbols(up_tcga)
down_tcga_symbols <- map_symbols(down_tcga)
up_cgga_symbols   <- map_symbols(up_cgga)
down_cgga_symbols <- map_symbols(down_cgga)
install.packages("ggVennDiagram")
library(ggVennDiagram)

# Upregulated Venn
venn_up <- list(
  TCGA = up_tcga_symbols,
  CGGA = up_cgga_symbols
)

p_up <- ggVennDiagram(venn_up, label_alpha=0, label_col="black") +
  scale_fill_gradient(low="white", high="firebrick") +
  ggtitle("Overlap of Upregulated DEGs (TCGA vs CGGA)")

ggsave("Fig2F_Venn_Upregulated_TCGA_CGGA.png", p_up, width=5, height=5)


# Downregulated Venn
venn_down <- list(
  TCGA = down_tcga_symbols,
  CGGA = down_cgga_symbols
)

p_down <- ggVennDiagram(venn_down, label_alpha=0, label_col="black") +
  scale_fill_gradient(low="white", high="royalblue") +
  ggtitle("Overlap of Downregulated DEGs (TCGA vs CGGA)")

ggsave("Fig2G_Venn_Downregulated_TCGA_CGGA.png", p_down, width=5, height=5)

