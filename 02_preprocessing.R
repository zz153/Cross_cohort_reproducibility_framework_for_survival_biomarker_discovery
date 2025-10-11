# ============================================================
# 02_preprocessing.R
# Downloading TCGA, GTEx, and CGGA Expression Data
# Harmonization of Gene Identifiers and Quality Control
# ============================================================
# --------------------------------------
# 1. Download TCGA GBM gene expression
# --------------------------------------
query_tcga <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(query_tcga)
saveRDS(query_tcga, file = "TCGA_GBM_query.rds")  # Save query for reproducibility

tcga_data <- GDCprepare(query_tcga)
saveRDS(tcga_data, file = "TCGA_GBM_data.rds")    # Save SummarizedExperiment object

# --------------------------------------
# 2. Download GTEx BRAIN expression data (via recount3)
# --------------------------------------
gtex_proj <- subset(available_projects("human"), project == "BRAIN" & file_source == "gtex")
rse_brain <- create_rse(gtex_proj)
saveRDS(rse_brain, file = "GTEx_Brain_rse.rds")

# --------------------------------------
# 3. Download GENCODE v38 annotation GTF (for gene types/info)
# --------------------------------------
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
gtf_file <- "gencode.v38.annotation.gtf.gz"
gtf_unzipped <- "gencode.v38.annotation.gtf"
options(timeout = 300)
if (!file.exists(gtf_file)) {
  tryCatch({
    download.file(gtf_url, destfile = gtf_file, method = "curl")
  }, error = function(e) {
    download.file(gtf_url, destfile = gtf_file)
  })
}
if (!file.exists(gtf_unzipped)) {
  R.utils::gunzip(gtf_file, overwrite = TRUE)
}
# You can now load this GTF in downstream steps

cat("All downloads complete!\n")

# ---------------------------------------------------
# 4. Load CGGA expression and clinical data
# ---------------------------------------------------

# CGGA mRNA expression (693 samples, raw read counts)
cgga_expr <- read.delim("CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt",
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        check.names = FALSE)

# Inspect
dim(cgga_expr)
head(cgga_expr[,1:5])

# CGGA clinical data
cgga_clin <- read.delim("CGGA.mRNAseq_693_clinical.20200506.txt",
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        check.names = FALSE)

# Inspect
dim(cgga_clin)
head(cgga_clin)

# Save for downstream scripts
saveRDS(cgga_expr, file = "CGGA_mRNA_expr.rds")
saveRDS(cgga_clin, file = "CGGA_clinical.rds")

cat("✅ CGGA data successfully loaded and saved!\n")

# ---------------------------------------------------
# Harmonize and extract mRNA/lncRNA matrices for TCGA, GTEx, and CGGA
# Save harmonized counts as CSV for downstream analysis
# ---------------------------------------------------

# --- Load libraries ---
libs <- c("SummarizedExperiment", "rtracklayer", "dplyr")
invisible(lapply(libs, require, character.only = TRUE))

# --- Load Data (from script 1 outputs) ---
tcga_data <- readRDS("TCGA_GBM_data.rds")
rse_brain <- readRDS("GTEx_Brain_rse.rds")
cgga_expr <- readRDS("CGGA_mRNA_expr.rds")
cgga_clin <- readRDS("CGGA_clinical.rds")

# --- Load & Process Annotation ---
gtf <- rtracklayer::import("gencode.v38.annotation.gtf")
genes_gtf <- gtf[gtf$type == "gene"]
gene_info <- data.frame(
  ensembl_gene_id    = gsub("\\..*", "", genes_gtf$gene_id),
  gene_biotype       = genes_gtf$gene_type,
  external_gene_name = genes_gtf$gene_name,
  stringsAsFactors   = FALSE
)
saveRDS(gene_info, file = "gene_info.rds")

# ---------------------------------------------------
# 1. Extract & Save Raw Counts
# ---------------------------------------------------

# TCGA counts (from SummarizedExperiment)
tcga_counts <- assay(tcga_data, "unstranded")
rownames(tcga_counts) <- sub("\\..*", "", rownames(tcga_counts))
head(tcga_counts)
tcga_counts_unique <- tcga_counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_gene_id") %>%
  group_by(ensembl_gene_id) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("ensembl_gene_id")

# GTEx counts (from recount3 RSE)
gtex_counts <- assay(rse_brain, "raw_counts")
rownames(gtex_counts) <- sub("\\..*", "", rownames(gtex_counts))
write.csv(gtex_counts, "GTEx_Brain_counts_by_ensg.csv")

# CGGA counts (already tab-delimited input)
## GBM ID extraction
head(cgga_clin)
# Extract GBM-only cases from CGGA clinical data
cgga_gbm_clin <- cgga_clin %>%
  dplyr::filter(Histology == "GBM")

# Inspect the subset
dim(cgga_gbm_clin)
head(cgga_gbm_clin)

# Save for downstream
write.csv(cgga_gbm_clin, "CGGA_GBM_clinical_only.csv", row.names = FALSE)


## GBM Only expression matrix

# 1. Get GBM sample IDs from clinical table
cgga_gbm_ids <- cgga_gbm_clin$CGGA_ID

# 2. Keep only the columns that match these IDs, plus the gene_name column
cgga_gbm_expr <- cgga_expr %>%
  dplyr::select(gene_name, all_of(cgga_gbm_ids))

# 3. Inspect the new table
dim(cgga_gbm_expr)
head(cgga_gbm_expr[, 1:6])

# 4. Save for downstream harmonization
write.csv(cgga_gbm_expr, "CGGA_GBM_counts.csv", row.names = FALSE)

# gene_info has: ensembl_gene_id, gene_biotype, external_gene_name

# 1. Rename gene_name to match gene_info
cgga_gbm_expr_ensg <- cgga_gbm_expr %>%
  rename(external_gene_name = gene_name) %>%
  inner_join(gene_info[, c("ensembl_gene_id", "external_gene_name")],
             by = "external_gene_name")

# 2. Put ENSG as first column
cgga_gbm_expr_ensg <- cgga_gbm_expr_ensg %>%
  relocate(ensembl_gene_id, .before = external_gene_name)

# 3. Inspect
head(cgga_gbm_expr_ensg[, 1:6])
dim(cgga_gbm_expr_ensg)
library(tibble)

# Aggregate duplicated Ensembl IDs by summing counts
cgga_gbm_expr_unique <- cgga_gbm_expr_ensg %>%
  dplyr::select(-external_gene_name) %>%   # explicitly use dplyr
  group_by(ensembl_gene_id) %>%
  summarise(across(everything(), sum), .groups = "drop")

cgga_gbm_expr_unique <- cgga_gbm_expr_ensg %>%
  select(-external_gene_name) %>%                   # drop gene symbol
  group_by(ensembl_gene_id) %>%                     # group by Ensembl ID
  summarise(across(everything(), sum)) %>%          # sum counts across duplicates
  ungroup()

# Embed Ensembl IDs into row names
cgga_gbm_expr_mat <- cgga_gbm_expr_unique %>%
  column_to_rownames("ensembl_gene_id")

# ✅ Now you have a clean Ensembl × samples matrix
dim(cgga_gbm_expr_mat)
head(cgga_gbm_expr_mat[, 1:5])

# Save output
write.csv(cgga_gbm_expr_mat, "CGGA_GBM_counts_by_ensg.csv")

cat("✅ Counts and clinical files exported for TCGA, GTEx, and CGGA!\n")

# ---------------------------------------------------
# 2. Harmonize gene IDs across TCGA, GTEx, and CGGA
# ---------------------------------------------------

# Step 1: Get rownames (Ensembl IDs) for each dataset
genes_tcga <- rownames(tcga_counts_unique)
genes_gtex <- rownames(gtex_counts)
genes_cgga <- rownames(cgga_gbm_expr_mat)

# Step 2: Find common genes across all three
common_genes <- Reduce(intersect, list(genes_tcga, genes_gtex, genes_cgga))
length(common_genes)

# Step 3: Subset datasets to common genes
tcga_h <- tcga_counts_unique[common_genes, , drop = FALSE]
gtex_h <- gtex_counts[common_genes, , drop = FALSE]
cgga_h <- cgga_gbm_expr_mat[common_genes, , drop = FALSE]

# Step 4: Save harmonized datasets
saveRDS(tcga_h, "TCGA_GBM_harmonized_counts.rds")
saveRDS(gtex_h, "GTEx_Brain_harmonized_counts.rds")
saveRDS(cgga_h, "CGGA_GBM_harmonized_counts.rds")

write.csv(tcga_h, "TCGA_GBM_harmonized_counts.csv")
write.csv(gtex_h, "GTEx_Brain_harmonized_counts.csv")
write.csv(cgga_h, "CGGA_GBM_harmonized_counts.csv")

cat("✅ Harmonization complete! Common genes across datasets:", length(common_genes), "\n")

# ---------------------------------------------------
# 3. QC and Dataset Summary
# ---------------------------------------------------

# CPM calculator
calc_cpm <- function(mat) {
  mat[is.na(mat)] <- 0
  log2(edgeR::cpm(mat, log = FALSE) + 1)
}
tcga_cpm <- calc_cpm(tcga_h)
gtex_cpm <- calc_cpm(gtex_h)
cgga_cpm <- calc_cpm(cgga_h)
dim(cgga_cpm)
# Identify TCGA tumor vs normal
tcga_meta <- colData(tcga_data)
tumor_samples  <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])
normal_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])

# (A) Median log2 CPM
meds <- data.frame(
  sample = c(colnames(tcga_h), colnames(gtex_h), colnames(cgga_h)),
  group  = c(
    ifelse(colnames(tcga_h) %in% tumor_samples, "TCGA_Tumor", "TCGA_Normal"),
    rep("GTEx_Normal", ncol(gtex_h)),
    rep("CGGA_Tumor", ncol(cgga_h))
  ),
  median_cpm = c(apply(tcga_cpm, 2, median),
                 apply(gtex_cpm, 2, median),
                 apply(cgga_cpm, 2, median))
)
pA <- ggplot(meds, aes(x = group, y = median_cpm, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  labs(title = "(A) Sample-wise median log2 CPM",
       x = "", y = "Sample median (log2 CPM)") +
  scale_fill_manual(values = c("TCGA_Normal"="#009E73",
                               "GTEx_Normal"="#56B4E9",
                               "TCGA_Tumor"="#D55E00",
                               "CGGA_Tumor"="#CC79A7")) +
  theme(legend.position="none")

# (B) IQR log2 CPM
iqrs <- data.frame(
  sample = meds$sample,
  group  = meds$group,
  iqr_cpm = c(apply(tcga_cpm, 2, IQR),
              apply(gtex_cpm, 2, IQR),
              apply(cgga_cpm, 2, IQR))
)
pB <- ggplot(iqrs, aes(x = group, y = iqr_cpm, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  labs(title = "(B) Sample-wise IQR of log2 CPM",
       x = "", y = "Sample IQR (log2 CPM)") +
  scale_fill_manual(values = c("TCGA_Normal"="#009E73",
                               "GTEx_Normal"="#56B4E9",
                               "TCGA_Tumor"="#D55E00",
                               "CGGA_Tumor"="#CC79A7")) +
  theme(legend.position="none")

# (C) Dataset Summary Table
summary_tbl <- data.frame(
  DataType = "Gene Expression",
  Subset   = c("All Samples","TCGA Tumor","TCGA Normal","GTEx Normal","CGGA Tumor"),
  Features = c(length(common_genes),
               length(common_genes),
               length(common_genes),
               length(common_genes),
               length(common_genes)),
  Samples  = c(ncol(tcga_h) + ncol(gtex_h) + ncol(cgga_h),
               length(tumor_samples),
               length(normal_samples),
               ncol(gtex_h),
               ncol(cgga_h))
)
pC <- gridExtra::tableGrob(summary_tbl, rows=NULL)

# Save outputs
qc_plot <- ggpubr::ggarrange(pA, pB, ncol = 2, labels = c("A","B"))
ggsave("Figure1_QC_panels_TCGA_GTEx_CGGA.pdf", qc_plot, width=9, height=4)
ggsave("Figure1C_SummaryTable_TCGA_GTEx_CGGA.pdf", pC, width=6, height=2.8)

cat("✅ QC plots and summary table saved!\n")

# ---------------------------------------------------
# 4. PCA for Batch Effect Assessment
# ---------------------------------------------------
install.packages("FactoMineR")
install.packages("factoextra")
suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
})

# --- Build combined logCPM matrix ---
all_cpm <- cbind(tcga_cpm, gtex_cpm, cgga_cpm)

# --- Metadata for PCA plot ---
pca_meta <- data.frame(
  sample = colnames(all_cpm),
  group = c(
    ifelse(colnames(tcga_cpm) %in% tumor_samples, "TCGA_Tumor", "TCGA_Normal"),
    rep("GTEx_Normal", ncol(gtex_cpm)),
    rep("CGGA_Tumor", ncol(cgga_cpm))
  ),
  batch = c(
    rep("TCGA", ncol(tcga_cpm)),
    rep("GTEx", ncol(gtex_cpm)),
    rep("CGGA", ncol(cgga_cpm))
  )
)

# --- Run PCA (remove zero-variance genes) ---
# Identify genes with non-zero variance
nzv_genes <- apply(all_cpm, 1, function(x) var(x, na.rm = TRUE) > 0)

all_cpm_nzv <- all_cpm[nzv_genes, ]

# Now run PCA
pca_res <- prcomp(t(all_cpm_nzv), scale. = TRUE)

# --- Extract variance explained ---
pvar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2], 1)

# --- Plot PCA ---
pPCA <- ggplot(data.frame(pca_res$x[,1:2], pca_meta),
               aes(x = PC1, y = PC2, color = group, shape = batch)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "(C) PCA of Samples (log2 CPM)",
    x = paste0("PC1 (", pvar[1], "%)"),
    y = paste0("PC2 (", pvar[2], "%)")
  ) +
  scale_color_manual(values = c("TCGA_Normal"="#009E73",
                                "GTEx_Normal"="#56B4E9",
                                "TCGA_Tumor"="#D55E00",
                                "CGGA_Tumor"="#CC79A7"))

ggsave("Figure1D_PCA_TCGA_GTEx_CGGA.pdf", pPCA, width = 6, height = 5)

cat("✅ PCA plot saved as Figure1D_PCA_TCGA_GTEx_CGGA.pdf\n")

# --- Run PCA (remove zero-variance genes) ---
nzv_genes <- apply(all_cpm, 1, function(x) var(x, na.rm = TRUE) > 0)
all_cpm_nzv <- all_cpm[nzv_genes, ]

# PCA
pca_res <- prcomp(t(all_cpm_nzv), scale. = TRUE)

# Variance explained
pvar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:2], 1)

# --- Remove outliers ---
# Define threshold: e.g. samples > ±4 SD on PC1 or PC2
pc_df <- data.frame(pca_res$x[,1:2], pca_meta)
z_pc1 <- scale(pc_df$PC1)
z_pc2 <- scale(pc_df$PC2)
pc_df$keep <- !(abs(z_pc1) > 4 | abs(z_pc2) > 4)

pc_df_filt <- pc_df[pc_df$keep, ]
pca_res$x <- pca_res$x[pc_df$keep, 1:2]

# --- Plot PCA (color only, no shapes) ---
pPCA <- ggplot(pc_df_filt, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw() +
  labs(
    title = "(C) PCA of Samples (log2 CPM, outliers removed)",
    x = paste0("PC1 (", pvar[1], "%)"),
    y = paste0("PC2 (", pvar[2], "%)")
  ) +
  scale_color_manual(values = c(
    "TCGA_Normal"="#009E73",
    "GTEx_Normal"="#56B4E9",
    "TCGA_Tumor"="#D55E00",
    "CGGA_Tumor"="#CC79A7"
  ))

ggsave("Figure1D_PCA_TCGA_GTEx_CGGA_clean.pdf", pPCA, width = 6, height = 5)

cat("✅ Clean PCA plot saved with outliers removed!\n")


