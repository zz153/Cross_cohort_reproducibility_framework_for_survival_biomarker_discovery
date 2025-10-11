#!/usr/bin/env Rscript
# ============================================================
# 01_install_and_setup.R
# Cross-Cohort Reproducibility Framework for Survival Biomarker Discovery
# Author: Dr. Zohaib Rana, University of Otago
# ============================================================

cat("\n🔧 Initializing package installation and environment setup...\n")

# ----------------------------
# Helper function to install/load packages
# ----------------------------
install_and_load <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("📦 Installing:", pkg))
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# ============================================================
# 🧬 Bioconductor Packages
# ============================================================
bioc_pkgs <- c(
  "TCGAbiolinks", "recount3", "rtracklayer", "SummarizedExperiment",
  "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db"
)
install_and_load(bioc_pkgs, bioc = TRUE)

# ============================================================
# 📦 CRAN Packages
# ============================================================
cran_pkgs <- c(
  "R.utils", "FactoMineR", "factoextra", "tidyverse", "ggrepel",
  "pheatmap", "RColorBrewer", "UpSetR", "VennDiagram", "ggVennDiagram",
  "gridExtra", "survival", "survminer", "metafor", "biomaRt",
  "qvalue", "purrr", "stringr", "tibble", "dplyr", "tidyr", "ggplot2", "readr"
)
install_and_load(cran_pkgs, bioc = FALSE)

# ============================================================
# ✅ Verify installation
# ============================================================
cat("\n✅ Package setup complete. Loaded libraries:\n")
loaded <- c(bioc_pkgs, cran_pkgs)
print(loaded[loaded %in% .packages()])
#!/usr/bin/env Rscript
# ============================================================
# 01_install_and_setup.R
# Cross-Cohort Reproducibility Framework for Survival Biomarker Discovery
# Author: Dr. Zohaib Rana, University of Otago
# ============================================================

cat("\n🔧 Initializing package installation and environment setup...\n")

# ----------------------------
# Helper function to install/load packages
# ----------------------------
install_and_load <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("📦 Installing:", pkg))
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
    }
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# ============================================================
# 🧬 Bioconductor Packages
# ============================================================
bioc_pkgs <- c(
  "TCGAbiolinks", "recount3", "rtracklayer", "SummarizedExperiment",
  "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db"
)
install_and_load(bioc_pkgs, bioc = TRUE)

# ============================================================
# 📦 CRAN Packages
# ============================================================
cran_pkgs <- c(
  "R.utils", "FactoMineR", "factoextra", "tidyverse", "ggrepel",
  "pheatmap", "RColorBrewer", "UpSetR", "VennDiagram", "ggVennDiagram",
  "gridExtra", "survival", "survminer", "metafor", "biomaRt",
  "qvalue", "purrr", "stringr", "tibble", "dplyr", "tidyr", "ggplot2", "readr"
)
install_and_load(cran_pkgs, bioc = FALSE)

# ============================================================
# ✅ Verify installation
# ============================================================
cat("\n✅ Package setup complete. Loaded libraries:\n")
loaded <- c(bioc_pkgs, cran_pkgs)
print(loaded[loaded %in% .packages()])
library(grid)
cat("\n🚀 Environment ready. You can now proceed to 02_preprocessing.R\n")
