# 🧬 Cross-Cohort Survival Biomarker Discovery

A reproducible, cross-cohort framework for survival biomarker discovery and validation using multi-omics data. Implements standardized preprocessing, differential expression, Cox and Kaplan–Meier survival modeling, meta-analysis, and permutation benchmarking across TCGA, GTEx, and CGGA datasets.

---

## 🔍 Overview

Reproducibility is a key challenge in biomarker discovery — many published signatures fail to replicate across independent cohorts due to differences in data processing, cohort composition, and modeling strategies.

This framework provides a unified and transparent workflow to:

- Harmonize multi-cohort expression and clinical data
- Identify reproducible survival-associated genes or signatures
- Quantify reproducibility through permutation-based benchmarking
- Visualize survival outcomes, hazard ratios, and cross-cohort concordance

---

## ⚙️ Features

- **Multi-cohort harmonization** — supports TCGA, GTEx, and CGGA datasets
- **Survival modeling** — univariate and multivariate Cox regression, Kaplan–Meier stratification
- **Permutation benchmarking** — empirical p-values to assess robustness
- **Meta-analysis** — integration of hazard ratios across cohorts
- **Visualization suite** — forest plots, risk score distributions, KM curves
- **Reproducibility-first design** — modular scripts and transparent code structure

---

## 🧠 Methods

| Step | Description |
|---|---|
| Preprocessing | Harmonization of expression matrices and clinical metadata (Ensembl-based, batch-corrected) |
| Differential Expression | Tumor vs normal comparison (TCGA + GTEx normal brain) |
| Survival Analysis | Cox proportional hazards modeling and Kaplan–Meier stratification |
| Permutation Testing | Randomized survival shuffling (1,000 permutations) for empirical significance |
| Meta-Analysis | Integration of survival statistics across TCGA and CGGA cohorts |

---

## 📁 Pipeline Structure

```
01_install_and_setup.R                        # Package installation and environment setup
02_preprocessing.R                            # Data loading, normalisation, harmonization
03_differential_expression.R                  # DEG analysis (tumor vs normal)
04_cross_cohort_survival_analysis.R           # Cox and KM survival modeling
05_visualization_and_survival_plots.R         # KM curves, forest plots, volcano plots
06_meta_cox_analysis_and_forest_plots.R       # Meta-analysis of hazard ratios
07_meta_permutation_and_candidate_validation.R # Permutation benchmarking
08_bootstrap_survival_reproducibility.R       # Bootstrap reproducibility scoring
```

---

## 📊 Example Outputs

- Volcano plots and survival forest plots
- Concordance of hazard ratios across cohorts
- Kaplan–Meier curves for top reproducible biomarkers
- Reproducibility distribution plots under the null model

---

## 🚀 Getting Started

```r
# Install dependencies
source("01_install_and_setup.R")

# Run pipeline sequentially
source("02_preprocessing.R")
source("03_differential_expression.R")
# ... continue through scripts 04-08
```

---

*University of Otago · Department of Biochemistry · Dunedin, NZ*
