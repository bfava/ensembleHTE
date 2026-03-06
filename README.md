# ensembleHTE

<!-- badges: start -->
[![R-CMD-check](https://github.com/bfava/ensembleHTE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bfava/ensembleHTE/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bfava/ensembleHTE/actions/workflows/pkgdown.yaml/badge.svg)](https://bfava.github.io/ensembleHTE/)
<!-- badges: end -->

## Overview

`ensembleHTE` detects and characterizes heterogeneous treatment effects (HTE) in randomized experiments, and supports standard prediction tasks with valid downstream inference. It uses **repeated cross-fitting**—where every observation is used for both training and evaluation—and **combines predictions from multiple ML algorithms** to produce stable estimates and avoid the multiple-testing issues that arise from having to pick a single best algorithm.

Two main entry points:

- **`ensemble_hte()`** — estimate heterogeneous *treatment* effects.
- **`ensemble_pred()`** — predict an outcome (no treatment structure needed).

A suite of analysis functions—`blp()`, `gates()`, `gavs()`, `clan()`, and their restricted counterparts—lets you test for heterogeneity, characterize who is in each group, and compare targeting strategies.

For the statistical foundations, see [Fava (2025)](https://bfava.com/files/Bruno_Fava_JMP.pdf).

## Installation

You can install the development version of ensembleHTE from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bfava/ensembleHTE")
```

## Quick Start

```r
library(ensembleHTE)
data(microcredit)

# 1. Fit ensemble HTE model
covariates <- c("css_creditscorefinal", "lower_window",
                "own_anybus", "max_yearsinbusiness", "css_assetvalue")

fit <- ensemble_hte(
  Y = "exp_yrly_end",
  D = "treat",
  X = covariates,
  data       = microcredit,
  prop_score = "prop_score",
  algorithms = c("grf", "glmnet", "xgboost"),
  M = 20,
  K = 4
)

# 2. View results
print(fit)
summary(fit)

# 3. Analyze treatment effect heterogeneity
blp_results   <- blp(fit)                    # Best Linear Predictor
gates_results <- gates(fit, n_groups = 3)    # Group Average Treatment Effects
clan_results  <- clan(fit, n_groups = 3)     # Classification Analysis

# 4. Visualize
plot(gates_results)
plot(clan_results)

# --- Prediction task (no treatment structure) ---

# Predict bank profits (observed only for borrowers who took a loan)
has_loan <- microcredit$treat == 1 & microcredit$loan_size > 0

fit_pred <- ensemble_pred(
  Y    = "bank_profits_pp",
  X    = covariates,
  data = microcredit,
  train_idx  = has_loan,
  algorithms = c("grf", "glmnet", "xgboost"),
  M = 20,
  K = 4
)

summary(fit_pred)
gavs_results <- gavs(fit_pred, n_groups = 3)
plot(gavs_results)
```

## Main Functions

### Estimation
- `ensemble_hte()`: Fit ensemble heterogeneous treatment effect model
- `ensemble_pred()`: Fit ensemble prediction model (without treatment effects)

### Analysis
- `blp()` / `blp_pred()`: Best Linear Predictor (test for heterogeneity / calibration)
- `gates()`: Group Average Treatment Effects
- `gavs()`: Group Averages for prediction tasks
- `clan()`: Classification Analysis — characterize high/low effect groups

### Restricted comparisons
- `gates_restricted()`: Compare GATES under unrestricted vs. equity-constrained targeting
- `gavs_restricted()`: Compare GAVS under unrestricted vs. equity-constrained targeting

## Getting Help

If you encounter a bug or have a feature request, please [file an issue](https://github.com/bfava/ensembleHTE/issues).

## Citation

If you use this package in your research, please cite:

```bibtex
@article{fava2025training,
  title={Training and Testing with Multiple Splits: A Central Limit Theorem for Split-Sample Estimators},
  author={Fava, Bruno},
  journal={arXiv preprint arXiv:2511.04957},
  year={2025}
}
```

## License

Apache License 2.0 — see LICENSE file for details
