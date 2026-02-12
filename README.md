# ensembleHTE

<!-- badges: start -->
<!-- badges: end -->

## Overview

`ensembleHTE` implements an ensemble method for learning features of heterogeneous treatment effects in randomized controlled trials. The package focuses on estimating **Group Average Treatment Effects (GATES)** by sorting individuals into groups based on predicted treatment effects from multiple machine learning algorithms.

The key innovation is combining predictions from multiple ML algorithms through Best Linear Predictor (BLP) weights while using K-fold cross-fitting for model training and separate L-fold splits for ensemble calibration. This approach uses the **entire sample** for final GATES estimation, providing higher statistical power than existing split-sample methods while maintaining proper coverage.

## Installation

You can install the development version of ensembleHTE from GitHub:

```r
# install.packages("devtools")
devtools::install_github("bfava/ensembleHTE")
```

## Quick Start

```r
library(ensembleHTE)

# 1. Fit ensemble HTE model
fit <- ensemble_hte(
  formula = Y ~ X1 + X2 + X3,
  treatment = D,
  data = mydata,
  algorithms = c("lm", "grf"),
  M = 5,  # Number of repetitions
  K = 3   # Number of cross-fitting folds
)

# 2. View results
print(fit)
summary(fit)

# 3. Analyze treatment effect heterogeneity
gates_results <- gates(fit, n_groups = 3)  # Group Average Treatment Effects
blp_results <- blp(fit)                    # Best Linear Predictor
clan_results <- clan(fit)                  # Classification Analysis

# 4. Visualize results
plot(gates_results)
plot(clan_results)
```

## Key Features

- **GATES Estimation**: Sort individuals into groups by predicted treatment effects and estimate group-specific treatment effects
- **Ensemble Learning**: Combines multiple ML algorithms (Random Forest, XGBoost, Neural Nets, etc.) using Best Linear Predictor weights
- **K-fold Cross-Fitting**: Trains ML models on K-1 folds, generates out-of-sample predictions
- **L-fold Calibration**: Uses separate fold structure to estimate optimal ensemble weights
- **Full-Sample Inference**: Uses entire dataset for final GATES regression (no data wasted on hold-out sets)
- **Higher Power**: Demonstrated improvements over CDDF and sequential aggregation methods
- **Repeated Sample-Splitting**: Aggregates across M repetitions for stability

## Main Functions

### Estimation
- `ensemble_hte()`: Fit ensemble heterogeneous treatment effect model
- `ensemble_pred()`: Fit ensemble prediction model (without treatment effects)
- `combine_ensembles()`: Combine multiple ensemble fits from different sessions

### Analysis
- `gates()`: Compute Group Average Treatment Effects (GATES)
- `blp()`: Compute Best Linear Predictor of CATE
- `clan()`: Classification Analysis - characterize high/low effect groups
- `gavs()`: Compute Group Averages for prediction tasks
- `blp_pred()`: Best Linear Predictor for predictions

### Comparison
- `gates_compare()`: Compare GATES between unrestricted and restricted ranking
- `gavs_compare()`: Compare GAVS between unrestricted and restricted ranking

## Methodology

This package implements the ensemble GATES estimator developed in Fava (2025). The approach:

1. **K-fold Cross-Fitting**: Splits data into K folds; for each fold k, trains A machine learning algorithms on the remaining K-1 folds to predict individual treatment effects (ITEs)
2. **L-fold Calibration**: Uses a separate L-fold split to estimate Best Linear Predictor weights that optimally combine the A algorithms' predictions
3. **Ensemble Aggregation**: Combines predictions from multiple algorithms using estimated BLP weights: τ̂ᵢ = Σₐ β̂ₐ τ̂ᵢ,ₐ
4. **Group Formation**: Sorts individuals by predicted ITEs (τ̂ᵢ) into J quantile groups
5. **GATES Estimation**: Runs weighted regression on the **entire sample** to estimate group-specific average treatment effects
6. **Repeated Splitting**: Repeats process M times and averages estimates for stability

The key advantage over existing methods (Chernozhukov et al. 2025, Wager & Walther 2024) is using the full sample for inference rather than holding out data, which provides higher statistical power.

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

MIT License - see LICENSE file for details
