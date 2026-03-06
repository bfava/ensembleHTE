# Summary Method for ensemble_pred_fit Objects

Provides a comprehensive summary of an ensemble prediction fit,
including descriptive statistics, prediction accuracy metrics, BLP
calibration test, and GAVS group averages. When `train_idx` was used,
accuracy metrics are computed only on training observations (where Y is
observed).

## Usage

``` r
# S3 method for class 'ensemble_pred_fit'
summary(object, n_groups = 3, ...)
```

## Arguments

- object:

  An object of class `ensemble_pred_fit` from
  [`ensemble_pred`](https://bfava.github.io/ensembleHTE/reference/ensemble_pred.md).

- n_groups:

  Number of groups for GAVS analysis (default: 3).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns a `summary.ensemble_pred_fit` object containing:

- call: The original function call

- outcome: Name of the outcome variable

- n: Total number of observations

- n_train: Number of training observations

- M: Number of repetitions

- metrics: Prediction accuracy metrics (R-squared, RMSE, MAE,
  correlation)

- blp: BLP calibration test results

- gavs: GAVS group average results

## Examples

``` r
if (FALSE) { # \dontrun{
data(microcredit)
covars <- c("age", "gender", "education", "hhinc_yrly_base",
            "css_creditscorefinal")
dat <- microcredit[, c("bank_profits_pp", covars)]
fit <- ensemble_pred(
  bank_profits_pp ~ ., data = dat,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
summary(fit)
summary(fit, n_groups = 5)
} # }
```
