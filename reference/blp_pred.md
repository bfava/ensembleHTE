# Compute BLP (Best Linear Predictor) for Prediction

Computes the Best Linear Predictor (BLP) of the outcome using the
ensemble predictions (from `ensemble_pred`) or ITE predictions (from
`ensemble_hte`). This is a simple regression of Y on the predicted
values.

This function implements the multiple-split estimation strategy
developed in Fava (2025), which combines predictions from multiple
machine learning algorithms into an ensemble and averages BLP estimates
across M repetitions of K-fold cross-fitting to improve statistical
power.

## Usage

``` r
blp_pred(ensemble_fit, outcome = NULL, subset = NULL)
```

## Arguments

- ensemble_fit:

  An object of class \`ensemble_pred_fit\` from \`ensemble_pred()\` or
  \`ensemble_hte_fit\` from \`ensemble_hte()\`.

- outcome:

  Either:

  - NULL (default): uses the same outcome as in ensemble fitting

  - Character string: column name in the \`data\` used in ensemble
    fitting

  - Numeric vector: custom outcome variable (must have same length as
    data, or same length as subset if subset is provided)

  This allows computing BLP for a different outcome than the one used
  for prediction.

- subset:

  For \`ensemble_pred_fit\` with subset training (train_idx), controls
  which observations to use:

  - NULL (default): uses training observations if default outcome and
    subset training was used, otherwise uses all observations

  - "train": uses only training observations (requires train_idx in
    ensemble_fit)

  - "all": uses all observations

  For \`ensemble_hte_fit\`, this can be a logical or integer vector
  specifying which observations to include.

## Value

An object of class \`blp_pred_results\` containing:

- estimates: data.table with BLP estimates averaged across repetitions

- outcome: outcome variable used

- targeted_outcome: original outcome from ensemble fitting

- fit_type: "hte" or "pred" depending on input

- n_used: number of observations used

- M: number of repetitions

- call: the function call

## Estimation Procedure

For each repetition \\m = 1, \ldots, M\\:

1.  The ensemble predictions from repetition \\m\\ are used as the
    regressor. These predictions were generated via cross-fitting in
    [`ensemble_pred()`](https://bfava.github.io/ensembleHTE/reference/ensemble_pred.md)
    or
    [`ensemble_hte()`](https://bfava.github.io/ensembleHTE/reference/ensemble_hte.md),
    so each observation's prediction is out-of-sample.

2.  A single ordinary least squares regression is run: \$\$Y_i =
    \\alpha + \\beta \\, \\hat{Y}\_i + \\varepsilon_i\$\$ where
    \\\\hat{Y}\_i\\ is the predicted value. HC1 robust standard errors
    are computed (or cluster-robust SEs when `individual_id` was
    specified in the ensemble fit).

The final reported estimates and standard errors are the simple averages
of the per-repetition estimates and standard errors across all \\M\\
repetitions.

**Interpretation:** A coefficient (\\\\beta\\) close to 1 indicates good
calibration; significantly different from 1 suggests over- or
under-prediction.

When using `ensemble_hte_fit` objects, this allows testing whether
predicted treatment effects correlate with a different outcome variable
(e.g., an endline measure that may only be observed for a subset of the
data).

## References

Fava, B. (2025). Training and Testing with Multiple Splits: A Central
Limit Theorem for Split-Sample Estimators. *arXiv preprint
arXiv:2511.04957*.

## Examples

``` r
# \donttest{
data(microcredit)
covars <- c("age", "gender", "education", "hhinc_yrly_base",
            "css_creditscorefinal")
dat <- microcredit[, c("bank_profits_pp", covars)]

fit <- ensemble_pred(
  bank_profits_pp ~ ., data = dat,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
result <- blp_pred(fit)
print(result)
#> BLP Results (Best Linear Predictor - Prediction)
#> =================================================
#> 
#> Fit type: Prediction (ensemble_pred)
#> Outcome analyzed: bank_profits_pp
#> Observations used: 650
#> Repetitions: 3
#> 
#> Coefficients:
#>   intercept: Regression intercept (0 = well-calibrated)
#>   beta: Prediction loading (1 = well-calibrated)
#> 
#>         Term    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>    intercept       -0.00        0.60     -0.01       0.994 
#>         beta        0.98        0.07     13.59       0.000 ***
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
```
