# Compute GAVS (Group Averages)

Computes Group Averages (GAVS), which measures average outcomes for
groups defined by predicted value quantiles. This works for both
heterogeneous treatment effect estimation (groups by predicted ITE) and
prediction problems (groups by predicted Y).

This function implements the multiple-split estimation strategy
developed in Fava (2025), which combines predictions from multiple
machine learning algorithms into an ensemble and averages GAVS estimates
across M repetitions of K-fold cross-fitting to improve statistical
power.

## Usage

``` r
gavs(
  ensemble_fit,
  n_groups = 3,
  outcome = NULL,
  subset = NULL,
  restrict_by = NULL,
  group_on = c("auto", "all", "analysis")
)
```

## Arguments

- ensemble_fit:

  An object of class `ensemble_hte_fit` from
  [`ensemble_hte()`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md)
  or `ensemble_pred_fit` from
  [`ensemble_pred()`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md).

- n_groups:

  Number of groups to divide the sample into (default: 3)

- outcome:

  Either:

  - NULL (default): uses the same outcome as in the ensemble function

  - Character string: column name in the `data` used in the ensemble
    function

  - Numeric vector: custom outcome variable (must have appropriate
    length)

  This allows computing GAVS for a different outcome than the one used
  for prediction.

- subset:

  Controls which observations to use for evaluation:

  - `NULL` (default): For `ensemble_pred_fit` with `train_idx`, uses
    training observations when `outcome = NULL`; otherwise uses all. For
    `ensemble_hte_fit`, uses all observations.

  - `"train"`: Use only training observations (only valid for
    `ensemble_pred_fit` with `train_idx`).

  - `"all"`: Use all observations.

  - Logical vector: TRUE/FALSE for each observation (length must equal
    number of rows in data).

  - Integer vector: Indices of observations to use.

  See **Subsample Usage** section for guidance on when to use each
  option.

- restrict_by:

  Optional. Stratification variable for restricted ranking:

  - NULL (default): unrestricted ranking across full sample within folds

  - Character string: column name in the `data` for stratified ranking

  - Numeric/factor vector: group indicator (must have same length as
    data)

  When specified, predicted values are ranked within each stratum (and
  fold), rather than across the full sample.

- group_on:

  Character controlling which observations define the quantile cutoffs
  used to form groups. One of:

  - `"auto"` (default): Uses the ML training population. For
    `ensemble_hte_fit` this is all observations. For `ensemble_pred_fit`
    with `train_idx`, it is the training subset. This ensures an
    observation's group assignment does not change when you vary the
    analysis subset.

  - `"all"`: Always form groups using all observations.

  - `"analysis"`: Form groups within whatever observations are being
    analyzed (i.e. the `subset`).

  Has no effect when `subset = NULL` and all observations are used.

## Value

An object of class `gavs_results` containing:

- estimates: data.table with GAVS estimates averaged across repetitions.
  Columns: `group` (integer group index, 1 = lowest predicted values),
  `estimate` (group-specific mean outcome), `se` (standard error),
  `n_reps`, `t_value`, `p_value`

- top_bottom: data.table with the top-bottom difference test. Columns:
  `estimate`, `se`, `n_reps`, `t_value`, `p_value`

- all: data.table with the overall mean (weighted avg of all groups).
  Columns: `estimate`, `se`, `n_reps`, `t_value`, `p_value`

- top_all: data.table with the top minus average test. Columns:
  `estimate`, `se`, `n_reps`, `t_value`, `p_value`

- n_groups: number of groups used

- outcome: the outcome variable used for GAVS

- targeted_outcome: the outcome used for prediction

- fit_type: "hte" or "pred" depending on input

- restrict_by: the restrict_by variable used (if any)

- group_on: how groups are formed ("auto", "all", or "analysis")

- n_used: number of observations used

- M: number of repetitions

- call: the function call

## Estimation Procedure

For each repetition \\m = 1, \ldots, M\\:

1.  Observations are assigned to `n_groups` quantile-based groups by
    ranking the ensemble predictions from repetition \\m\\ *within each
    fold*. Group 1 contains the lowest predicted values and group
    `n_groups` the highest. Forming groups within folds ensures that
    group assignment is independent of the model used to generate
    predictions for that observation (since predictions are
    out-of-sample within each fold).

2.  A single ordinary least squares regression is run on all
    observations: \$\$Y_i = \sum\_{g=1}^{G} \mu_g \\ \mathbf{1}\\i \in
    g\\ + \varepsilon_i\$\$ This is a regression of Y on group dummies
    (no intercept), so each \\\mu_g\\ directly estimates the average
    outcome in group \\g\\.

3.  HC1 heteroskedasticity-robust standard errors are computed (or
    cluster-robust SEs when `individual_id` was specified).

The final reported estimates and standard errors are the simple averages
of the per-repetition estimates and standard errors across all \\M\\
repetitions.

Three summary tests are reported:

- **Top-Bottom**: difference between top and bottom group averages.

- **All**: weighted average of all groups (estimates the overall mean).

- **Top-All**: difference between the top group and the overall mean.

## Subsample Usage

The `subset` parameter controls which observations are used for
evaluation. This is useful when:

- The ML model was trained on a subset (e.g., using `train_idx` in
  [`ensemble_pred()`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md))
  and you want to evaluate on the same or different subset.

- You want to evaluate treatment effect targeting on an outcome that is
  only observed for a subset of observations.

A message is printed when either the ML model was trained on a subset or
the evaluation uses a subset, to help avoid unintended subsetting.

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
dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]

fit <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat,
  prop_score = microcredit$prop_score,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
#> Warning: Some propensity scores are below 0.20 or above 0.80. This package is designed for randomized controlled trials (RCTs), where propensity scores are typically well-balanced. Extreme propensity scores may indicate an observational study or a heavily unbalanced design. Please verify your experimental design.
result <- gavs(fit, n_groups = 3)
print(result)
#> GAVS Results (Group Averages)
#> =============================
#> 
#> Outcome analyzed: hhinc_yrly_end
#> Number of groups: 3
#> Repetitions: 3
#> 
#> Group Average Outcomes (groups by predicted ITE):
#> 
#>   Group    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>       1    11429.93     1680.35      6.80       0.000 ***
#>       2    12158.02     1500.99      8.10       0.000 ***
#>       3    12931.36     1298.40      9.96       0.000 ***
#> 
#> Heterogeneity Tests:
#>   ----------------------------------------------------
#>           Test    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>     Top-Bottom     1501.43     2132.22      0.70       0.481 
#>        Top-All      758.21     1152.55      0.66       0.511 
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
plot(result)

# }
```
