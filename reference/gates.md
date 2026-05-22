# Compute GATES (Group Average Treatment Effects)

Computes Group Average Treatment Effects (GATES), an estimand introduced
by Chernozhukov et al. (2025) to measure average treatment effects for
groups defined by predicted value quantiles.

This function works with both `ensemble_hte_fit` (groups by predicted
ITE) and `ensemble_pred_fit` (groups by predicted Y). For prediction
fits, a treatment variable must be specified.

This function implements the multiple-split estimation strategy
developed by Fava (2025), which combines predictions from multiple
machine learning algorithms into an ensemble and averages GATES
estimates across M repetitions of K-fold cross-fitting to improve
statistical power.

## Usage

``` r
gates(
  ensemble_fit,
  n_groups = 3,
  outcome = NULL,
  treatment = NULL,
  prop_score = NULL,
  controls = NULL,
  baseline_as_control = NULL,
  restrict_by = NULL,
  subset = NULL,
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

  - Numeric vector: custom outcome variable (must have same length as
    data)

  This allows computing GATES for a different outcome than the one used
  for prediction.

- treatment:

  For `ensemble_pred_fit` only. The treatment variable:

  - Character string: column name in the `data` used in
    [`ensemble_pred()`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md)

  - Numeric vector: binary treatment variable (must have same length as
    data)

  Ignored for `ensemble_hte_fit` (uses the treatment from the fit).

- prop_score:

  For `ensemble_pred_fit` only. Propensity score:

  - NULL (default): estimated as mean of treatment variable

  - Character string: column name in the `data` used in the ensemble
    function

  - Numeric value: constant propensity score for all observations

  - Numeric vector: observation-specific propensity scores

  For `ensemble_hte_fit`, uses the propensity score from the fit.

- controls:

  Character vector of control variable names to include in regression.
  These must be column names present in the `data` argument used when
  calling the ensemble function.

- baseline_as_control:

  Logical or NULL. Whether to include stored baseline predictions as
  control variables in the regression. Options:

  - NULL (default): include baseline if it was stored in the fit
    (`store_baseline != "none"`); omit otherwise

  - TRUE: always include; errors if no baseline was stored in the fit

  - FALSE: exclude even if baseline was stored

  Applies to `ensemble_hte_fit` objects only. For
  `store_baseline = "ensemble"` a single column named `"baseline"` is
  added per repetition. For `store_baseline = "all"` one column per
  algorithm is added.

- restrict_by:

  Optional. Stratification variable for restricted ranking:

  - NULL (default): unrestricted ranking across full sample within folds

  - Character string: column name in the `data` for stratified ranking

  - Numeric/factor vector: group indicator (must have same length as
    data)

  When specified, predicted values are ranked within each stratum (and
  fold), rather than across the full sample.

- subset:

  Which observations to use for the GATES analysis. Options:

  - NULL (default): uses all observations (or training obs for
    `ensemble_pred_fit` when using the default outcome with subset
    training)

  - `"train"`: uses only training observations (for `ensemble_pred_fit`
    only)

  - `"all"`: explicitly uses all observations

  - Logical vector: TRUE/FALSE for each observation (must have same
    length as data)

  - Integer vector: indices of observations to include (1-indexed)

  This allows evaluating GATES on a subset of observations. Note that
  GATES requires observations from both treatment and control groups in
  the subset.

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

An object of class `gates_results` containing:

- estimates: data.table with GATES estimates averaged across
  repetitions. Columns: `group` (integer group index, 1 = lowest
  predicted effects), `estimate` (group-specific treatment effect), `se`
  (standard error), `n_reps`, `t_value`, `p_value`

- top_bottom: data.table with the top-bottom difference test. Columns:
  `estimate`, `se`, `n_reps`, `t_value`, `p_value`

- all: data.table with the average treatment effect (weighted avg of all
  groups). Columns: `estimate`, `se`, `n_reps`, `t_value`, `p_value`

- top_all: data.table with the top minus average test. Columns:
  `estimate`, `se`, `n_reps`, `t_value`, `p_value`

- n_groups: number of groups used

- outcome: the outcome variable used for GATES

- targeted_outcome: the outcome used for prediction

- fit_type: "hte" or "pred" depending on input

- restrict_by: the restrict_by variable used (if any)

- controls: control variables used (if any)

- group_on: how groups are formed ("auto", "all", or "analysis")

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

2.  A single weighted least squares regression is run on all
    observations: \$\$Y_i = \alpha + \sum\_{g=1}^{G} \gamma_g \\
    \mathbf{1}\\i \in g\\(D_i - e(X_i)) + \varepsilon_i\$\$ where
    \\e(X_i)\\ is the propensity score. Weights are \\1 / (e(X_i)(1 -
    e(X_i)))\\. Each \\\gamma_g\\ estimates the average treatment effect
    for group \\g\\.

3.  HC1 heteroskedasticity-robust standard errors are computed (or
    cluster-robust SEs when `individual_id` was specified).

The final reported estimates and standard errors are the simple averages
of the per-repetition estimates and standard errors across all \\M\\
repetitions.

Three heterogeneity tests are reported:

- **Top-Bottom**: difference between the top and bottom groups.

- **All**: weighted average of all group effects (estimates the overall
  ATE).

- **Top-All**: difference between the top group and the overall ATE.

## References

Chernozhukov, V., Demirer, M., Duflo, E., & Fernández-Val, I. (2025).
Fisher–Schultz Lecture: Generic Machine Learning Inference on
Heterogeneous Treatment Effects in Randomized Experiments, with an
Application to Immunization in India. *Econometrica*, 93(4), 1121-1164.

Fava, B. (2025). Training and Testing with Multiple Splits: A Central
Limit Theorem for Split-Sample Estimators. *arXiv preprint
arXiv:2511.04957*.

## Examples

``` r
if (FALSE) { # \dontrun{
data(microcredit)
covars <- c("age", "gender", "education", "hhinc_yrly_base",
            "css_creditscorefinal")
dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]

fit <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat,
  prop_score = microcredit$prop_score,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
result <- gates(fit, n_groups = 3)
print(result)
plot(result)
} # }
```
