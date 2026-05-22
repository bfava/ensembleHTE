# Compute CLAN (Classification Analysis)

Computes Classification Analysis (CLAN), an analysis introduced by
Chernozhukov et al. (2025) to characterize which units have the highest
and lowest predicted treatment effects or predicted outcomes.

This function implements the multiple-split estimation strategy
developed by Fava (2025), which combines predictions from multiple
machine learning algorithms into an ensemble and averages CLAN estimates
across M repetitions of K-fold cross-fitting to improve statistical
power.

## Usage

``` r
clan(
  ensemble_fit,
  variables = NULL,
  n_groups = 3,
  na_rm = FALSE,
  scale = FALSE,
  subset = NULL,
  group_on = c("auto", "all", "analysis")
)
```

## Arguments

- ensemble_fit:

  An object of class \`ensemble_hte_fit\` from \`ensemble_hte()\` or
  \`ensemble_pred_fit\` from \`ensemble_pred()\`.

- variables:

  Either:

  - `NULL` (default): Uses all covariates from the ensemble fit (i.e.,
    `names(ensemble_fit$X)`)

  - Character vector of variable names present in the \`data\` used in
    the ensemble function

  - A data.frame/data.table with covariates to analyze (must have same
    number of rows as data)

- n_groups:

  Number of groups to divide the sample into (default: 3). CLAN compares
  the top group (highest predicted ITE or Y) to others.

- na_rm:

  Logical, whether to remove NA values when computing means and
  variances (default: FALSE). If FALSE and NAs are present, results will
  be NA.

- scale:

  Logical, whether to scale non-binary variables to have mean 0 and
  standard deviation 1 before computing differences (default: FALSE).
  This makes coefficients comparable across variables with different
  scales.

- subset:

  Which observations to use for the CLAN analysis. Options:

  - NULL (default): uses all observations (or training obs for
    `ensemble_pred_fit` when using the default outcome with subset
    training)

  - `"train"`: uses only training observations (for `ensemble_pred_fit`
    only)

  - `"all"`: explicitly uses all observations

  - Logical vector: TRUE/FALSE for each observation (must have same
    length as data)

  - Integer vector: indices of observations to include (1-indexed)

  This allows evaluating CLAN on a subset of observations.

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

An object of class \`clan_results\` containing:

- estimates: data.table with CLAN estimates averaged across repetitions,
  including group means (`mean_top`, `mean_bottom`, `mean_else`,
  `mean_all`), their standard errors (`se_top`, `se_bottom`, `se_else`,
  `se_all`), differences (`diff_top_bottom`, `diff_top_else`,
  `diff_top_all`), difference standard errors (`se_diff_top_bottom`,
  `se_diff_top_else`, `se_diff_top_all`), t-values (`t_diff_top_bottom`,
  etc.), and p-values (`p_diff_top_bottom`, etc.)

- variables: variables analyzed

- n_groups: number of groups used

- targeted_outcome: the outcome used for prediction

- fit_type: "hte" or "pred" depending on input

- M: number of repetitions

- scaled: whether variables were scaled

- variable_sds: named numeric vector of per-variable standard deviations
  (NA for binary variables). Used by the `plot` method for on-the-fly
  rescaling.

- group_on: how groups are formed (\\auto\\, \\all\\, or \\analysis\\)

- call: the function call

## Estimation Procedure

For each repetition \\m = 1, \ldots, M\\:

1.  Observations are assigned to `n_groups` quantile-based groups by
    ranking the ensemble predictions from repetition \\m\\ *within each
    fold*. Group 1 contains the lowest predicted values and group
    `n_groups` the highest. Forming groups within folds ensures
    independence between group assignment and out-of-sample predictions.

2.  For each covariate, the mean is computed within the top group, the
    bottom group, the "else" group (all groups except the top), and all
    observations. Standard errors are computed analytically from group
    variances (or via regression with cluster-robust SEs when
    `individual_id` was specified).

3.  Three differences are computed: top minus bottom, top minus else,
    and top minus all.

The final reported estimates and standard errors are the simple averages
of the per-repetition estimates and standard errors across all \\M\\
repetitions.

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
result <- clan(fit)
print(result)
plot(result)
} # }
```
