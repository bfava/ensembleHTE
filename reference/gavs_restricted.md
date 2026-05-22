# Compare GAVS: Unrestricted vs Restricted Ranking

Compares Group Averages (GAVS) between an unrestricted strategy (ranking
predictions across the full sample) and a restricted strategy (ranking
predictions within groups defined by a restrict_by variable). Both
strategies use the \*\*same\*\* ensemble fit — they differ only in how
observations are ranked into groups.

This function implements the analysis from Fava (2025) to test whether
ranking by predicted values within subgroups (e.g., income quintiles,
education levels) yields different targeting results than global
ranking.

## Usage

``` r
gavs_restricted(
  ensemble_fit,
  restrict_by,
  n_groups = 3,
  outcome = NULL,
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

- restrict_by:

  The variable defining groups for restricted ranking:

  - Character string: column name in the `data` used when fitting the
    ensemble model. The variable must exist in that data at the time of
    fitting; columns added to the data.frame after fitting are not
    available.

  - Numeric/factor vector: group indicator (must have same length as
    data)

- n_groups:

  Number of groups to divide the sample into (default: 3)

- outcome:

  Either:

  - NULL (default): uses the same outcome as in the ensemble function

  - Character string: column name in the `data` used in the ensemble
    function

  - Numeric vector: custom outcome variable (must have appropriate
    length)

- subset:

  Which observations to use for the GAVS comparison. Options:

  - NULL (default): uses all observations (or training obs for
    `ensemble_pred_fit` when using the default outcome with subset
    training)

  - `"train"`: uses only training observations (for `ensemble_pred_fit`
    only)

  - `"all"`: explicitly uses all observations

  - Logical vector: TRUE/FALSE for each observation (must have same
    length as data)

  - Integer vector: indices of observations to include (1-indexed)

  This allows evaluating GAVS comparison on a subset of observations.

- group_on:

  Character controlling which observations define the quantile cutoffs
  used to form groups. One of `"auto"` (default), `"all"`, or
  `"analysis"`. See
  [`gavs`](https://bfava.com/ensembleHTE/reference/gavs.md) for details.

## Value

An object of class `gavs_restricted_results` containing:

- unrestricted: `gavs_results` object for unrestricted strategy

- restricted: `gavs_results` object for restricted strategy

- difference: data.table with the difference (unrestricted - restricted)
  for each group, with properly computed standard errors

- top_bottom_diff: data.table with difference in top-bottom estimates

- strata_var: name of the stratification variable

- strata_levels: unique levels of the stratification variable

- n_groups: number of groups used

- outcome: the outcome variable used

- targeted_outcome: the outcome used for prediction

- fit_type: "hte" or "pred" depending on input

- n_used: number of observations used

- M: number of repetitions

- call: the function call

## Estimation Procedure

For each repetition \\m = 1, \ldots, M\\, two GAVS regressions are run
on the same data using the same ensemble predictions from repetition
\\m\\:

1.  **Unrestricted**: groups are formed by ranking predictions within
    each fold (as in
    [`gavs`](https://bfava.com/ensembleHTE/reference/gavs.md)).

2.  **Restricted**: groups are formed by ranking predictions within each
    fold **and** within each level of the `restrict_by` variable.

The difference between unrestricted and restricted estimates is computed
for each repetition, with standard errors accounting for the
cross-covariance between the two regressions (since they share the same
data). The final reported estimates and standard errors are the simple
averages across all \\M\\ repetitions.

## References

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

# Compare unrestricted vs restricted ranking by gender
comparison <- gavs_restricted(fit, restrict_by = "gender", n_groups = 3)
print(comparison)
plot(comparison)
} # }
```
