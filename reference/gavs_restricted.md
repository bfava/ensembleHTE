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
  [`ensemble_hte()`](https://bfava.github.io/ensembleHTE/reference/ensemble_hte.md)
  or `ensemble_pred_fit` from
  [`ensemble_pred()`](https://bfava.github.io/ensembleHTE/reference/ensemble_pred.md).

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
  [`gavs`](https://bfava.github.io/ensembleHTE/reference/gavs.md) for
  details.

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
    [`gavs`](https://bfava.github.io/ensembleHTE/reference/gavs.md)).

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

# Compare unrestricted vs restricted ranking by gender
comparison <- gavs_restricted(fit, restrict_by = "gender", n_groups = 3)
print(comparison)
#> 
#> GAVS Comparison: Unrestricted vs Restricted Ranking
#> ==================================================== 
#> 
#> Strategy comparison:
#>   - Unrestricted: Rank predictions across full sample within folds
#>   - Restricted: Rank predictions within groups ('gender')
#>   - Restrict_by levels: 0, 1
#> 
#> Groups (3) defined by: predicted ITE
#> Outcome: hhinc_yrly_end
#> Observations: 1113
#> Repetitions: 3
#> 
#> Unrestricted GAVS Estimates:
#> ---------------------------------------- 
#>  Group Estimate SE      t     p-value    
#>  1     13137.12 1174.74 11.18 0.000   ***
#>  2     10631.47 1941.46 5.48  0.000   ***
#>  3     12738.40 1173.68 10.85 0.000   ***
#> 
#> Top-Bottom: -398.73 (SE: 1730.94, p = 0.818) 
#> All: 12173.14 (SE: 869.57, p = 0.000) ***
#> Top-All: 565.26 (SE: 1112.41, p = 0.611) 
#> 
#> Restricted GAVS Estimates:
#> ---------------------------------------- 
#>  Group Estimate SE      t     p-value    
#>  1     12918.12 1152.64 11.21 0.000   ***
#>  2     10920.25 1955.19 5.59  0.000   ***
#>  3     12670.28 1176.65 10.77 0.000   ***
#> 
#> Top-Bottom: -247.84 (SE: 1720.47, p = 0.885) 
#> All: 12173.14 (SE: 870.15, p = 0.000) ***
#> Top-All: 497.14 (SE: 1115.24, p = 0.656) 
#> 
#> Difference (Unrestricted - Restricted):
#> ---------------------------------------- 
#>  Group Estimate SE     t     p-value 
#>  1     219.00   284.57 0.77  0.442   
#>  2     -288.78  341.96 -0.84 0.398   
#>  3     68.12    198.43 0.34  0.731   
#> 
#> Top-Bottom Diff: -150.89 (SE: 355.67, p = 0.671) 
#> All Diff: 0.00 (SE: 32.99, p = 1.000) 
#> Top-All Diff: 68.12 (SE: 199.20, p = 0.732) 
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
plot(comparison)

# }
```
