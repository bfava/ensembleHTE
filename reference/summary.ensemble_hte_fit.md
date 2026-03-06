# Summary Method for ensemble_hte_fit Objects

Provides a summary of the estimated individual treatment effects (ITEs)
from an ensemble HTE fit, including descriptive statistics, BLP
coefficients for average treatment effect and heterogeneity, and GATES
group estimates.

## Usage

``` r
# S3 method for class 'ensemble_hte_fit'
summary(object, n_groups = 3, group_on = c("auto", "all", "analysis"), ...)
```

## Arguments

- object:

  An object of class `ensemble_hte_fit` from
  [`ensemble_hte`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md).

- n_groups:

  Integer. Number of groups for GATES analysis (default: 3).

- group_on:

  Character. How to form groups when `train_idx` was used. Passed
  through to
  [`gates`](https://bfava.com/ensembleHTE/reference/gates.md). One of
  `"auto"` (default), `"all"`, or `"analysis"`.

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object.

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
summary(fit)
summary(fit, n_groups = 5)
} # }
```
