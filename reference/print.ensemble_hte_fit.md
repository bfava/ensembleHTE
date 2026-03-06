# Print Method for ensemble_hte_fit Objects

Displays a formatted summary of an ensemble HTE fit, including data
dimensions, model specification, and split-sample parameters.

## Usage

``` r
# S3 method for class 'ensemble_hte_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `ensemble_hte_fit` from
  [`ensemble_hte`](https://bfava.github.io/ensembleHTE/reference/ensemble_hte.md).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the input object `x`.

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
print(fit)
} # }
```
