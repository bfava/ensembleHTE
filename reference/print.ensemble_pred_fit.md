# Print Method for ensemble_pred_fit Objects

Displays a formatted summary of an ensemble prediction fit, including
data dimensions, model specification, and split-sample parameters.

## Usage

``` r
# S3 method for class 'ensemble_pred_fit'
print(x, ...)
```

## Arguments

- x:

  An object of class `ensemble_pred_fit` from
  [`ensemble_pred`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md).

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
dat <- microcredit[, c("bank_profits_pp", covars)]
fit <- ensemble_pred(
  bank_profits_pp ~ ., data = dat,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
print(fit)
} # }
```
