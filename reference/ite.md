# Extract ITE Predictions as a Matrix

Extracts the individual treatment effect (ITE) predictions from an
`ensemble_hte_fit` object and returns them as a plain numeric matrix.
Each column corresponds to one repetition (rep_1, rep_2, ..., rep_M).

This is a convenience function to avoid potential issues with
data.table-to-matrix coercion when working with `fit$ite` directly.

## Usage

``` r
ite(fit)
```

## Arguments

- fit:

  An object of class `ensemble_hte_fit` from
  [`ensemble_hte`](https://bfava.github.io/ensembleHTE/reference/ensemble_hte.md).

## Value

A numeric matrix with `n` rows and `M` columns, where each column
contains the ITE predictions from one repetition.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
                    algorithms = c("lm", "grf"), M = 5, K = 3)

ite_mat <- ite(fit)
} # }
```
