# Plot method for gates_restricted_results objects

Creates a comparison plot showing GATES estimates with confidence
intervals for both unrestricted and restricted strategies side by side.

## Usage

``` r
# S3 method for class 'gates_restricted_results'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  An object of class `gates_restricted_results` from
  [`gates_restricted()`](https://bfava.com/ensembleHTE/reference/gates_restricted.md)

- alpha:

  Significance level for confidence intervals (default 0.05)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
