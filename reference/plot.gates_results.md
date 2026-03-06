# Plot method for gates_results objects

Creates a coefficient plot showing GATES estimates with confidence
intervals for each group, including top-bottom heterogeneity test
results.

## Usage

``` r
# S3 method for class 'gates_results'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  An object of class `gates_results` from
  [`gates()`](https://bfava.github.io/ensembleHTE/reference/gates.md)

- alpha:

  Significance level for confidence intervals (default 0.05)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
