# Plot method for gavs_results objects

Creates a coefficient plot showing GAVS estimates with confidence
intervals for each group, including top-bottom heterogeneity test
results.

## Usage

``` r
# S3 method for class 'gavs_results'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  An object of class `gavs_results` from
  [`gavs()`](https://bfava.com/ensembleHTE/reference/gavs.md)

- alpha:

  Significance level for confidence intervals (default 0.05)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
