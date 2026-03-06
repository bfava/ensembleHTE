# Plot method for gavs_restricted_results objects

Creates a comparison plot showing GAVS estimates with confidence
intervals for both unrestricted and restricted strategies side by side.

## Usage

``` r
# S3 method for class 'gavs_restricted_results'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  An object of class `gavs_restricted_results` from
  [`gavs_restricted()`](https://bfava.com/ensembleHTE/reference/gavs_restricted.md)

- alpha:

  Significance level for confidence intervals (default 0.05)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
