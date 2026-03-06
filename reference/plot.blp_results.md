# Plot method for blp_results objects

Creates a coefficient plot showing BLP estimates (ATE and heterogeneity
loading) with confidence intervals. Each coefficient is shown in its own
panel with an independent y-axis scale. A dashed reference line marks
the calibration target: 0 for beta1 (ATE) and 1 for beta2 (heterogeneity
loading).

## Usage

``` r
# S3 method for class 'blp_results'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  An object of class `blp_results` from
  [`blp()`](https://bfava.github.io/ensembleHTE/reference/blp.md)

- alpha:

  Significance level for confidence intervals (default 0.05)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
