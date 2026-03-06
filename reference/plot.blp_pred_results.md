# Plot method for blp_pred_results objects

Creates a coefficient plot showing BLP prediction estimates (intercept
and prediction loading) with confidence intervals. Each coefficient is
shown in its own panel with an independent y-axis scale. A dashed
reference line marks the calibration target: 0 for the intercept and 1
for the prediction loading.

## Usage

``` r
# S3 method for class 'blp_pred_results'
plot(x, alpha = 0.05, ...)
```

## Arguments

- x:

  An object of class `blp_pred_results` from
  [`blp_pred()`](https://bfava.com/ensembleHTE/reference/blp_pred.md)

- alpha:

  Significance level for confidence intervals (default 0.05)

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
