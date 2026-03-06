# Plot method for clan_results objects

Creates a bar plot showing CLAN estimates (differences in covariate
means between top group and comparison group) with confidence intervals.
Bars are colored based on whether the difference is positive or
negative, and ordered from lowest to highest.

## Usage

``` r
# S3 method for class 'clan_results'
plot(
  x,
  comparison = c("top_bottom", "top_else", "top_all"),
  alpha = 0.05,
  scale = TRUE,
  ...
)
```

## Arguments

- x:

  An object of class `clan_results` from
  [`clan()`](https://bfava.github.io/ensembleHTE/reference/clan.md)

- comparison:

  Which comparison to plot: "top_bottom", "top_else", or "top_all"

- alpha:

  Significance level for confidence intervals (default 0.05)

- scale:

  Logical. If `TRUE` (default), non-binary variables are rescaled to
  standard-deviation units for plotting, making coefficients comparable
  across variables with different scales. This rescaling is applied on
  the fly regardless of whether `scale` was set in the original
  [`clan()`](https://bfava.github.io/ensembleHTE/reference/clan.md)
  call. Set `FALSE` to plot in the original units of each variable.

- ...:

  Additional arguments (currently unused)

## Value

A ggplot object (invisibly)
