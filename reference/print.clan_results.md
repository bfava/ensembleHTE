# Print method for clan_results objects

Displays CLAN results in two compact panels: group means (with standard
errors) and differences from the top group (with significance stars and
standard errors). When there are many variables, output is truncated to
`max_rows` variables; the full results are always accessible via
`$estimates`.

## Usage

``` r
# S3 method for class 'clan_results'
print(x, max_rows = 15, ...)
```

## Arguments

- x:

  An object of class `clan_results` from
  [`clan()`](https://bfava.github.io/ensembleHTE/reference/clan.md)

- max_rows:

  Maximum number of variables to display (default: 15). Set to `Inf` to
  show all.

- ...:

  Additional arguments (currently unused)
