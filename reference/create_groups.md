# Create Quantile Groups

Divides a numeric vector into g quantile-based groups. Attempts to use
\`cut_number()\` for equal-sized bins, falls back to rank-based groups
if needed.

## Usage

``` r
create_groups(x, n_groups)
```

## Arguments

- x:

  Numeric vector to be divided into groups

- n_groups:

  Integer, number of groups to create

## Value

Integer vector of group assignments (1 to n_groups)
