# Compute difference between two regression results

Computes the difference in group estimates and test statistics between
two regression results (e.g., unrestricted vs restricted ranking). Uses
cross-regression covariance for proper SE calculation when both
regressions use the same data.

## Usage

``` r
.compute_reg_diff(res1, res2, n_groups)
```

## Arguments

- res1:

  Result from .gates_single or .gavs_single (first regression)

- res2:

  Result from .gates_single or .gavs_single (second regression)

- n_groups:

  Number of groups

## Value

List with group_diff (differences in group estimates), top_bottom_diff,
all_diff, top_all_diff
