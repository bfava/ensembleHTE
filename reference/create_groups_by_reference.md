# Create Quantile Groups Using a Reference Population

Assigns all observations in `x` to quantile groups, where the group
boundaries (cutoffs) are determined by a reference subset of
observations indicated by `ref_mask`. Observations outside the reference
subset are assigned to groups based on where their values fall relative
to the reference cutoffs. When `ref_mask` is `NULL` or all `TRUE`, this
is equivalent to
[`create_groups`](https://bfava.com/ensembleHTE/reference/create_groups.md).

This is used when the grouping population differs from the analysis
population (e.g., groups defined by the ML training sample, then applied
to all data).

## Usage

``` r
create_groups_by_reference(x, n_groups, ref_mask = NULL)
```

## Arguments

- x:

  Numeric vector of all values to assign to groups

- n_groups:

  Integer, number of groups to create

- ref_mask:

  Logical vector (same length as `x`) indicating which observations form
  the reference population for computing group cutoffs. If `NULL` or all
  `TRUE`, equivalent to `create_groups(x, n_groups)`.

## Value

Integer vector of group assignments (1 to n_groups)
