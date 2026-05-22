# Extract baseline predictions for one repetition as a data.table

Extract baseline predictions for one repetition as a data.table

## Usage

``` r
.extract_baseline_rep(baseline, m, idx = NULL)
```

## Arguments

- baseline:

  The `$baseline` field of an `ensemble_hte_fit`: a data.table
  (`store_baseline = "ensemble"`) or 3D array (`store_baseline = "all"`)

- m:

  Integer repetition index (1-based)

- idx:

  Optional logical vector of length n to subset rows (NULL = all rows)

## Value

A data.table. Single column `"baseline"` for ensemble mode; one column
per algorithm for all mode.
