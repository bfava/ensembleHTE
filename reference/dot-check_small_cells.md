# Check for small fold (or fold x restrict_by) cells and warn

Checks whether any fold (or fold x restrict_by) cell has fewer
observations than the requested number of groups across all repetitions.
When this happens, some groups will be empty in those cells and
observations are assigned deterministically to the lowest groups. This
is not an error, but users should be aware of it.

## Usage

``` r
.check_small_cells(fold_list, n_groups, restrict_by = NULL, func_name = "")
```

## Arguments

- fold_list:

  A list of integer vectors of fold assignments (one per repetition), or
  a single integer vector. When a list is provided, all repetitions are
  scanned and the worst case is reported.

- n_groups:

  Number of groups requested

- restrict_by:

  Optional factor/integer vector for restricted ranking (same length as
  each fold vector). Assumed constant across repetitions (only fold
  assignments change).

- func_name:

  Character name of the calling function (for the warning message)

## Value

Invisible NULL. Emits a warning if any cell is too small.
