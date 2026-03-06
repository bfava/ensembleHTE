# Print a "Data Usage" block for analysis result objects

Displays how data was partitioned between ML training, the analysis
subset, and group formation. Only prints when there is something
non-trivial to report (i.e., when not all observations are used for
everything).

## Usage

``` r
.print_data_usage(n, n_train, n_used, group_on = NULL, n_groups = NULL)
```

## Arguments

- n:

  Total observations in the original data

- n_train:

  Observations used to train the ML model

- n_used:

  Observations used for this analysis

- group_on:

  group_on value (\\auto\\, \\all\\, \\analysis\\), or NULL for
  functions like BLP that do not form groups

- n_groups:

  Number of groups (NULL for BLP)
