# Internal function to compute CLAN for a single repetition

Internal function to compute CLAN for a single repetition

## Usage

``` r
.clan_single(
  predicted_values,
  fold,
  n_groups,
  variables,
  variable_data,
  na_rm = FALSE,
  group_ref_idx = NULL,
  analysis_idx = NULL,
  cluster_id = NULL
)
```

## Arguments

- predicted_values:

  Predicted values (ITE or Y) vector for this repetition (full length)

- fold:

  Fold assignment vector for this repetition (full length)

- n_groups:

  Number of groups for CLAN

- variables:

  Character vector of variable names to analyze

- variable_data:

  data.table/data.frame with the variables to analyze (full length)

- na_rm:

  Logical, whether to remove NA values in calculations (default: FALSE)

- group_ref_idx:

  Optional logical vector (same length as predicted_values). When
  provided, group cutoffs are computed from only the reference
  observations, then applied to all observations. If `NULL`, all
  observations are used.

- analysis_idx:

  Optional logical vector (same length as predicted_values). When
  provided, means and SEs are computed only for these observations
  (after group assignment). If `NULL`, all observations are used.

- cluster_id:

  Optional vector of cluster identifiers for cluster-robust SEs. When
  provided, uses a regression-based approach with clustered inference
  instead of the default analytical SE formulas.

## Value

data.table with CLAN estimates for each variable
