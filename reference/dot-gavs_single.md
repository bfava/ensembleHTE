# Internal function to compute GAVS for a single repetition

Internal function to compute GAVS for a single repetition

## Usage

``` r
.gavs_single(
  Y,
  predicted_values,
  fold,
  n_groups = 5,
  restrict_by = NULL,
  group_ref_idx = NULL,
  analysis_idx = NULL,
  cluster_id = NULL
)
```

## Arguments

- Y:

  Outcome vector (full length)

- predicted_values:

  Predicted values vector for this repetition (full length)

- fold:

  Fold assignment vector for this repetition (full length)

- n_groups:

  Number of groups for GAVS

- restrict_by:

  Optional restrict_by indicator for restricted ranking (factor/integer
  vector, full length)

- group_ref_idx:

  Optional logical vector (same length as Y). When provided, group
  cutoffs are computed from only the reference observations, then
  applied to all observations. If `NULL`, all observations are used for
  group formation.

- analysis_idx:

  Optional logical vector (same length as Y). When provided, the
  regression is run only on these observations (after group assignment).
  If `NULL`, all observations are used.

- cluster_id:

  Optional vector of cluster identifiers for cluster-robust SEs. When
  provided, passed to `reg_from_formula` for clustered inference.

## Value

List with group_estimates, tests (top_bottom, all, top_all), and reg
object
