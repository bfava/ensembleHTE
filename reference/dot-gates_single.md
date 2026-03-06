# Internal function to compute GATES for a single repetition

Internal function to compute GATES for a single repetition

## Usage

``` r
.gates_single(
  Y,
  D,
  prop_score,
  weight,
  predicted_values,
  fold,
  n_groups = 5,
  restrict_by = NULL,
  controls = NULL,
  control_data = NULL,
  group_ref_idx = NULL,
  analysis_idx = NULL,
  cluster_id = NULL
)
```

## Arguments

- Y:

  Outcome vector (full length)

- D:

  Treatment vector (full length)

- prop_score:

  Propensity score vector (full length)

- weight:

  Weight vector (full length)

- predicted_values:

  Predicted values (ITE or Y) vector for this repetition (full length)

- fold:

  Fold assignment vector for this repetition (full length)

- n_groups:

  Number of groups for GATES

- restrict_by:

  Optional restrict_by indicator for restricted ranking (factor/integer
  vector, full length)

- controls:

  Optional control variables (character vector of column names)

- control_data:

  Optional data.table/data.frame with control variables (full length)

- group_ref_idx:

  Optional logical vector (same length as Y). When provided, group
  cutoffs are computed from only the reference observations, then
  applied to all observations. If `NULL`, all observations are used for
  group formation.

- analysis_idx:

  Optional logical vector (same length as Y). When provided, the
  regression is run only on these observations (after group assignment).
  If `NULL`, all observations are used for the regression.

- cluster_id:

  Optional vector of cluster identifiers for cluster-robust SEs. When
  provided, passed to `reg_from_formula` for clustered inference.

## Value

List with group_estimates, tests (top_bottom, all, top_all), and reg
object
