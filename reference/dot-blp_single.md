# Internal function to compute BLP for a single repetition

Internal function to compute BLP for a single repetition

## Usage

``` r
.blp_single(
  Y,
  D,
  prop_score,
  weight,
  predicted_ite,
  controls = NULL,
  control_data = NULL,
  cluster_id = NULL
)
```

## Arguments

- Y:

  Outcome vector

- D:

  Treatment vector

- prop_score:

  Propensity score vector

- weight:

  Weight vector

- predicted_ite:

  Predicted ITE vector for this repetition

- controls:

  Optional control variables (character vector of column names)

- control_data:

  Optional data.table/data.frame with control variables

- cluster_id:

  Optional vector of cluster identifiers for cluster-robust SEs. When
  provided, passed to `reg_from_formula` for clustered inference.

## Value

data.table with BLP estimates (beta1 for ATE, beta2 for heterogeneity)
