# Internal function to compute BLP for prediction (single repetition)

Internal function to compute BLP for prediction (single repetition)

## Usage

``` r
.blp_pred_single(Y, predicted_y, cluster_id = NULL)
```

## Arguments

- Y:

  Outcome vector

- predicted_y:

  Predicted Y vector for this repetition

- cluster_id:

  Optional vector of cluster identifiers for cluster-robust SEs.

## Value

data.table with BLP estimates
