# Compute cross-regression covariance

Computes the covariance between coefficient estimates from two different
regressions using the same data. This is needed to correctly compute the
standard error of the difference between estimates from two regressions.

## Usage

``` r
.calc_cov_cross_reg(X1, X2, resid1, resid2, weights = NULL, cluster_id = NULL)
```

## Arguments

- X1:

  Model matrix from regression 1

- X2:

  Model matrix from regression 2

- resid1:

  Residuals from regression 1

- resid2:

  Residuals from regression 2

- weights:

  Optional WLS weights (same for both regressions). When provided, uses
  the WLS sandwich formula with weighted bread and meat. NULL for OLS.

- cluster_id:

  Optional vector of cluster identifiers for cluster-robust
  cross-covariance. When provided, computes the meat matrix by summing
  score contributions within clusters (analogous to
  [`sandwich::vcovCL`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)).

## Value

Covariance matrix between coefficients of regression 1 (rows) and
regression 2 (cols)
