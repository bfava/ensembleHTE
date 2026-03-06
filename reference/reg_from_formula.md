# Regression from Formula with Robust Standard Errors

Fits a linear regression model from a formula (or formula string) with
optional weights and computes HC1 robust standard errors. When
`cluster_id` is provided, computes cluster-robust (CR1) standard errors
using
[`sandwich::vcovCL`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
instead of heteroskedasticity-robust (HC1) standard errors.

## Usage

``` r
reg_from_formula(formula, data, weights = NULL, cluster_id = NULL)
```

## Arguments

- formula:

  Formula object or character string representing the regression formula

- data:

  Data frame containing the regression data

- weights:

  Optional numeric vector of weights for weighted regression. Must have
  the same length as the number of rows in data.

- cluster_id:

  Optional vector identifying clusters (e.g., individual IDs in panel
  data). When provided, cluster-robust standard errors are computed via
  [`sandwich::vcovCL`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html).
  When `NULL` (default), heteroskedasticity- robust (HC1) standard
  errors are computed.

## Value

List containing:

- coef: Coefficient test results with robust standard errors

- vcov: Robust variance-covariance matrix (HC1 or cluster-robust)

- model_matrix: Model matrix from the regression

- residuals: Regression residuals

- fitted: Fitted values
