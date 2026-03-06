# Compute difference between two coefficients with proper SE

Compute difference between two coefficients with proper SE

## Usage

``` r
.compute_coef_diff(coef, vcov, coef1, coef2)
```

## Arguments

- coef:

  Coefficient matrix from regression

- vcov:

  Variance-covariance matrix

- coef1:

  Name of first coefficient

- coef2:

  Name of second coefficient

## Value

data.table with estimate and se
