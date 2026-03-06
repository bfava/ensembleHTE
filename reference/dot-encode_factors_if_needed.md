# Encode factor/character columns as dummy variables

Used internally for algorithms that don't support factor features
natively.

## Usage

``` r
.encode_factors_if_needed(.X)
```

## Arguments

- .X:

  data.frame of covariates (may contain factor/character columns)

## Value

data.frame with factor/character columns replaced by dummy variables
