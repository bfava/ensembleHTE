# Validate Common Inputs

Validates common input parameters for ensemble functions. Internal
function.

## Usage

``` r
validate_common_inputs(M, K, algorithms, ensemble_folds, n_cores)
```

## Arguments

- M:

  Number of repetitions

- K:

  Number of folds

- algorithms:

  Vector of algorithm names

- ensemble_folds:

  Number of ensemble folds

- n_cores:

  Number of cores

## Value

The validated (and possibly deduplicated) algorithms vector, invisibly.
