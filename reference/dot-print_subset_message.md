# Internal function to print informative message about subset usage

Prints a message when the ensemble_fit was trained on a subset or when
the analysis function is evaluated on a subset. The goal is to help
users avoid mistakes when specifying subsets.

## Usage

``` r
.print_subset_message(
  func_name,
  n_total,
  n_fit,
  n_eval,
  fit_type,
  has_train_idx
)
```

## Arguments

- func_name:

  Name of the analysis function (e.g., "GAVS", "GATES", "BLP")

- n_total:

  Total number of observations in the data

- n_fit:

  Number of observations used to fit the ML model (n_train or n_total)

- n_eval:

  Number of observations used for evaluation

- fit_type:

  "hte" or "pred"

- has_train_idx:

  Whether the ensemble_fit used subset training

## Value

NULL (prints message as side effect)
