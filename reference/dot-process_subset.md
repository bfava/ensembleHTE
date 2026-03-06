# Internal function to process subset argument for analysis functions

Processes the subset argument and returns a logical vector indicating
which observations to use. Works for both ensemble_hte_fit and
ensemble_pred_fit.

## Usage

``` r
.process_subset(
  subset,
  n,
  fit_type,
  has_train_idx,
  train_idx,
  using_default_outcome,
  Y
)
```

## Arguments

- subset:

  The subset argument from the analysis function. Can be:

  - NULL: smart default based on fit type and outcome

  - "train": use training observations only (pred fits only)

  - "all": use all observations

  - Logical vector: TRUE/FALSE for each observation

  - Integer vector: indices of observations to use

- n:

  Total number of observations

- fit_type:

  "hte" or "pred"

- has_train_idx:

  Whether the ensemble_fit used subset training

- train_idx:

  The training indices from ensemble_fit (if available)

- using_default_outcome:

  Whether the default outcome is being used

- Y:

  The outcome vector (to check for NAs)

## Value

A logical vector of length n indicating which observations to use
