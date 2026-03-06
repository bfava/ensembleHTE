# Resolve group_on to a reference mask for group formation

Determines which observations form the reference population for
computing group cutoffs, based on the `group_on` argument and the fit
object.

- `"auto"`: Use the ML training population. For `ensemble_hte_fit` this
  is always all observations. For `ensemble_pred_fit` with `train_idx`,
  it is the training subset.

- `"all"`: Always use all observations for group formation.

- `"analysis"`: Use whatever observations are being analyzed (the
  `analysis_idx`).

## Usage

``` r
.resolve_group_ref_idx(group_on, analysis_idx, train_idx = NULL)
```

## Arguments

- group_on:

  Character: `"auto"`, `"all"`, or `"analysis"`

- analysis_idx:

  Logical vector (length n) indicating which observations are being
  analyzed.

- train_idx:

  Logical vector (length n) or NULL. The training index from an
  `ensemble_pred_fit`, if available.

## Value

A logical vector (length n) indicating reference observations for group
formation, or `NULL` if all observations should be used (equivalent to
all `TRUE`).
