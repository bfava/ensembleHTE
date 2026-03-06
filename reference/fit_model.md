# Fit a model using the appropriate framework

Dispatcher function that routes to mlr3 or grf based on algorithm

## Usage

``` r
fit_model(
  .Y,
  .X,
  algorithm,
  task_type = c("regr", "classif"),
  tune = FALSE,
  tune_params = list(time = 30, cv_folds = 3, stagnation_iters = 250,
    stagnation_threshold = 0.01, measure = NULL),
  weights = NULL,
  learner_params = NULL
)
```

## Arguments

- .Y:

  Numeric vector of outcomes

- .X:

  data.frame of covariates

- algorithm:

  Character string specifying the ML algorithm

- task_type:

  Character: "regr" or "classif"

- tune:

  Logical, whether to tune hyperparameters

- tune_params:

  List of tuning parameters (see details in `ensemble_hte`)

- weights:

  Optional numeric vector of observation weights

- learner_params:

  Optional named list of parameter-value pairs for this specific mlr3
  learner (e.g., `list(num.trees = 500)` for ranger). Ignored for grf.
  Dispatched per-algorithm from `ensemble_hte`/`ensemble_pred`.
