# Fit a model using mlr3

Fit a model using mlr3

## Usage

``` r
.fit_model_mlr3(
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

- learner_params:

  Optional named list of parameter-value pairs for this specific mlr3
  learner (e.g., `list(num.trees = 500)` for ranger). Applied after
  algorithm-specific defaults and overrides them if there is a conflict.
