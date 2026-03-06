# Predict Individual Treatment Effects

Internal function that predicts individual treatment effects (ITEs)
using the specified metalearner strategy and machine learning algorithm.

## Usage

``` r
predict_ite(
  Y,
  X,
  D,
  prop_score,
  W,
  train_idx,
  algorithm,
  metalearner = c("t", "s", "x", "r"),
  r_learner = "grf",
  task_type = c("regr", "classif"),
  tune = FALSE,
  tune_params = list(time = 30, cv_folds = 3, stagnation_iters = 250,
    stagnation_threshold = 0.01, measure = NULL),
  test_idx = NULL,
  learner_params = NULL
)
```

## Arguments

- Y:

  Numeric vector of outcomes

- X:

  data.frame of covariates

- D:

  Numeric vector of treatment indicators (0/1)

- prop_score:

  Numeric vector of propensity scores

- W:

  Numeric vector of inverse propensity weights

- train_idx:

  Logical vector indicating training observations

- algorithm:

  Character string specifying the ML algorithm

- metalearner:

  Character string specifying the metalearner strategy

- r_learner:

  Character string specifying the R-learner algorithm

- task_type:

  Character string: "regr" or "classif"

- tune:

  Logical, whether to tune hyperparameters

- tune_params:

  List of tuning parameters

- test_idx:

  Optional logical vector indicating which observations to predict for.
  If `NULL` (default), predicts for `!train_idx`.

- learner_params:

  Optional named list of parameter-value pairs for this specific mlr3
  learner (e.g., `list(num.trees = 500)` for ranger).

## Value

List with predicted_y0, predicted_y1, and predicted_ite
