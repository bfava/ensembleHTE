# Fit Ensemble Prediction Model

Predicts an outcome variable Y using an ensemble of machine learning
algorithms combined with multiple sample splitting. This function
implements the estimation strategy developed by Fava (2025), which
improves statistical power by averaging predictions across M repetitions
of K-fold cross-fitting.

Unlike
[`ensemble_hte`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md)
which estimates heterogeneous treatment effects, this function performs
standard prediction of Y given X without any treatment variable or
causal structure.

The function supports two interfaces:

- **Formula interface**: Specify `formula` and `data`

- **Matrix interface**: Specify `Y` and `X` directly

In the matrix interface, `Y` and `X` can be provided as column name(s)
instead of raw data. When `Y` is a single string and/or `X` is a
character vector, the corresponding columns are extracted from `data`.

## Usage

``` r
ensemble_pred(
  formula = NULL,
  data = NULL,
  Y = NULL,
  X = NULL,
  train_idx = NULL,
  M = 2,
  K = 3,
  algorithms = c("lm", "grf"),
  ensemble_folds = 5,
  task_type = NULL,
  scale_covariates = TRUE,
  tune = FALSE,
  tune_params = list(time = 30, cv_folds = 3, stagnation_iters = 250,
    stagnation_threshold = 0.01, measure = NULL),
  learner_params = NULL,
  ensemble_strategy = c("cv", "average"),
  individual_id = NULL,
  n_cores = 1
)
```

## Arguments

- formula:

  A formula specifying the outcome and covariates (e.g., `Y ~ X1 + X2`
  or `Y ~ .`). Use `~ . - Z` to exclude variables. Required if `Y` and
  `X` are not provided.

- data:

  A data.frame or data.table containing the variables. Required for the
  formula interface or when `Y`/`X` are passed as column name(s).

- Y:

  Numeric vector of outcomes, or a single string naming a column in
  `data`. Use this with `X` as an alternative to the formula interface.
  Can contain NA values for observations where `train_idx = FALSE`.

- X:

  Matrix, data.frame, or character vector of column names in `data`. Use
  this with `Y` as an alternative to the formula interface. When a
  character vector is supplied (e.g., `X = c("age", "income")` or
  `X = microcredit_covariates`), the corresponding columns are extracted
  from `data`.

- train_idx:

  Optional logical or integer vector indicating which observations to
  use for training. If `NULL` (default), all observations are used. If
  provided:

  - Logical vector: `TRUE` for training observations

  - Integer vector: indices of training observations

  Y must not have NA values for training observations. See **Training on
  a Subset** section for details on how this affects cross-fitting and
  predictions.

- M:

  Integer. Number of sample splitting repetitions (default: 2). Higher
  values improve stability but increase computation time.

- K:

  Integer. Number of cross-fitting folds within each repetition
  (default: 3). Each observation appears in exactly one test fold per
  repetition.

- algorithms:

  Character vector of ML algorithms to include in the ensemble. Default
  is `c("lm", "grf")`. Algorithms can come from two sources:

  - **grf package**: Use `"grf"` for generalized random forest (via
    [`grf::regression_forest`](https://rdrr.io/pkg/grf/man/regression_forest.html)
    or
    [`grf::probability_forest`](https://rdrr.io/pkg/grf/man/probability_forest.html))

  - **mlr3 learners**: Any algorithm available in mlr3 or its
    extensions. Specify just the algorithm name without the task prefix
    (e.g., use `"ranger"` not `"regr.ranger"`). The function will
    automatically add the appropriate prefix based on the task type.
    Common examples include:

    - `"lm"`: Linear regression

    - `"ranger"`: Random forest

    - `"glmnet"`: Elastic net regularization

    - `"xgboost"`: Gradient boosting

    - `"nnet"`: Neural network

    - `"kknn"`: K-nearest neighbors

    - `"svm"`: Support vector machine

    To see all available learners, run `mlr3::mlr_learners$keys()`.
    Additional learners may require installing `mlr3learners` or
    `mlr3extralearners` packages.

- ensemble_folds:

  Integer. Number of folds for cross-validated ensemble weight
  estimation (default: 5).

- task_type:

  Character. Type of prediction task: `"regr"` for continuous outcomes
  or `"classif"` for binary outcomes. If `NULL` (default), automatically
  detected from the outcome.

- scale_covariates:

  Logical. Whether to standardize non-binary numeric covariates to mean
  0 and standard deviation 1 before ML training (default: `TRUE`).
  Binary variables (0/1) are not scaled. The original data is preserved
  in the returned object.

- tune:

  Logical. Whether to perform hyperparameter tuning for ML algorithms
  (default: `FALSE`). When `TRUE`, uses random search with early
  stopping.

- tune_params:

  List of tuning parameters. Supports two modes:

  **Simple mode** (default)

  :   A list of scalar parameters that configure the built-in
      auto-tuner:

      - `time`: Maximum tuning time in seconds (default: 30)

      - `cv_folds`: Number of CV folds for tuning (default: 3)

      - `stagnation_iters`: Stop if no improvement for this many
        iterations (default: 250)

      - `stagnation_threshold`: Minimum improvement threshold (default:
        0.01)

      - `measure`: Performance measure string (default: R² for
        regression, AUC for classification)

  **Advanced mode**

  :   Pass mlr3tuning objects directly for full control:

      - `tuner`: A `Tuner` object (e.g.,
        `mlr3tuning::tnr("grid_search")`)

      - `terminator`: A `Terminator` object (e.g.,
        `mlr3tuning::trm("evals", n_evals = 50)`)

      - `resampling`: A `Resampling` object (e.g.,
        `mlr3::rsmp("holdout")`)

      - `search_space`: A `ParamSet` or tuning space object

      - `measure`: A `Measure` object or string

- learner_params:

  Optional named list of algorithm-specific parameters for mlr3
  learners. Each element name should match an algorithm in `algorithms`,
  and the value should be a list of parameter-value pairs. This only
  affects **mlr3-based** algorithms; it is ignored for `"grf"` (which
  uses its own internal defaults).

  Example:

      learner_params = list(
          ranger = list(num.trees = 1000, min.node.size = 5),
          glmnet = list(alpha = 0),        # ridge regression
          xgboost = list(nrounds = 200, max_depth = 6)
        )

  Parameters are applied after algorithm-specific defaults and override
  them when there is a conflict. To see which parameters are available
  for a given learner, run `mlr3::lrn("regr.<algorithm>")$param_set`.

- ensemble_strategy:

  Character. Strategy for combining algorithm predictions:

  - `"cv"` (default): Cross-validated OLS regression of Y on algorithm
    predictions. Learns optimal weights for each algorithm.

  - `"average"`: Simple average of all algorithm predictions. No weight
    learning; all algorithms contribute equally. Works with a single
    algorithm.

  See **Ensemble Strategy** section for more details.

- individual_id:

  Required when the dataset is a panel (e.g., individuals observed over
  multiple time periods). Specifies the column that identifies
  individuals so that (1) all observations for the same individual are
  placed in the same cross-fitting fold, and (2) cluster-robust standard
  errors are used in all downstream analyses.

  Example: for a panel of students observed across semesters, set
  `individual_id = student_id`.

  Can be an unquoted column name, a quoted string (`"student_id"`), or a
  vector of identifiers.

- n_cores:

  Integer. Number of cores for parallel processing of repetitions.
  Default is 1 (sequential). Set to higher values to parallelize the M
  repetitions. Uses the `future` framework, so users can also set up
  their own parallel backend via
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  before calling this function.

## Value

An object of class `ensemble_pred_fit` containing:

- predictions:

  data.table of predictions with M columns (one per repetition)

- call:

  The matched function call

- formula:

  The formula used (or constructed from Y/X)

- data:

  The original data (or constructed data.table from Y/X)

- Y:

  Vector of outcomes

- X:

  data.table of covariates (unscaled)

- train_idx:

  Logical vector indicating training observations

- splits:

  List of fold assignments for each repetition

- n:

  Number of observations

- n_train:

  Number of training observations

- M, K:

  Number of repetitions and folds

- algorithms:

  Algorithms used in ensemble

- ensemble_folds:

  Number of ensemble CV folds

- ensemble_strategy:

  Strategy used for combining predictions ("cv" or "average")

- task_type:

  Task type (regr or classif)

- scale_covariates:

  Whether covariates were scaled

- tune, tune_params:

  Tuning settings

- individual_id:

  Vector of individual identifiers (if panel data)

- n_cores:

  Number of cores used for parallel processing

## Cross-Fitting Procedure

The estimation proceeds as follows:

1.  The data is randomly split into \\K\\ folds. This random splitting
    is repeated \\M\\ times (each time with a fresh random partition).

2.  For each repetition \\m = 1, \ldots, M\\ and each fold \\k = 1,
    \ldots, K\\:

    - Each ML algorithm in `algorithms` is trained on the \\K - 1\\
      folds that exclude fold \\k\\.

    - Out-of-sample predictions of Y are generated for all observations
      in fold \\k\\.

3.  The per-algorithm predictions are combined into a single ensemble
    prediction using the `ensemble_strategy` (cross-validated OLS or
    simple average).

4.  This produces one complete vector of out-of-sample predictions per
    repetition (each observation appears in exactly one test fold per
    repetition).

The resulting \\M\\ prediction vectors are stored and used by downstream
analysis functions
([`blp_pred`](https://bfava.com/ensembleHTE/reference/blp_pred.md),
[`gavs`](https://bfava.com/ensembleHTE/reference/gavs.md),
[`gates`](https://bfava.com/ensembleHTE/reference/gates.md),
[`clan`](https://bfava.com/ensembleHTE/reference/clan.md)), which
compute their estimands separately for each repetition and then average
the estimates and standard errors across the \\M\\ repetitions.

## Training on a Subset

The `train_idx` parameter allows training models on a subset of
observations while generating predictions for all observations. This is
useful when the outcome Y is only observed for some units (e.g., only
treated units in an experiment) but you want predictions for everyone.

When `train_idx` is provided:

- Y can have NA values for observations where `train_idx = FALSE`

- Cross-fitting splits ALL observations into K folds, stratifying by
  `train_idx`

- For each fold k, models are trained on training observations NOT in
  fold k

- Predictions are generated for ALL observations in fold k (both
  training and non-training)

- Each observation gets exactly one prediction per repetition (from the
  model where they were in the test fold)

- Ensemble weights are estimated using only training observations

- Summary statistics are computed only on training observations

## Ensemble Strategy

The ensemble combines predictions from multiple ML algorithms. The
`ensemble_strategy` parameter controls how predictions are combined:

- `"cv"` (default): Uses cross-validated OLS regression of Y on the
  predicted values from each algorithm to learn optimal combination
  weights.

- `"average"`: Uses simple averaging across all algorithm predictions
  with equal weights. This is useful when you want a simpler ensemble or
  when using only a single algorithm.

## References

Fava, B. (2025). Training and Testing with Multiple Splits: A Central
Limit Theorem for Split-Sample Estimators. *arXiv preprint
arXiv:2511.04957*.

## See also

[`gavs`](https://bfava.com/ensembleHTE/reference/gavs.md) for Group
Averages analysis,
[`blp_pred`](https://bfava.com/ensembleHTE/reference/blp_pred.md) for
Best Linear Predictor analysis,
[`clan`](https://bfava.com/ensembleHTE/reference/clan.md) for
Classification Analysis,
[`ensemble_hte`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md)
for heterogeneous treatment effect estimation

## Examples

``` r
# \donttest{
# --- Predict bank profits using the microcredit data ---
# Bank profits are only observed for borrowers (loan_size > 0 & treat == 1).
# Use train_idx to train on borrowers, predict for the full sample.
# (Athey, Fava, Karlan, Osman & Zinman, 2025)
data(microcredit)

# Subset of covariates for speed (full set: microcredit_covariates)
covars <- c("age", "gender", "education", "hhinc_yrly_base",
            "css_creditscorefinal")

# Include hhinc_yrly_end and treat for downstream gates() analysis
dat <- microcredit[, c("bank_profits_pp", "hhinc_yrly_end", "treat",
                       "prop_score", covars)]
f <- as.formula(paste("bank_profits_pp ~", paste(covars, collapse = " + ")))

fit <- ensemble_pred(
  f, data = dat,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
print(fit)
#> Ensemble Prediction Fit
#> =======================
#> 
#> Call:
#> ensemble_pred(formula = f, data = dat, train_idx = microcredit$loan_size > 
#>     0 & microcredit$treat == 1, M = 3, K = 3, algorithms = c("lm", 
#>     "grf"))
#> 
#> Data:
#>   Observations:      1113
#>   Training obs:      650 (subset)
#>   Outcome:           bank_profits_pp
#>   Covariates:        5
#> 
#> Model specification:
#>   Algorithms:        lm, grf
#>   Task type:         regression (continuous outcome)
#> 
#> Split-sample parameters:
#>   Repetitions (M):   3
#>   Folds (K):         3
#>   Ensemble strategy: cross-validated OLS
#>   Ensemble folds:    5
#>   Covariate scaling: enabled
#>   Hyperparameter tuning: disabled
summary(fit)
#> Ensemble Prediction Summary
#> ===========================
#> 
#> Call:
#> ensemble_pred(formula = f, data = dat, train_idx = microcredit$loan_size > 
#>     0 & microcredit$treat == 1, M = 3, K = 3, algorithms = c("lm", 
#>     "grf"))
#> 
#> Outcome:     bank_profits_pp
#> Observations: 1113
#> Training obs: 650
#> Repetitions:  3
#> 
#> Prediction Accuracy (averaged across 3 repetitions):
#>   R-squared:         0.2
#>   RMSE:              15.33
#>   MAE:               11.03
#>   Correlation:       0.44
#> 
#> Best Linear Predictor (BLP):
#>   intercept:         -0.01 (SE: 0.60, p: 0.993) 
#>   slope:             0.97 (SE: 0.07, p: 0.000) ***
#>   -> Intercept close to 0 and slope close to 1 indicate good calibration
#> Note: ML model trained on 650 of 1113 observations; GAVS evaluated on 650 observations. Use subset = "all" to use all observations.
#> 
#> Group Averages (GAVS) with 3 groups:
#>   Group    Estimate   Std.Error    Pr(>|t|)
#>   --------------------------------------------
#>       1       -9.05        1.32       0.000 ***
#>       2       -0.81        0.91       0.370 
#>       3        7.64        1.03       0.000 ***
#> 
#>   Top - Bottom:  16.68 (SE: 1.67, p: 0.000) ***
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Can we predict bank profits? (blp_pred and gavs use training obs by default)
blp_pred(fit)
#> BLP Results (Best Linear Predictor - Prediction)
#> =================================================
#> 
#> Fit type: Prediction (ensemble_pred)
#> Outcome analyzed: bank_profits_pp
#> Observations used: 650
#> Repetitions: 3
#> 
#> Coefficients:
#>   intercept: Regression intercept (0 = well-calibrated)
#>   beta: Prediction loading (1 = well-calibrated)
#> 
#>         Term    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>    intercept       -0.01        0.60     -0.01       0.993 
#>         beta        0.97        0.07     13.46       0.000 ***
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
gavs(fit, n_groups = 3)
#> Note: ML model trained on 650 of 1113 observations; GAVS evaluated on 650 observations. Use subset = "all" to use all observations.
#> GAVS Results (Group Averages)
#> =============================
#> 
#> Outcome analyzed: bank_profits_pp
#> Number of groups: 3
#> Repetitions: 3
#> 
#> Data usage:
#>   ML trained on:         650 of 1113 obs
#>   Analysis evaluated on: 650 of 1113 obs
#>   Groups (3) formed on: all observations (1113 obs)
#> 
#> Group Average Outcomes (groups by predicted Y):
#> 
#>   Group    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>       1       -9.05        1.32     -6.85       0.000 ***
#>       2       -0.81        0.91     -0.90       0.370 
#>       3        7.64        1.03      7.43       0.000 ***
#> 
#> Heterogeneity Tests:
#>   ----------------------------------------------------
#>           Test    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>     Top-Bottom       16.68        1.67      9.96       0.000 ***
#>        Top-All        7.64        0.81      9.40       0.000 ***
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Do predicted profits correlate with treatment effects on household income?
gates(fit, outcome = "hhinc_yrly_end", treatment = "treat",
      prop_score = dat$prop_score, subset = "all", n_groups = 3)
#> Note: ML model trained on 650 of 1113 observations; GATES evaluated on 1113 observations. Use subset = "train" to restrict to training observations.
#> GATES Results
#> =============
#> 
#> Fit type: Prediction (ensemble_pred)
#> Outcome analyzed: hhinc_yrly_end
#>   (Groups based on: bank_profits_pp predicted Y)
#> Number of groups: 3
#> Repetitions: 3
#> 
#> Data usage:
#>   ML trained on:         650 of 1113 obs
#>   Analysis evaluated on: 1113 of 1113 obs
#>   Groups (3) formed on: all observations (1113 obs)
#> 
#> Group Average Treatment Effects:
#> 
#>   Group    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>       1      240.98     2212.75      0.11       0.913 
#>       2     -142.88     1560.51     -0.09       0.927 
#>       3     4554.22     4499.98      1.01       0.312 
#> 
#> Heterogeneity Tests:
#>   ----------------------------------------------------
#>           Test    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>     Top-Bottom     4313.24     5022.80      0.86       0.390 
#>        Top-All     3003.10     3134.93      0.96       0.338 
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
if (FALSE) { # \dontrun{
# --- Additional interface examples ---

# Matrix interface
fit <- ensemble_pred(
  Y = microcredit$bank_profits_pp,
  X = microcredit[, covars],
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "ranger"), M = 5, K = 5
)

# Using all covariates from Athey et al. (2025)
dat_full <- microcredit[, c("bank_profits_pp", microcredit_covariates)]
fit_full <- ensemble_pred(
  bank_profits_pp ~ ., data = dat_full,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "grf"), M = 5, K = 5
)

# Column-name interface (equivalent, no need to subset data)
fit_cov <- ensemble_pred(
  Y = "bank_profits_pp",
  X = microcredit_covariates,
  data = microcredit,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("lm", "grf"), M = 5, K = 5
)

# With parallel processing (4 cores)
fit <- ensemble_pred(
  bank_profits_pp ~ ., data = dat,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  M = 10, K = 5, n_cores = 4
)

# With algorithm-specific learner parameters (mlr3 algorithms only)
fit <- ensemble_pred(
  bank_profits_pp ~ ., data = dat,
  train_idx = microcredit$loan_size > 0 & microcredit$treat == 1,
  algorithms = c("ranger", "glmnet", "lm"),
  learner_params = list(
    ranger = list(num.trees = 1000, min.node.size = 5),
    glmnet = list(alpha = 0)   # ridge regression
  )
)

# Panel data: individuals observed across multiple time periods
panel_data <- data.frame(
  id = rep(1:100, each = 3),
  time = rep(1:3, 100),
  Y = rnorm(300),
  X1 = rnorm(300),
  X2 = rnorm(300)
)
fit_panel <- ensemble_pred(
  formula = Y ~ X1 + X2,
  data = panel_data,
  individual_id = id,
  algorithms = c("lm", "grf"),
  M = 5, K = 3
)
gavs(fit_panel, n_groups = 3)
blp_pred(fit_panel)
} # }
```
