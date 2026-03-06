# Fit Ensemble Heterogeneous Treatment Effects Model

Estimates heterogeneous treatment effects (HTEs) using an ensemble of
machine learning algorithms combined with multiple sample splitting.
This function implements the estimation strategy developed by Fava
(2025), which improves statistical power by averaging predictions across
M repetitions of K-fold cross-fitting.

By default, the function uses the R-learner metalearner strategy with
generalized random forest (`grf`) as the CATE estimator.

The function supports two interfaces:

- **Formula interface**: Specify `formula`, `treatment`, and `data`

- **Matrix interface**: Specify `Y`, `X`, and `D` directly (as
  vectors/matrices, or as column names to be looked up in `data`)

## Usage

``` r
ensemble_hte(
  formula = NULL,
  treatment = NULL,
  data = NULL,
  Y = NULL,
  X = NULL,
  D = NULL,
  prop_score = NULL,
  M = 2,
  K = 3,
  algorithms = c("lm", "grf"),
  metalearner = c("r", "t", "s", "x"),
  r_learner = "grf",
  ensemble_folds = 5,
  task_type = NULL,
  scale_covariates = TRUE,
  tune = FALSE,
  tune_params = list(time = 30, cv_folds = 3, stagnation_iters = 250,
    stagnation_threshold = 0.01, measure = NULL),
  learner_params = NULL,
  train_idx = NULL,
  ensemble_strategy = c("cv", "average"),
  individual_id = NULL,
  n_cores = 1
)
```

## Arguments

- formula:

  A formula specifying the outcome and covariates (e.g., `Y ~ X1 + X2`
  or `Y ~ .`). Use `~ . - Z` to exclude variables. Required if `Y`, `X`,
  `D` are not provided.

- treatment:

  The treatment variable. Must be coded as 0 (control) and 1 (treated).
  Can be specified as:

  - An unquoted variable name: `treatment = D`

  - A quoted string: `treatment = "D"`

  - A variable containing the column name:
    `treat_col <- "D"; treatment = treat_col`

  - Ignored when using matrix interface (use `D` parameter instead)

- data:

  A data.frame or data.table containing the variables referenced in the
  formula, or in `Y`, `X`, `D` when those are given as column names.

- Y:

  Numeric vector of outcomes, or a string with the column name in
  `data`. Use this with `X` and `D` as an alternative to the formula
  interface.

- X:

  Matrix, data.frame, or character vector of column names in `data`.
  When a character vector is provided, the columns are extracted from
  `data`. Use this with `Y` and `D` as an alternative to the formula
  interface.

  Example: `X = c("age", "gender", "income")` or
  `X = microcredit_covariates`.

- D:

  Numeric vector of treatment indicators (0/1), or a string with the
  column name in `data`. Must contain only 0s and 1s.

- prop_score:

  Numeric vector of propensity scores (probability of treatment given
  covariates), a single string naming a column in `data`, or a scalar
  constant. Must be strictly between 0 and 1 (exclusive). If `NULL`
  (default), assumes constant propensity equal to the sample treatment
  proportion (appropriate for randomized experiments).

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

- metalearner:

  Character (single value). The metalearner strategy for ITE estimation.
  Exactly one of:

  - `"r"` (default): R-learner with Robinson transformation

  - `"t"`: T-learner with separate models per treatment arm

  - `"s"`: S-learner with treatment as a feature

  - `"x"`: X-learner with imputed counterfactuals

  Only one metalearner can be used per call. See **Metalearners**
  section below for detailed descriptions.

- r_learner:

  Character. When `metalearner = "r"`, specifies the algorithm for
  estimating the conditional average treatment effect (CATE) in the
  final stage. Default is `"grf"`
  ([`grf::causal_forest`](https://rdrr.io/pkg/grf/man/causal_forest.html)).
  Can be:

  - `"grf"`: Uses
    [`grf::causal_forest`](https://rdrr.io/pkg/grf/man/causal_forest.html)
    (recommended)

  - Any mlr3 learner: e.g., `"ranger"`, `"xgboost"`, `"glmnet"`

  This does not need to be in the `algorithms` list. Only used when
  `metalearner = "r"`.

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

- train_idx:

  Optional logical or integer vector indicating which observations to
  use for training. If `NULL` (default), all observations are used. If
  provided:

  - Logical vector: `TRUE` for training observations

  - Integer vector: indices of training observations

  This is useful for multi-arm trials where you want to fit HTE using
  only one treatment-control pair but generate ITE predictions for all
  units. When `train_idx` is provided:

  - Cross-fitting splits ALL observations into K folds, stratifying by
    `train_idx`

  - Models are trained only on training observations

  - ITE predictions are generated for ALL observations in each test fold

  - Ensemble weights are estimated using only training observations

- ensemble_strategy:

  Character. Strategy for combining algorithm predictions. One of `"cv"`
  (default) or `"average"`. `"cv"` uses cross-validated BLP regression
  to learn optimal weights; `"average"` uses a simple unweighted average
  of algorithm predictions. See **Ensemble Strategy** section for
  details.

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

An object of class `ensemble_hte_fit` containing:

- ite:

  data.table of ITE predictions with M columns (one per repetition)

- call:

  The matched function call

- formula:

  The formula used (or constructed from Y/X/D)

- treatment:

  Name of the treatment variable

- data:

  The original data (or constructed data.table from Y/X/D)

- Y:

  Vector of outcomes

- X:

  data.table of covariates (unscaled)

- D:

  Vector of treatment indicators

- prop_score:

  Vector of propensity scores

- weights:

  Inverse propensity weights

- splits:

  List of fold assignments for each repetition

- n:

  Number of observations

- n_train:

  Number of training observations

- train_idx:

  Logical vector indicating training observations

- M, K:

  Number of repetitions and folds

- algorithms:

  Algorithms used in ensemble

- metalearner:

  Metalearner strategy used

- r_learner:

  R-learner algorithm (if applicable)

- ensemble_folds:

  Number of ensemble CV folds

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
      folds that exclude fold \\k\\, using the chosen `metalearner`
      strategy.

    - Out-of-sample ITE predictions are generated for all observations
      in fold \\k\\.

3.  The per-algorithm ITE predictions are combined into a single
    ensemble prediction using the `ensemble_strategy` (cross-validated
    BLP or simple average).

4.  This produces one complete vector of out-of-sample ITE predictions
    per repetition (each observation appears in exactly one test fold
    per repetition).

The resulting \\M\\ vectors of ITE predictions are stored and used by
the downstream analysis functions
([`blp`](https://bfava.com/ensembleHTE/reference/blp.md),
[`gates`](https://bfava.com/ensembleHTE/reference/gates.md),
[`clan`](https://bfava.com/ensembleHTE/reference/clan.md),
[`gavs`](https://bfava.com/ensembleHTE/reference/gavs.md)), which
compute their estimands separately for each repetition and then average
the estimates and standard errors across the \\M\\ repetitions.

## Metalearners

The function supports four metalearner strategies for estimating
individual treatment effects (ITEs):

- **R-learner** (default): Robinson transformation with
  residual-on-residual regression. Uses
  [`grf::causal_forest`](https://rdrr.io/pkg/grf/man/causal_forest.html)
  by default for the final CATE model.

- **T-learner**: Trains separate models for treated and control groups

- **S-learner**: Trains a single model with treatment as a feature

- **X-learner**: Two-stage approach that imputes counterfactual outcomes

See Nie & Wager (2021) for R-learner and Künzel et al. (2019) for
T/S/X-learners.

## Ensemble Strategy

The ensemble combines predictions from multiple ML algorithms into a
single ITE estimate per observation. Two strategies are available:

- `"cv"`:

  (Default) Uses a cross-validated Best Linear Predictor (BLP)
  regression to learn optimal weights for each algorithm. Weights are
  derived from a weighted least squares regression of outcomes on
  algorithm predictions, using only training observations. This is the
  recommended approach and the one described in the paper.

- `"average"`:

  Combines algorithm predictions using a simple (unweighted) average.
  This is faster and more robust with small samples or few algorithms,
  but does not adapt weights to algorithm performance.

## References

Fava, B. (2025). Training and Testing with Multiple Splits: A Central
Limit Theorem for Split-Sample Estimators. *arXiv preprint
arXiv:2511.04957*.

Nie, X., & Wager, S. (2021). Quasi-Oracle Estimation of Heterogeneous
Treatment Effects. *Biometrika*, 108(2), 299-319.

Künzel, S.R., Sekhon, J.S., Bickel, P.J., & Yu, B. (2019). Metalearners
for estimating heterogeneous treatment effects using machine learning.
*Proceedings of the National Academy of Sciences*, 116(10), 4156-4165.

## See also

[`gates`](https://bfava.com/ensembleHTE/reference/gates.md) for Group
Average Treatment Effects analysis,
[`blp`](https://bfava.com/ensembleHTE/reference/blp.md) for Best Linear
Predictor analysis,
[`clan`](https://bfava.com/ensembleHTE/reference/clan.md) for
Classification Analysis,
[`ensemble_pred`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md)
for standard prediction without treatment effects

## Examples

``` r
# \donttest{
# --- HTE estimation on the Philippine microcredit experiment ---
# Outcome: household income; Treatment: microloan offer
data(microcredit)

# Subset of covariates for speed (full set: object microcredit_covariates)
covars <- c("age", "gender", "education", "hhinc_yrly_base",
            "css_creditscorefinal")
dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]

fit <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat,
  prop_score = microcredit$prop_score,
  algorithms = c("lm", "grf"), M = 3, K = 3
)
#> Warning: Some propensity scores are below 0.20 or above 0.80. This package is designed for randomized controlled trials (RCTs), where propensity scores are typically well-balanced. Extreme propensity scores may indicate an observational study or a heavily unbalanced design. Please verify your experimental design.
print(fit)
#> Ensemble HTE Fit
#> ================
#> 
#> Call:
#> ensemble_hte(formula = hhinc_yrly_end ~ ., treatment = treat, 
#>     data = dat, prop_score = microcredit$prop_score, M = 3, K = 3, 
#>     algorithms = c("lm", "grf"))
#> 
#> Data:
#>   Observations:      1113
#>   Targeted outcome:  hhinc_yrly_end
#>   Treatment:         treat
#>   Covariates:        5
#> 
#> Model specification:
#>   Algorithms:        lm, grf
#>   Metalearner:       R-learner
#>   Task type:         regression (continuous outcome)
#>   R-learner method:  grf
#> 
#> Split-sample parameters:
#>   Repetitions (M):   3
#>   Folds (K):         3
#>   Ensemble strategy: cross-validated BLP
#>   Ensemble folds:    5
#>   Covariate scaling: enabled
#>   Hyperparameter tuning: disabled
summary(fit)
#> Ensemble HTE Summary
#> ====================
#> 
#> Call:
#> ensemble_hte(formula = hhinc_yrly_end ~ ., treatment = treat, 
#>     data = dat, prop_score = microcredit$prop_score, M = 3, K = 3, 
#>     algorithms = c("lm", "grf"))
#> 
#> Outcome:     hhinc_yrly_end
#> Treatment:   treat
#> Observations: 1113
#> Repetitions:  3
#> 
#> Best Linear Predictor (BLP):
#>   beta1 (ATE):       1664.22 (SE: 1683.04, p: 0.323) 
#>   beta2 (HET):       -0.63 (SE: 0.73, p: 0.386) 
#>   -> No significant heterogeneity detected (p >= 0.05)
#> 
#> Group Average Treatment Effects (GATES) with 3 groups:
#>   Group    Estimate   Std.Error    Pr(>|t|)
#>   --------------------------------------------
#>       1     2835.66     2565.05       0.269 
#>       2      831.96     1963.17       0.672 
#>       3     1125.47     4201.18       0.789 
#> 
#>   Top - Bottom:  -1710.19 (SE: 4959.58, p: 0.730) 
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Downstream analysis
gates(fit, n_groups = 3)
#> GATES Results
#> =============
#> 
#> Fit type: HTE (ensemble_hte)
#> Outcome analyzed: hhinc_yrly_end
#> Number of groups: 3
#> Repetitions: 3
#> 
#> Group Average Treatment Effects:
#> 
#>   Group    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>       1     2835.66     2565.05      1.11       0.269 
#>       2      831.96     1963.17      0.42       0.672 
#>       3     1125.47     4201.18      0.27       0.789 
#> 
#> Heterogeneity Tests:
#>   ----------------------------------------------------
#>           Test    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>     Top-Bottom    -1710.19     4959.58     -0.34       0.730 
#>        Top-All     -474.29     3008.83     -0.16       0.875 
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
blp(fit)
#> BLP Results (Best Linear Predictor of CATE)
#> ============================================
#> 
#> Fit type: HTE (ensemble_hte)
#> Outcome analyzed: hhinc_yrly_end
#> Repetitions: 3
#> 
#> Coefficients:
#>   beta1 (ATE): Average Treatment Effect
#>   beta2 (HET): Heterogeneity loading (significant = ML captures heterogeneity)
#> 
#>     Term    Estimate   Std.Error   t value    Pr(>|t|)
#>   ----------------------------------------------------
#>    beta1     1664.22     1683.04      0.99       0.323 
#>    beta2       -0.63        0.73     -0.87       0.386 
#> 
#> ---
#> Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# }
if (FALSE) { # \dontrun{
# --- Additional interface examples ---

# Matrix interface
covars <- c("age", "gender", "education", "hhinc_yrly_base",
            "css_creditscorefinal")
fit <- ensemble_hte(
  Y = microcredit$hhinc_yrly_end,
  X = microcredit[, covars],
  D = microcredit$treat,
  prop_score = microcredit$prop_score,
  algorithms = c("lm", "ranger"), M = 5, K = 5
)

# Using all covariates from the original paper
dat_full <- microcredit[, c("hhinc_yrly_end", "treat", microcredit_covariates)]
fit_full <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat_full,
  prop_score = microcredit$prop_score,
  algorithms = c("lm", "grf"), M = 5, K = 5
)

# Column-name interface (equivalent, no need to subset data)
fit_names <- ensemble_hte(
  Y = "hhinc_yrly_end",
  X = microcredit_covariates,
  D = "treat",
  data = microcredit,
  prop_score = "prop_score",
  algorithms = c("lm", "grf"), M = 5, K = 5
)

# With propensity scores and X-learner
fit <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat,
  prop_score = microcredit$prop_score,
  metalearner = "x"
)

# With parallel processing (4 cores)
fit <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat,
  prop_score = microcredit$prop_score,
  M = 10, K = 5, n_cores = 4
)

# With algorithm-specific learner parameters (mlr3 algorithms only)
fit <- ensemble_hte(
  hhinc_yrly_end ~ ., treatment = treat, data = dat,
  algorithms = c("ranger", "glmnet", "lm"),
  learner_params = list(
    ranger = list(num.trees = 1000, min.node.size = 5),
    glmnet = list(alpha = 0)   # ridge regression
  )
)

# Panel data: individuals observed across multiple time periods
# All observations for the same individual are kept in the same fold,
# and downstream analyses use cluster-robust standard errors.
panel_data <- data.frame(
  id = rep(1:100, each = 3),
  time = rep(1:3, 100),
  Y = rnorm(300),
  D = rbinom(300, 1, 0.5),
  X1 = rnorm(300),
  X2 = rnorm(300)
)
fit_panel <- ensemble_hte(
  formula = Y ~ X1 + X2,
  treatment = D,
  data = panel_data,
  individual_id = id,
  algorithms = c("lm", "grf"),
  M = 5, K = 3
)
gates(fit_panel, n_groups = 3)
blp(fit_panel)
clan(fit_panel, c("X1", "X2"))
} # }
```
