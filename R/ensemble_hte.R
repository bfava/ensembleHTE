#' Fit Ensemble Heterogeneous Treatment Effects Model
#'
#' @description
#' Estimates heterogeneous treatment effects (HTEs) using an ensemble of machine
#' learning algorithms combined with multiple sample splitting. This function
#' implements the estimation strategy developed by Fava (2025), which improves
#' statistical power by averaging predictions across M repetitions of K-fold
#' cross-fitting.
#'
#' By default, the function uses the R-learner metalearner strategy with
#' generalized random forest (\code{grf}) as the CATE estimator.
#'
#' The function supports two interfaces:
#' \itemize{
#'   \item \strong{Formula interface}: Specify \code{formula}, \code{treatment}, and \code{data}
#'   \item \strong{Matrix interface}: Specify \code{Y}, \code{X}, and \code{D} directly
#'     (as vectors/matrices, or as column names to be looked up in \code{data})
#' }
#'
#' @section Cross-Fitting Procedure:
#' The estimation proceeds as follows:
#' \enumerate{
#'   \item The data is randomly split into \eqn{K} folds. This random splitting
#'     is repeated \eqn{M} times (each time with a fresh random partition).
#'   \item For each repetition \eqn{m = 1, \ldots, M} and each fold
#'     \eqn{k = 1, \ldots, K}:
#'     \itemize{
#'       \item Each ML algorithm in \code{algorithms} is trained on the
#'         \eqn{K - 1} folds that exclude fold \eqn{k}, using the chosen
#'         \code{metalearner} strategy.
#'       \item Out-of-sample ITE predictions are generated for all observations
#'         in fold \eqn{k}.
#'     }
#'   \item The per-algorithm ITE predictions are combined into a single ensemble
#'     prediction using the \code{ensemble_strategy} (cross-validated BLP or
#'     simple average).
#'   \item This produces one complete vector of out-of-sample ITE predictions per
#'     repetition (each observation appears in exactly one test fold per repetition).
#' }
#' The resulting \eqn{M} vectors of ITE predictions are stored and used by the
#' downstream analysis functions (\code{\link{blp}}, \code{\link{gates}},
#' \code{\link{clan}}, \code{\link{gavs}}), which compute their estimands
#' separately for each repetition and then average the estimates and standard
#' errors across the \eqn{M} repetitions.
#'
#' @section Metalearners:
#' The function supports four metalearner strategies for estimating individual
#' treatment effects (ITEs):
#' \itemize{
#'   \item \strong{R-learner} (default): Robinson transformation with residual-on-residual
#'     regression. Uses \code{grf::causal_forest} by default for the final CATE model.
#'   \item \strong{T-learner}: Trains separate models for treated and control groups
#'   \item \strong{S-learner}: Trains a single model with treatment as a feature
#'   \item \strong{X-learner}: Two-stage approach that imputes counterfactual outcomes
#' }
#' See Nie & Wager (2021) for R-learner and Künzel et al. (2019) for T/S/X-learners.
#'
#' @section Ensemble Strategy:
#' The ensemble combines predictions from multiple ML algorithms into a single
#' ITE estimate per observation. Two strategies are available:
#' \describe{
#'   \item{\code{"cv"}}{(Default) Uses a cross-validated Best Linear
#'     Predictor (BLP) regression to learn optimal weights for each algorithm.
#'     Weights are derived from a weighted least squares regression of
#'     outcomes on algorithm predictions, using only training observations.
#'     This is the recommended approach and the one described in the paper.}
#'   \item{\code{"average"}}{Combines algorithm predictions using a simple
#'     (unweighted) average. This is faster and more robust with small samples
#'     or few algorithms, but does not adapt weights to algorithm performance.}
#' }
#'
#' @references
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#'
#' Nie, X., & Wager, S. (2021). Quasi-Oracle Estimation of Heterogeneous
#' Treatment Effects. \emph{Biometrika}, 108(2), 299-319.
#'
#' Künzel, S.R., Sekhon, J.S., Bickel, P.J., & Yu, B. (2019). Metalearners for
#' estimating heterogeneous treatment effects using machine learning.
#' \emph{Proceedings of the National Academy of Sciences}, 116(10), 4156-4165.
#'
#' @param formula A formula specifying the outcome and covariates (e.g.,
#'   \code{Y ~ X1 + X2} or \code{Y ~ .}). Use \code{~ . - Z} to exclude variables.
#'   Required if \code{Y}, \code{X}, \code{D} are not provided.
#' @param treatment The treatment variable. Must be coded as 0 (control) and 1
#'   (treated). Can be specified as:
#'   \itemize{
#'     \item An unquoted variable name: \code{treatment = D}
#'     \item A quoted string: \code{treatment = "D"}
#'     \item A variable containing the column name: \code{treat_col <- "D"; treatment = treat_col}
#'     \item Ignored when using matrix interface (use \code{D} parameter instead)
#'   }
#' @param data A data.frame or data.table containing the variables referenced in
#'   the formula, or in \code{Y}, \code{X}, \code{D} when those are given as
#'   column names.
#' @param Y Numeric vector of outcomes, or a string with the column name in
#'   \code{data}. Use this with \code{X} and \code{D} as an alternative to the
#'   formula interface.
#' @param X Matrix, data.frame, or character vector of column names in
#'   \code{data}. When a character vector is provided, the columns are extracted
#'   from \code{data}. Use this with \code{Y} and \code{D} as an alternative to
#'   the formula interface.
#'
#'   Example: \code{X = c("age", "gender", "income")} or
#'   \code{X = microcredit_covariates}.
#' @param D Numeric vector of treatment indicators (0/1), or a string with the
#'   column name in \code{data}. Must contain only 0s and 1s.
#' @param prop_score Numeric vector of propensity scores (probability of treatment
#'   given covariates), a single string naming a column in \code{data}, or a
#'   scalar constant. Must be strictly between 0 and 1 (exclusive). If \code{NULL}
#'   (default), assumes constant propensity equal to the sample treatment proportion
#'   (appropriate for randomized experiments).
#' @param M Integer. Number of sample splitting repetitions (default: 2). Higher
#'   values improve stability but increase computation time.
#' @param K Integer. Number of cross-fitting folds within each repetition
#'   (default: 3). Each observation appears in exactly one test fold per repetition.
#' @param algorithms Character vector of ML algorithms to include in the ensemble.
#'   Default is \code{c("lm", "grf")}. Algorithms can come from two sources:
#'   \itemize{
#'     \item \strong{grf package}: Use \code{"grf"} for generalized random forest
#'       (via \code{grf::regression_forest} or \code{grf::probability_forest})
#'     \item \strong{mlr3 learners}: Any algorithm available in mlr3 or its extensions.
#'       Specify just the algorithm name without the task prefix (e.g., use \code{"ranger"}
#'       not \code{"regr.ranger"}). The function will automatically add the appropriate
#'       prefix based on the task type. Common examples include:
#'       \itemize{
#'         \item \code{"lm"}: Linear regression
#'         \item \code{"ranger"}: Random forest
#'         \item \code{"glmnet"}: Elastic net regularization
#'         \item \code{"xgboost"}: Gradient boosting
#'         \item \code{"nnet"}: Neural network
#'         \item \code{"kknn"}: K-nearest neighbors
#'         \item \code{"svm"}: Support vector machine
#'       }
#'       To see all available learners, run \code{mlr3::mlr_learners$keys()}.
#'       Additional learners may require installing \code{mlr3learners} or
#'       \code{mlr3extralearners} packages.
#'   }
#' @param metalearner Character (single value). The metalearner strategy for ITE estimation.
#'   Exactly one of:
#'   \itemize{
#'     \item \code{"r"} (default): R-learner with Robinson transformation
#'     \item \code{"t"}: T-learner with separate models per treatment arm
#'     \item \code{"s"}: S-learner with treatment as a feature
#'     \item \code{"x"}: X-learner with imputed counterfactuals
#'   }
#'   Only one metalearner can be used per call.
#'   See \strong{Metalearners} section below for detailed descriptions.
#' @param r_learner Character. When \code{metalearner = "r"}, specifies the algorithm
#'   for estimating the conditional average treatment effect (CATE) in the final
#'   stage. Default is \code{"grf"} (\code{grf::causal_forest}). Can be:
#'   \itemize{
#'     \item \code{"grf"}: Uses \code{grf::causal_forest} (recommended)
#'     \item Any mlr3 learner: e.g., \code{"ranger"}, \code{"xgboost"}, \code{"glmnet"}
#'   }
#'   This does not need to be in the \code{algorithms} list. Only used when
#'   \code{metalearner = "r"}.
#' @param ensemble_folds Integer. Number of folds for cross-validated ensemble
#'   weight estimation (default: 5).
#' @param task_type Character. Type of prediction task: \code{"regr"} for
#'   continuous outcomes or \code{"classif"} for binary outcomes. If \code{NULL}
#'   (default), automatically detected from the outcome.
#' @param scale_covariates Logical. Whether to standardize non-binary numeric
#'   covariates to mean 0 and standard deviation 1 before ML training (default:
#'   \code{TRUE}). Binary variables (0/1) are not scaled. The original data is
#'   preserved in the returned object.
#' @param tune Logical. Whether to perform hyperparameter tuning for ML algorithms
#'   (default: \code{FALSE}). When \code{TRUE}, uses random search with early
#'   stopping.
#' @param tune_params List of tuning parameters. Supports two modes:
#'   \describe{
#'     \item{\strong{Simple mode} (default)}{A list of scalar parameters that
#'       configure the built-in auto-tuner:
#'       \itemize{
#'         \item \code{time}: Maximum tuning time in seconds (default: 30)
#'         \item \code{cv_folds}: Number of CV folds for tuning (default: 3)
#'         \item \code{stagnation_iters}: Stop if no improvement for this many iterations (default: 250)
#'         \item \code{stagnation_threshold}: Minimum improvement threshold (default: 0.01)
#'         \item \code{measure}: Performance measure string (default: R² for regression, AUC for classification)
#'       }
#'     }
#'     \item{\strong{Advanced mode}}{Pass mlr3tuning objects directly for full control:
#'       \itemize{
#'         \item \code{tuner}: A \code{Tuner} object (e.g., \code{mlr3tuning::tnr("grid_search")})
#'         \item \code{terminator}: A \code{Terminator} object (e.g., \code{mlr3tuning::trm("evals", n_evals = 50)})
#'         \item \code{resampling}: A \code{Resampling} object (e.g., \code{mlr3::rsmp("holdout")})
#'         \item \code{search_space}: A \code{ParamSet} or tuning space object
#'         \item \code{measure}: A \code{Measure} object or string
#'       }
#'     }
#'   }
#' @param learner_params Optional named list of algorithm-specific parameters
#'   for mlr3 learners. Each element name should match an algorithm in
#'   \code{algorithms}, and the value should be a list of parameter-value pairs.
#'   This only affects \strong{mlr3-based} algorithms; it is ignored for
#'   \code{"grf"} (which uses its own internal defaults).
#'
#'   Example:
#'   \preformatted{learner_params = list(
#'     ranger = list(num.trees = 1000, min.node.size = 5),
#'     glmnet = list(alpha = 0),        # ridge regression
#'     xgboost = list(nrounds = 200, max_depth = 6)
#'   )}
#'
#'   Parameters are applied after algorithm-specific defaults and override them
#'   when there is a conflict. To see which parameters are available for a given
#'   learner, run \code{mlr3::lrn("regr.<algorithm>")$param_set}.
#' @param train_idx Optional logical or integer vector indicating which observations
#'   to use for training. If \code{NULL} (default), all observations are used.
#'   If provided:
#'   \itemize{
#'     \item Logical vector: \code{TRUE} for training observations
#'     \item Integer vector: indices of training observations
#'   }
#'   This is useful for multi-arm trials where you want to fit HTE using only
#'   one treatment-control pair but generate ITE predictions for all units.
#'   When \code{train_idx} is provided:
#'   \itemize{
#'     \item Cross-fitting splits ALL observations into K folds, stratifying by \code{train_idx}
#'     \item Models are trained only on training observations
#'     \item ITE predictions are generated for ALL observations in each test fold
#'     \item Ensemble weights are estimated using only training observations
#'   }
#' @param ensemble_strategy Character. Strategy for combining algorithm predictions.
#'   One of \code{"cv"} (default) or \code{"average"}.
#'   \code{"cv"} uses cross-validated BLP regression to learn optimal weights;
#'   \code{"average"} uses a simple unweighted average of algorithm predictions.
#'   See \strong{Ensemble Strategy} section for details.
#' @param individual_id Required when the dataset is a panel (e.g., individuals
#'   observed over multiple time periods). Specifies the column that identifies
#'   individuals so that (1) all observations for the same individual are placed
#'   in the same cross-fitting fold, and (2) cluster-robust standard errors are
#'   used in all downstream analyses.
#'
#'   Example: for a panel of students observed across semesters, set
#'   \code{individual_id = student_id}.
#'
#'   Can be an unquoted column name, a quoted string (\code{"student_id"}),
#'   or a vector of identifiers.
#' @param n_cores Integer. Number of cores for parallel processing of repetitions.
#'   Default is 1 (sequential). Set to higher values to parallelize the M repetitions.
#'   Uses the \code{future} framework, so users can also set up their own parallel
#'   backend via \code{future::plan()} before calling this function.
#'
#' @return An object of class \code{ensemble_hte_fit} containing:
#' \describe{
#'   \item{ite}{data.table of ITE predictions with M columns (one per repetition)}
#'   \item{call}{The matched function call}
#'   \item{formula}{The formula used (or constructed from Y/X/D)}
#'   \item{treatment}{Name of the treatment variable}
#'   \item{data}{The original data (or constructed data.table from Y/X/D)}
#'   \item{Y}{Vector of outcomes}
#'   \item{X}{data.table of covariates (unscaled)}
#'   \item{D}{Vector of treatment indicators}
#'   \item{prop_score}{Vector of propensity scores}
#'   \item{weights}{Inverse propensity weights}
#'   \item{splits}{List of fold assignments for each repetition}
#'   \item{n}{Number of observations}
#'   \item{n_train}{Number of training observations}
#'   \item{train_idx}{Logical vector indicating training observations}
#'   \item{M, K}{Number of repetitions and folds}
#'   \item{algorithms}{Algorithms used in ensemble}
#'   \item{metalearner}{Metalearner strategy used}
#'   \item{r_learner}{R-learner algorithm (if applicable)}
#'   \item{ensemble_folds}{Number of ensemble CV folds}
#'   \item{task_type}{Task type (regr or classif)}
#'   \item{scale_covariates}{Whether covariates were scaled}
#'   \item{tune, tune_params}{Tuning settings}
#'   \item{individual_id}{Vector of individual identifiers (if panel data)}
#'   \item{n_cores}{Number of cores used for parallel processing}
#' }
#'
#' @examples
#' \donttest{
#' # --- HTE estimation on the Philippine microcredit experiment ---
#' # Outcome: household income; Treatment: microloan offer
#' data(microcredit)
#'
#' # Subset of covariates for speed (full set: object microcredit_covariates)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#'
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
#' print(fit)
#' summary(fit)
#'
#' # Downstream analysis
#' gates(fit, n_groups = 3)
#' blp(fit)
#' }
#' \dontrun{
#' # --- Additional interface examples ---
#'
#' # Matrix interface
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' fit <- ensemble_hte(
#'   Y = microcredit$hhinc_yrly_end,
#'   X = microcredit[, covars],
#'   D = microcredit$treat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "ranger"), M = 5, K = 5
#' )
#'
#' # Using all covariates from the original paper
#' dat_full <- microcredit[, c("hhinc_yrly_end", "treat", microcredit_covariates)]
#' fit_full <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat_full,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 5, K = 5
#' )
#'
#' # Column-name interface (equivalent, no need to subset data)
#' fit_names <- ensemble_hte(
#'   Y = "hhinc_yrly_end",
#'   X = microcredit_covariates,
#'   D = "treat",
#'   data = microcredit,
#'   prop_score = "prop_score",
#'   algorithms = c("lm", "grf"), M = 5, K = 5
#' )
#'
#' # With propensity scores and X-learner
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   metalearner = "x"
#' )
#'
#' # With parallel processing (4 cores)
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   M = 10, K = 5, n_cores = 4
#' )
#'
#' # With algorithm-specific learner parameters (mlr3 algorithms only)
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   algorithms = c("ranger", "glmnet", "lm"),
#'   learner_params = list(
#'     ranger = list(num.trees = 1000, min.node.size = 5),
#'     glmnet = list(alpha = 0)   # ridge regression
#'   )
#' )
#'
#' # Panel data: individuals observed across multiple time periods
#' # All observations for the same individual are kept in the same fold,
#' # and downstream analyses use cluster-robust standard errors.
#' panel_data <- data.frame(
#'   id = rep(1:100, each = 3),
#'   time = rep(1:3, 100),
#'   Y = rnorm(300),
#'   D = rbinom(300, 1, 0.5),
#'   X1 = rnorm(300),
#'   X2 = rnorm(300)
#' )
#' fit_panel <- ensemble_hte(
#'   formula = Y ~ X1 + X2,
#'   treatment = D,
#'   data = panel_data,
#'   individual_id = id,
#'   algorithms = c("lm", "grf"),
#'   M = 5, K = 3
#' )
#' gates(fit_panel, n_groups = 3)
#' blp(fit_panel)
#' clan(fit_panel, c("X1", "X2"))
#' }
#'
#' @seealso
#' \code{\link{gates}} for Group Average Treatment Effects analysis,
#' \code{\link{blp}} for Best Linear Predictor analysis,
#' \code{\link{clan}} for Classification Analysis,
#' \code{\link{ensemble_pred}} for standard prediction without treatment effects
#'
#' @export
ensemble_hte <- function(formula = NULL, treatment = NULL, data = NULL,
                         Y = NULL, X = NULL, D = NULL,
                         prop_score = NULL, M = 2, K = 3, 
                         algorithms = c("lm", "grf"), 
                         metalearner = c("r", "t", "s", "x"), 
                         r_learner = "grf", 
                         ensemble_folds = 5, 
                         task_type = NULL, 
                         scale_covariates = TRUE,
                         tune = FALSE, 
                         tune_params = list(
                           time = 30,
                           cv_folds = 3,
                           stagnation_iters = 250,
                           stagnation_threshold = 0.01,
                           measure = NULL
                         ),
                         learner_params = NULL,
                         train_idx = NULL,
                         ensemble_strategy = c("cv", "average"),
                         individual_id = NULL,
                         n_cores = 1) {
  
  # Capture the call
  cl <- match.call()
  
  # Validate metalearner argument and use first as default
  metalearner <- match.arg(metalearner)
  
  # Validate ensemble_strategy
  ensemble_strategy <- match.arg(ensemble_strategy)
  
  # Validate common inputs
  algorithms <- validate_common_inputs(M, K, algorithms, ensemble_folds, n_cores)
  
  # Resolve prop_score column name (if a single string referencing data)
  prop_score_var <- NULL
  if (!is.null(prop_score) && is.character(prop_score) && length(prop_score) == 1) {
    if (!is.null(data)) {
      if (prop_score %in% names(data)) {
        prop_score_var <- prop_score
        prop_score <- data[[prop_score_var]]
      } else {
        stop("Column '", prop_score, "' (passed as prop_score) not found in data.")
      }
    } else {
      stop("'data' must be provided when prop_score is a column name string.")
    }
  }
  
  # Validate prop_score if provided (after column-name resolution)
  if (!is.null(prop_score)) {
    if (!is.numeric(prop_score) || any(prop_score <= 0 | prop_score >= 1)) {
      stop("Propensity scores must be strictly between 0 and 1 (exclusive). ",
           "Values equal to 0 or 1 cause division by zero in inverse propensity weighting.")
    }
    if (any(prop_score < 0.20 | prop_score > 0.80)) {
      warning("Some propensity scores are below 0.20 or above 0.80. ",
              "This package is designed for randomized controlled trials (RCTs), ",
              "where propensity scores are typically well-balanced. ",
              "Extreme propensity scores may indicate an observational study or ",
              "a heavily unbalanced design. Please verify your experimental design.")
    }
  }
  
  # Resolve Y, X, D when given as column name(s) referencing data
  # Y: scalar string -> column name for outcome
  # X: character vector -> column names for covariates
  # D: scalar string -> column name for treatment
  Y_is_name <- is.character(Y) && length(Y) == 1 && !is.null(data)
  X_is_names <- is.character(X) && !is.null(data)
  D_is_name <- is.character(D) && length(D) == 1 && !is.null(data)

  if (Y_is_name || X_is_names || D_is_name) {
    if (is.null(data)) {
      stop("'data' must be provided when Y, X, or D are column name strings.")
    }
    if (!is.data.table(data)) {
      data <- as.data.table(data)
    }

    if (Y_is_name) {
      outcome_var <- Y
      if (!outcome_var %in% names(data)) {
        stop("Column '", outcome_var, "' (passed as Y) not found in data.")
      }
      Y <- data[[outcome_var]]
    }
    if (D_is_name) {
      treatment_var <- D
      if (!treatment_var %in% names(data)) {
        stop("Column '", treatment_var, "' (passed as D) not found in data.")
      }
      D <- data[[treatment_var]]
    }
    if (X_is_names) {
      missing_cols <- setdiff(X, names(data))
      if (length(missing_cols) > 0) {
        stop("The following columns (passed as X) are not in data: ",
             paste(missing_cols, collapse = ", "))
      }
      covariate_vars <- X
      X <- data[, ..covariate_vars]
    }
  }

  # Determine which interface is being used
  use_matrix_interface <- !is.null(Y) && !is.null(X) && !is.null(D)
  use_formula_interface <- !is.null(formula) && !is.null(data)
  
  # Warn if both interfaces are provided simultaneously
  if (use_matrix_interface && use_formula_interface) {
    warning("Both formula/data and Y/X/D arguments were provided. ",
            "Using the matrix interface (Y, X, D). ",
            "Please use one interface or the other to avoid confusion.")
  }
  

  if (use_matrix_interface) {
    # Matrix interface: Y, X, D provided directly (or resolved from column names)
    
    # Validate inputs
    if (is.factor(Y)) {
      stop("Y is a factor. Please convert to numeric using as.numeric() or as.numeric(as.character())")
    }
    if (!is.numeric(Y)) {
      stop("Y must be a numeric vector")
    }
    if (!is.data.frame(X) && !is.matrix(X)) {
      stop("X must be a data.frame or matrix")
    }
    if (is.factor(D)) {
      stop("D is a factor. Please convert to numeric (0/1) using as.numeric() or as.numeric(as.character())")
    }
    if (!is.numeric(D)) {
      stop("D must be a numeric vector")
    }
    
    n <- length(Y)
    if (nrow(X) != n) {
      stop(paste0("X has ", nrow(X), " rows but Y has ", n, " elements"))
    }
    if (length(D) != n) {
      stop(paste0("D has ", length(D), " elements but Y has ", n, " elements"))
    }
    
    # Generate default column names if needed
    if (is.matrix(X) && is.null(colnames(X))) {
      colnames(X) <- paste0("X", seq_len(ncol(X)))
    }
    
    # Convert X to data.table
    if (!is.data.table(X)) {
      X <- as.data.table(X)
    }
    
    # Generate default column names if needed (fallback for data.frames)
    if (is.null(names(X)) || any(names(X) == "")) {
      names(X) <- paste0("X", seq_len(ncol(X)))
    }
    
    # Set variable names if not already set (from column-name resolution)
    if (!exists("outcome_var", inherits = FALSE)) outcome_var <- "Y"
    if (!exists("treatment_var", inherits = FALSE)) treatment_var <- "D"
    if (!exists("covariate_vars", inherits = FALSE)) covariate_vars <- names(X)
    
    # Create combined data for storage (if not already provided)
    if (is.null(data)) {
      data <- cbind(data.table(Y = Y, D = D), X)
    }
    
    # Create a formula for consistency
    formula <- as.formula(paste(outcome_var, "~",
                                paste(covariate_vars, collapse = " + ")))
    
  } else if (use_formula_interface) {
    # Formula interface
    
    # Convert data to data.table if needed
    if (!is.data.table(data)) {
      data <- as.data.table(data)
    }
    
    # Handle treatment variable - can be quoted, unquoted, or a variable containing a string
    treatment_expr <- substitute(treatment)
    treatment_var <- parse_column_name(treatment_expr, parent.frame(), "treatment", data)
    
    if (is.null(treatment_var)) {
      stop("treatment must be specified when using formula interface")
    }
    
    # Determine outcome variable
    outcome_var <- all.vars(formula)[1]
    
    # Determine covariates from formula RHS
    if (length(formula) == 3) {
      # Has RHS (e.g., Y ~ X1 + X2 or Y ~ . or Y ~ . - Z)
      rhs <- formula[[3]]
      rhs_text <- paste(deparse(rhs), collapse = "")
      
      # Check if formula contains "." (meaning "all other variables")
      if (grepl("\\.", rhs_text)) {
        # User specified Y ~ . or Y ~ . - Z
        # Get all columns except outcome and treatment
        all_cols <- names(data)
        covariate_vars <- setdiff(all_cols, c(outcome_var, treatment_var))
        if (!is.null(prop_score_var)) {
          covariate_vars <- setdiff(covariate_vars, prop_score_var)
        }
        
        # Handle exclusions (e.g., Y ~ . - Z - W)
        rhs_vars <- all.vars(rhs)
        rhs_vars <- setdiff(rhs_vars, ".")
        if (length(rhs_vars) > 0) {
          covariate_vars <- setdiff(covariate_vars, rhs_vars)
        }
      } else {
        # User specified specific variables (e.g., Y ~ X1 + X2)
        covariate_vars <- all.vars(rhs)
        covariate_vars <- setdiff(covariate_vars, treatment_var)
      }
    } else {
      stop("Formula must include covariates (e.g., Y ~ . or Y ~ X1 + X2)")
    }
    
    # Rebuild formula to reflect actual covariates used
    formula <- as.formula(paste(outcome_var, "~",
                                paste(covariate_vars, collapse = " + ")))
    
    # Extract data
    Y <- data[[outcome_var]]
    X <- data[, ..covariate_vars]
    D <- data[[treatment_var]]
    n <- nrow(data)
    
    # Validate inputs for formula interface
    if (is.factor(Y)) {
      stop("Outcome variable '", outcome_var, "' is a factor. For HTE estimation, the outcome ",
           "should be numeric. If you have a binary outcome, convert it to 0/1 numeric. ",
           "Use as.numeric(as.character(Y)) or as.numeric(Y) - 1 for factor levels.")
    }
    if (!is.numeric(Y)) {
      stop("Outcome variable '", outcome_var, "' must be numeric, not ", class(Y)[1])
    }
    if (is.factor(D)) {
      stop("Treatment variable '", treatment_var, "' is a factor. Treatment should be ",
           "a numeric 0/1 indicator. Use as.numeric(as.character(D)) or as.numeric(D) - 1.")
    }
    if (!is.numeric(D)) {
      stop("Treatment variable '", treatment_var, "' must be numeric (0/1), not ", class(D)[1])
    }
    
  } else {
    # Check if the user tried the matrix interface but some args are NULL
    # (e.g., D = df$nonexistent silently returns NULL)
    null_args <- c(
      if ("Y" %in% names(cl) && is.null(Y)) "Y",
      if ("X" %in% names(cl) && is.null(X)) "X",
      if ("D" %in% names(cl) && is.null(D)) "D"
    )
    if (length(null_args) > 0) {
      stop(paste0(paste(null_args, collapse = ", "),
                  if (length(null_args) == 1) " is NULL. " else " are NULL. ",
                  "Did you pass a non-existent column name (e.g., df$wrong_name)?"))
    }
    stop("Must provide either (formula, treatment, data) or (Y, X, D)")
  }
  
  # Check that we have at least one covariate
  if (length(covariate_vars) == 0 || ncol(X) == 0) {
    stop("Formula must include at least one covariate (e.g., Y ~ X1 + X2). ",
         "Y ~ 1 is not supported.")
  }
  
  # Check for duplicate column names
  if (any(duplicated(names(data)))) {
    dup_names <- unique(names(data)[duplicated(names(data))])
    warning("Data contains duplicate column names: ",
            paste(dup_names, collapse = ", "),
            ". This may cause unexpected behavior. ",
            "Please rename or remove duplicate columns.")
  }
  
  # Check for NAs in covariates
  na_cols <- names(X)[vapply(X, function(col) anyNA(col), logical(1))]
  if (length(na_cols) > 0) {
    stop("Covariates contain NA values in column(s): ",
         paste(na_cols, collapse = ", "),
         ". Please remove or impute missing values before calling ensemble_hte().")
  }
  
  # Check for NAs in outcome
  if (any(is.na(Y))) {
    stop("Outcome variable contains NA values. Please remove or impute missing ",
         "values before calling ensemble_hte().")
  }
  
  # Warn if outcome has no variation
  if (stats::var(Y) < .Machine$double.eps) {
    warning("Outcome variable has zero or near-zero variance. ",
            "Treatment effect estimation may not be meaningful.")
  }
  
  # Validate treatment is binary {0,1}
  unique_D <- unique(D[!is.na(D)])
  if (!all(unique_D %in% c(0, 1))) {
    stop("Treatment variable must be coded as 0/1. Found values: ",
         paste(sort(unique_D), collapse = ", "), ". Please recode.")
  }
  if (length(unique_D) < 2) {
    stop("Treatment variable must contain both treated (1) and control (0) units. ",
         "Found only value(s): ", paste(sort(unique_D), collapse = ", "), ".")
  }
  if (any(is.na(D))) {
    stop("Treatment variable contains NA values. Please remove or impute missing ",
         "values before calling ensemble_hte().")
  }
  
  # Handle train_idx
  if (is.null(train_idx)) {
    train_idx <- rep(TRUE, n)
  } else if (is.logical(train_idx)) {
    if (length(train_idx) != n) {
      stop(paste0("train_idx has length ", length(train_idx), " but data has ", n, " rows"))
    }
  } else if (is.numeric(train_idx)) {
    # Convert integer indices to logical
    if (any(train_idx < 1) || any(train_idx > n)) {
      stop("train_idx contains indices outside the range 1 to n")
    }
    train_idx_logical <- rep(FALSE, n)
    train_idx_logical[train_idx] <- TRUE
    train_idx <- train_idx_logical
  } else {
    stop("train_idx must be NULL, a logical vector, or an integer vector of indices")
  }
  
  n_train <- sum(train_idx)
  if (n_train < K) {
    stop(paste0("Number of training observations (", n_train,
                ") must be at least K (", K, ")"))
  }
  
  # Validate Y for training observations (when train_idx is used)
  if (any(is.na(Y[train_idx]))) {
    stop("Y must not have NA values for training observations (where train_idx = TRUE)")
  }
  
  # Handle individual_id for panel data
  individual_id_vec <- NULL
  if (!is.null(individual_id)) {
    # Try to resolve individual_id as column name (like treatment)
    id_expr <- substitute(individual_id)
    id_resolved <- tryCatch(
      parse_column_name(id_expr, parent.frame(), "individual_id", data),
      error = function(e) NULL
    )
    
    if (!is.null(id_resolved) && id_resolved %in% names(data)) {
      individual_id_vec <- data[[id_resolved]]
    } else if (is.character(individual_id) && length(individual_id) == 1 && individual_id %in% names(data)) {
      individual_id_vec <- data[[individual_id]]
    } else if (length(individual_id) == n) {
      individual_id_vec <- individual_id
    } else {
      stop("individual_id must be a column name in the data or a vector of length n (", n, ")")
    }
    
    if (anyNA(individual_id_vec)) {
      stop("individual_id contains NA values. All observations must have an individual identifier.")
    }
    
    n_individuals <- length(unique(individual_id_vec))
    if (n_individuals == n) {
      warning("Every observation has a unique individual_id (", n_individuals,
              " unique IDs for ", n, " observations). ",
              "This is equivalent to no panel structure. ",
              "If your data is not panel data, you can omit individual_id.")
    }
    if (n_individuals < K) {
      stop("Number of unique individuals (", n_individuals,
           ") must be at least K (", K, ")")
    }
    
  }
  
  # Propensity score handling
  if (is.null(prop_score)) {
    prop_score <- rep(mean(D), n)
    if (any(prop_score < 0.20 | prop_score > 0.80)) {
      warning("The estimated propensity score (sample treatment proportion) is below 0.20 or above 0.80. ",
              "This package is designed for randomized controlled trials (RCTs), ",
              "where propensity scores are typically well-balanced. ",
              "Extreme propensity scores may indicate an observational study or ",
              "a heavily unbalanced design. Please verify your experimental design.")
    }
  } else if (length(prop_score) == 1) {
    # Expand scalar propensity score to a vector of length n
    prop_score <- rep(prop_score, n)
  } else if (length(prop_score) != n) {
    stop("prop_score has length ", length(prop_score), " but data has ", n,
         " rows. prop_score must be a single value or a vector of length n.")
  }
  W <- (prop_score * (1 - prop_score))^(-1)
  
  # Convert character columns to factor (ML backends require factor, not character)
  char_cols <- names(X)[sapply(X, is.character)]
  if (length(char_cols) > 0) {
    for (col in char_cols) {
      X[, (col) := as.factor(get(col))]
    }
  }
  
  # Scale non-binary numeric covariates for ML training if requested
  # Use only training observations for computing scaling parameters
  # Factor columns are left as-is (encoding is handled per-algorithm in fit_model)
  X_scaled <- copy(X)
  if (scale_covariates) {
    for (var in names(X_scaled)) {
      col_vals <- X_scaled[[var]]
      # Only scale numeric columns
      if (is.numeric(col_vals)) {
        unique_vals <- unique(col_vals[train_idx & !is.na(col_vals)])
        # Check if binary (only 0 and 1)
        is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
        if (!is_binary) {
          col_mean <- mean(col_vals[train_idx], na.rm = FALSE)
          col_sd <- sd(col_vals[train_idx], na.rm = FALSE)
          if (col_sd > 0) {
            X_scaled[, (var) := (col_vals - col_mean) / col_sd]
          }
        }
      }
    }
  }
  
  # Auto-detect task type if not provided
  if (is.null(task_type)) {
    Y_for_detection <- Y[train_idx]
    if (is.factor(Y_for_detection) || is.character(Y_for_detection) || length(unique(Y_for_detection)) == 2) {
      task_type <- "classif"
    } else {
      task_type <- "regr"
    }
  } else {
    task_type <- match.arg(task_type, c("regr", "classif"))
    # Warn if classification is specified but outcome looks continuous
    if (task_type == "classif") {
      n_unique <- length(unique(Y[train_idx]))
      if (n_unique > 2) {
        warning("task_type = 'classif' was specified but the outcome has ", n_unique,
                " unique values. If the outcome is continuous, use task_type = 'regr' instead.")
      }
    }
  }

  # Split the sample - stratify by D and train_idx
  stratify_var <- interaction(D, train_idx, drop = TRUE)
  splits <- create_folds(n, M = M, K = K, stratify_var = stratify_var,
                         cluster_id = individual_id_vec)

  # Train learners and compute ensemble predictions
  ite_cols <- paste0("ite_", algorithms)
  y0_cols <- paste0("y0_", algorithms)
  
  # Define function to process a single repetition
  process_repetition <- function(m) {
    # Temporary storage for this repetition
    predictions_m <- data.frame(matrix(NA_real_, nrow = n, ncol = 2 * length(algorithms)))
    colnames(predictions_m) <- c(ite_cols, y0_cols)
    
    for (k in 1:K) {
      # Test fold: all observations in fold k
      test_fold_idx <- which(splits[[m]] == k)
      
      # Training for ML: training observations NOT in fold k
      ml_train_idx <- which(train_idx & splits[[m]] != k)
      ml_train_logical <- rep(FALSE, n)
      ml_train_logical[ml_train_idx] <- TRUE
      
      # Test fold as logical vector (for predict_ite)
      test_fold_logical <- rep(FALSE, n)
      test_fold_logical[test_fold_idx] <- TRUE
      
      for (algorithm in algorithms) { 
        algo_params <- learner_params[[algorithm]]
        result <- predict_ite(
          Y, X_scaled, D, prop_score, W, ml_train_logical, algorithm,
          metalearner, r_learner, task_type, tune, tune_params,
          test_idx = test_fold_logical,
          learner_params = algo_params
        )
        predictions_m[test_fold_idx, paste0("ite_", algorithm)] <- result$predicted_ite
        predictions_m[test_fold_idx, paste0("y0_", algorithm)] <- result$predicted_y0
      }
    }
    
    # Ensemble at repetition level
    if (ensemble_strategy == "average") {
      # Simple unweighted average of ITE predictions across algorithms
      ite_rep <- rowMeans(predictions_m[, ite_cols, drop = FALSE], na.rm = TRUE)
      return(ite_rep)
    }
    
    # Cross-validated BLP ensemble (default "cv" strategy)
    ens_splits <- create_folds(n, M = 1, K = ensemble_folds, stratify_var = splits[[m]])[[1]]
    dt_ens <- data.table(
      Y = Y, 
      D = D,
      prop_score = prop_score,
      weight = W, 
      ens_fold = ens_splits, 
      fold = splits[[m]], 
      W1 = D - prop_score,
      predictions_m
    )
    
    # Create W2 columns (interaction terms) by fold
    w2_cols <- paste0("W2_ite_", algorithms)
    dt_ens[, (w2_cols) := lapply(.SD, function(x) (D - prop_score) * (x - mean(x))), 
          .SDcols = ite_cols, by = fold]
    
    # Check if any algorithm has zero W2 variance (constant predictions)
    # and exclude them from the ensemble formula
    valid_w2_cols <- character(0)
    for (col in w2_cols) {
      if (var(dt_ens[[col]]) < .Machine$double.eps) {
        warning(paste0("Algorithm predictions have zero variance for term '", col, 
                       "'. This algorithm will be excluded from the ensemble. ",
                       "This may happen with small samples or when the base learner ",
                       "cannot detect heterogeneity. ",
                       "(This is expected when using S-learner with linear models. ",
                       "Consider using flexible learners such as ranger or xgboost.)"))
      } else {
        valid_w2_cols <- c(valid_w2_cols, col)
      }
    }
    
    # If all algorithms have zero variance, fall back to simple average
    if (length(valid_w2_cols) == 0) {
      warning("All algorithms produced constant predictions. Using simple average of ITE predictions.")
      ite_rep <- rowMeans(predictions_m[, ite_cols, drop = FALSE], na.rm = TRUE)
      return(ite_rep)
    }
    
    # Ensemble regression formula
    # Only include y0_cols if we have valid Y0 predictions (not for R-learner)
    has_y0 <- !all(is.na(predictions_m[[y0_cols[1]]]))
    if (has_y0) {
      formula_ens <- as.formula(paste("Y ~ W1 +", paste(c(valid_w2_cols, y0_cols), collapse = " + ")))
    } else {
      formula_ens <- as.formula(paste("Y ~ W1 +", paste(valid_w2_cols, collapse = " + ")))
    }
    
    # Cross-validated ensemble predictions
    ite_rep <- rep(NA_real_, n)
    for (ell in 1:ensemble_folds) {
      # Fit on training obs NOT in fold ell, predict for ALL obs in fold ell
      fit_idx <- which(train_idx & ens_splits != ell)
      test_idx_ens <- which(ens_splits == ell)
      fit_ens <- lm(formula_ens, data = dt_ens[fit_idx], weights = weight)
      coefs <- coef(fit_ens)[!is.na(coef(fit_ens))]
      w2_coef_names <- names(coefs)[startsWith(names(coefs), "W2_")]
      ite_pred_names <- sub("^W2_", "", w2_coef_names)
      ite_rep[test_idx_ens] <- as.matrix(predictions_m[test_idx_ens, ite_pred_names, drop = FALSE]) %*% coefs[w2_coef_names]
    }
    
    return(ite_rep)
  }
  
  # Run repetitions (parallel or sequential)
  if (n_cores > 1) {
    # Set up parallel backend
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = n_cores)
    
    # Run in parallel
    results_list <- future.apply::future_lapply(
      1:M, 
      process_repetition,
      future.seed = TRUE
    )
  } else {
    # Run sequentially
    results_list <- lapply(1:M, process_repetition)
  }
  
  # Combine results into matrix
  ite_blp <- do.call(cbind, results_list)
  colnames(ite_blp) <- paste0("rep_", 1:M)

  # Construct the ensemble_fit object
  structure(
    list(
      ite = as.data.table(ite_blp),
      call = cl,
      formula = formula,
      treatment = treatment_var,
      data = data,
      Y = Y,
      X = X,
      D = D,
      prop_score = prop_score,
      weights = W,
      train_idx = train_idx,
      individual_id = individual_id_vec,
      splits = splits,
      n = n,
      n_train = n_train,
      M = M,
      K = K,
      algorithms = algorithms,
      metalearner = metalearner,
      r_learner = r_learner,
      ensemble_folds = ensemble_folds,
      ensemble_strategy = ensemble_strategy,
      task_type = task_type,
      scale_covariates = scale_covariates,
      tune = tune,
      tune_params = tune_params,
      learner_params = learner_params,
      n_cores = n_cores
    ),
    class = "ensemble_hte_fit"
  )
}


#' Extract ITE Predictions as a Matrix
#'
#' @description
#' Extracts the individual treatment effect (ITE) predictions from an
#' \code{ensemble_hte_fit} object and returns them as a plain numeric matrix.
#' Each column corresponds to one repetition (rep_1, rep_2, ..., rep_M).
#'
#' This is a convenience function to avoid potential issues with
#' data.table-to-matrix coercion when working with \code{fit$ite} directly.
#'
#' @param fit An object of class \code{ensemble_hte_fit} from \code{\link{ensemble_hte}}.
#'
#' @return A numeric matrix with \code{n} rows and \code{M} columns, where
#'   each column contains the ITE predictions from one repetition.
#'
#' @examples
#' \dontrun{
#' fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
#'                     algorithms = c("lm", "grf"), M = 5, K = 3)
#'
#' ite_mat <- ite(fit)
#' }
#'
#' @export
ite <- function(fit) {
  if (!inherits(fit, "ensemble_hte_fit")) {
    stop("fit must be an object of class 'ensemble_hte_fit' from ensemble_hte()")
  }
  as.matrix(fit$ite)
}


#' Print Method for ensemble_hte_fit Objects
#'
#' @description
#' Displays a formatted summary of an ensemble HTE fit, including data dimensions,
#' model specification, and split-sample parameters.
#'
#' @param x An object of class \code{ensemble_hte_fit} from \code{\link{ensemble_hte}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' \dontrun{
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
#' print(fit)
#' }
#'
#' @export
print.ensemble_hte_fit <- function(x, ...) {
  # Format metalearner name
  metalearner_full <- switch(x$metalearner,
    "t" = "T-learner",
    "s" = "S-learner",
    "x" = "X-learner",
    "r" = "R-learner",
    x$metalearner
  )
  
  # Format task type
  task_type_full <- switch(x$task_type,
    "regr" = "regression (continuous outcome)",
    "classif" = "classification (binary outcome)",
    x$task_type
  )
  
  # Format status indicators
  tune_status <- if (x$tune) "enabled" else "disabled"
  scale_status <- if (x$scale_covariates) "enabled" else "disabled"
  
  # Get outcome name
  outcome_name <- if (!is.null(x$formula)) all.vars(x$formula)[1] else "Y"
  
  # Check if subset training was used
  subset_training <- !is.null(x$n_train) && x$n_train < x$n
  
  cat("Ensemble HTE Fit\n")
  cat("================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Data:\n")
  cat("  Observations:      ", x$n, "\n", sep = "")
  if (subset_training) {
    cat("  Training obs:      ", x$n_train, " (subset)\n", sep = "")
  }
  cat("  Targeted outcome:  ", outcome_name, "\n", sep = "")
  cat("  Treatment:         ", x$treatment, "\n", sep = "")
  cat("  Covariates:        ", ncol(x$X), "\n", sep = "")
  if (!is.null(x$individual_id)) {
    n_individuals <- length(unique(x$individual_id))
    cat("  Panel data:        ", n_individuals, " individuals\n", sep = "")
  }
  cat("\n")
  cat("Model specification:\n")
  cat("  Algorithms:        ", paste(x$algorithms, collapse = ", "), "\n", sep = "")
  cat("  Metalearner:       ", metalearner_full, "\n", sep = "")
  cat("  Task type:         ", task_type_full, "\n", sep = "")
  if (x$metalearner == "r") {
    cat("  R-learner method:  ", x$r_learner, "\n", sep = "")
  }
  cat("\n")
  ensemble_strat_desc <- if (x$ensemble_strategy == "average") "simple average" else "cross-validated BLP"
  cat("Split-sample parameters:\n")
  cat("  Repetitions (M):   ", x$M, "\n", sep = "")
  cat("  Folds (K):         ", x$K, "\n", sep = "")
  cat("  Ensemble strategy: ", ensemble_strat_desc, "\n", sep = "")
  if (x$ensemble_strategy == "cv") {
    cat("  Ensemble folds:    ", x$ensemble_folds, "\n", sep = "")
  }
  cat("  Covariate scaling: ", scale_status, "\n", sep = "")
  cat("  Hyperparameter tuning: ", tune_status, "\n", sep = "")
  if (!is.null(x$individual_id)) {
    cat("  Standard errors:   cluster-robust (at individual level)\n", sep = "")
  }
  invisible(x)
}


#' Summary Method for ensemble_hte_fit Objects
#'
#' @description
#' Provides a summary of the estimated individual treatment effects (ITEs) from
#' an ensemble HTE fit, including descriptive statistics, BLP coefficients for
#' average treatment effect and heterogeneity, and GATES group estimates.
#'
#' @param object An object of class \code{ensemble_hte_fit} from \code{\link{ensemble_hte}}.
#' @param n_groups Integer. Number of groups for GATES analysis (default: 3).
#' @param group_on Character. How to form groups when \code{train_idx} was used.
#'   Passed through to \code{\link{gates}}. One of \code{"auto"} (default),
#'   \code{"all"}, or \code{"analysis"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' \dontrun{
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
#' summary(fit)
#' summary(fit, n_groups = 5)
#' }
#'
#' @export
summary.ensemble_hte_fit <- function(object, n_groups = 3, group_on = c("auto", "all", "analysis"), ...) {
  
  group_on <- match.arg(group_on)
  
  # Get outcome name
  outcome_name <- if (!is.null(object$formula)) all.vars(object$formula)[1] else "Y"
  
  cat("Ensemble HTE Summary\n")
  cat("====================\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  cat("Outcome:     ", outcome_name, "\n", sep = "")
  cat("Treatment:   ", object$treatment, "\n", sep = "")
  cat("Observations: ", object$n, "\n", sep = "")
  cat("Repetitions:  ", object$M, "\n", sep = "")
  
  # BLP Analysis
  blp_result <- blp(object)
  cat("\nBest Linear Predictor (BLP):\n")
  cat("  beta1 (ATE):       ", sprintf("%.2f", blp_result$estimates[term == "beta1", estimate]),
      " (SE: ", sprintf("%.2f", blp_result$estimates[term == "beta1", se]), 
      ", p: ", sprintf("%.3f", blp_result$estimates[term == "beta1", p_value]), ") ",
      get_stars(blp_result$estimates[term == "beta1", p_value]), "\n", sep = "")
  cat("  beta2 (HET):       ", sprintf("%.2f", blp_result$estimates[term == "beta2", estimate]),
      " (SE: ", sprintf("%.2f", blp_result$estimates[term == "beta2", se]),
      ", p: ", sprintf("%.3f", blp_result$estimates[term == "beta2", p_value]), ") ",
      get_stars(blp_result$estimates[term == "beta2", p_value]), "\n", sep = "")
  
  # Interpretation
  het_pval <- blp_result$estimates[term == "beta2", p_value]
  if (het_pval < 0.05) {
    cat("  -> Significant heterogeneity detected (p < 0.05)\n")
  } else {
    cat("  -> No significant heterogeneity detected (p >= 0.05)\n")
  }
  
  # GATES Analysis
  gates_result <- gates(object, n_groups = n_groups, group_on = group_on)
  cat("\nGroup Average Treatment Effects (GATES) with ", n_groups, " groups:\n", sep = "")
  
  cat(sprintf("  %5s  %10s  %10s  %10s\n", "Group", "Estimate", "Std.Error", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 44), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(gates_result$estimates)) {
    stars <- get_stars(gates_result$estimates$p_value[i])
    cat(sprintf("  %5d  %10.2f  %10.2f  %10.3f %s\n",
                gates_result$estimates$group[i],
                gates_result$estimates$estimate[i],
                gates_result$estimates$se[i],
                gates_result$estimates$p_value[i],
                stars))
  }
  
  # Top-bottom test
  tb <- gates_result$top_bottom
  tb_stars <- get_stars(tb$p_value)
  cat("\n  Top - Bottom:  ", sprintf("%.2f", tb$estimate),
      " (SE: ", sprintf("%.2f", tb$se),
      ", p: ", sprintf("%.3f", tb$p_value), ") ", tb_stars, "\n", sep = "")
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  invisible(object)
}
