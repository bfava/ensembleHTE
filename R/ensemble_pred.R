#' Fit Ensemble Prediction Model
#'
#' @description
#' Predicts an outcome variable Y using an ensemble of machine learning
#' algorithms combined with multiple sample splitting. This function
#' implements the estimation strategy developed by Fava (2025), which improves
#' statistical power by averaging predictions across M repetitions of K-fold
#' cross-fitting.
#'
#' Unlike \code{\link{ensemble_hte}} which estimates heterogeneous treatment effects,
#' this function performs standard prediction of Y given X without any treatment
#' variable or causal structure.
#'
#' The function supports two interfaces:
#' \itemize{
#'   \item \strong{Formula interface}: Specify \code{formula} and \code{data}
#'   \item \strong{Matrix interface}: Specify \code{Y} and \code{X} directly
#' }
#'
#' @section Training on a Subset:
#' The \code{train_idx} parameter allows training models on a subset of observations
#' while generating predictions for all observations. This is useful when the outcome
#' Y is only observed for some units (e.g., only treated units in an experiment) but
#' you want predictions for everyone.
#'
#' When \code{train_idx} is provided:
#' \itemize{
#'   \item Y can have NA values for observations where \code{train_idx = FALSE}
#'   \item Cross-fitting splits ALL observations into K folds, stratifying by \code{train_idx}
#'   \item For each fold k, models are trained on training observations NOT in fold k
#'   \item Predictions are generated for ALL observations in fold k (both training and non-training)
#'   \item Each observation gets exactly one prediction per repetition (from the model where they were in the test fold)
#'   \item Ensemble weights are estimated using only training observations
#'   \item Summary statistics are computed only on training observations
#' }
#'
#' @section Ensemble Strategy:
#' The ensemble combines predictions from multiple ML algorithms. The
#' \code{ensemble_strategy} parameter controls how predictions are combined:
#' \itemize{
#'   \item \code{"cv"} (default): Uses cross-validated OLS regression of Y on the
#'     predicted values from each algorithm to learn optimal combination weights.
#'   \item \code{"average"}: Uses simple averaging across all algorithm predictions
#'     with equal weights. This is useful when you want a simpler ensemble or when
#'     using only a single algorithm.
#' }
#'
#' @references
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#'
#' @param formula A formula specifying the outcome and covariates (e.g.,
#'   \code{Y ~ X1 + X2} or \code{Y ~ .}). Use \code{~ . - Z} to exclude variables.
#'   Required if \code{Y} and \code{X} are not provided.
#' @param data A data.frame or data.table containing the variables in the formula.
#'   Required if using formula interface; ignored if \code{Y} and \code{X}
#'   are provided.
#' @param Y Numeric vector of outcomes. Use this with \code{X} as an
#'   alternative to the formula interface. Can contain NA values for observations
#'   where \code{train_idx = FALSE}.
#' @param X Matrix or data.frame of covariates. Use this with \code{Y}
#'   as an alternative to the formula interface.
#' @param train_idx Optional logical or integer vector indicating which observations
#'   to use for training. If \code{NULL} (default), all observations are used.
#'   If provided:
#'   \itemize{
#'     \item Logical vector: \code{TRUE} for training observations
#'     \item Integer vector: indices of training observations
#'   }
#'   Y must not have NA values for training observations.
#'   See \strong{Training on a Subset} section for details on how this affects
#'   cross-fitting and predictions.
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
#' @param tune_params List of tuning parameters:
#'   \itemize{
#'     \item \code{time}: Maximum tuning time in seconds (default: 30)
#'     \item \code{cv_folds}: Number of CV folds for tuning (default: 3)
#'     \item \code{stagnation_iters}: Stop if no improvement for this many iterations (default: 250)
#'     \item \code{stagnation_threshold}: Minimum improvement threshold (default: 0.01)
#'     \item \code{measure}: Performance measure (default: R² for regression, AUC for classification)
#'   }
#' @param ensemble_strategy Character. Strategy for combining algorithm predictions:
#'   \itemize{
#'     \item \code{"cv"} (default): Cross-validated OLS regression of Y on algorithm
#'       predictions. Learns optimal weights for each algorithm.
#'     \item \code{"average"}: Simple average of all algorithm predictions. No weight
#'       learning; all algorithms contribute equally. Works with a single algorithm.
#'   }
#'   See \strong{Ensemble Strategy} section for more details.
#' @param n_cores Integer. Number of cores for parallel processing of repetitions.
#'   Default is 1 (sequential). Set to higher values to parallelize the M repetitions.
#'   Uses the \code{future} framework, so users can also set up their own parallel
#'   backend via \code{future::plan()} before calling this function.
#'
#' @return An object of class \code{ensemble_pred_fit} containing:
#' \describe{
#'   \item{predictions}{data.table of predictions with M columns (one per repetition)}
#'   \item{call}{The matched function call}
#'   \item{formula}{The formula used (or constructed from Y/X)}
#'   \item{data}{The original data (or constructed data.table from Y/X)}
#'   \item{Y}{Vector of outcomes}
#'   \item{X}{data.table of covariates (unscaled)}
#'   \item{train_idx}{Logical vector indicating training observations}
#'   \item{splits}{List of fold assignments for each repetition}
#'   \item{n}{Number of observations}
#'   \item{n_train}{Number of training observations}
#'   \item{M, K}{Number of repetitions and folds}
#'   \item{algorithms}{Algorithms used in ensemble}
#'   \item{ensemble_folds}{Number of ensemble CV folds}
#'   \item{ensemble_strategy}{Strategy used for combining predictions ("cv" or "average")}
#'   \item{task_type}{Task type (regr or classif)}
#'   \item{scale_covariates}{Whether covariates were scaled}
#'   \item{tune, tune_params}{Tuning settings}
#'   \item{n_cores}{Number of cores used for parallel processing}
#' }
#'
#' @examples
#' \dontrun{
#' # Formula interface
#' fit <- ensemble_pred(
#'   formula = outcome ~ age + income + education,
#'   data = mydata,
#'   algorithms = c("lm", "ranger", "grf"),
#'   M = 5, K = 5
#' )
#'
#' # Matrix interface
#' fit <- ensemble_pred(
#'   Y = mydata$outcome,
#'   X = mydata[, c("age", "income", "education")],
#'   algorithms = c("lm", "ranger"),
#'   M = 5, K = 5
#' )
#'
#' # Training on a subset (e.g., outcome only observed for treated units)
#' # Y has NA for control units
#' fit <- ensemble_pred(
#'   formula = outcome ~ age + income,
#'   data = mydata,
#'   train_idx = mydata$treated == 1,  # Only train on treated
#'   M = 5, K = 5
#' )
#' # Predictions are available for ALL observations
#' # Use clan() or gavs() to analyze predictions
#'
#' # With parallel processing (4 cores)
#' fit <- ensemble_pred(
#'   Y ~ .,
#'   data = mydata,
#'   M = 10, K = 5,
#'   n_cores = 4
#' )
#'
#' # Print and summarize results
#' print(fit)
#' summary(fit)
#'
#' # Use for downstream analysis
#' gavs_results <- gavs(fit)
#' blp_results <- blp_pred(fit)
#' }
#'
#' @seealso
#' \code{\link{gavs}} for Group Averages analysis,
#' \code{\link{blp_pred}} for Best Linear Predictor analysis,
#' \code{\link{clan}} for Classification Analysis,
#' \code{\link{ensemble_hte}} for heterogeneous treatment effect estimation
#'
#' @export
ensemble_pred <- function(formula = NULL, data = NULL,
                          Y = NULL, X = NULL,
                          train_idx = NULL,
                          M = 2, K = 3, 
                          algorithms = c("lm", "grf"), 
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
                          ensemble_strategy = c("cv", "average"),
                          n_cores = 1) {
  
  # Capture the call
  cl <- match.call()
  
  # Validate ensemble_strategy
  ensemble_strategy <- match.arg(ensemble_strategy)
  
  # Validate common inputs
  validate_common_inputs(M, K, algorithms, ensemble_folds, n_cores)
  
  # Determine which interface is being used
  use_matrix_interface <- !is.null(Y) && !is.null(X)
  use_formula_interface <- !is.null(formula) && !is.null(data)
  
  if (use_matrix_interface) {
    # Matrix interface: Y, X provided directly
    
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
    
    n <- length(Y)
    if (nrow(X) != n) {
      stop(paste0("X has ", nrow(X), " rows but Y has ", n, " elements"))
    }
    
    # Convert X to data.table
    X <- as.data.table(X)
    
    # Generate default column names if needed
    if (is.null(names(X)) || any(names(X) == "")) {
      names(X) <- paste0("X", seq_len(ncol(X)))
    }
    
    # Create combined data for storage
    data <- cbind(data.table(Y = Y), X)
    
    # Create a formula for consistency
    formula <- as.formula(paste("Y ~", paste(names(X), collapse = " + ")))
    outcome_var <- "Y"
    covariate_vars <- names(X)
    
  } else if (use_formula_interface) {
    # Formula interface
    
    # Convert data to data.table if needed
    if (!is.data.table(data)) {
      data <- as.data.table(data)
    }
    
    outcome_var <- all.vars(formula)[1]
    
    # Handle formula processing
    if (length(formula) == 3) {
      # Has RHS (e.g., Y ~ X1 + X2 or Y ~ . or Y ~ . - Z)
      rhs <- formula[[3]]
      rhs_text <- paste(deparse(rhs), collapse = "")
      
      # Check if formula contains "." (meaning "all other variables")
      if (grepl("\\.", rhs_text)) {
        # User specified Y ~ . or Y ~ . - Z
        # Get all columns except outcome
        all_cols <- names(data)
        covariate_vars <- setdiff(all_cols, outcome_var)
        
        # Handle exclusions (e.g., Y ~ . - Z - W)
        rhs_vars <- all.vars(rhs)
        rhs_vars <- setdiff(rhs_vars, ".")
        if (length(rhs_vars) > 0) {
          covariate_vars <- setdiff(covariate_vars, rhs_vars)
        }
      } else {
        # User specified specific variables (e.g., Y ~ X1 + X2)
        covariate_vars <- all.vars(rhs)
      }
    } else {
      stop("Formula must include covariates (e.g., Y ~ . or Y ~ X1 + X2)")
    }
    
    # Extract data
    Y <- data[[outcome_var]]
    X <- data[, ..covariate_vars]
    n <- nrow(data)
    
    # Validate inputs for formula interface
    if (is.factor(Y)) {
      stop("Outcome variable '", outcome_var, "' is a factor. For prediction, the outcome ",
           "should be numeric. If you have a binary outcome, convert it to 0/1 numeric. ",
           "Use as.numeric(as.character(Y)) or as.numeric(Y) - 1 for factor levels.")
    }
    if (!is.numeric(Y)) {
      stop("Outcome variable '", outcome_var, "' must be numeric, not ", class(Y)[1])
    }
    
  } else {
    stop("Must provide either (formula, data) or (Y, X)")
  }
  
  # Handle train_idx
  if (is.null(train_idx)) {
    # Default: use all observations for training
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
    stop(paste0("Number of training observations (", n_train, ") must be at least K (", K, ")"))
  }
  
  # Validate Y for training observations
  if (any(is.na(Y[train_idx]))) {
    stop("Y must not have NA values for training observations (where train_idx = TRUE)")
  }
  
  # Scale non-binary numeric covariates for ML training if requested
  # Use only training observations for computing scaling parameters
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
  
  # Auto-detect task type if not provided (use only training observations)
  if (is.null(task_type)) {
    Y_train <- Y[train_idx]
    if (is.factor(Y_train) || is.character(Y_train) || length(unique(Y_train)) == 2) {
      task_type <- "classif"
    } else {
      task_type <- "regr"
    }
  }

  # Split the sample - stratify by train_idx to ensure each fold has both types
  splits <- create_folds(n, M = M, K = K, stratify_var = train_idx)

  # Train learners and compute ensemble predictions
  pred_cols <- paste0("pred_", algorithms)
  
  # Define function to process a single repetition
  process_repetition <- function(m) {
    # Temporary storage for this repetition
    predictions_m <- data.frame(matrix(NA_real_, nrow = n, ncol = length(algorithms)))
    colnames(predictions_m) <- pred_cols
    
    for (k in 1:K) {
      # Test fold: all observations in fold k
      test_fold_idx <- which(splits[[m]] == k)
      
      # Training for ML: training observations NOT in fold k
      ml_train_idx <- which(train_idx & splits[[m]] != k)
      
      for (algorithm in algorithms) { 
        # Fit model on training subset (training obs not in fold k)
        model <- fit_model(
          Y[ml_train_idx], 
          X_scaled[ml_train_idx, , drop = FALSE], 
          algorithm, 
          task_type, 
          tune, 
          tune_params
        )
        
        # Predict on ALL observations in test fold (both training and non-training)
        predicted_y <- .predict_model(model, X_scaled[test_fold_idx, , drop = FALSE])$response
        predictions_m[test_fold_idx, paste0("pred_", algorithm)] <- predicted_y
      }
    }
    
    # Ensemble at repetition level
    pred_rep <- rep(NA_real_, n)
    
    if (ensemble_strategy == "average") {
      # Simple average of all algorithm predictions
      pred_rep <- rowMeans(predictions_m, na.rm = TRUE)
    } else {
      # Cross-validated OLS ensemble (default "cv" strategy)
      # Create ensemble folds for ALL observations, stratified by original K-fold
      ens_splits <- create_folds(n, M = 1, K = ensemble_folds, 
                                 stratify_var = splits[[m]])[[1]]
      
      # Check for zero-variance predictors and exclude them
      # Use only training observations to compute variance
      valid_pred_cols <- character(0)
      for (col in pred_cols) {
        col_var <- var(predictions_m[[col]][train_idx], na.rm = TRUE)
        if (is.na(col_var) || col_var < .Machine$double.eps) {
          warning(paste0("Algorithm predictions have zero variance for '", col, 
                         "'. This algorithm will be excluded from the ensemble. ",
                         "This may happen with small samples or when the learner ",
                         "cannot capture the signal."))
        } else {
          valid_pred_cols <- c(valid_pred_cols, col)
        }
      }
      
      # If all algorithms have zero variance, fall back to simple average
      if (length(valid_pred_cols) == 0) {
        warning("All algorithms produced constant predictions. Using simple average.")
        pred_rep <- rowMeans(predictions_m, na.rm = TRUE)
        return(pred_rep)
      }
      
      # Prepare ensemble data
      dt_ens <- data.table(Y = Y, predictions_m)
      formula_ens <- as.formula(paste("Y ~", paste(valid_pred_cols, collapse = " + ")))
      
      # Cross-validated ensemble predictions
      for (ell in 1:ensemble_folds) {
        # Fit on training obs NOT in fold ell, predict for ALL obs in fold ell
        fit_idx <- which(train_idx & ens_splits != ell)
        test_idx <- which(ens_splits == ell)
        fit_ens <- lm(formula_ens, data = dt_ens[fit_idx])
        pred_rep[test_idx] <- predict(fit_ens, newdata = dt_ens[test_idx])
      }
    }
    
    return(pred_rep)
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
  predictions <- do.call(cbind, results_list)
  colnames(predictions) <- paste0("rep_", 1:M)

  # Construct the ensemble_pred_fit object
  structure(
    list(
      predictions = as.data.table(predictions),
      call = cl,
      formula = formula,
      data = data,
      Y = Y,
      X = X,
      train_idx = train_idx,
      splits = splits,
      n = n,
      n_train = n_train,
      M = M,
      K = K,
      algorithms = algorithms,
      ensemble_folds = ensemble_folds,
      ensemble_strategy = ensemble_strategy,
      task_type = task_type,
      scale_covariates = scale_covariates,
      tune = tune,
      tune_params = tune_params,
      n_cores = n_cores
    ),
    class = "ensemble_pred_fit"
  )
}


#' Print Method for ensemble_pred_fit Objects
#'
#' @description
#' Displays a formatted summary of an ensemble prediction fit, including data
#' dimensions, model specification, and split-sample parameters.
#'
#' @param x An object of class \code{ensemble_pred_fit} from \code{\link{ensemble_pred}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object \code{x}.
#'
#' @examples
#' \dontrun{
#' fit <- ensemble_pred(Y ~ ., data = mydata)
#' print(fit)
#' }
#'
#' @export
print.ensemble_pred_fit <- function(x, ...) {
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
  
  cat("Ensemble Prediction Fit\n")
  cat("=======================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Data:\n")
  cat("  Observations:      ", x$n, "\n", sep = "")
  if (subset_training) {
    cat("  Training obs:      ", x$n_train, " (subset)\n", sep = "")
  }
  cat("  Outcome:           ", outcome_name, "\n", sep = "")
  cat("  Covariates:        ", ncol(x$X), "\n", sep = "")
  cat("\n")
  cat("Model specification:\n")
  cat("  Algorithms:        ", paste(x$algorithms, collapse = ", "), "\n", sep = "")
  cat("  Task type:         ", task_type_full, "\n", sep = "")
  cat("\n")
  cat("Split-sample parameters:\n")
  cat("  Repetitions (M):   ", x$M, "\n", sep = "")
  cat("  Folds (K):         ", x$K, "\n", sep = "")
  ensemble_strat_desc <- if (x$ensemble_strategy == "average") "simple average" else "cross-validated OLS"
  cat("  Ensemble strategy: ", ensemble_strat_desc, "\n", sep = "")
  if (x$ensemble_strategy == "cv") {
    cat("  Ensemble folds:    ", x$ensemble_folds, "\n", sep = "")
  }
  cat("  Covariate scaling: ", scale_status, "\n", sep = "")
  cat("  Hyperparameter tuning: ", tune_status, "\n", sep = "")
  invisible(x)
}


#' Summary Method for ensemble_pred_fit Objects
#'
#' @description
#' Provides a comprehensive summary of an ensemble prediction fit, including
#' descriptive statistics, prediction accuracy metrics, BLP calibration test,
#' and GAVS group averages.
#' When \code{train_idx} was used, accuracy metrics are computed only on
#' training observations (where Y is observed).
#'
#' @param object An object of class \code{ensemble_pred_fit} from \code{\link{ensemble_pred}}.
#' @param n_groups Number of groups for GAVS analysis (default: 5).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a \code{summary.ensemble_pred_fit} object containing:
#' \itemize{
#'   \item call: The original function call
#'   \item outcome: Name of the outcome variable
#'   \item n: Total number of observations
#'   \item n_train: Number of training observations
#'   \item M: Number of repetitions
#'   \item prediction_summary: Summary statistics for predictions
#'   \item metrics: Prediction accuracy metrics (R-squared, RMSE, MAE, correlation)
#'   \item blp: BLP calibration test results
#'   \item gavs: GAVS group average results
#' }
#'
#' @examples
#' \dontrun{
#' fit <- ensemble_pred(Y ~ ., data = mydata)
#' summary(fit)
#' summary(fit, n_groups = 3)
#' }
#'
#' @export
summary.ensemble_pred_fit <- function(object, n_groups = 5, ...) {
  pred_mean <- rowMeans(object$predictions)
  
  # Get outcome name
  outcome_name <- if (!is.null(object$formula)) all.vars(object$formula)[1] else "Y"
  
  # Check if subset training was used
  subset_training <- !is.null(object$n_train) && object$n_train < object$n
  train_idx <- if (!is.null(object$train_idx)) object$train_idx else rep(TRUE, object$n)
  
  cat("Ensemble Prediction Summary\n")
  cat("===========================\n\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  cat("Outcome:     ", outcome_name, "\n", sep = "")
  cat("Observations: ", object$n, "\n", sep = "")
  if (subset_training) {
    cat("Training obs: ", object$n_train, " (accuracy computed on these)\n", sep = "")
  }
  cat("Repetitions:  ", object$M, "\n", sep = "")
  
  cat("\nPrediction Summary (averaged across ", object$M, " repetitions):\n", sep = "")
  cat("  All observations:\n")
  print(summary(pred_mean))
  
  if (subset_training) {
    cat("  Training observations only:\n")
    print(summary(pred_mean[train_idx]))
  }
  
  # Compute prediction accuracy metrics per repetition, then average
  # Only use training observations for accuracy computation
  Y_train <- object$Y[train_idx]
  var_Y <- var(Y_train)
  
  # Calculate metrics for each repetition
  metrics_by_rep <- lapply(1:object$M, function(m) {
    pred_m <- object$predictions[[m]][train_idx]
    residuals_m <- Y_train - pred_m
    data.frame(
      r_squared = 1 - var(residuals_m) / var_Y,
      rmse = sqrt(mean(residuals_m^2)),
      mae = mean(abs(residuals_m)),
      correlation = cor(Y_train, pred_m)
    )
  })
  
  # Average across repetitions
  metrics_df <- do.call(rbind, metrics_by_rep)
  
  metric_note <- if (subset_training) " (training obs only)" else ""
  cat("\nPrediction Accuracy", metric_note, " (averaged across ", object$M, " repetitions):\n", sep = "")
  cat("  R-squared:         ", round(mean(metrics_df$r_squared), 4), "\n", sep = "")
  cat("  RMSE:              ", round(mean(metrics_df$rmse), 4), "\n", sep = "")
  cat("  MAE:               ", round(mean(metrics_df$mae), 4), "\n", sep = "")
  cat("  Correlation:       ", round(mean(metrics_df$correlation), 4), "\n", sep = "")
  
  # BLP Analysis
  blp_result <- blp_pred(object)
  cat("\nBest Linear Predictor (BLP):\n")
  cat("  intercept:         ", sprintf("%.4f", blp_result$estimates[blp_result$estimates$term == "intercept", "estimate"]),
      " (SE: ", sprintf("%.4f", blp_result$estimates[blp_result$estimates$term == "intercept", "se"]), 
      ", p: ", sprintf("%.4f", blp_result$estimates[blp_result$estimates$term == "intercept", "p_value"]), ")\n", sep = "")
  cat("  slope:             ", sprintf("%.4f", blp_result$estimates[blp_result$estimates$term == "beta", "estimate"]),
      " (SE: ", sprintf("%.4f", blp_result$estimates[blp_result$estimates$term == "beta", "se"]),
      ", p: ", sprintf("%.4f", blp_result$estimates[blp_result$estimates$term == "beta", "p_value"]), ")\n", sep = "")
  
  # Interpretation
  cat("  -> Intercept close to 0 and slope close to 1 indicate good calibration\n")
  
  # GAVS Analysis
  gavs_result <- gavs(object, n_groups = n_groups)
  cat("\nGroup Averages (GAVS) with ", n_groups, " groups:\n", sep = "")
  
  cat(sprintf("  %5s  %10s  %10s  %10s\n", "Group", "Estimate", "Std.Error", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 44), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(gavs_result$estimates)) {
    stars <- get_stars(gavs_result$estimates$p_value[i])
    cat(sprintf("  %5d  %10.4f  %10.4f  %10.4f %s\n",
                gavs_result$estimates$group[i],
                gavs_result$estimates$estimate[i],
                gavs_result$estimates$se[i],
                gavs_result$estimates$p_value[i],
                stars))
  }
  
  # Top-bottom test
  tb <- gavs_result$top_bottom
  tb_stars <- get_stars(tb$p_value)
  cat("\n  Top - Bottom:  ", sprintf("%.4f", tb$estimate),
      " (SE: ", sprintf("%.4f", tb$se),
      ", p: ", sprintf("%.4f", tb$p_value), ") ", tb_stars, "\n", sep = "")
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  # Create and return summary object
  result <- structure(
    list(
      call = object$call,
      outcome = outcome_name,
      n = object$n,
      n_train = if (subset_training) object$n_train else object$n,
      M = object$M,
      prediction_summary = summary(pred_mean),
      metrics = colMeans(metrics_df),
      blp = blp_result,
      gavs = gavs_result
    ),
    class = "summary.ensemble_pred_fit"
  )
  
  invisible(result)
}
