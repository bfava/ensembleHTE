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
#' }
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
#' The ensemble combines predictions from multiple ML algorithms using a Best
#' Linear Predictor (BLP) approach. For each repetition, algorithm predictions
#' are combined via weighted least squares where weights are derived from a
#' cross-validated BLP regression. The \code{ensemble_strategy} parameter is
#' reserved for future ensemble methods (currently only "cv" is implemented).
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
#' @param treatment The treatment variable. Can be specified as:
#'   \itemize{
#'     \item An unquoted variable name: \code{treatment = D}
#'     \item A quoted string: \code{treatment = "D"}
#'     \item A variable containing the column name: \code{treat_col <- "D"; treatment = treat_col}
#'     \item Ignored when using matrix interface (use \code{D} parameter instead)
#'   }
#' @param data A data.frame or data.table containing the variables in the formula.
#'   Required if using formula interface; ignored if \code{Y}, \code{X}, \code{D}
#'   are provided.
#' @param Y Numeric vector of outcomes. Use this with \code{X} and \code{D} as an
#'   alternative to the formula interface.
#' @param X Matrix or data.frame of covariates. Use this with \code{Y} and \code{D}
#'   as an alternative to the formula interface.
#' @param D Numeric vector of treatment indicators (0/1). Use this with \code{Y}
#'   and \code{X} as an alternative to the formula interface.
#' @param prop_score Numeric vector of propensity scores (probability of treatment
#'   given covariates). If \code{NULL} (default), assumes constant propensity equal
#'   to the sample treatment proportion (appropriate for randomized experiments).
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
#' @param metalearner Character. The metalearner strategy for ITE estimation.
#'   One of:
#'   \itemize{
#'     \item \code{"r"} (default): R-learner with Robinson transformation
#'     \item \code{"t"}: T-learner with separate models per treatment arm
#'     \item \code{"s"}: S-learner with treatment as a feature
#'     \item \code{"x"}: X-learner with imputed counterfactuals
#'   }
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
#' @param tune_params List of tuning parameters:
#'   \itemize{
#'     \item \code{time}: Maximum tuning time in seconds (default: 30)
#'     \item \code{cv_folds}: Number of CV folds for tuning (default: 3)
#'     \item \code{stagnation_iters}: Stop if no improvement for this many iterations (default: 250)
#'     \item \code{stagnation_threshold}: Minimum improvement threshold (default: 0.01)
#'     \item \code{measure}: Performance measure (default: R² for regression, AUC for classification)
#'   }
#' @param ensemble_strategy Character. Strategy for combining algorithm predictions.
#'   Currently only \code{"cv"} is available, which uses cross-validated Best Linear
#'   Predictor (BLP) regression to learn optimal weights for each algorithm.
#'   See \strong{Ensemble Strategy} section for details.
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
#'   \item{M, K}{Number of repetitions and folds}
#'   \item{algorithms}{Algorithms used in ensemble}
#'   \item{metalearner}{Metalearner strategy used}
#'   \item{r_learner}{R-learner algorithm (if applicable)}
#'   \item{ensemble_folds}{Number of ensemble CV folds}
#'   \item{task_type}{Task type (regr or classif)}
#'   \item{scale_covariates}{Whether covariates were scaled}
#'   \item{tune, tune_params}{Tuning settings}
#'   \item{n_cores}{Number of cores used for parallel processing}
#' }
#'
#' @examples
#' \dontrun{
#' # Formula interface with unquoted treatment
#' fit <- ensemble_hte(
#'   formula = outcome ~ age + income + education,
#'   treatment = treated,
#'   data = mydata,
#'   algorithms = c("lm", "ranger", "grf"),
#'   M = 5, K = 5
#' )
#'
#' # Formula interface with quoted treatment
#' fit <- ensemble_hte(
#'   formula = outcome ~ age + income + education,
#'   treatment = "treated",
#'   data = mydata,
#'   algorithms = c("lm", "grf"),
#'   M = 5, K = 5
#' )
#'
#' # Treatment column name stored in a variable
#' treat_col <- "treated"
#' fit <- ensemble_hte(
#'   formula = outcome ~ .,
#'   treatment = treat_col,
#'   data = mydata
#' )
#'
#' # Matrix interface
#' fit <- ensemble_hte(
#'   Y = mydata$outcome,
#'   X = mydata[, c("age", "income", "education")],
#'   D = mydata$treated,
#'   algorithms = c("lm", "ranger"),
#'   M = 5, K = 5
#' )
#'
#' # With propensity scores
#' fit <- ensemble_hte(
#'   Y ~ . - other_outcome,
#'   treatment = treat,
#'   data = mydata,
#'   prop_score = mydata$pscore,
#'   metalearner = "x"
#' )
#'
#' # With parallel processing (4 cores)
#' fit <- ensemble_hte(
#'   Y ~ .,
#'   treatment = treat,
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
#' gates_results <- gates(fit)
#' blp_results <- blp(fit)
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
                         ensemble_strategy = "cv",
                         n_cores = 1) {
  
  # Capture the call
  cl <- match.call()
  
  # Validate metalearner argument and use first as default
  metalearner <- match.arg(metalearner)
  
  # Validate common inputs
  validate_common_inputs(M, K, algorithms, ensemble_folds, n_cores)
  
  # Validate prop_score if provided
  if (!is.null(prop_score)) {
    if (!is.numeric(prop_score) || any(prop_score < 0 | prop_score > 1)) {
      stop("prop_score must be between 0 and 1")
    }
  }
  
  # Determine which interface is being used
  use_matrix_interface <- !is.null(Y) && !is.null(X) && !is.null(D)
  use_formula_interface <- !is.null(formula) && !is.null(data)
  

  if (use_matrix_interface) {
    # Matrix interface: Y, X, D provided directly
    
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
    
    # Convert X to data.table
    X <- as.data.table(X)
    
    # Generate default column names if needed
    if (is.null(names(X)) || any(names(X) == "")) {
      names(X) <- paste0("X", seq_len(ncol(X)))
    }
    
    # Create combined data for storage
    data <- cbind(data.table(Y = Y, D = D), X)
    
    # Create a formula for consistency
    formula <- as.formula(paste("Y ~", paste(names(X), collapse = " + ")))
    treatment_var <- "D"
    outcome_var <- "Y"
    covariate_vars <- names(X)
    
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
    
    outcome_var <- all.vars(formula)[1]
    
    # Handle formula processing
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
        if (!is.null(prop_score) && is.character(prop_score)) {
          covariate_vars <- setdiff(covariate_vars, prop_score)
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
    stop("Must provide either (formula, treatment, data) or (Y, X, D)")
  }
  
  # Propensity score handling
  if (is.null(prop_score)) {
    prop_score <- rep(mean(D), n)
  }
  W <- (prop_score * (1 - prop_score))^(-1)
  
  # Scale non-binary numeric covariates for ML training if requested
  X_scaled <- copy(X)
  if (scale_covariates) {
    for (var in names(X_scaled)) {
      col_vals <- X_scaled[[var]]
      # Only scale numeric columns
      if (is.numeric(col_vals)) {
        unique_vals <- unique(col_vals[!is.na(col_vals)])
        # Check if binary (only 0 and 1)
        is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
        if (!is_binary) {
          col_mean <- mean(col_vals, na.rm = FALSE)
          col_sd <- sd(col_vals, na.rm = FALSE)
          if (col_sd > 0) {
            X_scaled[, (var) := (col_vals - col_mean) / col_sd]
          }
        }
      }
    }
  }
  
  # Auto-detect task type if not provided
  if (is.null(task_type)) {
    if (is.factor(Y) || is.character(Y) || length(unique(Y)) == 2) {
      task_type <- "classif"
    } else {
      task_type <- "regr"
    }
  }

  # Split the sample
  splits <- create_folds(n, M = M, K = K, stratify_var = D)

  # Train learners and compute ensemble predictions
  ite_cols <- paste0("ite_", algorithms)
  y0_cols <- paste0("y0_", algorithms)
  
  # Define function to process a single repetition
  process_repetition <- function(m) {
    # Temporary storage for this repetition
    predictions_m <- data.frame(matrix(NA_real_, nrow = n, ncol = 2 * length(algorithms)))
    colnames(predictions_m) <- c(ite_cols, y0_cols)
    
    for (k in 1:K) {
      train_idx <- splits[[m]] != k
      test_idx <- which(!train_idx)
      for (algorithm in algorithms) { 
        result <- predict_ite(
          Y, X_scaled, D, prop_score, W, train_idx, algorithm,
          metalearner, r_learner, task_type, tune, tune_params
        )
        predictions_m[test_idx, paste0("ite_", algorithm)] <- result$predicted_ite
        predictions_m[test_idx, paste0("y0_", algorithm)] <- result$predicted_y0
      }
    }
    
    # Ensemble at repetition level
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
                       "cannot detect heterogeneity."))
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
      train_idx_ens <- ens_splits != ell
      test_idx_ens <- which(!train_idx_ens)
      fit_ens <- lm(formula_ens, data = dt_ens[train_idx_ens], weights = weight)
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
      splits = splits,
      n = n,
      M = M,
      K = K,
      algorithms = algorithms,
      metalearner = metalearner,
      r_learner = r_learner,
      ensemble_folds = ensemble_folds,
      task_type = task_type,
      scale_covariates = scale_covariates,
      tune = tune,
      tune_params = tune_params,
      n_cores = n_cores
    ),
    class = "ensemble_hte_fit"
  )
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
#' fit <- ensemble_hte(Y ~ ., treatment = D, data = mydata)
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
  
  cat("Ensemble HTE Fit\n")
  cat("================\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Data:\n")
  cat("  Observations:      ", x$n, "\n", sep = "")
  cat("  Targeted outcome:  ", outcome_name, "\n", sep = "")
  cat("  Treatment:         ", x$treatment, "\n", sep = "")
  cat("  Covariates:        ", ncol(x$X), "\n", sep = "")
  cat("\n")
  cat("Model specification:\n")
  cat("  Algorithms:        ", paste(x$algorithms, collapse = ", "), "\n", sep = "")
  cat("  Metalearner:       ", metalearner_full, "\n", sep = "")
  cat("  Task type:         ", task_type_full, "\n", sep = "")
  if (x$metalearner == "r") {
    cat("  R-learner method:  ", x$r_learner, "\n", sep = "")
  }
  cat("\n")
  cat("Split-sample parameters:\n")
  cat("  Repetitions (M):   ", x$M, "\n", sep = "")
  cat("  Folds (K):         ", x$K, "\n", sep = "")
  cat("  Ensemble folds:    ", x$ensemble_folds, "\n", sep = "")
  cat("  Covariate scaling: ", scale_status, "\n", sep = "")
  cat("  Hyperparameter tuning: ", tune_status, "\n", sep = "")
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
#' @param n_groups Integer. Number of groups for GATES analysis (default: 5).
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns the input object.
#'
#' @examples
#' \dontrun{
#' fit <- ensemble_hte(Y ~ ., treatment = D, data = mydata)
#' summary(fit)
#' summary(fit, n_groups = 3)
#' }
#'
#' @export
summary.ensemble_hte_fit <- function(object, n_groups = 5, ...) {
  
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
  
  # ITE Distribution - calculate per repetition, then average
  ite_stats_by_rep <- lapply(1:object$M, function(m) {
    ite_m <- object$ite[[m]]
    data.frame(
      min = min(ite_m),
      q1 = quantile(ite_m, 0.25),
      median = median(ite_m),
      mean = mean(ite_m),
      q3 = quantile(ite_m, 0.75),
      max = max(ite_m),
      sd = sd(ite_m),
      pct_positive = mean(ite_m > 0) * 100
    )
  })
  ite_stats <- do.call(rbind, ite_stats_by_rep)
  
  cat("\nITE Distribution (averaged across ", object$M, " repetitions):\n", sep = "")
  cat("  Min:               ", round(mean(ite_stats$min), 4), "\n", sep = "")
  cat("  1st Qu:            ", round(mean(ite_stats$q1), 4), "\n", sep = "")
  cat("  Median:            ", round(mean(ite_stats$median), 4), "\n", sep = "")
  cat("  Mean:              ", round(mean(ite_stats$mean), 4), "\n", sep = "")
  cat("  3rd Qu:            ", round(mean(ite_stats$q3), 4), "\n", sep = "")
  cat("  Max:               ", round(mean(ite_stats$max), 4), "\n", sep = "")
  cat("  Std. Dev:          ", round(mean(ite_stats$sd), 4), "\n", sep = "")
  cat("  % positive:        ", round(mean(ite_stats$pct_positive), 1), "%\n", sep = "")
  
  # BLP Analysis
  blp_result <- blp(object)
  cat("\nBest Linear Predictor (BLP):\n")
  cat("  beta1 (ATE):       ", sprintf("%.4f", blp_result$estimates[term == "beta1", estimate]),
      " (SE: ", sprintf("%.4f", blp_result$estimates[term == "beta1", se]), 
      ", p: ", sprintf("%.4f", blp_result$estimates[term == "beta1", p_value]), ")\n", sep = "")
  cat("  beta2 (HET):       ", sprintf("%.4f", blp_result$estimates[term == "beta2", estimate]),
      " (SE: ", sprintf("%.4f", blp_result$estimates[term == "beta2", se]),
      ", p: ", sprintf("%.4f", blp_result$estimates[term == "beta2", p_value]), ")\n", sep = "")
  
  # Interpretation
  het_pval <- blp_result$estimates[term == "beta2", p_value]
  if (het_pval < 0.05) {
    cat("  -> Significant heterogeneity detected (p < 0.05)\n")
  } else {
    cat("  -> No significant heterogeneity detected (p >= 0.05)\n")
  }
  
  # GATES Analysis
  gates_result <- gates(object, n_groups = n_groups)
  cat("\nGroup Average Treatment Effects (GATES) with ", n_groups, " groups:\n", sep = "")
  
  cat(sprintf("  %5s  %10s  %10s  %10s\n", "Group", "Estimate", "Std.Error", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 44), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(gates_result$estimates)) {
    stars <- get_stars(gates_result$estimates$p_value[i])
    cat(sprintf("  %5d  %10.4f  %10.4f  %10.4f %s\n",
                gates_result$estimates$group[i],
                gates_result$estimates$estimate[i],
                gates_result$estimates$se[i],
                gates_result$estimates$p_value[i],
                stars))
  }
  
  # Top-bottom test
  tb <- gates_result$top_bottom
  tb_stars <- get_stars(tb$p_value)
  cat("\n  Top - Bottom:  ", sprintf("%.4f", tb$estimate),
      " (SE: ", sprintf("%.4f", tb$se),
      ", p: ", sprintf("%.4f", tb$p_value), ") ", tb_stars, "\n", sep = "")
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  invisible(object)
}
