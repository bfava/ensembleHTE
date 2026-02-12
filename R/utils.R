

# K-fold splitting
create_folds <- function(n, K = 5, M = 1, stratify_var = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate M repetitions, each returning a vector of fold assignments
  lapply(1:M, function(rep) {
    # Create fold assignments for this repetition
    if (!is.null(stratify_var)) {
      # Stratified sampling: ensure balanced folds by stratify_var
      folds <- rep(NA, n)
      for (level in unique(stratify_var)) {
        idx <- which(stratify_var == level)
        n_level <- length(idx)
        folds[idx] <- sample(rep(1:K, length.out = n_level))
      }
    } else {
      # Simple random assignment
      folds <- sample(rep(1:K, length.out = n))
    }
    
    # Return fold assignment vector (values 1 to K)
    folds
  })
}



#' Create Quantile Groups
#'
#' @description
#' Divides a numeric vector into g quantile-based groups. Attempts to use
#' `cut_number()` for equal-sized bins, falls back to rank-based groups if needed.
#'
#' @param x Numeric vector to be divided into groups
#' @param n_groups Integer, number of groups to create
#'
#' @return Integer vector of group assignments (1 to n_groups)
#'
#' @importFrom ggplot2 cut_number
#'
#' @keywords internal
create_groups <- function(x, n_groups) {
  # Try cut_number for equal-sized bins (better when possible)

  groups <- tryCatch(
    {
      result <- as.numeric(ggplot2::cut_number(x, n_groups))
      # Verify we got the requested number of groups
      if (length(unique(result)) == n_groups) {
        return(result)
      }
      NULL  # Signal to use fallback
    },
    error = function(e) NULL
  )
  
  # Fallback to ntile-like behavior if cut_number failed or didn't produce n_groups groups
  # This is equivalent to dplyr::ntile without the dependency
  n <- length(x)
  as.integer(floor((n_groups * (rank(x, ties.method = "first") - 1) / n) + 1))
}

#' Regression from Formula with Robust Standard Errors
#'
#' @description
#' Fits a linear regression model from a formula (or formula string) with 
#' optional weights and computes HC1 robust standard errors.
#'
#' @param formula Formula object or character string representing the regression formula
#' @param data Data frame containing the regression data
#' @param weights Optional numeric vector of weights for weighted regression.
#'   Must have the same length as the number of rows in data.
#'
#' @return List containing:
#' \itemize{
#'   \item coef: Coefficient test results with robust standard errors
#'   \item vcov: Robust variance-covariance matrix (HC1)
#'   \item model_matrix: Model matrix from the regression
#'   \item residuals: Regression residuals
#'   \item fitted: Fitted values
#' }
#'
#' @importFrom stats lm as.formula model.matrix residuals fitted
#' @importFrom sandwich vcovHC
#' @importFrom lmtest coeftest
#'
#' @keywords internal
reg_from_formula <- function(formula, data, weights = NULL) {
  # Convert to formula if character string
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }
  
  # Fit model with or without weights
  if (is.null(weights)) {
    fit <- stats::lm(formula, data = data)
  } else {
    fit <- stats::lm(formula, data = data, weights = weights)
  }
  
  # Compute robust standard errors
  vcov_robust <- sandwich::vcovHC(fit, type = "HC1")
  coef_robust <- lmtest::coeftest(fit, vcov = vcov_robust)
  
  # Extract model components
  model_matrix <- stats::model.matrix(formula, data = data)
  reg_residuals <- stats::residuals(fit)
  reg_fitted <- stats::fitted(fit)
  
  list(
    coef = coef_robust,
    vcov = vcov_robust,
    model_matrix = model_matrix,
    residuals = reg_residuals,
    fitted = reg_fitted,
    weights = weights
  )
}


#' Get Significance Stars
#'
#' @description
#' Converts p-values to significance stars following standard convention:
#' \itemize{
#'   \item \code{***}: p < 0.001
#'   \item \code{**}: p < 0.01
#'   \item \code{*}: p < 0.05
#'   \item \code{.}: p < 0.1
#'   \item (empty): p >= 0.1
#' }
#'
#' @param p Numeric vector of p-values
#'
#' @return Character vector of significance stars
#'
#' @examples
#' \dontrun{
#' get_stars(c(0.0001, 0.005, 0.03, 0.08, 0.5))
#' # Returns: "***" "**" "*" "." ""
#' }
#'
#' @keywords internal
get_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*",
                       ifelse(p < 0.1, ".", ""))))
}


#' Combine Multiple Ensemble Objects
#'
#' @description
#' Combines multiple ensemble fit objects (either \code{ensemble_hte_fit} or
#' \code{ensemble_pred_fit}) that were run on the same data with the same
#' specification. This allows users to run multiple estimation sessions
#' (e.g., with different numbers of repetitions) and combine them to increase
#' the total number of repetitions M.
#'
#' This is particularly useful when running computations on servers or clusters,
#' where you can distribute the work across multiple sessions or nodes. For example,
#' you can run M=50 repetitions on 5 different servers with M=10 each, save the
#' results, and then combine them into a single object with M=50 total repetitions.
#'
#' The function validates that all objects have the same:
#' \itemize{
#'   \item Class (all HTE or all pred)
#'   \item Sample size (n)
#'   \item Outcome variable name
#'   \item Covariate names (must match exactly)
#'   \item Number of folds (K)
#'   \item Algorithms
#'   \item Ensemble folds
#'   \item Task type
#'   \item Metalearner (for HTE objects)
#'   \item Treatment variable (for HTE objects)
#' }
#'
#' @param ... Two or more ensemble fit objects of the same type
#'   (\code{ensemble_hte_fit} or \code{ensemble_pred_fit}).
#'
#' @return A combined ensemble fit object of the same class as the inputs,
#'   with the predictions/ITEs and splits concatenated across all inputs.
#'   The \code{M} value is updated to reflect the total number of repetitions.
#'
#' @examples
#' \dontrun{
#' # Run ensemble_hte in two separate sessions (e.g., on different servers)
#' # Session 1:
#' fit1 <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = mydata, M = 5)
#' saveRDS(fit1, "fit_session1.rds")
#'
#' # Session 2:
#' fit2 <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = mydata, M = 5)
#' saveRDS(fit2, "fit_session2.rds")
#'
#' # Later, combine the results:
#' fit1 <- readRDS("fit_session1.rds")
#' fit2 <- readRDS("fit_session2.rds")
#' fit_combined <- combine_ensembles(fit1, fit2)
#' print(fit_combined)  # Shows M = 10
#'
#' # Works the same for ensemble_pred
#' pred1 <- ensemble_pred(Y ~ X1 + X2, data = mydata, M = 3)
#' pred2 <- ensemble_pred(Y ~ X1 + X2, data = mydata, M = 3)
#' pred_combined <- combine_ensembles(pred1, pred2)
#' }
#'
#' @export
combine_ensembles <- function(...) {
  fits <- list(...)
  
  if (length(fits) < 2) {
    stop("At least two ensemble fit objects are required")
  }
  
  # Determine the class of the first object
  first_class <- class(fits[[1]])[1]
  
  if (!first_class %in% c("ensemble_hte_fit", "ensemble_pred_fit")) {
    stop("All objects must be of class 'ensemble_hte_fit' or 'ensemble_pred_fit'")
  }
  
  is_hte <- first_class == "ensemble_hte_fit"
  
  # Validate all objects have the same class
  for (i in 2:length(fits)) {
    if (class(fits[[i]])[1] != first_class) {
      stop(paste0("All objects must be of the same class. Object 1 is '", 
                  first_class, "' but object ", i, " is '", class(fits[[i]])[1], "'"))
    }
  }
  
  # Extract the first fit as reference
  ref <- fits[[1]]
  
  # Define which parameters must match
  params_to_check <- c("n", "K", "algorithms", "ensemble_folds", "task_type", 
                       "scale_covariates")
  if (is_hte) {
    params_to_check <- c(params_to_check, "metalearner", "treatment", "r_learner")
  }
  
  # Validate all parameters match
  for (i in 2:length(fits)) {
    fit_i <- fits[[i]]
    
    for (param in params_to_check) {
      ref_val <- ref[[param]]
      fit_val <- fit_i[[param]]
      
      # Handle vector comparisons
      if (length(ref_val) != length(fit_val) || !all(ref_val == fit_val)) {
        stop(paste0("Parameter '", param, "' differs between object 1 (", 
                    paste(ref_val, collapse = ", "), ") and object ", i, " (",
                    paste(fit_val, collapse = ", "), ")"))
      }
    }
    
    # Check formula matches (compare as character to handle attribute differences)
    # Use paste with collapse to handle long formulas that deparse to multiple lines
    if (paste(deparse(ref$formula), collapse = "") != paste(deparse(fit_i$formula), collapse = "")) {
      stop(paste0("Formula differs between object 1 and object ", i))
    }
    
    # Check outcome variable matches
    ref_outcome <- all.vars(ref$formula)[1]
    fit_outcome <- all.vars(fit_i$formula)[1]
    if (ref_outcome != fit_outcome) {
      stop(paste0("Outcome variable differs between object 1 ('", ref_outcome, 
                  "') and object ", i, " ('", fit_outcome, "')"))
    }
    
    # Check covariate names match
    ref_covars <- names(ref$X)
    fit_covars <- names(fit_i$X)
    if (length(ref_covars) != length(fit_covars) || !all(ref_covars == fit_covars)) {
      stop(paste0("Covariate names differ between object 1 (", 
                  paste(ref_covars, collapse = ", "), ") and object ", i, " (",
                  paste(fit_covars, collapse = ", "), ")"))
    }
  }
  
  # Combine predictions/ITEs and splits
  if (is_hte) {
    # Combine ITE data.tables
    all_ites <- lapply(fits, function(f) f$ite)
    combined_ite <- do.call(cbind, all_ites)
    
    # Rename columns to be sequential
    total_M <- sum(sapply(fits, function(f) f$M))
    colnames(combined_ite) <- paste0("rep_", 1:total_M)
    
    # Combine splits
    all_splits <- unlist(lapply(fits, function(f) f$splits), recursive = FALSE)
    
    # Create combined object
    result <- ref
    result$ite <- combined_ite
    result$splits <- all_splits
    result$M <- total_M
    result$call <- match.call()
    
  } else {
    # Combine prediction data.tables
    all_preds <- lapply(fits, function(f) f$predictions)
    combined_pred <- do.call(cbind, all_preds)
    
    # Rename columns to be sequential
    total_M <- sum(sapply(fits, function(f) f$M))
    colnames(combined_pred) <- paste0("rep_", 1:total_M)
    
    # Combine splits
    all_splits <- unlist(lapply(fits, function(f) f$splits), recursive = FALSE)
    
    # Create combined object
    result <- ref
    result$predictions <- combined_pred
    result$splits <- all_splits
    result$M <- total_M
    result$call <- match.call()
  }
  
  result
}

#' Validate Common Inputs
#'
#' @description
#' Validates common input parameters for ensemble functions.
#' Internal function.
#'
#' @param M Number of repetitions
#' @param K Number of folds
#' @param algorithms Vector of algorithm names
#' @param ensemble_folds Number of ensemble folds
#' @param n_cores Number of cores
#'
#' @return NULL if valid, raises error otherwise
#'
#' @keywords internal
validate_common_inputs <- function(M, K, algorithms, ensemble_folds, n_cores) {
  if (!is.numeric(M) || length(M) != 1 || M < 1 || M %% 1 != 0) {
    stop("M must be a positive integer")
  }
  if (!is.numeric(K) || length(K) != 1 || K < 2 || K %% 1 != 0) {
    stop("K must be an integer >= 2")
  }
  if (!is.character(algorithms) || length(algorithms) == 0) {
    stop("algorithms must be a non-empty character vector")
  }
  if (!is.numeric(ensemble_folds) || length(ensemble_folds) != 1 || ensemble_folds < 2 || ensemble_folds %% 1 != 0) {
    stop("ensemble_folds must be an integer >= 2")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1 || n_cores %% 1 != 0) {
    stop("n_cores must be a positive integer")
  }
}


#' Parse Column Name from Various Input Types
#'
#' @description
#' Internal helper to parse a column name specification that can be provided as:
#' \itemize{
#'   \item An unquoted symbol (e.g., \code{treatment = D})
#'   \item A quoted string (e.g., \code{treatment = "D"})
#'   \item A variable containing a string (e.g., \code{my_var <- "D"; treatment = my_var})
#' }
#'
#' @param expr The captured expression (from \code{substitute()})
#' @param env The environment to evaluate in (typically \code{parent.frame()})
#' @param arg_name Name of the argument (for error messages)
#' @param data Optional data.frame to validate column existence
#'
#' @return Character string with the column name
#'
#' @keywords internal
parse_column_name <- function(expr, env, arg_name = "variable", data = NULL) {
  # Handle NULL input
  if (is.null(expr)) {
    return(NULL)
  }
  
  # If it's already a string literal, use it directly
  if (is.character(expr)) {
    col_name <- expr
  } else if (is.symbol(expr)) {
    # It's a bare symbol - could be a column name or a variable containing one
    sym_name <- as.character(expr)
    
    # Check if this symbol exists in the environment and contains a string
    if (exists(sym_name, envir = env, inherits = TRUE)) {
      val <- get(sym_name, envir = env, inherits = TRUE)
      if (is.character(val) && length(val) == 1) {
        # It's a variable containing a column name string
        col_name <- val
      } else {
        # It's a symbol that should be interpreted as the column name itself
        col_name <- sym_name
      }
    } else {
      # Symbol doesn't exist as a variable, treat as column name
      col_name <- sym_name
    }
  } else if (is.call(expr)) {
    # It's an expression/call - try to evaluate it
    val <- tryCatch(
      eval(expr, envir = env),
      error = function(e) {
        stop(arg_name, " must be a column name (quoted or unquoted), not an expression: ", 
             deparse(expr))
      }
    )
    if (!is.character(val) || length(val) != 1) {
      stop(arg_name, " must evaluate to a single column name string")
    }
    col_name <- val
  } else {
    stop(arg_name, " must be a column name (quoted or unquoted)")
  }
  
  # Validate column exists in data if provided
  if (!is.null(data) && !col_name %in% names(data)) {
    stop(arg_name, " '", col_name, "' not found in data. ",
         "Available columns: ", paste(head(names(data), 10), collapse = ", "),
         if (length(names(data)) > 10) ", ..." else "")
  }
  
  col_name
}
