

# K-fold splitting
create_folds <- function(n, K = 5, M = 1, stratify_var = NULL, cluster_id = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Panel data: assign folds at the cluster (individual) level so all observations

  # from the same individual end up in the same fold.
  if (!is.null(cluster_id)) {
    if (length(cluster_id) != n) {
      stop("cluster_id must have the same length as the number of observations (", n, ")")
    }
    unique_ids <- unique(cluster_id)
    n_clusters <- length(unique_ids)
    if (n_clusters < K) {
      stop("Number of unique individuals (", n_clusters, ") must be at least K (", K, ")")
    }
    
    # Build a stratify variable at the cluster level if stratify_var is provided
    # Take the first value per cluster (stratify_var should be constant within cluster)
    if (!is.null(stratify_var)) {
      cluster_strat <- vapply(unique_ids, function(id) {
        stratify_var[which(cluster_id == id)[1]]
      }, stratify_var[1])
    }
    
    return(lapply(1:M, function(rep) {
      # Assign folds at the cluster level
      if (!is.null(stratify_var)) {
        cluster_folds <- rep(NA_integer_, n_clusters)
        for (level in unique(cluster_strat)) {
          idx <- which(cluster_strat == level)
          n_level <- length(idx)
          cluster_folds[idx] <- sample(rep(1:K, length.out = n_level))
        }
      } else {
        cluster_folds <- sample(rep(1:K, length.out = n_clusters))
      }
      # Map cluster-level folds back to observation level
      names(cluster_folds) <- unique_ids
      as.integer(cluster_folds[as.character(cluster_id)])
    }))
  }
  
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


#' Create Quantile Groups Using a Reference Population
#'
#' @description
#' Assigns all observations in \code{x} to quantile groups, where the group
#' boundaries (cutoffs) are determined by a reference subset of observations
#' indicated by \code{ref_mask}. Observations outside the reference subset are
#' assigned to groups based on where their values fall relative to the reference
#' cutoffs. When \code{ref_mask} is \code{NULL} or all \code{TRUE}, this is
#' equivalent to \code{\link{create_groups}}.
#'
#' This is used when the grouping population differs from the analysis population
#' (e.g., groups defined by the ML training sample, then applied to all data).
#'
#' @param x Numeric vector of all values to assign to groups
#' @param n_groups Integer, number of groups to create
#' @param ref_mask Logical vector (same length as \code{x}) indicating which
#'   observations form the reference population for computing group cutoffs.
#'   If \code{NULL} or all \code{TRUE}, equivalent to \code{create_groups(x, n_groups)}.
#'
#' @return Integer vector of group assignments (1 to n_groups)
#'
#' @keywords internal
create_groups_by_reference <- function(x, n_groups, ref_mask = NULL) {
  # If no reference mask or all are reference, delegate to standard create_groups
  if (is.null(ref_mask) || all(ref_mask)) {
    return(create_groups(x, n_groups))
  }

  x_ref <- x[ref_mask]
  n_ref <- length(x_ref)

  # Compute quantile breakpoints from the reference population
  probs <- seq(0, 1, length.out = n_groups + 1)
  breaks <- as.numeric(quantile(x_ref, probs = probs))

  # If the reference population is too small or has too many ties,
  # the quantile breaks will not be unique and cut() would fail.
  # Fall back to rank-based grouping on all observations.
  # (No warning here — .check_small_cells already warns once before the loop.)
  if (length(unique(breaks)) < n_groups + 1) {
    return(create_groups(x, n_groups))
  }

  # Extend endpoints so values outside the reference range are captured
  breaks[1] <- min(breaks[1], min(x)) - .Machine$double.eps
  breaks[n_groups + 1] <- max(breaks[n_groups + 1], max(x)) + .Machine$double.eps

  # Assign all observations using the reference-derived breakpoints
  groups <- as.integer(cut(x, breaks = breaks, include.lowest = TRUE, labels = FALSE))

  if (anyNA(groups)) {
    stop("Internal error in create_groups_by_reference: some observations ",
         "could not be assigned to a group. Please report this bug.")
  }

  groups
}


#' Regression from Formula with Robust Standard Errors
#'
#' @description
#' Fits a linear regression model from a formula (or formula string) with 
#' optional weights and computes HC1 robust standard errors. When
#' \code{cluster_id} is provided, computes cluster-robust (CR1) standard
#' errors using \code{sandwich::vcovCL} instead of heteroskedasticity-robust
#' (HC1) standard errors.
#'
#' @param formula Formula object or character string representing the regression formula
#' @param data Data frame containing the regression data
#' @param weights Optional numeric vector of weights for weighted regression.
#'   Must have the same length as the number of rows in data.
#' @param cluster_id Optional vector identifying clusters (e.g., individual IDs
#'   in panel data). When provided, cluster-robust standard errors are computed
#'   via \code{sandwich::vcovCL}. When \code{NULL} (default), heteroskedasticity-
#'   robust (HC1) standard errors are computed.
#'
#' @return List containing:
#' \itemize{
#'   \item coef: Coefficient test results with robust standard errors
#'   \item vcov: Robust variance-covariance matrix (HC1 or cluster-robust)
#'   \item model_matrix: Model matrix from the regression
#'   \item residuals: Regression residuals
#'   \item fitted: Fitted values
#' }
#'
#' @importFrom stats lm as.formula model.matrix residuals fitted
#' @importFrom sandwich vcovHC vcovCL
#' @importFrom lmtest coeftest
#'
#' @keywords internal
reg_from_formula <- function(formula, data, weights = NULL, cluster_id = NULL) {
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
  if (!is.null(cluster_id)) {
    # Cluster-robust SEs (CR1, analogous to HC1 with clustering)
    vcov_robust <- sandwich::vcovCL(fit, cluster = cluster_id, type = "HC1")
  } else {
    # Heteroskedasticity-robust SEs (HC1)
    vcov_robust <- sandwich::vcovHC(fit, type = "HC1")
  }
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
    weights = weights,
    cluster_id = cluster_id
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
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#'
#' # Run ensemble_hte in two separate sessions (e.g., on different servers)
#' # Session 1:
#' fit1 <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 5
#' )
#' saveRDS(fit1, "fit_session1.rds")
#'
#' # Session 2:
#' fit2 <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 5
#' )
#' saveRDS(fit2, "fit_session2.rds")
#'
#' # Later, combine the results:
#' fit1 <- readRDS("fit_session1.rds")
#' fit2 <- readRDS("fit_session2.rds")
#' fit_combined <- combine_ensembles(fit1, fit2)
#' print(fit_combined)  # Shows M = 10
#'
#' # Works the same for ensemble_pred
#' dat_pred <- microcredit[, c("bank_profits_pp", covars)]
#' pred1 <- ensemble_pred(
#'   bank_profits_pp ~ ., data = dat_pred,
#'   train_idx = microcredit$loan_size > 0 & microcredit$treat == 1, M = 3
#' )
#' pred2 <- ensemble_pred(
#'   bank_profits_pp ~ ., data = dat_pred,
#'   train_idx = microcredit$loan_size > 0 & microcredit$treat == 1, M = 3
#' )
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
    
    # Check individual_id consistency
    ref_has_id <- !is.null(ref$individual_id)
    fit_has_id <- !is.null(fit_i$individual_id)
    if (ref_has_id != fit_has_id) {
      stop(paste0("individual_id usage differs between object 1 (", 
                  if (ref_has_id) "present" else "absent", ") and object ", i, " (",
                  if (fit_has_id) "present" else "absent", ")"))
    }
    if (ref_has_id && fit_has_id) {
      if (length(ref$individual_id) != length(fit_i$individual_id) ||
          !all(ref$individual_id == fit_i$individual_id)) {
        stop(paste0("individual_id differs between object 1 and object ", i))
      }
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
#' @return The validated (and possibly deduplicated) algorithms vector, invisibly.
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
  if (anyDuplicated(algorithms) > 0) {
    dup <- unique(algorithms[duplicated(algorithms)])
    warning("Duplicate algorithms detected: ", paste(dup, collapse = ", "),
            ". Removing duplicates.")
    algorithms <- unique(algorithms)
  }
  if (!is.numeric(ensemble_folds) || length(ensemble_folds) != 1 || ensemble_folds < 2 || ensemble_folds %% 1 != 0) {
    stop("ensemble_folds must be an integer >= 2")
  }
  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1 || n_cores %% 1 != 0) {
    stop("n_cores must be a positive integer")
  }
  invisible(algorithms)
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


#' Check for small fold (or fold x restrict_by) cells and warn
#'
#' @description
#' Checks whether any fold (or fold x restrict_by) cell has fewer observations than
#' the requested number of groups across all repetitions. When this happens,
#' some groups will be empty in those cells and observations are assigned
#' deterministically to the lowest groups. This is not an error, but users
#' should be aware of it.
#'
#' @param fold_list A list of integer vectors of fold assignments (one per
#'   repetition), or a single integer vector. When a list is provided, all
#'   repetitions are scanned and the worst case is reported.
#' @param n_groups Number of groups requested
#' @param restrict_by Optional factor/integer vector for restricted ranking
#'   (same length as each fold vector). Assumed constant across repetitions
#'   (only fold assignments change).
#' @param func_name Character name of the calling function (for the warning
#'   message)
#'
#' @return Invisible NULL. Emits a warning if any cell is too small.
#' @keywords internal
.check_small_cells <- function(fold_list, n_groups, restrict_by = NULL, func_name = "") {
  # Accept a single vector for backward compatibility
  if (!is.list(fold_list)) fold_list <- list(fold_list)

  cell_type <- if (is.null(restrict_by)) "fold" else "fold \u00d7 restrict_by"
  worst_min <- Inf
  worst_n_small <- 0L
  worst_n_total <- 0L
  n_reps_affected <- 0L

  for (fold in fold_list) {
    cell_sizes <- if (is.null(restrict_by)) {
      as.integer(table(fold))
    } else {
      as.integer(table(fold, restrict_by))
    }
    small <- cell_sizes[cell_sizes > 0 & cell_sizes < n_groups]
    if (length(small) > 0) {
      n_reps_affected <- n_reps_affected + 1L
      this_min <- min(small)
      if (this_min < worst_min) {
        worst_min <- this_min
        worst_n_small <- length(small)
        worst_n_total <- sum(cell_sizes > 0)
      }
    }
  }

  if (n_reps_affected > 0) {
    M <- length(fold_list)

    # --- Build a compact, readable warning ---
    # Line 1: what happened
    line1 <- paste0(
      func_name, ": some ", cell_type,
      " cells are too small for ", n_groups, " groups."
    )

    # Line 2: scope — how many reps affected
    if (M > 1) {
      pct <- round(100 * n_reps_affected / M)
      line2 <- paste0(
        "  Affected repetitions: ", n_reps_affected, "/", M,
        " (", pct, "%)"
      )
    } else {
      line2 <- NULL
    }

    # Line 3: severity — smallest cell
    line3 <- paste0(
      "  Smallest cell: ", worst_min,
      " obs (need >= ", n_groups, " for ", n_groups, " groups)"
    )

    # Line 4: what to do
    fixes <- "reduce n_groups or use fewer folds (K)"
    if (!is.null(restrict_by)) {
      fixes <- paste0(fixes, " or coarser restrict_by groups")
    }
    line4 <- paste0("  To fix: ", fixes, ".")

    warning(
      paste(c(line1, line2, line3, line4), collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' Resolve group_on to a reference mask for group formation
#'
#' @description
#' Determines which observations form the reference population for computing
#' group cutoffs, based on the \code{group_on} argument and the fit object.
#'
#' \itemize{
#'   \item \code{"auto"}: Use the ML training population. For
#'     \code{ensemble_hte_fit} this is always all observations. For
#'     \code{ensemble_pred_fit} with \code{train_idx}, it is the training subset.
#'   \item \code{"all"}: Always use all observations for group formation.
#'   \item \code{"analysis"}: Use whatever observations are being analyzed (the
#'     \code{analysis_idx}).
#' }
#'
#' @param group_on Character: \code{"auto"}, \code{"all"}, or \code{"analysis"}
#' @param analysis_idx Logical vector (length n) indicating which observations
#'   are being analyzed.
#' @param train_idx Logical vector (length n) or NULL. The training index from
#'   an \code{ensemble_pred_fit}, if available.
#'
#' @return A logical vector (length n) indicating reference observations for
#'   group formation, or \code{NULL} if all observations should be used
#'   (equivalent to all \code{TRUE}).
#' @keywords internal
.resolve_group_ref_idx <- function(group_on, analysis_idx, train_idx = NULL) {
  # When all observations are being analyzed, reference mask doesn't matter
  if (all(analysis_idx)) return(NULL)

  switch(group_on,
    "auto" = {
      # Default: use all observations for group formation so that group
      # definitions remain consistent across different analyses (e.g., gates
      # on one outcome and gavs on another will have identical groups).
      NULL
    },
    "all" = {
      # All observations — no reference mask needed
      NULL
    },
    "analysis" = {
      # Same as the analysis subset
      analysis_idx
    },
    stop("group_on must be 'auto', 'all', or 'analysis'")
  )
}
