
#' Internal function to process subset argument for analysis functions
#' 
#' @description
#' Processes the subset argument and returns a logical vector indicating which
#' observations to use. Works for both ensemble_hte_fit and ensemble_pred_fit.
#' 
#' @param subset The subset argument from the analysis function. Can be:
#'   \itemize{
#'     \item NULL: smart default based on fit type and outcome
#'     \item "train": use training observations only (pred fits only)
#'     \item "all": use all observations
#'     \item Logical vector: TRUE/FALSE for each observation
#'     \item Integer vector: indices of observations to use
#'   }
#' @param n Total number of observations
#' @param fit_type "hte" or "pred"
#' @param has_train_idx Whether the ensemble_fit used subset training
#' @param train_idx The training indices from ensemble_fit (if available)
#' @param using_default_outcome Whether the default outcome is being used
#' @param Y The outcome vector (to check for NAs)
#' 
#' @return A logical vector of length n indicating which observations to use
#' @keywords internal
.process_subset <- function(subset, n, fit_type, has_train_idx, train_idx,
                            using_default_outcome, Y) {
  
  if (is.null(subset)) {
    # Smart default
    if (fit_type == "pred" && using_default_outcome && has_train_idx) {
      # For pred fits with default outcome and subset training, use training obs
      use_idx <- train_idx
    } else {
      # Otherwise use all observations
      use_idx <- rep(TRUE, n)
    }
  } else if (is.character(subset) && length(subset) == 1) {
    if (subset == "train") {
      if (fit_type == "hte") {
        stop("subset = 'train' is not valid for ensemble_hte_fit. ",
             "For HTE fits, use a logical/integer vector to specify a subset.")
      }
      if (!has_train_idx) {
        stop("subset = 'train' specified but no train_idx was used in ensemble_pred()")
      }
      use_idx <- train_idx
    } else if (subset == "all") {
      use_idx <- rep(TRUE, n)
      # Warn if using default outcome with NAs
      if (fit_type == "pred" && using_default_outcome && has_train_idx && any(is.na(Y))) {
        warning("Using all observations but outcome has NA values for non-training observations")
      }
    } else {
      stop("subset must be NULL, 'train', 'all', a logical vector, or an integer vector")
    }
  } else if (is.logical(subset)) {
    if (length(subset) != n) {
      stop(paste0("subset logical vector has length ", length(subset), 
                  " but data has ", n, " rows"))
    }
    use_idx <- subset
  } else if (is.numeric(subset) && all(subset == floor(subset))) {
    # Integer indices
    if (any(subset < 1) || any(subset > n)) {
      stop(paste0("subset indices must be between 1 and ", n))
    }
    use_idx <- rep(FALSE, n)
    use_idx[subset] <- TRUE
  } else {
    stop("subset must be NULL, 'train', 'all', a logical vector, or an integer vector of indices")
  }
  
  use_idx
}


#' Internal function to print informative message about subset usage
#' 
#' @description
#' Prints a message when the ensemble_fit was trained on a subset or when
#' the analysis function is evaluated on a subset. The goal is to help users
#' avoid mistakes when specifying subsets.
#' 
#' @param func_name Name of the analysis function (e.g., "GAVS", "GATES", "BLP")
#' @param n_total Total number of observations in the data
#' @param n_fit Number of observations used to fit the ML model (n_train or n_total)
#' @param n_eval Number of observations used for evaluation
#' @param fit_type "hte" or "pred"
#' @param has_train_idx Whether the ensemble_fit used subset training
#' 
#' @return NULL (prints message as side effect)
#' @keywords internal
.print_subset_message <- function(func_name, n_total, n_fit, n_eval, fit_type, has_train_idx) {
  
  # Only print message if there's a discrepancy worth noting
  # Don't print if all observations are used for both fitting and evaluation
  if (n_fit == n_total && n_eval == n_total) {
    return(invisible(NULL))
  }
  
  msg_parts <- character(0)
  
  if (n_fit < n_total) {
    msg_parts <- c(msg_parts, 
      sprintf("ML model trained on %d of %d observations", n_fit, n_total))
  }
  
  if (n_eval < n_total || n_eval != n_fit) {
    msg_parts <- c(msg_parts,
      sprintf("%s evaluated on %d observations", func_name, n_eval))
  }
  
  if (length(msg_parts) > 0) {
    message(paste0("Note: ", paste(msg_parts, collapse = "; "), "."))
  }
  
  invisible(NULL)
}


#' Internal function to compute GATES for a single repetition
#' 
#' @param Y Outcome vector
#' @param D Treatment vector
#' @param prop_score Propensity score vector
#' @param weight Weight vector
#' @param predicted_values Predicted values (ITE or Y) vector for this repetition
#' @param fold Fold assignment vector for this repetition
#' @param n_groups Number of groups for GATES
#' @param strata Optional strata indicator for restricted ranking (factor/integer vector)
#' @param controls Optional control variables (character vector of column names)
#' @param control_data Optional data.table/data.frame with control variables
#' 
#' @return List with group_estimates, tests (top_bottom, all, top_all), and reg object
#' @keywords internal
.gates_single <- function(Y, D, prop_score, weight, predicted_values, fold, 
                          n_groups = 5, strata = NULL, controls = NULL, 
                          control_data = NULL) {
  
  # Build data.table for GATES regression
  dt_gates <- data.table(
    Y = Y,
    D = D,
    prop_score = prop_score,
    weight = weight,
    predicted_values = predicted_values,
    fold = fold
  )
  
  # Create groups within each fold (and optionally within strata for restricted ranking)
  if (is.null(strata)) {
    dt_gates[, group := as.integer(create_groups(predicted_values, n_groups)), by = fold]
  } else {
    dt_gates[, strata := strata]
    dt_gates[, group := as.integer(create_groups(predicted_values, n_groups)), by = .(fold, strata)]
  }
  
  # Create dummy variables for groups and interact with (D - prop_score)
  group_cols <- paste0("group_", 1:n_groups)
  for (g_idx in 1:n_groups) {
    col_name <- paste0("group_", g_idx)
    dt_gates[, (col_name) := as.integer(group == g_idx) * (D - prop_score)]
  }
  
  # Add control variables if provided
  if (!is.null(controls) && !is.null(control_data)) {
    dt_gates <- cbind(dt_gates, control_data[, ..controls])
  }
  
  # Build regression formula for group estimates
  # GATES regression: Y ~ 1 + group_1 + group_2 + ... + group_n (+ controls)
  formula_str <- paste("Y ~", paste(group_cols, collapse = " + "))
  if (!is.null(controls)) {
    formula_str <- paste(formula_str, "+", paste(controls, collapse = " + "))
  }
  
  # Run regression with robust SEs
  reg <- reg_from_formula(formula_str, dt_gates, weights = dt_gates$weight)
  coef_gates <- reg$coef
  vcov_gates <- reg$vcov
  
  # Extract estimates and SEs for each group
  group_estimates <- data.table(
    group = 1:n_groups,
    estimate = coef_gates[group_cols, 1],
    se = coef_gates[group_cols, 2]
  )
  
  # Define column names
  top_col <- paste0("group_", n_groups)
  bottom_col <- "group_1"
  
  # 1. Top-Bottom difference (group n - group 1)
  top_bottom <- .compute_coef_diff(coef_gates, vcov_gates, top_col, bottom_col)
  
  # 2. All: weighted average of all groups (= ATE)
  # Use group proportions as weights (for equal-sized groups, this is 1/n_groups each)
  n_obs <- nrow(dt_gates)
  group_props <- vapply(1:n_groups, function(g) sum(dt_gates$group == g) / n_obs, numeric(1))
  all_estimate <- sum(group_props * coef_gates[group_cols, 1])
  # SE: sqrt(w' * Vcov * w) where w is the weight vector
  all_se <- sqrt(as.numeric(t(group_props) %*% vcov_gates[group_cols, group_cols] %*% group_props))
  all_test <- data.table(estimate = all_estimate, se = all_se, weights = list(group_props))
  
  # 3. Top-All: top group minus weighted average of all groups
  # This is equivalent to top - sum(w_i * group_i)
  # = (1 - w_n) * top - sum_{i<n}(w_i * group_i)
  top_all_weights <- -group_props
  top_all_weights[n_groups] <- 1 - group_props[n_groups]
  top_all_estimate <- sum(top_all_weights * coef_gates[group_cols, 1])
  top_all_se <- sqrt(as.numeric(t(top_all_weights) %*% vcov_gates[group_cols, group_cols] %*% top_all_weights))
  top_all <- data.table(estimate = top_all_estimate, se = top_all_se)
  
  list(
    group_estimates = group_estimates,
    top_bottom = top_bottom,
    all = all_test,
    top_all = top_all,
    reg = reg  # Include full regression object for cross-covariance calculations
  )
}


#' Compute difference between two coefficients with proper SE
#' 
#' @param coef Coefficient matrix from regression
#' @param vcov Variance-covariance matrix
#' @param coef1 Name of first coefficient
#' @param coef2 Name of second coefficient
#' 
#' @return data.table with estimate and se
#' @keywords internal
.compute_coef_diff <- function(coef, vcov, coef1, coef2) {
  diff_estimate <- coef[coef1, 1] - coef[coef2, 1]
  diff_se <- sqrt(vcov[coef1, coef1] + vcov[coef2, coef2] - 2 * vcov[coef1, coef2])
  data.table(estimate = diff_estimate, se = diff_se)
}


#' Compute GATES (Group Average Treatment Effects)
#' 
#' @description
#' Computes Group Average Treatment Effects (GATES), an estimand introduced by
#' Chernozhukov et al. (2025) to measure average treatment effects for groups
#' defined by predicted value quantiles.
#' 
#' This function works with both \code{ensemble_hte_fit} (groups by predicted ITE)
#' and \code{ensemble_pred_fit} (groups by predicted Y). For prediction fits, a
#' treatment variable must be specified.
#' 
#' This function implements the multiple-split estimation strategy developed by
#' Fava (2025), which combines predictions from multiple machine learning algorithms
#' into an ensemble and averages GATES estimates across M repetitions of K-fold
#' cross-fitting to improve statistical power.
#' 
#' The method uses weighted least squares regression of Y on group indicators
#' interacted with (D - propensity_score), with weights equal to
#' 1 / (propensity_score * (1 - propensity_score)). Groups are formed by
#' ranking predicted values within each fold.
#' Robust (HC1) standard errors are computed.
#' 
#' @references
#' Chernozhukov, V., Demirer, M., Duflo, E., & Fernández-Val, I. (2025).
#' Fisher–Schultz Lecture: Generic Machine Learning Inference on Heterogeneous
#' Treatment Effects in Randomized Experiments, with an Application to Immunization
#' in India. \emph{Econometrica}, 93(4), 1121-1164.
#' 
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class \code{ensemble_hte_fit} from \code{ensemble_hte()}
#'   or \code{ensemble_pred_fit} from \code{ensemble_pred()}.
#' @param n_groups Number of groups to divide the sample into (default: 3)
#' @param outcome Either:
#'   \itemize{
#'     \item NULL (default): uses the same outcome as in the ensemble function
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric vector: custom outcome variable (must have same length as data)
#'   }
#'   This allows computing GATES for a different outcome than the one used for prediction.
#' @param treatment For \code{ensemble_pred_fit} only. The treatment variable:
#'   \itemize{
#'     \item Character string: column name in the \code{data} used in \code{ensemble_pred()}
#'     \item Numeric vector: binary treatment variable (must have same length as data)
#'   }
#'   Ignored for \code{ensemble_hte_fit} (uses the treatment from the fit).
#' @param prop_score For \code{ensemble_pred_fit} only. Propensity score:
#'   \itemize{
#'     \item NULL (default): estimated as mean of treatment variable
#'     \item Numeric value: constant propensity score for all observations
#'     \item Numeric vector: observation-specific propensity scores
#'   }
#'   For \code{ensemble_hte_fit}, uses the propensity score from the fit.
#' @param controls Character vector of control variable names to include in regression.
#'   These must be column names present in the \code{data} argument used when calling
#'   the ensemble function.
#' @param strata Optional. Stratification variable for restricted ranking:
#'   \itemize{
#'     \item NULL (default): unrestricted ranking across full sample within folds
#'     \item Character string: column name in the \code{data} for stratified ranking
#'     \item Numeric/factor vector: strata indicator (must have same length as data)
#'   }
#'   When specified, predicted values are ranked within each stratum (and fold),
#'   rather than across the full sample.
#' @param subset Which observations to use for the GATES analysis. Options:
#'   \itemize{
#'     \item NULL (default): uses all observations (or training obs for \code{ensemble_pred_fit}
#'       when using the default outcome with subset training)
#'     \item \code{"train"}: uses only training observations (for \code{ensemble_pred_fit} only)
#'     \item \code{"all"}: explicitly uses all observations
#'     \item Logical vector: TRUE/FALSE for each observation (must have same length as data)
#'     \item Integer vector: indices of observations to include (1-indexed)
#'   }
#'   This allows evaluating GATES on a subset of observations. Note that GATES 
#'   requires observations from both treatment and control groups in the subset.
#' 
#' @return An object of class \code{gates_results} containing:
#' \itemize{
#'   \item estimates: data.table with GATES estimates averaged across repetitions
#'   \item top_bottom: data.table with the top-bottom difference test
#'   \item all: data.table with the average treatment effect (weighted avg of all groups)
#'   \item top_all: data.table with the top minus average test
#'   \item n_groups: number of groups used
#'   \item outcome: the outcome variable used for GATES
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item strata: the stratification variable used (if any)
#'   \item controls: control variables used (if any)
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#' 
#' @examples
#' \donttest{
#' # Simulate data
#' set.seed(123)
#' n <- 500
#' X1 <- rnorm(n); X2 <- rnorm(n)
#' D <- rbinom(n, 1, 0.5)
#' Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n)
#' data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
#'
#' fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
#'                     algorithms = c("lm", "grf"), M = 3, K = 3)
#' result <- gates(fit, n_groups = 3)
#' print(result)
#' plot(result)
#' }
#'
#' @export
gates <- function(ensemble_fit, n_groups = 3, outcome = NULL, treatment = NULL,
                  prop_score = NULL, controls = NULL, strata = NULL, subset = NULL) {
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$ite[[m]])
    fit_type <- "hte"
    # For HTE fits, treatment info comes from the fit
    D <- ensemble_fit$D
    ps <- ensemble_fit$prop_score
    weight <- ensemble_fit$weights
  } else if (inherits(ensemble_fit, "ensemble_pred_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$predictions[[m]])
    fit_type <- "pred"
    # For pred fits, treatment must be specified
    if (is.null(treatment)) {
      stop("treatment must be specified when using ensemble_pred_fit")
    }
  } else {
    stop("ensemble_fit must be of class 'ensemble_hte_fit' or 'ensemble_pred_fit'")
  }
  
  cl <- match.call()
  n <- ensemble_fit$n
  
  # Get targeted outcome from ensemble_fit
  targeted_outcome <- all.vars(ensemble_fit$formula)[1]
  
  # Determine which outcome to use for GATES
  if (is.null(outcome)) {
    Y <- ensemble_fit$Y
    outcome_var <- targeted_outcome
  } else if (is.character(outcome)) {
    if (length(outcome) != 1) {
      stop("outcome must be a single column name or a numeric vector")
    }
    if (!outcome %in% names(ensemble_fit$data)) {
      stop(paste0("outcome '", outcome, "' not found in the data"))
    }
    Y <- ensemble_fit$data[[outcome]]
    outcome_var <- outcome
  } else if (is.numeric(outcome)) {
    if (length(outcome) != n) {
      stop(paste0("outcome vector has length ", length(outcome), 
                  " but data has ", n, " rows"))
    }
    Y <- outcome
    outcome_var <- "custom_outcome"
  } else {
    stop("outcome must be NULL, a character string, or a numeric vector")
  }
  
  # For pred fits, process treatment and propensity score
  if (fit_type == "pred") {
    # Process treatment variable
    if (is.character(treatment)) {
      if (length(treatment) != 1) {
        stop("treatment must be a single column name or a numeric vector")
      }
      if (!treatment %in% names(ensemble_fit$data)) {
        stop(paste0("treatment '", treatment, "' not found in the data"))
      }
      D <- ensemble_fit$data[[treatment]]
    } else if (is.numeric(treatment)) {
      if (length(treatment) != n) {
        stop(paste0("treatment vector has length ", length(treatment), 
                    " but data has ", n, " rows"))
      }
      D <- treatment
    } else {
      stop("treatment must be a character string or a numeric vector")
    }
    
    # Validate treatment is binary
    unique_D <- unique(D[!is.na(D)])
    if (!all(unique_D %in% c(0, 1))) {
      stop("treatment must be binary (0/1)")
    }
    
    # Process propensity score
    if (is.null(prop_score)) {
      # Estimate as mean of treatment
      ps <- rep(mean(D), n)
    } else if (length(prop_score) == 1) {
      ps <- rep(prop_score, n)
    } else if (length(prop_score) == n) {
      ps <- prop_score
    } else {
      stop("prop_score must be NULL, a single value, or a vector of length n")
    }
    
    # Compute weights
    weight <- 1 / (ps * (1 - ps))
  }
  
  # Process strata argument
  strata_vec <- NULL
  strata_var <- NULL
  if (!is.null(strata)) {
    if (is.character(strata)) {
      if (length(strata) != 1) {
        stop("strata must be a single column name or a vector")
      }
      if (!strata %in% names(ensemble_fit$data)) {
        stop(paste0("strata '", strata, "' not found in the data"))
      }
      strata_vec <- ensemble_fit$data[[strata]]
      strata_var <- strata
    } else if (is.numeric(strata) || is.factor(strata)) {
      if (length(strata) != n) {
        stop(paste0("strata vector has length ", length(strata), 
                    " but data has ", n, " rows"))
      }
      strata_vec <- strata
      strata_var <- "custom_strata"
    } else {
      stop("strata must be NULL, a character string, or a numeric/factor vector")
    }
    # Convert to factor if not already
    if (!is.factor(strata_vec)) {
      strata_vec <- as.factor(strata_vec)
    }
  }
  
  # Check if outcome is the default (targeted) outcome
  using_default_outcome <- is.null(outcome)
  
  # Check for train_idx in ensemble_fit
  has_train_idx <- !is.null(ensemble_fit$train_idx)
  
  # Determine which observations to use
  train_idx <- if (has_train_idx) ensemble_fit$train_idx else NULL
  use_idx <- .process_subset(
    subset = subset,
    n = n,
    fit_type = fit_type,
    has_train_idx = has_train_idx,
    train_idx = train_idx,
    using_default_outcome = using_default_outcome,
    Y = Y
  )
  
  # Validate no NAs in outcome for used observations
  if (any(is.na(Y[use_idx]))) {
    stop("outcome has NA values for the observations being used")
  }
  
  # Validate no NAs in treatment for used observations
  if (any(is.na(D[use_idx]))) {
    stop("treatment has NA values for the observations being used")
  }
  
  # Validate no NAs in strata for used observations (if strata provided)
  if (!is.null(strata_vec) && any(is.na(strata_vec[use_idx]))) {
    stop("strata has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (fit_type == "pred" && has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("GATES", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Get control data if controls specified (subset to used observations)
  control_data <- if (!is.null(controls)) ensemble_fit$data[use_idx, ] else NULL
  
  # Subset strata if provided
  strata_used <- if (!is.null(strata_vec)) strata_vec[use_idx] else NULL
  
  # Compute GATES for each repetition
  gates_by_rep <- lapply(1:M, function(m) {
    .gates_single(
      Y = Y[use_idx],
      D = D[use_idx],
      prop_score = ps[use_idx],
      weight = weight[use_idx],
      predicted_values = predictions_list[[m]][use_idx],
      fold = splits[[m]][use_idx],
      n_groups = n_groups,
      strata = strata_used,
      controls = controls,
      control_data = control_data
    )
  })
  
  # Combine group estimates across repetitions (mean method from Fava et al.)
  all_group_estimates <- rbindlist(lapply(gates_by_rep, `[[`, "group_estimates"), idcol = "repetition")
  
  # Aggregate across repetitions: mean of estimates, mean of SE
  combined <- all_group_estimates[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  
  # Add t-values and p-values for combined estimates
  combined[, `:=`(
    t_value = estimate / se,
    p_value = 2 * pnorm(-abs(estimate / se))
  )]
  
  # Helper function to aggregate test results
  .aggregate_test <- function(test_list) {
    all_tests <- rbindlist(test_list, idcol = "repetition")
    result <- all_tests[, .(estimate = mean(estimate), se = mean(se))]
    result[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
    result
  }
  
  # Combine test estimates across repetitions
  top_bottom_combined <- .aggregate_test(lapply(gates_by_rep, `[[`, "top_bottom"))
  all_combined <- .aggregate_test(lapply(gates_by_rep, `[[`, "all"))
  top_all_combined <- .aggregate_test(lapply(gates_by_rep, `[[`, "top_all"))
  
  # Construct gates_results object
  structure(
    list(
      estimates = combined,
      top_bottom = top_bottom_combined,
      all = all_combined,
      top_all = top_all_combined,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      strata = strata_var,
      controls = controls,
      M = M,
      call = cl
    ),
    class = "gates_results"
  )
}


#' Internal function to compute BLP for a single repetition
#' 
#' @param Y Outcome vector
#' @param D Treatment vector
#' @param prop_score Propensity score vector
#' @param weight Weight vector
#' @param predicted_ite Predicted ITE vector for this repetition
#' @param controls Optional control variables (character vector of column names)
#' @param control_data Optional data.table/data.frame with control variables
#' 
#' @return data.table with BLP estimates (beta1 for ATE, beta2 for heterogeneity)
#' @keywords internal
.blp_single <- function(Y, D, prop_score, weight, predicted_ite, 
                        controls = NULL, control_data = NULL) {
  
  # Build data.table for BLP regression
  dt_blp <- data.table(
    Y = Y,
    D = D,
    prop_score = prop_score,
    weight = weight,
    predicted_ite = predicted_ite
  )
  
  # Create BLP regressors
  # W1 = (D - prop_score): coefficient is ATE
  # W2 = (D - prop_score) * (ite - mean(ite)): coefficient tests heterogeneity
  dt_blp[, W1 := D - prop_score]
  dt_blp[, W2 := (D - prop_score) * (predicted_ite - mean(predicted_ite))]
  
  # Add control variables if provided
  if (!is.null(controls) && !is.null(control_data)) {
    dt_blp <- cbind(dt_blp, control_data[, ..controls])
  }
  
  # Build regression formula
  # BLP regression: Y ~ 1 + W1 + W2 (+ controls)
  # The intercept captures the baseline outcome
  formula_str <- "Y ~ W1 + W2"
  if (!is.null(controls)) {
    formula_str <- paste(formula_str, "+", paste(controls, collapse = " + "))
  }
  
  # Run regression with robust SEs
  reg <- reg_from_formula(formula_str, dt_blp, weights = dt_blp$weight)
  coef_blp <- reg$coef
  vcov_blp <- reg$vcov
  
  # Extract estimates and SEs
  blp_estimates <- data.table(
    term = c("beta1", "beta2"),
    estimate = coef_blp[c("W1", "W2"), 1],
    se = coef_blp[c("W1", "W2"), 2]
  )
  
  list(
    estimates = blp_estimates,
    vcov = vcov_blp
  )
}


#' Compute BLP (Best Linear Predictor) of CATE
#' 
#' @description
#' Computes the Best Linear Predictor (BLP) of the Conditional Average Treatment
#' Effect (CATE), an estimand introduced by Chernozhukov et al. (2025) to test
#' whether machine learning predictions capture meaningful treatment effect heterogeneity.
#' 
#' This function works with both \code{ensemble_hte_fit} (uses predicted ITE)
#' and \code{ensemble_pred_fit} (uses predicted Y). For prediction fits, a
#' treatment variable must be specified.
#' 
#' This function implements the multiple-split estimation strategy developed by
#' Fava (2025), which combines predictions from multiple machine learning algorithms
#' into an ensemble and averages BLP estimates across M repetitions of K-fold
#' cross-validation to improve statistical power.
#' 
#' The method uses weighted least squares regression of Y on:
#' \itemize{
#'   \item W1 = (D - propensity_score): coefficient (beta1) estimates the ATE
#'   \item W2 = (D - propensity_score) * (predicted - mean(predicted)): coefficient (beta2) tests heterogeneity
#' }
#' Weights are set to 1 / (propensity_score * (1 - propensity_score)).
#' Robust (HC1) standard errors are computed.
#' 
#' A significant beta2 indicates that the ML predictions capture meaningful
#' treatment effect heterogeneity.
#' 
#' @references
#' Chernozhukov, V., Demirer, M., Duflo, E., & Fernández-Val, I. (2025).
#' Fisher–Schultz Lecture: Generic Machine Learning Inference on Heterogeneous
#' Treatment Effects in Randomized Experiments, with an Application to Immunization
#' in India. \emph{Econometrica}, 93(4), 1121-1164.
#' 
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class \code{ensemble_hte_fit} from \code{ensemble_hte()}
#'   or \code{ensemble_pred_fit} from \code{ensemble_pred()}.
#' @param outcome Either:
#'   \itemize{
#'     \item NULL (default): uses the same outcome as in the ensemble function
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric vector: custom outcome variable (must have same length as data)
#'   }
#'   This allows computing BLP for a different outcome than the one used for prediction.
#' @param treatment For \code{ensemble_pred_fit} only. The treatment variable:
#'   \itemize{
#'     \item Character string: column name in the \code{data} used in \code{ensemble_pred()}
#'     \item Numeric vector: binary treatment variable (must have same length as data)
#'   }
#'   Ignored for \code{ensemble_hte_fit} (uses the treatment from the fit).
#' @param prop_score For \code{ensemble_pred_fit} only. Propensity score:
#'   \itemize{
#'     \item NULL (default): estimated as mean of treatment variable
#'     \item Numeric value: constant propensity score for all observations
#'     \item Numeric vector: observation-specific propensity scores
#'   }
#'   For \code{ensemble_hte_fit}, uses the propensity score from the fit.
#' @param controls Character vector of control variable names to include in regression.
#'   These must be column names present in the \code{data} argument used when calling
#'   the ensemble function.
#' @param subset Which observations to use for the BLP analysis. Options:
#'   \itemize{
#'     \item NULL (default): uses all observations (or training obs for \code{ensemble_pred_fit}
#'       when using the default outcome with subset training)
#'     \item \code{"train"}: uses only training observations (for \code{ensemble_pred_fit} only)
#'     \item \code{"all"}: explicitly uses all observations
#'     \item Logical vector: TRUE/FALSE for each observation (must have same length as data)
#'     \item Integer vector: indices of observations to include (1-indexed)
#'   }
#'   This allows evaluating BLP on a subset of observations. Note that BLP 
#'   requires observations from both treatment and control groups in the subset.
#' 
#' @return An object of class `blp_results` containing:
#' \itemize{
#'   \item estimates: data.table with BLP estimates averaged across repetitions
#'   \item outcome: outcome variable used
#'   \item targeted_outcome: original outcome from ensemble fitting
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item controls: control variables used (if any)
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 500
#' X1 <- rnorm(n); X2 <- rnorm(n)
#' D <- rbinom(n, 1, 0.5)
#' Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n)
#' data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
#'
#' fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
#'                     algorithms = c("lm", "grf"), M = 3, K = 3)
#' result <- blp(fit)
#' print(result)
#' }
#'
#' @export
blp <- function(ensemble_fit, outcome = NULL, treatment = NULL, 
                prop_score = NULL, controls = NULL, subset = NULL) {
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$ite[[m]])
    fit_type <- "hte"
    # For HTE fits, treatment info comes from the fit
    D <- ensemble_fit$D
    ps <- ensemble_fit$prop_score
    weight <- ensemble_fit$weights
  } else if (inherits(ensemble_fit, "ensemble_pred_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$predictions[[m]])
    fit_type <- "pred"
    # For pred fits, treatment must be specified
    if (is.null(treatment)) {
      stop("treatment must be specified when using ensemble_pred_fit")
    }
  } else {
    stop("ensemble_fit must be of class 'ensemble_hte_fit' or 'ensemble_pred_fit'")
  }
  
  cl <- match.call()
  n <- ensemble_fit$n
  
  # Get targeted outcome from ensemble_fit
  targeted_outcome <- all.vars(ensemble_fit$formula)[1]
  
  # Determine which outcome to use for BLP
  if (is.null(outcome)) {
    Y <- ensemble_fit$Y
    outcome_var <- targeted_outcome
  } else if (is.character(outcome)) {
    if (length(outcome) != 1) {
      stop("outcome must be a single column name or a numeric vector")
    }
    if (!outcome %in% names(ensemble_fit$data)) {
      stop(paste0("outcome '", outcome, "' not found in the data"))
    }
    Y <- ensemble_fit$data[[outcome]]
    outcome_var <- outcome
  } else if (is.numeric(outcome)) {
    if (length(outcome) != n) {
      stop(paste0("outcome vector has length ", length(outcome), 
                  " but data has ", n, " rows"))
    }
    Y <- outcome
    outcome_var <- "custom_outcome"
  } else {
    stop("outcome must be NULL, a character string, or a numeric vector")
  }
  
  # For pred fits, process treatment and propensity score
  if (fit_type == "pred") {
    # Process treatment variable
    if (is.character(treatment)) {
      if (length(treatment) != 1) {
        stop("treatment must be a single column name or a numeric vector")
      }
      if (!treatment %in% names(ensemble_fit$data)) {
        stop(paste0("treatment '", treatment, "' not found in the data"))
      }
      D <- ensemble_fit$data[[treatment]]
    } else if (is.numeric(treatment)) {
      if (length(treatment) != n) {
        stop(paste0("treatment vector has length ", length(treatment), 
                    " but data has ", n, " rows"))
      }
      D <- treatment
    } else {
      stop("treatment must be a character string or a numeric vector")
    }
    
    # Validate treatment is binary
    unique_D <- unique(D[!is.na(D)])
    if (!all(unique_D %in% c(0, 1))) {
      stop("treatment must be binary (0/1)")
    }
    
    # Process propensity score
    if (is.null(prop_score)) {
      # Estimate as mean of treatment
      ps <- rep(mean(D), n)
    } else if (length(prop_score) == 1) {
      ps <- rep(prop_score, n)
    } else if (length(prop_score) == n) {
      ps <- prop_score
    } else {
      stop("prop_score must be NULL, a single value, or a vector of length n")
    }
    
    # Compute weights
    weight <- 1 / (ps * (1 - ps))
  }
  
  # Check if outcome is the default (targeted) outcome
  using_default_outcome <- is.null(outcome)
  
  # Check for train_idx in ensemble_fit
  has_train_idx <- !is.null(ensemble_fit$train_idx)
  
  # Determine which observations to use
  train_idx <- if (has_train_idx) ensemble_fit$train_idx else NULL
  use_idx <- .process_subset(
    subset = subset,
    n = n,
    fit_type = fit_type,
    has_train_idx = has_train_idx,
    train_idx = train_idx,
    using_default_outcome = using_default_outcome,
    Y = Y
  )
  
  # Validate no NAs in outcome for used observations
  if (any(is.na(Y[use_idx]))) {
    stop("outcome has NA values for the observations being used")
  }
  
  # Validate no NAs in treatment for used observations
  if (any(is.na(D[use_idx]))) {
    stop("treatment has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (fit_type == "pred" && has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("BLP", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  M <- ensemble_fit$M
  
  # Get control data if controls specified (subset to used observations)
  control_data <- if (!is.null(controls)) ensemble_fit$data[use_idx, ] else NULL
  
  # Compute BLP for each repetition
  blp_by_rep <- lapply(1:M, function(m) {
    .blp_single(
      Y = Y[use_idx],
      D = D[use_idx],
      prop_score = ps[use_idx],
      weight = weight[use_idx],
      predicted_ite = predictions_list[[m]][use_idx],
      controls = controls,
      control_data = control_data
    )
  })
  
  # Combine results across repetitions (average estimates and SEs)
  all_estimates <- rbindlist(lapply(blp_by_rep, `[[`, "estimates"), idcol = "repetition")
  
  combined <- all_estimates[, .(
    estimate = mean(estimate),
    se = mean(se)
  ), by = term]
  
  # Add t-values and p-values
  combined[, `:=`(
    t_value = estimate / se,
    p_value = 2 * pnorm(-abs(estimate / se))
  )]
  
  # Construct blp_results object
  structure(
    list(
      estimates = combined,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      controls = controls,
      M = M,
      call = cl
    ),
    class = "blp_results"
  )
}


#' Print method for blp_results objects
#' @param x An object of class \code{blp_results} from \code{blp()}
#' @param ... Additional arguments (currently unused)
#' @export
print.blp_results <- function(x, ...) {
  cat("BLP Results (Best Linear Predictor of CATE)\n")
  cat("============================================\n\n")
  
  # Show fit type
  fit_label <- if (!is.null(x$fit_type) && x$fit_type == "hte") {
    "HTE (ensemble_hte)"
  } else if (!is.null(x$fit_type) && x$fit_type == "pred") {
    "Prediction (ensemble_pred)"
  } else {
    "HTE (ensemble_hte)"  # Default for backwards compatibility
  }
  cat("Fit type: ", fit_label, "\n", sep = "")
  
  # Show outcome information
  cat("Outcome analyzed: ", x$outcome, "\n", sep = "")
  if (x$outcome != x$targeted_outcome) {
    pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "ITE"
    cat("  (Predictions based on: ", x$targeted_outcome, " ", pred_label, ")\n", sep = "")
  }
  cat("Repetitions: ", x$M, "\n", sep = "")
  if (!is.null(x$controls)) {
    cat("Controls: ", paste(x$controls, collapse = ", "), "\n", sep = "")
  }
  
  cat("\nCoefficients:\n")
  cat("  beta1 (ATE): Average Treatment Effect\n")
  cat("  beta2 (HET): Heterogeneity loading (significant = ML captures heterogeneity)\n\n")
  
  # Format output table
  out <- copy(x$estimates)
  out[, stars := get_stars(p_value)]
  
  # Print as formatted text
  cat(sprintf("  %6s  %10s  %10s  %8s  %10s\n", 
              "Term", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %6s  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
                out$term[i], out$estimate[i], out$se[i], 
                out$t_value[i], out$p_value[i], out$stars[i]))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}


#' Print method for gates_results objects
#' @param x An object of class \code{gates_results} from \code{gates()}
#' @param ... Additional arguments (currently unused)
#' @export
print.gates_results <- function(x, ...) {
  cat("GATES Results\n")
  cat("=============\n\n")
  
  # Show fit type
  fit_label <- if (x$fit_type == "hte") "HTE (ensemble_hte)" else "Prediction (ensemble_pred)"
  cat("Fit type: ", fit_label, "\n", sep = "")
  
  # Show outcome information
  cat("Outcome analyzed: ", x$outcome, "\n", sep = "")
  if (x$outcome != x$targeted_outcome) {
    group_label <- if (x$fit_type == "hte") "ITE" else "predicted Y"
    cat("  (Groups based on: ", x$targeted_outcome, " ", group_label, ")\n", sep = "")
  }
  cat("Number of groups: ", x$n_groups, "\n", sep = "")
  cat("Repetitions: ", x$M, "\n", sep = "")
  if (!is.null(x$controls)) {
    cat("Controls: ", paste(x$controls, collapse = ", "), "\n", sep = "")
  }
  cat("\nGroup Average Treatment Effects:\n\n")
  
  # Format output table for group estimates
  out <- copy(x$estimates)
  out[, stars := get_stars(p_value)]
  
  # Print as formatted text to avoid data.table column type headers
  cat(sprintf("  %5s  %10s  %10s  %8s  %10s\n", 
              "Group", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %5d  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
                out$group[i], out$estimate[i], out$se[i], 
                out$t_value[i], out$p_value[i], out$stars[i]))
  }
  
  # Print heterogeneity tests
  cat("\nHeterogeneity Tests:\n")
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  cat(sprintf("  %12s  %10s  %10s  %8s  %10s\n", 
              "Test", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  
  # Top-Bottom test
  tb <- x$top_bottom
  tb_stars <- get_stars(tb$p_value)
  cat(sprintf("  %12s  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
              "Top-Bottom", tb$estimate, tb$se, tb$t_value, tb$p_value, tb_stars))
  
  # Top-All test (if available)
  if (!is.null(x$top_all)) {
    ta <- x$top_all
    ta_stars <- get_stars(ta$p_value)
    cat(sprintf("  %12s  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
                "Top-All", ta$estimate, ta$se, ta$t_value, ta$p_value, ta_stars))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}


#' Internal function to compute CLAN for a single repetition
#' 
#' @param predicted_values Predicted values (ITE or Y) vector for this repetition
#' @param fold Fold assignment vector for this repetition
#' @param n_groups Number of groups for CLAN
#' @param variables Character vector of variable names to analyze
#' @param variable_data data.table/data.frame with the variables to analyze
#' @param na_rm Logical, whether to remove NA values in calculations (default: FALSE)
#' 
#' @return data.table with CLAN estimates for each variable
#' @keywords internal
.clan_single <- function(predicted_values, fold, n_groups, variables, variable_data, na_rm = FALSE) {
  
  # Build data.table for CLAN analysis
  dt_clan <- data.table(
    predicted_values = predicted_values,
    fold = fold
  )
  
  # Add variables to analyze
  dt_clan <- cbind(dt_clan, variable_data[, ..variables])
  
  # Create groups within each fold (to ensure proper ranking)
  dt_clan[, group := as.integer(create_groups(predicted_values, n_groups)), by = fold]
  
  # Identify top, bottom, and else groups
  dt_clan[, `:=`(
    top = as.integer(group == n_groups),
    bottom = as.integer(group == 1),
    else_group = as.integer(group != n_groups)
  )]
  
  # Compute means and differences for each variable
  results_list <- lapply(variables, function(var) {
    var_values <- dt_clan[[var]]
    
    # Values by group
    top_values <- var_values[dt_clan$top == 1]
    bottom_values <- var_values[dt_clan$bottom == 1]
    else_values <- var_values[dt_clan$else_group == 1]
    
    # Sample sizes (account for NAs if na_rm = TRUE)
    if (na_rm) {
      n_top <- sum(!is.na(top_values))
      n_bottom <- sum(!is.na(bottom_values))
      n_else <- sum(!is.na(else_values))
      n_all <- sum(!is.na(var_values))
    } else {
      n_top <- length(top_values)
      n_bottom <- length(bottom_values)
      n_else <- length(else_values)
      n_all <- length(var_values)
    }
    
    # Means by group
    mean_top <- mean(top_values, na.rm = na_rm)
    mean_bottom <- mean(bottom_values, na.rm = na_rm)
    mean_else <- mean(else_values, na.rm = na_rm)
    mean_all <- mean(var_values, na.rm = na_rm)
    
    # Variances for SE calculation
    var_top <- var(top_values, na.rm = na_rm) / n_top
    var_bottom <- var(bottom_values, na.rm = na_rm) / n_bottom
    var_else <- var(else_values, na.rm = na_rm) / n_else
    var_all <- var(var_values, na.rm = na_rm) / n_all
    
    # SE for top - bottom difference
    se_top_bottom <- sqrt(var_top + var_bottom)
    
    # SE for top - else difference
    se_top_else <- sqrt(var_top + var_else)
    
    # SE for top - all difference (accounting for overlap)
    # Var(top - all) = Var(top) + Var(all) - 2*Cov(top, all)
    # Since top is subset of all: Cov = Var(top) * n_top / n_all
    se_top_all <- sqrt(var_top + var_all - 2 * var_top * n_top / n_all)
    
    data.table(
      variable = var,
      mean_top = mean_top,
      mean_bottom = mean_bottom,
      mean_else = mean_else,
      mean_all = mean_all,
      diff_top_bottom = mean_top - mean_bottom,
      se_top_bottom = se_top_bottom,
      diff_top_else = mean_top - mean_else,
      se_top_else = se_top_else,
      diff_top_all = mean_top - mean_all,
      se_top_all = se_top_all
    )
  })
  
  rbindlist(results_list)
}


#' Compute CLAN (Classification Analysis)
#' 
#' @description
#' Computes Classification Analysis (CLAN), an analysis introduced by
#' Chernozhukov et al. (2025) to characterize which units have the highest
#' and lowest predicted treatment effects or predicted outcomes.
#' 
#' This function implements the multiple-split estimation strategy developed by
#' Fava (2025), which combines predictions from multiple machine learning algorithms
#' into an ensemble and averages CLAN estimates across M repetitions of K-fold
#' cross-fitting to improve statistical power.
#' 
#' CLAN compares the mean of covariates between units in the top group of predicted
#' values versus: (1) bottom group, (2) all other units (else), and (3) all units.
#' Groups are formed by ranking predictions within each fold.
#' 
#' @references
#' Chernozhukov, V., Demirer, M., Duflo, E., & Fernández-Val, I. (2025).
#' Fisher–Schultz Lecture: Generic Machine Learning Inference on Heterogeneous
#' Treatment Effects in Randomized Experiments, with an Application to Immunization
#' in India. \emph{Econometrica}, 93(4), 1121-1164.
#' 
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class `ensemble_hte_fit` from `ensemble_hte()`
#'   or `ensemble_pred_fit` from `ensemble_pred()`.
#' @param variables Either:
#'   \itemize{
#'     \item \code{NULL} (default): Uses all covariates from the ensemble fit (i.e., \code{names(ensemble_fit$X)})
#'     \item Character vector of variable names present in the `data` used in the ensemble function
#'     \item A data.frame/data.table with covariates to analyze (must have same number of rows as data)
#'   }
#' @param n_groups Number of groups to divide the sample into (default: 3).
#'   CLAN compares the top group (highest predicted ITE or Y) to others.
#' @param na_rm Logical, whether to remove NA values when computing means and
#'   variances (default: FALSE). If FALSE and NAs are present, results will be NA.
#' @param scale Logical, whether to scale non-binary variables to have mean 0
#'   and standard deviation 1 before computing differences (default: FALSE).
#'   This makes coefficients comparable across variables with different scales.
#' @param subset Which observations to use for the CLAN analysis. Options:
#'   \itemize{
#'     \item NULL (default): uses all observations (or training obs for \code{ensemble_pred_fit}
#'       when using the default outcome with subset training)
#'     \item \code{"train"}: uses only training observations (for \code{ensemble_pred_fit} only)
#'     \item \code{"all"}: explicitly uses all observations
#'     \item Logical vector: TRUE/FALSE for each observation (must have same length as data)
#'     \item Integer vector: indices of observations to include (1-indexed)
#'   }
#'   This allows evaluating CLAN on a subset of observations.
#' 
#' @return An object of class `clan_results` containing:
#' \itemize{
#'   \item estimates: data.table with CLAN estimates averaged across repetitions
#'   \item variables: variables analyzed
#'   \item n_groups: number of groups used
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item M: number of repetitions
#'   \item scaled: whether variables were scaled
#'   \item call: the function call
#' }
#' 
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 500
#' X1 <- rnorm(n); X2 <- rnorm(n)
#' D <- rbinom(n, 1, 0.5)
#' Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n)
#' data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
#'
#' fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
#'                     algorithms = c("lm", "grf"), M = 3, K = 3)
#' result <- clan(fit)
#' print(result)
#' plot(result)
#' }
#'
#' @export
clan <- function(ensemble_fit, variables = NULL, n_groups = 3, na_rm = FALSE, scale = FALSE, subset = NULL) {
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- ensemble_fit$ite
    fit_type <- "hte"
  } else if (inherits(ensemble_fit, "ensemble_pred_fit")) {
    predictions_list <- ensemble_fit$predictions
    fit_type <- "pred"
  } else {
    stop("ensemble_fit must be of class 'ensemble_hte_fit' or 'ensemble_pred_fit'")
  }
  
  cl <- match.call()
  
  # Default: use all covariates from the ensemble fit
  if (is.null(variables)) {
    variables <- names(ensemble_fit$X)
  }
  
  # Handle variables input: can be character vector or data.frame
  if (is.data.frame(variables)) {
    # variables is a data.frame/data.table
    if (nrow(variables) != nrow(ensemble_fit$data)) {
      stop(paste0("variables data.frame has ", nrow(variables), 
                  " rows but ensemble_fit data has ", nrow(ensemble_fit$data), " rows"))
    }
    variable_data <- as.data.table(copy(variables))
    var_names <- names(variable_data)
  } else if (is.character(variables)) {
    # variables is a character vector of column names
    missing_vars <- setdiff(variables, names(ensemble_fit$data))
    if (length(missing_vars) > 0) {
      stop(paste0("Variables not found in data: ", paste(missing_vars, collapse = ", ")))
    }
    variable_data <- copy(ensemble_fit$data[, ..variables])
    var_names <- variables
  } else {
    stop("variables must be a character vector or a data.frame")
  }
  
  # Scale non-binary variables if requested
  if (scale) {
    for (var in var_names) {
      col_vals <- variable_data[[var]]
      unique_vals <- unique(col_vals[!is.na(col_vals)])
      # Check if binary (only 0 and 1, or only two values)
      is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
      if (!is_binary) {
        col_mean <- mean(col_vals, na.rm = FALSE)
        col_sd <- sd(col_vals, na.rm = FALSE)
        if (col_sd > 0) {
          variable_data[, (var) := (col_vals - col_mean) / col_sd]
        }
      }
    }
  }
  
  # Get targeted outcome from ensemble_fit
  targeted_outcome <- all.vars(ensemble_fit$formula)[1]
  n <- ensemble_fit$n
  
  # Check for train_idx in ensemble_fit
  has_train_idx <- !is.null(ensemble_fit$train_idx)
  
  # For CLAN, we don't have an outcome variable to check for NAs
  # So using_default_outcome is effectively FALSE for CLAN's subset logic
  using_default_outcome <- FALSE
  
  # Determine which observations to use
  train_idx <- if (has_train_idx) ensemble_fit$train_idx else NULL
  use_idx <- .process_subset(
    subset = subset,
    n = n,
    fit_type = fit_type,
    has_train_idx = has_train_idx,
    train_idx = train_idx,
    using_default_outcome = using_default_outcome,
    Y = NULL
  )
  
  n_used <- sum(use_idx)
  n_fit <- if (fit_type == "pred" && has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("CLAN", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Subset variable data
  variable_data_used <- variable_data[use_idx, ]
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Compute CLAN for each repetition
  clan_by_rep <- lapply(1:M, function(m) {
    preds <- predictions_list[[m]][use_idx]
    fold <- splits[[m]][use_idx]
    .clan_single(
      predicted_values = preds,
      fold = fold,
      n_groups = n_groups,
      variables = var_names,
      variable_data = variable_data_used,
      na_rm = na_rm
    )
  })
  
  # Combine results across repetitions
  all_estimates <- rbindlist(clan_by_rep, idcol = "repetition")
  
  # Aggregate across repetitions: mean of estimates, mean of SE
  combined <- all_estimates[, .(
    mean_top = mean(mean_top),
    mean_bottom = mean(mean_bottom),
    mean_else = mean(mean_else),
    mean_all = mean(mean_all),
    diff_top_bottom = mean(diff_top_bottom),
    se_top_bottom = mean(se_top_bottom),
    diff_top_else = mean(diff_top_else),
    se_top_else = mean(se_top_else),
    diff_top_all = mean(diff_top_all),
    se_top_all = mean(se_top_all)
  ), by = variable]
  
  # Add t-values and p-values for top-bottom difference
  combined[, `:=`(
    t_value_top_bottom = diff_top_bottom / se_top_bottom,
    p_value_top_bottom = 2 * pnorm(-abs(diff_top_bottom / se_top_bottom)),
    t_value_top_else = diff_top_else / se_top_else,
    p_value_top_else = 2 * pnorm(-abs(diff_top_else / se_top_else)),
    t_value_top_all = diff_top_all / se_top_all,
    p_value_top_all = 2 * pnorm(-abs(diff_top_all / se_top_all))
  )]
  
  # Construct clan_results object
  structure(
    list(
      estimates = combined,
      variables = var_names,
      n_groups = n_groups,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      M = M,
      scaled = scale,
      call = cl
    ),
    class = "clan_results"
  )
}


#' Internal function to compute GAVS for a single repetition
#' 
#' @param Y Outcome vector
#' @param predicted_values Predicted values vector for this repetition
#' @param fold Fold assignment vector for this repetition
#' @param n_groups Number of groups for GAVS
#' @param strata Optional strata indicator for restricted ranking (factor/integer vector)
#' 
#' @return List with group_estimates, tests (top_bottom, all, top_all), and reg object
#' @keywords internal
.gavs_single <- function(Y, predicted_values, fold, n_groups = 5, strata = NULL) {
  
  # Build data.table for GAVS regression
  dt_gavs <- data.table(
    Y = Y,
    predicted_values = predicted_values,
    fold = fold
  )
  
  # Create groups within each fold (and optionally within strata for restricted ranking)
  if (is.null(strata)) {
    dt_gavs[, group := as.integer(create_groups(predicted_values, n_groups)), by = fold]
  } else {
    dt_gavs[, strata := strata]
    dt_gavs[, group := as.integer(create_groups(predicted_values, n_groups)), by = .(fold, strata)]
  }
  
  # Create dummy variables for groups
  group_cols <- paste0("group_", 1:n_groups)
  for (g_idx in 1:n_groups) {
    col_name <- paste0("group_", g_idx)
    dt_gavs[, (col_name) := as.integer(group == g_idx)]
  }
  
  # Build regression formula
  # GAVS regression: Y ~ group_1 + group_2 + ... + group_n - 1
  formula_str <- paste("Y ~", paste(group_cols, collapse = " + "), "- 1")
  
  # Run regression with robust SEs
  reg <- reg_from_formula(formula_str, dt_gavs)
  coef_gavs <- reg$coef
  vcov_gavs <- reg$vcov
  
  # Extract estimates and SEs for each group
  group_estimates <- data.table(
    group = 1:n_groups,
    estimate = coef_gavs[group_cols, 1],
    se = coef_gavs[group_cols, 2]
  )
  
  # Define column names
  top_col <- paste0("group_", n_groups)
  bottom_col <- "group_1"
  
  # 1. Top-Bottom difference (group n - group 1)
  top_bottom <- .compute_coef_diff(coef_gavs, vcov_gavs, top_col, bottom_col)
  
  # 2. All: weighted average of all groups (= overall mean)
  # Use group proportions as weights (for equal-sized groups, this is 1/n_groups each)
  n_obs <- nrow(dt_gavs)
  group_props <- vapply(1:n_groups, function(g) sum(dt_gavs$group == g) / n_obs, numeric(1))
  all_estimate <- sum(group_props * coef_gavs[group_cols, 1])
  # SE: sqrt(w' * Vcov * w) where w is the weight vector
  all_se <- sqrt(as.numeric(t(group_props) %*% vcov_gavs[group_cols, group_cols] %*% group_props))
  all_test <- data.table(estimate = all_estimate, se = all_se, weights = list(group_props))
  
  # 3. Top-All: top group minus weighted average of all groups
  top_all_weights <- -group_props
  top_all_weights[n_groups] <- 1 - group_props[n_groups]
  top_all_estimate <- sum(top_all_weights * coef_gavs[group_cols, 1])
  top_all_se <- sqrt(as.numeric(t(top_all_weights) %*% vcov_gavs[group_cols, group_cols] %*% top_all_weights))
  top_all <- data.table(estimate = top_all_estimate, se = top_all_se)
  
  list(
    group_estimates = group_estimates,
    top_bottom = top_bottom,
    all = all_test,
    top_all = top_all,
    reg = reg  # Include full regression object for cross-covariance calculations
  )
}


#' Compute GAVS (Group Averages)
#' 
#' @description
#' Computes Group Averages (GAVS), which measures average outcomes for groups
#' defined by predicted value quantiles. This works for both heterogeneous 
#' treatment effect estimation (groups by predicted ITE) and prediction problems
#' (groups by predicted Y).
#' 
#' This function implements the multiple-split estimation strategy developed in
#' Fava (2025), which combines predictions from multiple machine learning algorithms
#' into an ensemble and averages GAVS estimates across M repetitions of K-fold
#' cross-fitting to improve statistical power.
#' 
#' The method uses ordinary least squares regression of Y on group indicators.
#' Groups are formed by ranking predicted values within each fold.
#' Robust (HC1) standard errors are computed.
#' 
#' @section Subsample Usage:
#' The \code{subset} parameter controls which observations are used for evaluation.
#' This is useful when:
#' \itemize{
#'   \item The ML model was trained on a subset (e.g., using \code{train_idx} in
#'     \code{ensemble_pred()}) and you want to evaluate on the same or different subset.
#'   \item You want to evaluate treatment effect targeting on an outcome that is
#'     only observed for a subset of observations.
#' }
#' 
#' A message is printed when either the ML model was trained on a subset or
#' the evaluation uses a subset, to help avoid unintended subsetting.
#' 
#' @references
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class \code{ensemble_hte_fit} from 
#'   \code{ensemble_hte()} or \code{ensemble_pred_fit} from \code{ensemble_pred()}.
#' @param n_groups Number of groups to divide the sample into (default: 3)
#' @param outcome Either:
#'   \itemize{
#'     \item NULL (default): uses the same outcome as in the ensemble function
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric vector: custom outcome variable (must have appropriate length)
#'   }
#'   This allows computing GAVS for a different outcome than the one used for prediction.
#' @param subset Controls which observations to use for evaluation:
#'   \itemize{
#'     \item \code{NULL} (default): For \code{ensemble_pred_fit} with \code{train_idx},
#'       uses training observations when \code{outcome = NULL}; otherwise uses all.
#'       For \code{ensemble_hte_fit}, uses all observations.
#'     \item \code{"train"}: Use only training observations (only valid for
#'       \code{ensemble_pred_fit} with \code{train_idx}).
#'     \item \code{"all"}: Use all observations.
#'     \item Logical vector: TRUE/FALSE for each observation (length must equal
#'       number of rows in data).
#'     \item Integer vector: Indices of observations to use.
#'   }
#'   See \strong{Subsample Usage} section for guidance on when to use each option.
#' @param strata Optional. Stratification variable for restricted ranking:
#'   \itemize{
#'     \item NULL (default): unrestricted ranking across full sample within folds
#'     \item Character string: column name in the \code{data} for stratified ranking
#'     \item Numeric/factor vector: strata indicator (must have same length as data)
#'   }
#'   When specified, predicted values are ranked within each stratum (and fold),
#'   rather than across the full sample.
#' 
#' @return An object of class \code{gavs_results} containing:
#' \itemize{
#'   \item estimates: data.table with GAVS estimates averaged across repetitions
#'   \item top_bottom: data.table with the top-bottom difference test
#'   \item all: data.table with the overall mean (weighted avg of all groups)
#'   \item top_all: data.table with the top minus average test
#'   \item n_groups: number of groups used
#'   \item outcome: the outcome variable used for GAVS
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item strata: the stratification variable used (if any)
#'   \item n_used: number of observations used
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#' 
#' @examples
#' \donttest{
#' set.seed(123)
#' n <- 500
#' X1 <- rnorm(n); X2 <- rnorm(n)
#' D <- rbinom(n, 1, 0.5)
#' Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n)
#' data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
#'
#' fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
#'                     algorithms = c("lm", "grf"), M = 3, K = 3)
#' result <- gavs(fit, n_groups = 3)
#' print(result)
#' plot(result)
#' }
#'
#' @export
gavs <- function(ensemble_fit, n_groups = 3, outcome = NULL, subset = NULL, 
                 strata = NULL) {
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$ite[[m]])
    fit_type <- "hte"
    has_train_idx <- FALSE
  } else if (inherits(ensemble_fit, "ensemble_pred_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$predictions[[m]])
    fit_type <- "pred"
    # Check if subset training was used
    has_train_idx <- !is.null(ensemble_fit$train_idx) && 
                     !is.null(ensemble_fit$n_train) && 
                     ensemble_fit$n_train < ensemble_fit$n
  } else {
    stop("ensemble_fit must be of class 'ensemble_hte_fit' or 'ensemble_pred_fit'")
  }
  
  cl <- match.call()
  n <- ensemble_fit$n
  
  # Get targeted outcome from ensemble_fit
  targeted_outcome <- all.vars(ensemble_fit$formula)[1]
  
  # Determine which outcome to use for GAVS
  using_default_outcome <- is.null(outcome)
  
  if (is.null(outcome)) {
    Y <- ensemble_fit$Y
    outcome_var <- targeted_outcome
  } else if (is.character(outcome)) {
    if (length(outcome) != 1) {
      stop("outcome must be a single column name or a numeric vector")
    }
    if (!outcome %in% names(ensemble_fit$data)) {
      stop(paste0("outcome '", outcome, "' not found in the data"))
    }
    Y <- ensemble_fit$data[[outcome]]
    outcome_var <- outcome
  } else if (is.numeric(outcome)) {
    if (length(outcome) != n) {
      stop(paste0("outcome vector has length ", length(outcome), 
                  " but data has ", n, " rows"))
    }
    Y <- outcome
    outcome_var <- "custom_outcome"
  } else {
    stop("outcome must be NULL, a character string, or a numeric vector")
  }
  
  # Process strata argument
  strata_vec <- NULL
  strata_var <- NULL
  if (!is.null(strata)) {
    if (is.character(strata)) {
      if (length(strata) != 1) {
        stop("strata must be a single column name or a vector")
      }
      if (!strata %in% names(ensemble_fit$data)) {
        stop(paste0("strata '", strata, "' not found in the data"))
      }
      strata_vec <- ensemble_fit$data[[strata]]
      strata_var <- strata
    } else if (is.numeric(strata) || is.factor(strata)) {
      if (length(strata) != n) {
        stop(paste0("strata vector has length ", length(strata), 
                    " but data has ", n, " rows"))
      }
      strata_vec <- strata
      strata_var <- "custom_strata"
    } else {
      stop("strata must be NULL, a character string, or a numeric/factor vector")
    }
    # Convert to factor if not already
    if (!is.factor(strata_vec)) {
      strata_vec <- as.factor(strata_vec)
    }
  }
  
  # Determine which observations to use
  train_idx <- if (has_train_idx) ensemble_fit$train_idx else NULL
  use_idx <- .process_subset(
    subset = subset,
    n = n,
    fit_type = fit_type,
    has_train_idx = has_train_idx,
    train_idx = train_idx,
    using_default_outcome = using_default_outcome,
    Y = Y
  )
  
  # Validate no NAs in outcome for used observations
  if (any(is.na(Y[use_idx]))) {
    stop("outcome has NA values for the observations being used")
  }
  
  # Validate no NAs in strata for used observations (if strata provided)
  if (!is.null(strata_vec) && any(is.na(strata_vec[use_idx]))) {
    stop("strata has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (fit_type == "pred" && has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("GAVS", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Subset strata if provided
  strata_used <- if (!is.null(strata_vec)) strata_vec[use_idx] else NULL
  
  # Compute GAVS for each repetition
  gavs_by_rep <- lapply(1:M, function(m) {
    .gavs_single(
      Y = Y[use_idx],
      predicted_values = predictions_list[[m]][use_idx],
      fold = splits[[m]][use_idx],
      n_groups = n_groups,
      strata = strata_used
    )
  })
  
  # Combine group estimates across repetitions (mean method from Fava (2025))
  all_group_estimates <- rbindlist(lapply(gavs_by_rep, `[[`, "group_estimates"), idcol = "repetition")
  
  # Aggregate across repetitions: mean of estimates, mean of SE
  combined <- all_group_estimates[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  
  # Add t-values and p-values for combined estimates
  combined[, `:=`(
    t_value = estimate / se,
    p_value = 2 * pnorm(-abs(estimate / se))
  )]
  
  # Helper function to aggregate test results
  .aggregate_test <- function(test_list) {
    all_tests <- rbindlist(test_list, idcol = "repetition")
    result <- all_tests[, .(estimate = mean(estimate), se = mean(se))]
    result[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
    result
  }
  
  # Combine test estimates across repetitions
  top_bottom_combined <- .aggregate_test(lapply(gavs_by_rep, `[[`, "top_bottom"))
  all_combined <- .aggregate_test(lapply(gavs_by_rep, `[[`, "all"))
  top_all_combined <- .aggregate_test(lapply(gavs_by_rep, `[[`, "top_all"))
  
  # Construct gavs_results object
  structure(
    list(
      estimates = combined,
      top_bottom = top_bottom_combined,
      all = all_combined,
      top_all = top_all_combined,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      strata = strata_var,
      n_used = n_used,
      M = M,
      call = cl
    ),
    class = "gavs_results"
  )
}


#' Print method for clan_results objects
#' @param x An object of class \code{clan_results} from \code{clan()}
#' @param ... Additional arguments (currently unused)
#' @export
print.clan_results <- function(x, ...) {
  cat("CLAN Results (Classification Analysis)\n")
  cat("======================================\n\n")
  
cat("Targeted outcome: ", x$targeted_outcome, "\n", sep = "")
  cat("Number of groups: ", x$n_groups, "\n", sep = "")
  cat("Repetitions: ", x$M, "\n", sep = "")
  
  # Determine label based on fit type
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "Y" else "ITE"
  
  out <- copy(x$estimates)
  
  # Top vs Bottom
  cat(sprintf("\nTop vs Bottom (highest vs lowest predicted %s):\n\n", pred_label))
  out[, stars_tb := get_stars(p_value_top_bottom)]
  cat(sprintf("  %-20s  %10s  %10s  %10s  %10s  %10s\n", 
              "Variable", "Mean Top", "Mean Bot", "Diff", "Std.Error", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 78), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %-20s  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f %s\n",
                substr(out$variable[i], 1, 20),
                out$mean_top[i], out$mean_bottom[i],
                out$diff_top_bottom[i], out$se_top_bottom[i],
                out$p_value_top_bottom[i], out$stars_tb[i]))
  }
  
  # Top vs Else
  cat(sprintf("\nTop vs Else (highest predicted %s vs all others):\n\n", pred_label))
  out[, stars_te := get_stars(p_value_top_else)]
  cat(sprintf("  %-20s  %10s  %10s  %10s  %10s  %10s\n", 
              "Variable", "Mean Top", "Mean Else", "Diff", "Std.Error", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 78), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %-20s  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f %s\n",
                substr(out$variable[i], 1, 20),
                out$mean_top[i], out$mean_else[i],
                out$diff_top_else[i], out$se_top_else[i],
                out$p_value_top_else[i], out$stars_te[i]))
  }
  
  # Top vs All
  cat(sprintf("\nTop vs All (highest predicted %s vs full sample):\n\n", pred_label))
  out[, stars_ta := get_stars(p_value_top_all)]
  cat(sprintf("  %-20s  %10s  %10s  %10s  %10s  %10s\n", 
              "Variable", "Mean Top", "Mean All", "Diff", "Std.Error", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 78), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %-20s  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f %s\n",
                substr(out$variable[i], 1, 20),
                out$mean_top[i], out$mean_all[i],
                out$diff_top_all[i], out$se_top_all[i],
                out$p_value_top_all[i], out$stars_ta[i]))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}


#' Plot method for gates_results objects
#'
#' @description
#' Creates a coefficient plot showing GATES estimates with confidence intervals
#' for each group, including top-bottom heterogeneity test results.
#'
#' @param x An object of class \code{gates_results} from \code{gates()}
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   labs theme_minimal theme element_text scale_x_continuous
#' @export
plot.gates_results <- function(x, alpha = 0.05, ...) {
  
  # Compute critical value
  z_crit <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  # Prepare data for plotting
  plot_data <- copy(x$estimates)
  plot_data[, `:=`(
    ci_lower = estimate - z_crit * se,
    ci_upper = estimate + z_crit * se
  )]
  
  # Compute global ATE (average of group estimates - equal-sized groups)
  global_ate <- mean(plot_data$estimate)
  
  # Prepare top-bottom test annotation
  tb <- x$top_bottom
  tb_text <- sprintf(
    "Top-Bottom Test: Est = %.3f, SE = %.3f, p = %.3f%s",
    tb$estimate, tb$se, tb$p_value, get_stars(tb$p_value)
  )
  
  # Determine label based on fit type
  group_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "predicted ITE"
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = global_ate, linetype = "solid", color = "blue", linewidth = 0.8) +
    ggplot2::annotate("text", x = x$n_groups + 0.3, y = global_ate, 
                      label = sprintf("ATE = %.3f", global_ate), 
                      hjust = 0, vjust = -0.5, color = "blue", size = 3.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_continuous(breaks = 1:x$n_groups, 
                                limits = c(0.5, x$n_groups + 1)) +
    ggplot2::labs(
      title = "Group Average Treatment Effects (GATES)",
      subtitle = paste0("Outcome: ", x$outcome,
                        if (x$outcome != x$targeted_outcome) 
                          paste0(" (groups based on: ", x$targeted_outcome, ")") else "",
                        "\n", tb_text),
      x = paste0("Group (by ", group_label, " quantile)"),
      y = "Treatment Effect",
      caption = paste0(conf_level, "% confidence intervals. ", x$M, " repetitions. Blue line = global ATE.")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  invisible(p)
}


#' Print method for gavs_results objects
#' @param x An object of class \code{gavs_results} from \code{gavs()}
#' @param ... Additional arguments (currently unused)
#' @export
print.gavs_results <- function(x, ...) {
  cat("GAVS Results (Group Averages)\n")
  cat("=============================\n\n")
  
  # Determine label based on fit type
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "predicted ITE"
  
  # Show outcome information
  cat("Outcome analyzed: ", x$outcome, "\n", sep = "")
  if (x$outcome != x$targeted_outcome) {
    cat("Targeted outcome: ", x$targeted_outcome, "\n", sep = "")
  }
  cat("Number of groups: ", x$n_groups, "\n", sep = "")
  if (!is.null(x$n_used)) {
    cat("Observations used: ", x$n_used, "\n", sep = "")
  }
  cat("Repetitions: ", x$M, "\n", sep = "")
  cat("\nGroup Average Outcomes (groups by ", pred_label, "):\n\n", sep = "")
  
  # Format output table for group estimates
  out <- copy(x$estimates)
  out[, stars := get_stars(p_value)]
  
  # Print as formatted text
  cat(sprintf("  %5s  %10s  %10s  %8s  %10s\n", 
              "Group", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %5d  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
                out$group[i], out$estimate[i], out$se[i], 
                out$t_value[i], out$p_value[i], out$stars[i]))
  }
  
  # Print heterogeneity tests
  cat("\nHeterogeneity Tests:\n")
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  cat(sprintf("  %12s  %10s  %10s  %8s  %10s\n", 
              "Test", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  
  # Top-Bottom test
  tb <- x$top_bottom
  tb_stars <- get_stars(tb$p_value)
  cat(sprintf("  %12s  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
              "Top-Bottom", tb$estimate, tb$se, tb$t_value, tb$p_value, tb_stars))
  
  # Top-All test (if available)
  if (!is.null(x$top_all)) {
    ta <- x$top_all
    ta_stars <- get_stars(ta$p_value)
    cat(sprintf("  %12s  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
                "Top-All", ta$estimate, ta$se, ta$t_value, ta$p_value, ta_stars))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}


#' Plot method for gavs_results objects
#'
#' @description
#' Creates a coefficient plot showing GAVS estimates with confidence intervals
#' for each group, including top-bottom heterogeneity test results.
#'
#' @param x An object of class \code{gavs_results} from \code{gavs()}
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   labs theme_minimal theme element_text scale_x_continuous
#' @export
plot.gavs_results <- function(x, alpha = 0.05, ...) {
  
  # Compute critical value
  z_crit <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  # Determine label based on fit type
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "predicted ITE"
  
  # Prepare data for plotting
  plot_data <- copy(x$estimates)
  plot_data[, `:=`(
    ci_lower = estimate - z_crit * se,
    ci_upper = estimate + z_crit * se
  )]
  
  # Compute global average (average of group estimates - equal-sized groups)
  global_avg <- mean(plot_data$estimate)
  
  # Prepare top-bottom test annotation
  tb <- x$top_bottom
  tb_text <- sprintf(
    "Top-Bottom: %.3f (SE: %.3f, p: %.3f)",
    tb$estimate, tb$se, tb$p_value
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = estimate)) +
    ggplot2::geom_hline(yintercept = global_avg, linetype = "solid", color = "blue", linewidth = 0.8) +
    ggplot2::annotate("text", x = x$n_groups + 0.3, y = global_avg, 
                      label = sprintf("Mean = %.3f", global_avg), 
                      hjust = 0, vjust = -0.5, color = "blue", size = 3.5) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2
    ) +
    ggplot2::scale_x_continuous(breaks = 1:x$n_groups,
                                limits = c(0.5, x$n_groups + 1)) +
    ggplot2::labs(
      title = "Group Averages (GAVS)",
      subtitle = tb_text,
      x = paste0("Group (by ", pred_label, ")"),
      y = paste0("Mean ", x$outcome),
      caption = paste0(conf_level, "% confidence intervals. ", x$M, " repetitions. Blue line = global mean.")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  invisible(p)
}

#' Plot method for clan_results objects
#'
#' @description
#' Creates a bar plot showing CLAN estimates (differences in covariate means
#' between top group and comparison group) with confidence intervals. Bars are
#' colored based on whether the difference is positive or negative, and ordered
#' from lowest to highest.
#'
#' @param x An object of class \code{clan_results} from \code{clan()}
#' @param comparison Which comparison to plot: "top_bottom", "top_else", or "top_all"
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar geom_hline
#'   labs theme_minimal theme element_text coord_flip scale_fill_manual
#' @export
plot.clan_results <- function(x, comparison = c("top_bottom", "top_else", "top_all"), alpha = 0.05, ...) {
  
  comparison <- match.arg(comparison)
  
  # Compute critical value
  z_crit <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  # Determine label based on fit type
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "predicted ITE"
  
  # Select appropriate columns based on comparison
  plot_data <- copy(x$estimates)
  
  if (comparison == "top_bottom") {
    plot_data[, `:=`(
      diff = diff_top_bottom,
      se_diff = se_top_bottom,
      comparison_label = "Top vs Bottom"
    )]
  } else if (comparison == "top_else") {
    plot_data[, `:=`(
      diff = diff_top_else,
      se_diff = se_top_else,
      comparison_label = "Top vs Else"
    )]
  } else {
    plot_data[, `:=`(
      diff = diff_top_all,
      se_diff = se_top_all,
      comparison_label = "Top vs All"
    )]
  }
  
  plot_data[, `:=`(
    ci_lower = diff - z_crit * se_diff,
    ci_upper = diff + z_crit * se_diff,
    sign = ifelse(diff >= 0, "Positive", "Negative")
  )]
  
  # Reorder variables by effect size (lowest to highest)
  plot_data[, variable := factor(variable, levels = variable[order(diff)])]
  
  # Create bar plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = diff, fill = sign)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.2, linewidth = 0.5, color = "black"
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c("Negative" = "#E74C3C", "Positive" = "#27AE60"),
      guide = "none"
    ) +
    ggplot2::labs(
      title = "Classification Analysis (CLAN)",
      subtitle = paste0(plot_data$comparison_label[1], " - Groups by ", pred_label),
      x = NULL,
      y = "Difference in Means",
      caption = paste0(conf_level, "% confidence intervals. ", x$M, " repetitions. ", x$n_groups, " groups.",
                       if (x$scaled) " Variables scaled." else "")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  invisible(p)
}

