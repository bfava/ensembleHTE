
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
    if (using_default_outcome && has_train_idx) {
      # For fits with subset training and default outcome, use training obs
      use_idx <- train_idx
    } else {
      # Otherwise use all observations
      use_idx <- rep(TRUE, n)
    }
  } else if (is.character(subset) && length(subset) == 1) {
    if (subset == "train") {
      if (!has_train_idx) {
        stop("subset = 'train' specified but no train_idx was used in the ensemble fit")
      }
      use_idx <- train_idx
    } else if (subset == "all") {
      use_idx <- rep(TRUE, n)
      # Warn if using default outcome with NAs
      if (using_default_outcome && has_train_idx && any(is.na(Y))) {
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
    hint <- ""
    if (has_train_idx && n_eval == n_total) {
      hint <- ' Use subset = "train" to restrict to training observations.'
    } else if (has_train_idx && n_eval < n_total && n_eval == n_fit) {
      hint <- ' Use subset = "all" to use all observations.'
    }
    message(paste0("Note: ", paste(msg_parts, collapse = "; "), ".", hint))
  }
  
  invisible(NULL)
}


#' Internal function to compute GATES for a single repetition
#' 
#' @param Y Outcome vector (full length)
#' @param D Treatment vector (full length)
#' @param prop_score Propensity score vector (full length)
#' @param weight Weight vector (full length)
#' @param predicted_values Predicted values (ITE or Y) vector for this repetition (full length)
#' @param fold Fold assignment vector for this repetition (full length)
#' @param n_groups Number of groups for GATES
#' @param restrict_by Optional restrict_by indicator for restricted ranking (factor/integer vector, full length)
#' @param controls Optional control variables (character vector of column names)
#' @param control_data Optional data.table/data.frame with control variables (full length)
#' @param group_ref_idx Optional logical vector (same length as Y). When provided,
#'   group cutoffs are computed from only the reference observations, then applied
#'   to all observations. If \code{NULL}, all observations are used for group formation.
#' @param analysis_idx Optional logical vector (same length as Y). When provided,
#'   the regression is run only on these observations (after group assignment).
#'   If \code{NULL}, all observations are used for the regression.
#' @param cluster_id Optional vector of cluster identifiers for cluster-robust SEs.
#'   When provided, passed to \code{reg_from_formula} for clustered inference.
#' 
#' @return List with group_estimates, tests (top_bottom, all, top_all), and reg object
#' @keywords internal
.gates_single <- function(Y, D, prop_score, weight, predicted_values, fold, 
                          n_groups = 5, restrict_by = NULL, controls = NULL, 
                          control_data = NULL,
                          group_ref_idx = NULL, analysis_idx = NULL,
                          cluster_id = NULL) {
  
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
  # When group_ref_idx is provided, cutoffs come from the reference observations only
  if (is.null(restrict_by)) {
    dt_gates[, group := {
      ref <- if (!is.null(group_ref_idx)) group_ref_idx[.I] else NULL
      as.integer(create_groups_by_reference(predicted_values, n_groups, ref_mask = ref))
    }, by = fold]
  } else {
    dt_gates[, strata := restrict_by]
    dt_gates[, group := {
      ref <- if (!is.null(group_ref_idx)) group_ref_idx[.I] else NULL
      as.integer(create_groups_by_reference(predicted_values, n_groups, ref_mask = ref))
    }, by = .(fold, strata)]
  }
  
  # Subset to analysis observations (after group formation)
  if (!is.null(analysis_idx)) {
    dt_gates <- dt_gates[analysis_idx, ]
    if (!is.null(control_data)) {
      control_data <- control_data[analysis_idx, ]
    }
    if (!is.null(cluster_id)) {
      cluster_id <- cluster_id[analysis_idx]
    }
    # Check for empty groups after subsetting
    observed_groups <- sort(unique(dt_gates$group))
    empty_groups <- setdiff(1:n_groups, observed_groups)
    if (length(empty_groups) > 0) {
      stop("After subsetting, group(s) ", paste(empty_groups, collapse = ", "),
           " have no observations. This happens when groups are formed on all ",
           "observations but the analysis subset excludes entire groups. ",
           'To form groups within the subset instead, use group_on = "analysis".')
    }
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
  reg <- reg_from_formula(formula_str, dt_gates, weights = dt_gates$weight,
                          cluster_id = cluster_id)
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
#' @section Estimation Procedure:
#' For each repetition \eqn{m = 1, \ldots, M}:
#' \enumerate{
#'   \item Observations are assigned to \code{n_groups} quantile-based groups by
#'     ranking the ensemble predictions from repetition \eqn{m} \emph{within each
#'     fold}. Group 1 contains the lowest predicted values and group
#'     \code{n_groups} the highest. Forming groups within folds ensures that
#'     group assignment is independent of the model used to generate predictions
#'     for that observation (since predictions are out-of-sample within each fold).
#'   \item A single weighted least squares regression is run on all observations:
#'     \deqn{Y_i = \alpha + \sum_{g=1}^{G} \gamma_g \, \mathbf{1}\{i \in g\}(D_i - e(X_i)) + \varepsilon_i}
#'     where \eqn{e(X_i)} is the propensity score. Weights are
#'     \eqn{1 / (e(X_i)(1 - e(X_i)))}. Each \eqn{\gamma_g} estimates the
#'     average treatment effect for group \eqn{g}.
#'   \item HC1 heteroskedasticity-robust standard errors are computed
#'     (or cluster-robust SEs when \code{individual_id} was specified).
#' }
#' The final reported estimates and standard errors are the simple averages of
#' the per-repetition estimates and standard errors across all \eqn{M} repetitions.
#' 
#' Three heterogeneity tests are reported:
#' \itemize{
#'   \item \strong{Top-Bottom}: difference between the top and bottom groups.
#'   \item \strong{All}: weighted average of all group effects (estimates the overall ATE).
#'   \item \strong{Top-All}: difference between the top group and the overall ATE.
#' }
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
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric value: constant propensity score for all observations
#'     \item Numeric vector: observation-specific propensity scores
#'   }
#'   For \code{ensemble_hte_fit}, uses the propensity score from the fit.
#' @param controls Character vector of control variable names to include in regression.
#'   These must be column names present in the \code{data} argument used when calling
#'   the ensemble function.
#' @param restrict_by Optional. Stratification variable for restricted ranking:
#'   \itemize{
#'     \item NULL (default): unrestricted ranking across full sample within folds
#'     \item Character string: column name in the \code{data} for stratified ranking
#'     \item Numeric/factor vector: group indicator (must have same length as data)
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
#' @param group_on Character controlling which observations define the quantile
#'   cutoffs used to form groups. One of:
#'   \itemize{
#'     \item \code{"auto"} (default): Uses the ML training population. For
#'       \code{ensemble_hte_fit} this is all observations. For
#'       \code{ensemble_pred_fit} with \code{train_idx}, it is the training
#'       subset. This ensures an observation's group assignment does not change
#'       when you vary the analysis subset.
#'     \item \code{"all"}: Always form groups using all observations.
#'     \item \code{"analysis"}: Form groups within whatever observations are
#'       being analyzed (i.e. the \code{subset}).
#'   }
#'   Has no effect when \code{subset = NULL} and all observations are used.
#' 
#' @return An object of class \code{gates_results} containing:
#' \itemize{
#'   \item estimates: data.table with GATES estimates averaged across repetitions.
#'     Columns: \code{group} (integer group index, 1 = lowest predicted effects),
#'     \code{estimate} (group-specific treatment effect), \code{se} (standard error),
#'     \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item top_bottom: data.table with the top-bottom difference test.
#'     Columns: \code{estimate}, \code{se}, \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item all: data.table with the average treatment effect (weighted avg of all groups).
#'     Columns: \code{estimate}, \code{se}, \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item top_all: data.table with the top minus average test.
#'     Columns: \code{estimate}, \code{se}, \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item n_groups: number of groups used
#'   \item outcome: the outcome variable used for GATES
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item restrict_by: the restrict_by variable used (if any)
#'   \item controls: control variables used (if any)
#'   \item group_on: how groups are formed ("auto", "all", or "analysis")
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#' 
#' @examples
#' \donttest{
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#'
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
#' result <- gates(fit, n_groups = 3)
#' print(result)
#' plot(result)
#' }
#'
#' @export
gates <- function(ensemble_fit, n_groups = 3, outcome = NULL, treatment = NULL,
                  prop_score = NULL, controls = NULL, restrict_by = NULL, subset = NULL,
                  group_on = c("auto", "all", "analysis")) {
  
  # Match group_on argument
  group_on <- match.arg(group_on)
  
  # Validate n_groups
  if (!is.numeric(n_groups) || length(n_groups) != 1 || n_groups < 2 || n_groups %% 1 != 0) {
    stop("n_groups must be an integer >= 2")
  }
  if (n_groups >= ensemble_fit$n) {
    stop("n_groups (", n_groups, ") must be smaller than the number of ",
         "observations in the data (n = ", ensemble_fit$n, "). ",
         "Each group needs at least a few observations for valid inference.")
  }
  
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
        stop(paste0("treatment '", treatment, "' not found in the data. ",
                    "The fit object stores a snapshot of the data at fit time. ",
                    "If the treatment column was added after fitting, pass it as a ",
                    "numeric vector instead: treatment = your_data$", treatment))
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
    
    # Resolve prop_score column name (if a single string referencing data)
    if (!is.null(prop_score) && is.character(prop_score) && length(prop_score) == 1) {
      if (prop_score %in% names(ensemble_fit$data)) {
        prop_score <- ensemble_fit$data[[prop_score]]
      } else {
        stop("Column '", prop_score, "' (passed as prop_score) not found in the data")
      }
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
  
  # Process restrict_by argument
  strata_vec <- NULL
  strata_var <- NULL
  if (!is.null(restrict_by)) {
    if (is.character(restrict_by) && length(restrict_by) == 1) {
      # Single string: treat as column name
      if (!restrict_by %in% names(ensemble_fit$data)) {
        stop(paste0("restrict_by '", restrict_by, "' not found in the data"))
      }
      strata_vec <- ensemble_fit$data[[restrict_by]]
      strata_var <- restrict_by
    } else if (is.numeric(restrict_by) || is.factor(restrict_by) || is.character(restrict_by)) {
      # Vector passed directly (numeric, factor, or character)
      if (length(restrict_by) != n) {
        stop(paste0("restrict_by vector has length ", length(restrict_by), 
                    " but data has ", n, " rows"))
      }
      strata_vec <- restrict_by
      strata_var <- "custom_strata"
    } else {
      stop("restrict_by must be NULL, a column name (string), or a vector (numeric/factor/character)")
    }
    # Convert to factor if not already
    if (!is.factor(strata_vec)) {
      strata_vec <- as.factor(strata_vec)
    }
  }
  
  # Check if outcome is the default (targeted) outcome
  using_default_outcome <- is.null(outcome)
  
  # Check for train_idx in ensemble_fit
  has_train_idx <- !is.null(ensemble_fit$train_idx) &&
                   !is.null(ensemble_fit$n_train) &&
                   ensemble_fit$n_train < ensemble_fit$n
  
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
    n_na <- sum(is.na(Y[use_idx]))
    stop("The outcome variable has ", n_na, " NA value(s) among the ", sum(use_idx),
         " observations being used. ",
         if (has_train_idx) "If you used 'train_idx' during fitting, try subset = 'train' or subset = NULL to restrict to observations with observed outcomes." else "Remove or impute missing values before analysis.")
  }
  
  # Validate no NAs in treatment for used observations
  if (any(is.na(D[use_idx]))) {
    stop("treatment has NA values for the observations being used")
  }
  
  # Validate no NAs in strata for used observations (if restrict_by provided)
  if (!is.null(strata_vec) && any(is.na(strata_vec[use_idx]))) {
    stop("restrict_by has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("GATES", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Get control data if controls specified
  control_data <- if (!is.null(controls)) ensemble_fit$data else NULL
  
  # Resolve which observations define group cutoffs
  group_ref_idx <- .resolve_group_ref_idx(group_on, use_idx, train_idx)
  
  # For small-cell check, use the reference population (or full data if NULL)
  check_ref <- if (!is.null(group_ref_idx)) group_ref_idx else rep(TRUE, n)
  check_folds <- lapply(splits, function(s) s[check_ref])
  check_strata <- if (!is.null(strata_vec)) strata_vec[check_ref] else NULL
  .check_small_cells(check_folds, n_groups, restrict_by = check_strata, func_name = "GATES")
  
  # Extract individual_id for cluster-robust SEs (if panel data)
  cluster_id <- ensemble_fit$individual_id
  
  gates_by_rep <- lapply(1:M, function(m) {
    .gates_single(
      Y = Y,
      D = D,
      prop_score = ps,
      weight = weight,
      predicted_values = predictions_list[[m]],
      fold = splits[[m]],
      n_groups = n_groups,
      restrict_by = strata_vec,
      controls = controls,
      control_data = control_data,
      group_ref_idx = group_ref_idx,
      analysis_idx = if (!all(use_idx)) use_idx else NULL,
      cluster_id = cluster_id
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
    result <- all_tests[, .(estimate = mean(estimate), se = mean(se), n_reps = .N)]
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
      restrict_by = strata_var,
      controls = controls,
      group_on = group_on,
      n = n,
      n_train = n_fit,
      n_used = n_used,
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
#' @param cluster_id Optional vector of cluster identifiers for cluster-robust SEs.
#'   When provided, passed to \code{reg_from_formula} for clustered inference.
#' 
#' @return data.table with BLP estimates (beta1 for ATE, beta2 for heterogeneity)
#' @keywords internal
.blp_single <- function(Y, D, prop_score, weight, predicted_ite, 
                        controls = NULL, control_data = NULL,
                        cluster_id = NULL) {
  
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
  reg <- reg_from_formula(formula_str, dt_blp, weights = dt_blp$weight,
                          cluster_id = cluster_id)
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
#' cross-fitting to improve statistical power.
#' 
#' @section Estimation Procedure:
#' For each repetition \eqn{m = 1, \ldots, M}:
#' \enumerate{
#'   \item The ensemble predictions (ITEs or predicted Y) from repetition \eqn{m} are
#'     used as regressors. These predictions were generated by cross-fitting in
#'     \code{ensemble_hte()} or \code{ensemble_pred()}: for each fold \eqn{k}, the
#'     model was trained on all folds except \eqn{k} and predicted on fold \eqn{k},
#'     so each observation's prediction is out-of-sample.
#'   \item A single weighted least squares regression is run on all observations:
#'     \deqn{Y_i = \alpha + \beta_1 W_{1i} + \beta_2 W_{2i} + \varepsilon_i}
#'     where \eqn{W_{1i} = D_i - e(X_i)} and
#'     \eqn{W_{2i} = (D_i - e(X_i))(\hat{s}(X_i) - \bar{\hat{s}})},
#'     with \eqn{e(X_i)} the propensity score and \eqn{\hat{s}(X_i)} the predicted
#'     ITE (or predicted Y). Weights are \eqn{1 / (e(X_i)(1 - e(X_i)))}.
#'   \item HC1 heteroskedasticity-robust standard errors are computed
#'     (or cluster-robust SEs when \code{individual_id} was specified in the
#'     ensemble fit).
#' }
#' The final reported estimates and standard errors are the simple averages of
#' the per-repetition estimates and standard errors across all \eqn{M} repetitions.
#' 
#' \strong{Interpretation:}
#' \itemize{
#'   \item \eqn{\beta_1} (ATE): estimates the Average Treatment Effect.
#'   \item \eqn{\beta_2} (HET): a significant \eqn{\beta_2} indicates that the
#'     ML predictions capture meaningful treatment effect heterogeneity.
#' }
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
#'     \item Character string: column name in the \code{data} used in the ensemble function
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
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#'
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
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
        stop(paste0("treatment '", treatment, "' not found in the data. ",
                    "The fit object stores a snapshot of the data at fit time. ",
                    "If the treatment column was added after fitting, pass it as a ",
                    "numeric vector instead: treatment = your_data$", treatment))
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
    
    # Resolve prop_score column name (if a single string referencing data)
    if (!is.null(prop_score) && is.character(prop_score) && length(prop_score) == 1) {
      if (prop_score %in% names(ensemble_fit$data)) {
        prop_score <- ensemble_fit$data[[prop_score]]
      } else {
        stop("Column '", prop_score, "' (passed as prop_score) not found in the data")
      }
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
  has_train_idx <- !is.null(ensemble_fit$train_idx) &&
                   !is.null(ensemble_fit$n_train) &&
                   ensemble_fit$n_train < ensemble_fit$n
  
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
    n_na <- sum(is.na(Y[use_idx]))
    stop("The outcome variable has ", n_na, " NA value(s) among the ", sum(use_idx),
         " observations being used. ",
         if (has_train_idx) "If you used 'train_idx' during fitting, try subset = 'train' or subset = NULL to restrict to observations with observed outcomes." else "Remove or impute missing values before analysis.")
  }
  
  # Validate no NAs in treatment for used observations
  if (any(is.na(D[use_idx]))) {
    stop("treatment has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("BLP", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  M <- ensemble_fit$M
  
  # Get control data if controls specified (subset to used observations)
  control_data <- if (!is.null(controls)) ensemble_fit$data[use_idx, ] else NULL
  
  # Extract individual_id for cluster-robust SEs (if panel data)
  cluster_id <- if (!is.null(ensemble_fit$individual_id)) ensemble_fit$individual_id[use_idx] else NULL
  
  # Compute BLP for each repetition
  blp_by_rep <- lapply(1:M, function(m) {
    .blp_single(
      Y = Y[use_idx],
      D = D[use_idx],
      prop_score = ps[use_idx],
      weight = weight[use_idx],
      predicted_ite = predictions_list[[m]][use_idx],
      controls = controls,
      control_data = control_data,
      cluster_id = cluster_id
    )
  })
  
  # Combine results across repetitions (average estimates and SEs)
  all_estimates <- rbindlist(lapply(blp_by_rep, `[[`, "estimates"), idcol = "repetition")
  
  combined <- all_estimates[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
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
      n = n,
      n_train = n_fit,
      n_used = n_used,
      M = M,
      call = cl
    ),
    class = "blp_results"
  )
}


#' Print a "Data Usage" block for analysis result objects
#'
#' Displays how data was partitioned between ML training, the analysis subset,
#' and group formation. Only prints when there is something non-trivial to
#' report (i.e., when not all observations are used for everything).
#'
#' @param n Total observations in the original data
#' @param n_train Observations used to train the ML model
#' @param n_used Observations used for this analysis
#' @param group_on group_on value (\"auto\", \"all\", \"analysis\"), or NULL for
#'   functions like BLP that do not form groups
#' @param n_groups Number of groups (NULL for BLP)
#' @keywords internal
.print_data_usage <- function(n, n_train, n_used, group_on = NULL, n_groups = NULL) {
  # Nothing to report when everything uses the full data
  if (n_train == n && n_used == n) return(invisible(NULL))

  cat("\nData usage:\n")
  cat("  ML trained on:         ", n_train, " of ", n, " obs\n", sep = "")
  cat("  Analysis evaluated on: ", n_used, " of ", n, " obs\n", sep = "")
  if (!is.null(group_on) && !is.null(n_groups)) {
    group_desc <- switch(group_on,
      "auto" = paste0("all observations (", n, " obs)"),
      "all" = paste0("all observations (", n, " obs)"),
      "analysis" = paste0("analysis subset (", n_used, " obs)"),
      paste0(group_on, " (", n, " obs)")
    )
    cat("  Groups (", n_groups, ") formed on: ", group_desc, "\n", sep = "")
  }
  invisible(NULL)
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
  
  # Data usage block
  if (!is.null(x$n)) {
    .print_data_usage(x$n, x$n_train, x$n_used)
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
    cat(sprintf("  %6s  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
                out$term[i], out$estimate[i], out$se[i], 
                out$t_value[i], out$p_value[i], out$stars[i]))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}


#' Plot method for blp_results objects
#'
#' @description
#' Creates a coefficient plot showing BLP estimates (ATE and heterogeneity loading)
#' with confidence intervals. Each coefficient is shown in its own panel with an
#' independent y-axis scale. A dashed reference line marks the calibration target:
#' 0 for beta1 (ATE) and 1 for beta2 (heterogeneity loading).
#'
#' @param x An object of class \code{blp_results} from \code{blp()}
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   labs theme_minimal theme element_text facet_wrap
#' @export
plot.blp_results <- function(x, alpha = 0.05, ...) {
  
  z_crit <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  plot_data <- copy(x$estimates)
  plot_data[, `:=`(
    ci_lower  = estimate - z_crit * se,
    ci_upper  = estimate + z_crit * se,
    label     = ifelse(term == "beta1", "beta1 (ATE)", "beta2 (HET)"),
    ref_value = ifelse(term == "beta1", 0, 1)
  )]
  # Fix panel order: ATE first, HET second
  plot_data[, label := factor(label, levels = c("beta1 (ATE)", "beta2 (HET)"))]
  
  # Per-panel reference data (used by geom_hline inside facets)
  ref_df <- unique(plot_data[, .(label, ref_value)])
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = "", y = estimate)) +
    ggplot2::geom_hline(
      data = ref_df,
      ggplot2::aes(yintercept = ref_value),
      linetype = "dashed", color = "gray50"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15, linewidth = 0.8
    ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::facet_wrap(~ label, scales = "free_y") +
    ggplot2::labs(
      title    = "Best Linear Predictor (BLP)",
      subtitle = paste0("Outcome: ", x$outcome),
      x        = NULL,
      y        = "Estimate",
      caption  = paste0(conf_level, "% confidence intervals. ", x$M,
                        " repetitions. Dashed line: calibration target (0 for ATE, 1 for HET).")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold"),
      axis.title   = ggplot2::element_text(face = "bold"),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  
  print(p)
  invisible(p)
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
  
  # Data usage block
  if (!is.null(x$n)) {
    .print_data_usage(x$n, x$n_train, x$n_used, x$group_on, x$n_groups)
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
    cat(sprintf("  %5d  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
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
  cat(sprintf("  %12s  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
              "Top-Bottom", tb$estimate, tb$se, tb$t_value, tb$p_value, tb_stars))
  
  # Top-All test (if available)
  if (!is.null(x$top_all)) {
    ta <- x$top_all
    ta_stars <- get_stars(ta$p_value)
    cat(sprintf("  %12s  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
                "Top-All", ta$estimate, ta$se, ta$t_value, ta$p_value, ta_stars))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}


#' Internal function to compute CLAN for a single repetition
#' 
#' @param predicted_values Predicted values (ITE or Y) vector for this repetition (full length)
#' @param fold Fold assignment vector for this repetition (full length)
#' @param n_groups Number of groups for CLAN
#' @param variables Character vector of variable names to analyze
#' @param variable_data data.table/data.frame with the variables to analyze (full length)
#' @param na_rm Logical, whether to remove NA values in calculations (default: FALSE)
#' @param group_ref_idx Optional logical vector (same length as predicted_values).
#'   When provided, group cutoffs are computed from only the reference observations,
#'   then applied to all observations. If \code{NULL}, all observations are used.
#' @param analysis_idx Optional logical vector (same length as predicted_values).
#'   When provided, means and SEs are computed only for these observations (after
#'   group assignment). If \code{NULL}, all observations are used.
#' @param cluster_id Optional vector of cluster identifiers for cluster-robust SEs.
#'   When provided, uses a regression-based approach with clustered inference
#'   instead of the default analytical SE formulas.
#' 
#' @return data.table with CLAN estimates for each variable
#' @keywords internal
.clan_single <- function(predicted_values, fold, n_groups, variables, variable_data,
                         na_rm = FALSE, group_ref_idx = NULL, analysis_idx = NULL,
                         cluster_id = NULL) {
  
  # Build data.table for CLAN analysis
  dt_clan <- data.table(
    predicted_values = predicted_values,
    fold = fold
  )
  
  # Add variables to analyze
  dt_clan <- cbind(dt_clan, variable_data[, ..variables])
  
  # Create groups within each fold (using reference observations for cutoffs)
  dt_clan[, group := {
    ref <- if (!is.null(group_ref_idx)) group_ref_idx[.I] else NULL
    as.integer(create_groups_by_reference(predicted_values, n_groups, ref_mask = ref))
  }, by = fold]
  
  # Identify top, bottom, and else groups
  dt_clan[, `:=`(
    top = as.integer(group == n_groups),
    bottom = as.integer(group == 1),
    else_group = as.integer(group != n_groups)
  )]
  
  # Subset to analysis observations (after group formation)
  if (!is.null(analysis_idx)) {
    dt_clan <- dt_clan[analysis_idx, ]
    if (!is.null(cluster_id)) {
      cluster_id <- cluster_id[analysis_idx]
    }
    # Check for empty groups after subsetting
    observed_groups <- sort(unique(dt_clan$group))
    empty_groups <- setdiff(1:n_groups, observed_groups)
    if (length(empty_groups) > 0) {
      stop("After subsetting, group(s) ", paste(empty_groups, collapse = ", "),
           " have no observations. This happens when groups are formed on all ",
           "observations but the analysis subset excludes entire groups. ",
           'To form groups within the subset instead, use group_on = "analysis".')
    }
  }
  
  # Compute means and differences for each variable
  results_list <- lapply(variables, function(var) {
    var_values <- dt_clan[[var]]
    
    if (!is.null(cluster_id)) {
      # --- Cluster-robust approach: use regressions for proper clustered SEs ---
      # Regress variable on top + bottom + else_group (no intercept) for group means
      dt_reg <- data.table(
        var_value = var_values,
        top = dt_clan$top,
        bottom = dt_clan$bottom,
        else_group = dt_clan$else_group
      )
      
      # Handle NAs if na_rm = TRUE
      if (na_rm) {
        keep <- !is.na(dt_reg$var_value)
        dt_reg <- dt_reg[keep, ]
        cluster_id_var <- cluster_id[keep]
      } else {
        cluster_id_var <- cluster_id
      }
      
      # Regression for top, bottom, else means with cluster-robust SEs
      reg_groups <- reg_from_formula("var_value ~ top + bottom + else_group - 1", 
                                     dt_reg, cluster_id = cluster_id_var)
      
      mean_top <- reg_groups$coef["top", "Estimate"]
      mean_bottom <- reg_groups$coef["bottom", "Estimate"]
      mean_else <- reg_groups$coef["else_group", "Estimate"]
      se_top_raw <- reg_groups$coef["top", "Std. Error"]
      se_bottom_raw <- reg_groups$coef["bottom", "Std. Error"]
      se_else_raw <- reg_groups$coef["else_group", "Std. Error"]
      vcov_groups <- reg_groups$vcov
      
      # Overall mean from intercept-only regression
      reg_all <- reg_from_formula("var_value ~ 1", dt_reg, cluster_id = cluster_id_var)
      mean_all <- reg_all$coef["(Intercept)", "Estimate"]
      se_all_raw <- reg_all$coef["(Intercept)", "Std. Error"]
      
      # top - bottom difference using vcov
      diff_tb_var <- vcov_groups["top", "top"] + vcov_groups["bottom", "bottom"] - 
                     2 * vcov_groups["top", "bottom"]
      se_diff_top_bottom <- sqrt(max(diff_tb_var, 0))
      
      # top - else difference using vcov
      diff_te_var <- vcov_groups["top", "top"] + vcov_groups["else_group", "else_group"] - 
                     2 * vcov_groups["top", "else_group"]
      se_diff_top_else <- sqrt(max(diff_te_var, 0))
      
      # top - all: cross-covariance between group regression and intercept-only regression
      cross_cov <- .calc_cov_cross_reg(
        reg_groups$model_matrix, reg_all$model_matrix,
        reg_groups$residuals, reg_all$residuals,
        weights = NULL, cluster_id = cluster_id_var
      )
      # Var(top - all) = Var(top) + Var(all) - 2*Cov(top, all)
      top_idx <- which(rownames(reg_groups$coef) == "top")
      diff_ta_var <- se_top_raw^2 + se_all_raw^2 - 2 * cross_cov[top_idx, 1]
      se_diff_top_all <- sqrt(max(diff_ta_var, 0))
      
      data.table(
        variable = var,
        mean_top = mean_top,
        se_top = se_top_raw,
        mean_bottom = mean_bottom,
        se_bottom = se_bottom_raw,
        mean_else = mean_else,
        se_else = se_else_raw,
        mean_all = mean_all,
        se_all = se_all_raw,
        diff_top_bottom = mean_top - mean_bottom,
        se_diff_top_bottom = se_diff_top_bottom,
        diff_top_else = mean_top - mean_else,
        se_diff_top_else = se_diff_top_else,
        diff_top_all = mean_top - mean_all,
        se_diff_top_all = se_diff_top_all
      )
    } else {
      # --- Standard approach: analytical SEs (no clustering) ---
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
      se_diff_top_bottom <- sqrt(var_top + var_bottom)
      
      # SE for top - else difference
      se_diff_top_else <- sqrt(var_top + var_else)
      
      # SE for top - all difference (accounting for overlap)
      # Var(top - all) = Var(top) + Var(all) - 2*Cov(top, all)
      # Since top is subset of all: Cov = Var(top) * n_top / n_all
      se_diff_top_all <- sqrt(var_top + var_all - 2 * var_top * n_top / n_all)
      
      data.table(
        variable = var,
        mean_top = mean_top,
        se_top = sqrt(var_top),
        mean_bottom = mean_bottom,
        se_bottom = sqrt(var_bottom),
        mean_else = mean_else,
        se_else = sqrt(var_else),
        mean_all = mean_all,
        se_all = sqrt(var_all),
        diff_top_bottom = mean_top - mean_bottom,
        se_diff_top_bottom = se_diff_top_bottom,
        diff_top_else = mean_top - mean_else,
        se_diff_top_else = se_diff_top_else,
        diff_top_all = mean_top - mean_all,
        se_diff_top_all = se_diff_top_all
      )
    }
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
#' @section Estimation Procedure:
#' For each repetition \eqn{m = 1, \ldots, M}:
#' \enumerate{
#'   \item Observations are assigned to \code{n_groups} quantile-based groups
#'     by ranking the ensemble predictions from repetition \eqn{m}
#'     \emph{within each fold}. Group 1 contains the lowest predicted values
#'     and group \code{n_groups} the highest. Forming groups within folds
#'     ensures independence between group assignment and out-of-sample
#'     predictions.
#'   \item For each covariate, the mean is computed within the top group, the
#'     bottom group, the "else" group (all groups except the top), and all
#'     observations. Standard errors are computed analytically from group
#'     variances (or via regression with cluster-robust SEs when
#'     \code{individual_id} was specified).
#'   \item Three differences are computed: top minus bottom, top minus else,
#'     and top minus all.
#' }
#' The final reported estimates and standard errors are the simple averages of
#' the per-repetition estimates and standard errors across all \eqn{M}
#' repetitions.
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
#' @param group_on Character controlling which observations define the quantile
#'   cutoffs used to form groups. One of:
#'   \itemize{
#'     \item \code{"auto"} (default): Uses the ML training population. For
#'       \code{ensemble_hte_fit} this is all observations. For
#'       \code{ensemble_pred_fit} with \code{train_idx}, it is the training
#'       subset. This ensures an observation's group assignment does not change
#'       when you vary the analysis subset.
#'     \item \code{"all"}: Always form groups using all observations.
#'     \item \code{"analysis"}: Form groups within whatever observations are
#'       being analyzed (i.e. the \code{subset}).
#'   }
#'   Has no effect when \code{subset = NULL} and all observations are used.
#' 
#' @return An object of class `clan_results` containing:
#' \itemize{
#'   \item estimates: data.table with CLAN estimates averaged across repetitions,
#'     including group means (\code{mean_top}, \code{mean_bottom}, \code{mean_else},
#'     \code{mean_all}), their standard errors (\code{se_top}, \code{se_bottom},
#'     \code{se_else}, \code{se_all}), differences (\code{diff_top_bottom},
#'     \code{diff_top_else}, \code{diff_top_all}), difference standard errors
#'     (\code{se_diff_top_bottom}, \code{se_diff_top_else}, \code{se_diff_top_all}),
#'     t-values (\code{t_diff_top_bottom}, etc.), and p-values (\code{p_diff_top_bottom}, etc.)
#'   \item variables: variables analyzed
#'   \item n_groups: number of groups used
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item M: number of repetitions
#'   \item scaled: whether variables were scaled
#'   \item variable_sds: named numeric vector of per-variable standard deviations
#'     (NA for binary variables). Used by the \code{plot} method for on-the-fly rescaling.
#'   \item group_on: how groups are formed (\"auto\", \"all\", or \"analysis\")
#'   \item call: the function call
#' }
#' 
#' @examples
#' \donttest{
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#'
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
#' result <- clan(fit)
#' print(result)
#' plot(result)
#' }
#'
#' @export
clan <- function(ensemble_fit, variables = NULL, n_groups = 3, na_rm = FALSE, scale = FALSE, subset = NULL, group_on = c("auto", "all", "analysis")) {
  
  # Match group_on argument
  group_on <- match.arg(group_on)
  
  # Validate n_groups
  if (!is.numeric(n_groups) || length(n_groups) != 1 || n_groups < 2 || n_groups %% 1 != 0) {
    stop("n_groups must be an integer >= 2")
  }
  if (n_groups >= ensemble_fit$n) {
    stop("n_groups (", n_groups, ") must be smaller than the number of ",
         "observations in the data (n = ", ensemble_fit$n, "). ",
         "Each group needs at least a few observations for valid inference.")
  }
  
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
    variable_data <- as.data.table(as.data.frame(ensemble_fit$data)[, variables, drop = FALSE])
    var_names <- variables
  } else {
    stop("variables must be a character vector or a data.frame")
  }
  
  # Check for non-numeric covariates and skip them with a warning
  non_numeric_vars <- var_names[!vapply(variable_data[, ..var_names], is.numeric, logical(1))]
  if (length(non_numeric_vars) > 0) {
    warning("Non-numeric variable(s) skipped: ",
            paste(non_numeric_vars, collapse = ", "),
            ". CLAN computes group means and requires numeric variables.")
    var_names <- setdiff(var_names, non_numeric_vars)
    if (length(var_names) == 0) {
      stop("No numeric variables remaining for CLAN analysis after removing non-numeric variables.")
    }
  }
  
  # Compute per-variable SDs (always, for use by plot even when scale = FALSE)
  variable_sds <- vapply(var_names, function(var) {
    col_vals <- variable_data[[var]]
    unique_vals <- unique(col_vals[!is.na(col_vals)])
    is_binary <- length(unique_vals) == 2 && all(unique_vals %in% c(0, 1))
    if (is_binary) return(NA_real_)
    col_sd <- sd(col_vals, na.rm = FALSE)
    if (is.na(col_sd) || col_sd == 0) return(NA_real_)
    col_sd
  }, numeric(1))

  # Scale non-binary variables if requested
  if (scale) {
    for (var in var_names) {
      s <- variable_sds[[var]]
      if (!is.na(s)) {
        col_vals <- variable_data[[var]]
        col_mean <- mean(col_vals, na.rm = FALSE)
        variable_data[, (var) := (col_vals - col_mean) / s]
      }
    }
  }
  
  # Get targeted outcome from ensemble_fit
  targeted_outcome <- all.vars(ensemble_fit$formula)[1]
  n <- ensemble_fit$n
  
  # Check for train_idx in ensemble_fit
  has_train_idx <- !is.null(ensemble_fit$train_idx) &&
                   !is.null(ensemble_fit$n_train) &&
                   ensemble_fit$n_train < ensemble_fit$n
  
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
  n_fit <- if (has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("CLAN", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Resolve which observations define group cutoffs
  group_ref_idx <- .resolve_group_ref_idx(group_on, use_idx, train_idx)
  
  # For small-cell check, use the reference population (or full data if NULL)
  check_ref <- if (!is.null(group_ref_idx)) group_ref_idx else rep(TRUE, n)
  check_folds <- lapply(splits, function(s) s[check_ref])
  .check_small_cells(check_folds, n_groups, func_name = "CLAN")
  
  # Extract individual_id for cluster-robust SEs (if panel data)
  cluster_id <- ensemble_fit$individual_id
  
  # Check for NAs in analysis variables and warn if na_rm = FALSE
  if (!na_rm) {
    na_vars <- var_names[vapply(variable_data[use_idx, ..var_names], 
                                function(col) anyNA(col), logical(1))]
    if (length(na_vars) > 0) {
      warning("NA values detected in variable(s): ",
              paste(na_vars, collapse = ", "),
              ". Results will contain NAs for these variables. ",
              "Set na_rm = TRUE to exclude NA values from calculations.")
    }
  }
  
  clan_by_rep <- lapply(1:M, function(m) {
    .clan_single(
      predicted_values = predictions_list[[m]],
      fold = splits[[m]],
      n_groups = n_groups,
      variables = var_names,
      variable_data = variable_data,
      na_rm = na_rm,
      group_ref_idx = group_ref_idx,
      analysis_idx = if (!all(use_idx)) use_idx else NULL,
      cluster_id = cluster_id
    )
  })
  
  # Combine results across repetitions
  all_estimates <- rbindlist(clan_by_rep, idcol = "repetition")
  
  # Aggregate across repetitions: mean of estimates, mean of SE
  combined <- all_estimates[, .(
    mean_top = mean(mean_top),
    se_top = mean(se_top),
    mean_bottom = mean(mean_bottom),
    se_bottom = mean(se_bottom),
    mean_else = mean(mean_else),
    se_else = mean(se_else),
    mean_all = mean(mean_all),
    se_all = mean(se_all),
    diff_top_bottom = mean(diff_top_bottom),
    se_diff_top_bottom = mean(se_diff_top_bottom),
    diff_top_else = mean(diff_top_else),
    se_diff_top_else = mean(se_diff_top_else),
    diff_top_all = mean(diff_top_all),
    se_diff_top_all = mean(se_diff_top_all),
    n_reps = .N
  ), by = variable]
  
  # Add t-values and p-values for top-bottom difference
  combined[, `:=`(
    t_diff_top_bottom = diff_top_bottom / se_diff_top_bottom,
    p_diff_top_bottom = 2 * pnorm(-abs(diff_top_bottom / se_diff_top_bottom)),
    t_diff_top_else = diff_top_else / se_diff_top_else,
    p_diff_top_else = 2 * pnorm(-abs(diff_top_else / se_diff_top_else)),
    t_diff_top_all = diff_top_all / se_diff_top_all,
    p_diff_top_all = 2 * pnorm(-abs(diff_top_all / se_diff_top_all))
  )]
  
  # Construct clan_results object
  structure(
    list(
      estimates = combined,
      variables = var_names,
      n_groups = n_groups,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      n = n,
      n_train = n_fit,
      n_used = n_used,
      M = M,
      scaled = scale,
      variable_sds = variable_sds,
      group_on = group_on,
      call = cl
    ),
    class = "clan_results"
  )
}


#' Internal function to compute GAVS for a single repetition
#' 
#' @param Y Outcome vector (full length)
#' @param predicted_values Predicted values vector for this repetition (full length)
#' @param fold Fold assignment vector for this repetition (full length)
#' @param n_groups Number of groups for GAVS
#' @param restrict_by Optional restrict_by indicator for restricted ranking (factor/integer vector, full length)
#' @param group_ref_idx Optional logical vector (same length as Y). When provided,
#'   group cutoffs are computed from only the reference observations, then applied
#'   to all observations. If \code{NULL}, all observations are used for group formation.
#' @param analysis_idx Optional logical vector (same length as Y). When provided,
#'   the regression is run only on these observations (after group assignment).
#'   If \code{NULL}, all observations are used.
#' @param cluster_id Optional vector of cluster identifiers for cluster-robust SEs.
#'   When provided, passed to \code{reg_from_formula} for clustered inference.
#' 
#' @return List with group_estimates, tests (top_bottom, all, top_all), and reg object
#' @keywords internal
.gavs_single <- function(Y, predicted_values, fold, n_groups = 5, restrict_by = NULL,
                         group_ref_idx = NULL, analysis_idx = NULL,
                         cluster_id = NULL) {
  
  # Build data.table for GAVS regression
  dt_gavs <- data.table(
    Y = Y,
    predicted_values = predicted_values,
    fold = fold
  )
  
  # Create groups within each fold (and optionally within strata for restricted ranking)
  if (is.null(restrict_by)) {
    dt_gavs[, group := {
      ref <- if (!is.null(group_ref_idx)) group_ref_idx[.I] else NULL
      as.integer(create_groups_by_reference(predicted_values, n_groups, ref_mask = ref))
    }, by = fold]
  } else {
    dt_gavs[, strata := restrict_by]
    dt_gavs[, group := {
      ref <- if (!is.null(group_ref_idx)) group_ref_idx[.I] else NULL
      as.integer(create_groups_by_reference(predicted_values, n_groups, ref_mask = ref))
    }, by = .(fold, strata)]
  }
  
  # Subset to analysis observations (after group formation)
  if (!is.null(analysis_idx)) {
    dt_gavs <- dt_gavs[analysis_idx, ]
    if (!is.null(cluster_id)) {
      cluster_id <- cluster_id[analysis_idx]
    }
    # Check for empty groups after subsetting
    observed_groups <- sort(unique(dt_gavs$group))
    empty_groups <- setdiff(1:n_groups, observed_groups)
    if (length(empty_groups) > 0) {
      stop("After subsetting, group(s) ", paste(empty_groups, collapse = ", "),
           " have no observations. This happens when groups are formed on all ",
           "observations but the analysis subset excludes entire groups. ",
           'To form groups within the subset instead, use group_on = "analysis".')
    }
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
  reg <- reg_from_formula(formula_str, dt_gavs, cluster_id = cluster_id)
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
#' @section Estimation Procedure:
#' For each repetition \eqn{m = 1, \ldots, M}:
#' \enumerate{
#'   \item Observations are assigned to \code{n_groups} quantile-based groups by
#'     ranking the ensemble predictions from repetition \eqn{m} \emph{within each
#'     fold}. Group 1 contains the lowest predicted values and group
#'     \code{n_groups} the highest. Forming groups within folds ensures that
#'     group assignment is independent of the model used to generate predictions
#'     for that observation (since predictions are out-of-sample within each fold).
#'   \item A single ordinary least squares regression is run on all observations:
#'     \deqn{Y_i = \sum_{g=1}^{G} \mu_g \, \mathbf{1}\{i \in g\} + \varepsilon_i}
#'     This is a regression of Y on group dummies (no intercept), so each
#'     \eqn{\mu_g} directly estimates the average outcome in group \eqn{g}.
#'   \item HC1 heteroskedasticity-robust standard errors are computed
#'     (or cluster-robust SEs when \code{individual_id} was specified).
#' }
#' The final reported estimates and standard errors are the simple averages of
#' the per-repetition estimates and standard errors across all \eqn{M} repetitions.
#' 
#' Three summary tests are reported:
#' \itemize{
#'   \item \strong{Top-Bottom}: difference between top and bottom group averages.
#'   \item \strong{All}: weighted average of all groups (estimates the overall mean).
#'   \item \strong{Top-All}: difference between the top group and the overall mean.
#' }
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
#' @param restrict_by Optional. Stratification variable for restricted ranking:
#'   \itemize{
#'     \item NULL (default): unrestricted ranking across full sample within folds
#'     \item Character string: column name in the \code{data} for stratified ranking
#'     \item Numeric/factor vector: group indicator (must have same length as data)
#'   }
#'   When specified, predicted values are ranked within each stratum (and fold),
#'   rather than across the full sample.
#' @param group_on Character controlling which observations define the quantile
#'   cutoffs used to form groups. One of:
#'   \itemize{
#'     \item \code{"auto"} (default): Uses the ML training population. For
#'       \code{ensemble_hte_fit} this is all observations. For
#'       \code{ensemble_pred_fit} with \code{train_idx}, it is the training
#'       subset. This ensures an observation's group assignment does not change
#'       when you vary the analysis subset.
#'     \item \code{"all"}: Always form groups using all observations.
#'     \item \code{"analysis"}: Form groups within whatever observations are
#'       being analyzed (i.e. the \code{subset}).
#'   }
#'   Has no effect when \code{subset = NULL} and all observations are used.
#' 
#' @return An object of class \code{gavs_results} containing:
#' \itemize{
#'   \item estimates: data.table with GAVS estimates averaged across repetitions.
#'     Columns: \code{group} (integer group index, 1 = lowest predicted values),
#'     \code{estimate} (group-specific mean outcome), \code{se} (standard error),
#'     \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item top_bottom: data.table with the top-bottom difference test.
#'     Columns: \code{estimate}, \code{se}, \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item all: data.table with the overall mean (weighted avg of all groups).
#'     Columns: \code{estimate}, \code{se}, \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item top_all: data.table with the top minus average test.
#'     Columns: \code{estimate}, \code{se}, \code{n_reps}, \code{t_value}, \code{p_value}
#'   \item n_groups: number of groups used
#'   \item outcome: the outcome variable used for GAVS
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item restrict_by: the restrict_by variable used (if any)
#'   \item group_on: how groups are formed ("auto", "all", or "analysis")
#'   \item n_used: number of observations used
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#' 
#' @examples
#' \donttest{
#' data(microcredit)
#' covars <- c("age", "gender", "education", "hhinc_yrly_base",
#'             "css_creditscorefinal")
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", covars)]
#'
#' fit <- ensemble_hte(
#'   hhinc_yrly_end ~ ., treatment = treat, data = dat,
#'   prop_score = microcredit$prop_score,
#'   algorithms = c("lm", "grf"), M = 3, K = 3
#' )
#' result <- gavs(fit, n_groups = 3)
#' print(result)
#' plot(result)
#' }
#'
#' @export
gavs <- function(ensemble_fit, n_groups = 3, outcome = NULL, subset = NULL, 
                 restrict_by = NULL, group_on = c("auto", "all", "analysis")) {
  
  # Match group_on argument
  group_on <- match.arg(group_on)
  
  # Validate n_groups
  if (!is.numeric(n_groups) || length(n_groups) != 1 || n_groups < 2 || n_groups %% 1 != 0) {
    stop("n_groups must be an integer >= 2")
  }
  if (n_groups >= ensemble_fit$n) {
    stop("n_groups (", n_groups, ") must be smaller than the number of ",
         "observations in the data (n = ", ensemble_fit$n, "). ",
         "Each group needs at least a few observations for valid inference.")
  }
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$ite[[m]])
    fit_type <- "hte"
    has_train_idx <- !is.null(ensemble_fit$train_idx) &&
                     !is.null(ensemble_fit$n_train) &&
                     ensemble_fit$n_train < ensemble_fit$n
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
  
  # Process restrict_by argument
  strata_vec <- NULL
  strata_var <- NULL
  if (!is.null(restrict_by)) {
    if (is.character(restrict_by) && length(restrict_by) == 1) {
      # Single string: treat as column name
      if (!restrict_by %in% names(ensemble_fit$data)) {
        stop(paste0("restrict_by '", restrict_by, "' not found in the data"))
      }
      strata_vec <- ensemble_fit$data[[restrict_by]]
      strata_var <- restrict_by
    } else if (is.numeric(restrict_by) || is.factor(restrict_by) || is.character(restrict_by)) {
      # Vector passed directly (numeric, factor, or character)
      if (length(restrict_by) != n) {
        stop(paste0("restrict_by vector has length ", length(restrict_by), 
                    " but data has ", n, " rows"))
      }
      strata_vec <- restrict_by
      strata_var <- "custom_strata"
    } else {
      stop("restrict_by must be NULL, a column name (string), or a vector (numeric/factor/character)")
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
    n_na <- sum(is.na(Y[use_idx]))
    stop("The outcome variable has ", n_na, " NA value(s) among the ", sum(use_idx),
         " observations being used. ",
         if (has_train_idx) "If you used 'train_idx' during fitting, try subset = 'train' or subset = NULL to restrict to observations with observed outcomes." else "Remove or impute missing values before analysis.")
  }
  
  # Validate no NAs in strata for used observations (if restrict_by provided)
  if (!is.null(strata_vec) && any(is.na(strata_vec[use_idx]))) {
    stop("restrict_by has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("GAVS", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Resolve which observations define group cutoffs
  group_ref_idx <- .resolve_group_ref_idx(group_on, use_idx, train_idx)
  
  # For small-cell check, use the reference population (or full data if NULL)
  check_ref <- if (!is.null(group_ref_idx)) group_ref_idx else rep(TRUE, n)
  check_folds <- lapply(splits, function(s) s[check_ref])
  check_strata <- if (!is.null(strata_vec)) strata_vec[check_ref] else NULL
  .check_small_cells(check_folds, n_groups, restrict_by = check_strata, func_name = "GAVS")
  
  # Extract individual_id for cluster-robust SEs (if panel data)
  cluster_id <- ensemble_fit$individual_id
  
  gavs_by_rep <- lapply(1:M, function(m) {
    .gavs_single(
      Y = Y,
      predicted_values = predictions_list[[m]],
      fold = splits[[m]],
      n_groups = n_groups,
      restrict_by = strata_vec,
      group_ref_idx = group_ref_idx,
      analysis_idx = if (!all(use_idx)) use_idx else NULL,
      cluster_id = cluster_id
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
    result <- all_tests[, .(estimate = mean(estimate), se = mean(se), n_reps = .N)]
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
      restrict_by = strata_var,
      group_on = group_on,
      n = n,
      n_train = n_fit,
      n_used = n_used,
      M = M,
      call = cl
    ),
    class = "gavs_results"
  )
}


#' Print method for clan_results objects
#'
#' @description
#' Displays CLAN results in two compact panels: group means (with standard
#' errors) and differences from the top group (with significance stars and
#' standard errors). When there are many variables, output is truncated to
#' \code{max_rows} variables; the full results are always accessible via
#' \code{$estimates}.
#'
#' @param x An object of class \code{clan_results} from \code{clan()}
#' @param max_rows Maximum number of variables to display (default: 15).
#'   Set to \code{Inf} to show all.
#' @param ... Additional arguments (currently unused)
#' @export
print.clan_results <- function(x, max_rows = 15, ...) {
  cat("CLAN Results (Classification Analysis)\n")
  cat("=======================================\n\n")
  
  fit_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") " (prediction)" else " (treatment effects)"
  cat("Outcome: ", x$targeted_outcome, fit_label, " | ", sep = "")
  cat("Groups: ", x$n_groups, " | ", sep = "")
  cat("Reps: ", x$M, "\n", sep = "")
  if (!is.null(x$scaled) && x$scaled) cat("(Variables scaled to SD = 1)\n")
  
  # Data usage block
  if (!is.null(x$n)) {
    .print_data_usage(x$n, x$n_train, x$n_used, x$group_on, x$n_groups)
  }
  
  # Determine label based on fit type
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "Y" else "ITE"
  
  out <- copy(x$estimates)
  n_vars <- nrow(out)
  
  # Limit number of rows shown
  show_all <- n_vars <= max_rows
  if (!show_all) {
    out <- out[1:max_rows]
  }
  
  # --- Compute column widths dynamically from content ---
  max_name_len <- max(nchar(out$variable), 8)
  nw <- min(max_name_len, 20)  # name width, capped at 20
  
  # Panel 1 column width: must fit headers, formatted values, and formatted SEs
  p1_val_strs <- sprintf("%.2f", c(out$mean_top, out$mean_bottom,
                                    out$mean_else, out$mean_all))
  p1_se_strs  <- sprintf("(%.2f)", c(out$se_top, out$se_bottom,
                                      out$se_else, out$se_all))
  cw <- max(nchar(p1_val_strs), nchar(p1_se_strs),
            nchar(c("Top", "Bottom", "Else", "All")))
  
  # Panel 2 column width: must fit value+stars, SEs, and headers
  star_w <- 3
  p2_num_strs <- sprintf("%.2f", c(out$diff_top_bottom, out$diff_top_else,
                                    out$diff_top_all))
  p2_se_strs  <- sprintf("(%.2f)", c(out$se_diff_top_bottom, out$se_diff_top_else,
                                      out$se_diff_top_all))
  p2_headers  <- c("Top-Bot", "Top-Else", "Top-All")
  num_w <- max(nchar(p2_num_strs))
  dw <- max(num_w + star_w, max(nchar(p2_se_strs)), max(nchar(p2_headers)))
  num_w <- dw - star_w
  
  # --- Panel 1: Group Means ---
  p1_total <- nw + 4 * (2 + cw)
  cat(sprintf("\nGroup Means (by predicted %s):\n", pred_label))
  cat(sprintf("  %-*s  %*s  %*s  %*s  %*s\n",
              nw, "", cw, "Top", cw, "Bottom", cw, "Else", cw, "All"))
  cat("  ", strrep("-", p1_total), "\n", sep = "")
  for (i in 1:nrow(out)) {
    vname <- substr(out$variable[i], 1, nw)
    cat(sprintf("  %-*s  %*s  %*s  %*s  %*s\n",
                nw, vname,
                cw, sprintf("%.2f", out$mean_top[i]),
                cw, sprintf("%.2f", out$mean_bottom[i]),
                cw, sprintf("%.2f", out$mean_else[i]),
                cw, sprintf("%.2f", out$mean_all[i])))
    cat(sprintf("  %-*s  %*s  %*s  %*s  %*s\n",
                nw, "",
                cw, sprintf("(%.2f)", out$se_top[i]),
                cw, sprintf("(%.2f)", out$se_bottom[i]),
                cw, sprintf("(%.2f)", out$se_else[i]),
                cw, sprintf("(%.2f)", out$se_all[i])))
  }
  
  # --- Panel 2: Differences from Top Group ---
  p2_total <- nw + 3 * (2 + dw)
  cat(sprintf("\nDifferences from Top Group:\n"))
  cat(sprintf("  %-*s  %*s  %*s  %*s\n",
              nw, "", dw, "Top-Bot", dw, "Top-Else", dw, "Top-All"))
  cat("  ", strrep("-", p2_total), "\n", sep = "")
  for (i in 1:nrow(out)) {
    vname <- substr(out$variable[i], 1, nw)
    s_tb <- get_stars(out$p_diff_top_bottom[i])
    s_te <- get_stars(out$p_diff_top_else[i])
    s_ta <- get_stars(out$p_diff_top_all[i])
    cat(sprintf("  %-*s  %*.2f%-*s  %*.2f%-*s  %*.2f%-*s\n",
                nw, vname,
                num_w, out$diff_top_bottom[i], star_w, s_tb,
                num_w, out$diff_top_else[i], star_w, s_te,
                num_w, out$diff_top_all[i], star_w, s_ta))
    cat(sprintf("  %-*s  %*s  %*s  %*s\n",
                nw, "",
                dw, sprintf("(%.2f)", out$se_diff_top_bottom[i]),
                dw, sprintf("(%.2f)", out$se_diff_top_else[i]),
                dw, sprintf("(%.2f)", out$se_diff_top_all[i])))
  }
  
  if (!show_all) {
    cat(sprintf("\n  ... %d of %d variables shown. Access all via $estimates\n",
                max_rows, n_vars))
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
  cat("Repetitions: ", x$M, "\n", sep = "")
  
  # Data usage block
  if (!is.null(x$n)) {
    .print_data_usage(x$n, x$n_train, x$n_used, x$group_on, x$n_groups)
  }
  cat("\nGroup Average Outcomes (groups by ", pred_label, "):\n\n", sep = "")
  
  # Format output table for group estimates
  out <- copy(x$estimates)
  out[, stars := get_stars(p_value)]
  
  # Print as formatted text
  cat(sprintf("  %5s  %10s  %10s  %8s  %10s\n", 
              "Group", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %5d  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
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
  cat(sprintf("  %12s  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
              "Top-Bottom", tb$estimate, tb$se, tb$t_value, tb$p_value, tb_stars))
  
  # Top-All test (if available)
  if (!is.null(x$top_all)) {
    ta <- x$top_all
    ta_stars <- get_stars(ta$p_value)
    cat(sprintf("  %12s  %10.2f  %10.2f  %8.2f  %10.3f %s\n",
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
#' @param scale Logical. If \code{TRUE} (default), non-binary variables are
#'   rescaled to standard-deviation units for plotting, making coefficients
#'   comparable across variables with different scales. This rescaling is
#'   applied on the fly regardless of whether \code{scale} was set in the
#'   original \code{clan()} call. Set \code{FALSE} to plot in the original
#'   units of each variable.
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar geom_hline
#'   labs theme_minimal theme element_text coord_flip scale_fill_manual
#' @export
plot.clan_results <- function(x, comparison = c("top_bottom", "top_else", "top_all"), alpha = 0.05, scale = TRUE, ...) {
  
  comparison <- match.arg(comparison)
  
  # Compute critical value
  z_crit <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  # Determine label based on fit type
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "predicted ITE"
  
  # Select appropriate columns based on comparison
  plot_data <- copy(x$estimates)
  
  # Rescale on the fly if needed
  # If scale = TRUE and clan was not scaled, divide diffs/SEs by SD
  # If scale = FALSE and clan was scaled, multiply diffs/SEs by SD
  was_scaled <- isTRUE(x$scaled)
  sds <- x$variable_sds
  if (!is.null(sds) && scale != was_scaled) {
    diff_cols <- c("diff_top_bottom", "diff_top_else", "diff_top_all")
    se_cols   <- c("se_diff_top_bottom", "se_diff_top_else", "se_diff_top_all")
    mean_cols <- c("mean_top", "mean_bottom", "mean_else", "mean_all")
    se_mean_cols <- c("se_top", "se_bottom", "se_else", "se_all")
    all_cols <- c(diff_cols, se_cols, mean_cols, se_mean_cols)
    for (i in seq_len(nrow(plot_data))) {
      v <- plot_data$variable[i]
      s <- sds[[v]]
      if (!is.na(s)) {
        if (scale && !was_scaled) {
          # Convert from raw to SD units: divide by SD
          for (col in all_cols) {
            set(plot_data, i, col, plot_data[[col]][i] / s)
          }
        } else {
          # Convert from SD units to raw: multiply by SD
          for (col in all_cols) {
            set(plot_data, i, col, plot_data[[col]][i] * s)
          }
        }
      }
    }
  }
  
  # Determine the effective scaling state for the caption
  plot_scaled <- if (!is.null(sds)) scale else was_scaled
  
  if (comparison == "top_bottom") {
    plot_data[, `:=`(
      diff = diff_top_bottom,
      se_diff = se_diff_top_bottom,
      comparison_label = "Top vs Bottom"
    )]
  } else if (comparison == "top_else") {
    plot_data[, `:=`(
      diff = diff_top_else,
      se_diff = se_diff_top_else,
      comparison_label = "Top vs Else"
    )]
  } else {
    plot_data[, `:=`(
      diff = diff_top_all,
      se_diff = se_diff_top_all,
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
      y = if (plot_scaled) "Difference in Means (SD units)" else "Difference in Means",
      caption = paste0(conf_level, "% confidence intervals. ", x$M, " repetitions. ", x$n_groups, " groups.",
                       if (plot_scaled) " Variables scaled." else "")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
  
  print(p)
  invisible(p)
}

