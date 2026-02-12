#' Internal function to compute BLP for prediction (single repetition)
#' 
#' @param Y Outcome vector
#' @param predicted_y Predicted Y vector for this repetition
#' 
#' @return data.table with BLP estimates
#' @keywords internal
.blp_pred_single <- function(Y, predicted_y) {
  
  # Build data.table for BLP regression
  dt_blp <- data.table(
    Y = Y,
    predicted_y = predicted_y
  )
  
  # Build regression formula
  # BLP regression: Y ~ predicted_y
  formula_str <- "Y ~ predicted_y"
  
  # Run regression with robust SEs
  reg <- reg_from_formula(formula_str, dt_blp)
  coef_blp <- reg$coef
  
  # Extract estimates and SEs
  blp_estimates <- data.table(
    term = c("intercept", "beta"),
    estimate = coef_blp[c("(Intercept)", "predicted_y"), 1],
    se = coef_blp[c("(Intercept)", "predicted_y"), 2]
  )
  
  list(
    estimates = blp_estimates
  )
}


#' Compute BLP (Best Linear Predictor) for Prediction
#' 
#' @description
#' Computes the Best Linear Predictor (BLP) of the outcome using the ensemble
#' predictions (from `ensemble_pred`) or ITE predictions (from `ensemble_hte`).
#' This is a simple regression of Y on the predicted values.
#' 
#' This function implements the multiple-split estimation strategy developed in
#' Fava (2025), which combines predictions from multiple machine learning algorithms
#' into an ensemble and averages BLP estimates across M repetitions of K-fold
#' cross-fitting to improve statistical power.
#' 
#' The method uses ordinary least squares regression of Y on predicted values.
#' A coefficient (beta) close to 1 indicates good calibration; significantly
#' different from 1 suggests over- or under-prediction.
#' 
#' When using `ensemble_hte_fit` objects, this allows testing whether predicted
#' treatment effects correlate with a different outcome variable (e.g., an endline
#' measure that may only be observed for a subset of the data).
#' 
#' @references
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class `ensemble_pred_fit` from `ensemble_pred()`
#'   or `ensemble_hte_fit` from `ensemble_hte()`.
#' @param outcome Either:
#'   \itemize{
#'     \item NULL (default): uses the same outcome as in ensemble fitting
#'     \item Character string: column name in the `data` used in ensemble fitting
#'     \item Numeric vector: custom outcome variable (must have same length as data,
#'           or same length as subset if subset is provided)
#'   }
#'   This allows computing BLP for a different outcome than the one used for prediction.
#' @param subset For `ensemble_pred_fit` with subset training (train_idx), controls 
#'   which observations to use:
#'   \itemize{
#'     \item NULL (default): uses training observations if default outcome and subset 
#'           training was used, otherwise uses all observations
#'     \item "train": uses only training observations (requires train_idx in ensemble_fit)
#'     \item "all": uses all observations
#'   }
#'   For `ensemble_hte_fit`, this can be a logical or integer vector specifying which
#'   observations to include.
#' 
#' @return An object of class `blp_pred_results` containing:
#' \itemize{
#'   \item estimates: data.table with BLP estimates averaged across repetitions
#'   \item outcome: outcome variable used
#'   \item targeted_outcome: original outcome from ensemble fitting
#'   \item fit_type: "hte" or "pred" depending on input
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
#' Y <- 2 * X1 + X2 + rnorm(n)
#' data <- data.frame(Y = Y, X1 = X1, X2 = X2)
#'
#' fit <- ensemble_pred(Y ~ X1 + X2, data = data,
#'                      algorithms = c("lm", "ranger"), M = 3, K = 3)
#' result <- blp_pred(fit)
#' print(result)
#' }
#'
#' @export
blp_pred <- function(ensemble_fit, outcome = NULL, subset = NULL) {
  
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
  
  # Determine which observations to use
  using_default_outcome <- is.null(outcome)
  
  if (fit_type == "hte") {
    # For HTE: subset can be a logical or integer vector
    if (is.null(subset)) {
      use_idx <- rep(TRUE, n)
    } else if (is.logical(subset)) {
      if (length(subset) != n) {
        stop(paste0("subset logical vector has length ", length(subset), 
                    " but data has ", n, " rows"))
      }
      use_idx <- subset
    } else if (is.numeric(subset)) {
      # Integer indices
      if (any(subset < 1) || any(subset > n)) {
        stop("subset indices out of range")
      }
      use_idx <- rep(FALSE, n)
      use_idx[subset] <- TRUE
    } else {
      stop("For ensemble_hte_fit, subset must be NULL, a logical vector, or integer indices")
    }
  } else {
    # For pred: handle subset argument like gavs
    if (is.null(subset)) {
      # Smart default: use training obs if default outcome and subset training was used
      if (using_default_outcome && has_train_idx) {
        use_idx <- ensemble_fit$train_idx
      } else {
        use_idx <- rep(TRUE, n)
      }
    } else if (subset == "train") {
      if (!has_train_idx) {
        stop("subset = 'train' specified but no train_idx was used in ensemble_pred()")
      }
      use_idx <- ensemble_fit$train_idx
    } else if (subset == "all") {
      use_idx <- rep(TRUE, n)
    } else {
      stop("For ensemble_pred_fit, subset must be NULL, 'train', or 'all'")
    }
  }
  
  n_used <- sum(use_idx)
  
  # Determine which outcome to use for BLP
  if (is.null(outcome)) {
    Y <- ensemble_fit$Y[use_idx]
    outcome_var <- targeted_outcome
  } else if (is.character(outcome)) {
    if (length(outcome) != 1) {
      stop("outcome must be a single column name or a numeric vector")
    }
    if (!outcome %in% names(ensemble_fit$data)) {
      stop(paste0("outcome '", outcome, "' not found in the data"))
    }
    Y <- ensemble_fit$data[[outcome]][use_idx]
    outcome_var <- outcome
  } else if (is.numeric(outcome)) {
    # outcome can be full length (will be subsetted) or already subsetted
    if (length(outcome) == n) {
      Y <- outcome[use_idx]
    } else if (length(outcome) == n_used) {
      Y <- outcome
    } else {
      stop(paste0("outcome vector has length ", length(outcome), 
                  " but expected ", n, " (full data) or ", n_used, " (subset)"))
    }
    outcome_var <- "custom_outcome"
  } else {
    stop("outcome must be NULL, a character string, or a numeric vector")
  }
  
  # Validate no NAs in outcome for used observations
  if (any(is.na(Y))) {
    stop("outcome has NA values for the observations being used")
  }
  
  # Extract components from ensemble_fit
  M <- ensemble_fit$M
  
  # Compute BLP for each repetition
  blp_by_rep <- lapply(1:M, function(m) {
    .blp_pred_single(
      Y = Y,
      predicted_y = predictions_list[[m]][use_idx]
    )
  })
  
  # Combine estimates across repetitions
  all_estimates <- rbindlist(lapply(blp_by_rep, `[[`, "estimates"), idcol = "repetition")
  
  # Aggregate across repetitions: mean of estimates, mean of SE
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
  
  # Construct blp_pred_results object
  structure(
    list(
      estimates = combined,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      n_used = n_used,
      M = M,
      call = cl
    ),
    class = "blp_pred_results"
  )
}




# ==============================================================================
# Print methods
# ==============================================================================

#' Print method for blp_pred_results objects
#' @param x An object of class \code{blp_pred_results} from \code{blp_pred()}
#' @param ... Additional arguments (currently unused)
#' @export
print.blp_pred_results <- function(x, ...) {
  cat("BLP Results (Best Linear Predictor - Prediction)\n")
  cat("=================================================\n\n")
  
  # Show fit type
  fit_label <- if (!is.null(x$fit_type) && x$fit_type == "hte") {
    "HTE (ensemble_hte)"
  } else {
    "Prediction (ensemble_pred)"
  }
  cat("Fit type: ", fit_label, "\n", sep = "")
  
  # Show outcome information
  cat("Outcome analyzed: ", x$outcome, "\n", sep = "")
  if (x$outcome != x$targeted_outcome) {
    pred_label <- if (!is.null(x$fit_type) && x$fit_type == "hte") "ITE" else "predicted Y"
    cat("  (Predictions based on: ", x$targeted_outcome, " ", pred_label, ")\n", sep = "")
  }
  cat("Observations used: ", x$n_used, "\n", sep = "")
  cat("Repetitions: ", x$M, "\n", sep = "")
  
  cat("\nCoefficients:\n")
  cat("  intercept: Regression intercept (0 = well-calibrated)\n")
  cat("  beta: Prediction loading (1 = well-calibrated)\n\n")
  
  # Format output table
  out <- copy(x$estimates)
  out[, stars := get_stars(p_value)]
  
  # Print as formatted text
  cat(sprintf("  %10s  %10s  %10s  %8s  %10s\n", 
              "Term", "Estimate", "Std.Error", "t value", "Pr(>|t|)"))
  cat("  ", paste(rep("-", 52), collapse = ""), "\n", sep = "")
  for (i in 1:nrow(out)) {
    cat(sprintf("  %10s  %10.4f  %10.4f  %8.3f  %10.4f %s\n",
                out$term[i], out$estimate[i], out$se[i], 
                out$t_value[i], out$p_value[i], out$stars[i]))
  }
  
  cat("\n---\nSignif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  invisible(x)
}
