
#' Compute cross-regression covariance
#' 
#' @description
#' Computes the covariance between coefficient estimates from two different
#' regressions using the same data. This is needed to correctly compute the
#' standard error of the difference between estimates from two regressions.
#' 
#' @param X1 Model matrix from regression 1
#' @param X2 Model matrix from regression 2
#' @param resid1 Residuals from regression 1
#' @param resid2 Residuals from regression 2
#' @param weights Optional WLS weights (same for both regressions). When provided,
#'   uses the WLS sandwich formula with weighted bread and meat. NULL for OLS.
#' 
#' @return Covariance matrix between coefficients of regression 1 (rows) and regression 2 (cols)
#' @keywords internal
.calc_cov_cross_reg <- function(X1, X2, resid1, resid2, weights = NULL) {
  n <- nrow(X1)
  k <- ncol(X1)
  
  if (is.null(weights)) {
    # OLS: bread = (X'X)^{-1}, meat uses e1*e2
    bread1 <- solve(crossprod(X1))
    bread2 <- solve(crossprod(X2))
    meat <- crossprod(X1, X2 * (resid1 * resid2))
  } else {
    # WLS: bread = (X'WX)^{-1}, meat uses w^2 * e1 * e2
    bread1 <- solve(crossprod(X1, X1 * weights))
    bread2 <- solve(crossprod(X2, X2 * weights))
    meat <- crossprod(X1, X2 * (weights^2 * resid1 * resid2))
  }
  
  # HC1 correction factor (consistent with sandwich::vcovHC type="HC1")
  hc1 <- n / (n - k)
  
  cov <- hc1 * bread1 %*% meat %*% bread2
  return(cov)
}


#' Compute SE of difference between two regression coefficients
#' 
#' @description
#' Computes the estimate and standard error of the difference between
#' corresponding coefficients from two regressions (reg1 - reg2).
#' 
#' @param reg1 Result from reg_from_formula (list with coef, model_matrix, residuals)
#' @param reg2 Result from reg_from_formula (list with coef, model_matrix, residuals)
#' @param coef_name Name of the coefficient to compare
#' 
#' @return Named vector with 'estimate' and 'se'
#' @keywords internal
.calc_se_diff_regs <- function(reg1, reg2, coef_name) {
  # Get coefficient estimates and SEs
  est1 <- reg1$coef[coef_name, "Estimate"]
  se1 <- reg1$coef[coef_name, "Std. Error"]
  est2 <- reg2$coef[coef_name, "Estimate"]
  se2 <- reg2$coef[coef_name, "Std. Error"]
  
  # Compute cross-covariance (passing weights for WLS consistency)
  cross_cov <- .calc_cov_cross_reg(
    reg1$model_matrix, reg2$model_matrix,
    reg1$residuals, reg2$residuals,
    reg1$weights
  )
  
  # Get index for the coefficient in each regression
  idx1 <- which(rownames(reg1$coef) == coef_name)
  idx2 <- which(rownames(reg2$coef) == coef_name)
  
  # Estimate of difference
  diff_estimate <- est1 - est2
  
  # Variance of difference: var(A - B) = var(A) + var(B) - 2*cov(A,B)
  diff_var <- se1^2 + se2^2 - 2 * cross_cov[idx1, idx2]
  diff_se <- sqrt(max(diff_var, 0))  # Ensure non-negative
  
  c(estimate = diff_estimate, se = diff_se)
}


#' Compute difference between two regression results
#' 
#' @description
#' Computes the difference in group estimates and test statistics between two
#' regression results (e.g., unrestricted vs restricted ranking). Uses cross-regression
#' covariance for proper SE calculation when both regressions use the same data.
#' 
#' @param res1 Result from .gates_single or .gavs_single (first regression)
#' @param res2 Result from .gates_single or .gavs_single (second regression)
#' @param n_groups Number of groups
#' 
#' @return List with group_diff (differences in group estimates), top_bottom_diff, all_diff, top_all_diff
#' @keywords internal
.compute_reg_diff <- function(res1, res2, n_groups) {
  
  reg1 <- res1$reg
  reg2 <- res2$reg
  
  # Compute cross-covariance matrix (passing weights for WLS consistency)
  cross_cov <- .calc_cov_cross_reg(
    reg1$model_matrix, reg2$model_matrix,
    reg1$residuals, reg2$residuals,
    reg1$weights
  )
  
  group_cols <- paste0("group_", 1:n_groups)
  
  # Compute difference in group estimates
  group_diff <- data.table(group = 1:n_groups)
  for (g_idx in 1:n_groups) {
    col <- paste0("group_", g_idx)
    est1 <- reg1$coef[col, "Estimate"]
    se1 <- reg1$coef[col, "Std. Error"]
    est2 <- reg2$coef[col, "Estimate"]
    se2 <- reg2$coef[col, "Std. Error"]
    
    # Get indices
    idx1 <- which(rownames(reg1$coef) == col)
    idx2 <- which(rownames(reg2$coef) == col)
    
    # Variance of difference
    diff_var <- se1^2 + se2^2 - 2 * cross_cov[idx1, idx2]
    
    group_diff[g_idx, `:=`(
      estimate = est1 - est2,
      se = sqrt(max(diff_var, 0))
    )]
  }
  
  # Compute difference of difference tests (e.g., diff in top-bottom between two regressions)
  # For top-bottom: (G_n1 - G_11) - (G_n2 - G_12) = (G_n1 - G_n2) - (G_11 - G_12)
  # Need the full covariance structure
  
  top_col <- paste0("group_", n_groups)
  bottom_col <- "group_1"
  
  top_idx1 <- which(rownames(reg1$coef) == top_col)
  bottom_idx1 <- which(rownames(reg1$coef) == bottom_col)
  top_idx2 <- which(rownames(reg2$coef) == top_col)
  bottom_idx2 <- which(rownames(reg2$coef) == bottom_col)
  
  # Estimates of top-bottom for each regression
  tb1 <- reg1$coef[top_col, 1] - reg1$coef[bottom_col, 1]
  tb2 <- reg2$coef[top_col, 1] - reg2$coef[bottom_col, 1]
  diff_tb_est <- tb1 - tb2
  
  # Build covariance structure for [top1, bottom1, top2, bottom2]
  vcov1 <- reg1$vcov
  vcov2 <- reg2$vcov
  
  full_cov <- matrix(0, 4, 4)
  # Variances and covariances within reg1
  full_cov[1, 1] <- vcov1[top_col, top_col]
  full_cov[2, 2] <- vcov1[bottom_col, bottom_col]
  full_cov[1, 2] <- full_cov[2, 1] <- vcov1[top_col, bottom_col]
  # Variances and covariances within reg2
  full_cov[3, 3] <- vcov2[top_col, top_col]
  full_cov[4, 4] <- vcov2[bottom_col, bottom_col]
  full_cov[3, 4] <- full_cov[4, 3] <- vcov2[top_col, bottom_col]
  # Cross-covariances
  full_cov[1, 3] <- full_cov[3, 1] <- cross_cov[top_idx1, top_idx2]
  full_cov[1, 4] <- full_cov[4, 1] <- cross_cov[top_idx1, bottom_idx2]
  full_cov[2, 3] <- full_cov[3, 2] <- cross_cov[bottom_idx1, top_idx2]
  full_cov[2, 4] <- full_cov[4, 2] <- cross_cov[bottom_idx1, bottom_idx2]
  
  # Contrast for (top1 - bottom1) - (top2 - bottom2) = top1 - bottom1 - top2 + bottom2
  contrast_tb <- c(1, -1, -1, 1)
  diff_tb_var <- as.numeric(t(contrast_tb) %*% full_cov %*% contrast_tb)
  diff_tb_se <- sqrt(max(diff_tb_var, 0))
  
  top_bottom_diff <- data.table(
    estimate = diff_tb_est,
    se = diff_tb_se
  )
  
  # --- Difference in "all" test ---
  # all = weighted average of groups, each side uses its own group proportions
  weights1 <- res1$all$weights[[1]]  # Group proportions for reg1
  weights2 <- res2$all$weights[[1]]  # Group proportions for reg2
  # Difference: sum(w1_g * G_g1) - sum(w2_g * G_g2)
  all_diff_est <- sum(weights1 * res1$reg$coef[group_cols, 1]) - 
                  sum(weights2 * res2$reg$coef[group_cols, 1])
  
  # Build full covariance for all groups from both regressions
  # (g1_1, g2_1, ..., gn_1, g1_2, g2_2, ..., gn_2)
  full_cov_all <- matrix(0, 2 * n_groups, 2 * n_groups)
  # Reg1 block
  for (i in 1:n_groups) {
    for (j in 1:n_groups) {
      full_cov_all[i, j] <- vcov1[group_cols[i], group_cols[j]]
    }
  }
  # Reg2 block
  for (i in 1:n_groups) {
    for (j in 1:n_groups) {
      full_cov_all[n_groups + i, n_groups + j] <- vcov2[group_cols[i], group_cols[j]]
    }
  }
  # Cross block (use named indices to handle intercept/controls correctly)
  group_idx1 <- sapply(group_cols, function(g) which(rownames(reg1$coef) == g))
  group_idx2 <- sapply(group_cols, function(g) which(rownames(reg2$coef) == g))
  for (i in 1:n_groups) {
    for (j in 1:n_groups) {
      full_cov_all[i, n_groups + j] <- cross_cov[group_idx1[i], group_idx2[j]]
      full_cov_all[n_groups + j, i] <- cross_cov[group_idx1[i], group_idx2[j]]
    }
  }
  
  # Contrast: w1_1, ..., w1_n, -w2_1, ..., -w2_n (each side uses its own weights)
  contrast_all <- c(weights1, -weights2)
  all_diff_var <- as.numeric(t(contrast_all) %*% full_cov_all %*% contrast_all)
  all_diff_se <- sqrt(max(all_diff_var, 0))
  
  all_diff <- data.table(
    estimate = all_diff_est,
    se = all_diff_se
  )
  
  # --- Difference in "top-all" test ---
  # top-all = G_n - weighted_avg = G_n - sum(w_g * G_g)
  # Difference: (G_n1 - all1) - (G_n2 - all2)
  # Each side uses its own group proportions as weights
  top_all1 <- res1$top_all$estimate
  top_all2 <- res2$top_all$estimate
  top_all_diff_est <- top_all1 - top_all2
  
  # Contrast for top-all in combined vector, using each side's own weights
  # For reg1: G_n - sum(w1*G_g) -> coef on G_g is (1-w1_n if g=n, else -w1_g)
  contrast_ta1 <- rep(0, n_groups)
  for (g_idx in 1:n_groups) {
    if (g_idx == n_groups) {
      contrast_ta1[g_idx] <- 1 - weights1[g_idx]
    } else {
      contrast_ta1[g_idx] <- -weights1[g_idx]
    }
  }
  # For reg2: -(G_n - sum(w2*G_g)) -> coef on G_g is (-(1-w2_n) if g=n, else w2_g)
  contrast_ta2 <- rep(0, n_groups)
  for (g_idx in 1:n_groups) {
    if (g_idx == n_groups) {
      contrast_ta2[g_idx] <- -(1 - weights2[g_idx])
    } else {
      contrast_ta2[g_idx] <- weights2[g_idx]
    }
  }
  contrast_ta <- c(contrast_ta1, contrast_ta2)
  top_all_diff_var <- as.numeric(t(contrast_ta) %*% full_cov_all %*% contrast_ta)
  top_all_diff_se <- sqrt(max(top_all_diff_var, 0))
  
  top_all_diff <- data.table(
    estimate = top_all_diff_est,
    se = top_all_diff_se
  )
  
  list(
    group_diff = group_diff,
    top_bottom_diff = top_bottom_diff,
    all_diff = all_diff,
    top_all_diff = top_all_diff
  )
}


#' Compare GAVS: Unrestricted vs Restricted Ranking
#' 
#' @description
#' Compares Group Averages (GAVS) between an unrestricted strategy (ranking
#' predictions across the full sample) and a restricted strategy (ranking
#' predictions within strata of a stratification variable).
#' 
#' This function implements the analysis from Fava (2025) to test whether
#' ranking by predicted values within subgroups (e.g., income quintiles,
#' education levels) yields different targeting results than global ranking.
#' 
#' The key statistical challenge is computing correct standard errors for the
#' difference between the two strategies, since they use the same data. This
#' is handled by computing the cross-covariance between regression estimates.
#' 
#' @references
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class \code{ensemble_hte_fit} from 
#'   \code{ensemble_hte()} or \code{ensemble_pred_fit} from \code{ensemble_pred()}.
#' @param strata The stratification variable defining groups for restricted ranking:
#'   \itemize{
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric/factor vector: strata indicator (must have same length as data)
#'   }
#' @param n_groups Number of groups to divide the sample into (default: 3)
#' @param outcome Either:
#'   \itemize{
#'     \item NULL (default): uses the same outcome as in the ensemble function
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric vector: custom outcome variable (must have appropriate length)
#'   }
#' @param subset Which observations to use for the GAVS comparison. Options:
#'   \itemize{
#'     \item NULL (default): uses all observations (or training obs for \code{ensemble_pred_fit}
#'       when using the default outcome with subset training)
#'     \item \code{"train"}: uses only training observations (for \code{ensemble_pred_fit} only)
#'     \item \code{"all"}: explicitly uses all observations
#'     \item Logical vector: TRUE/FALSE for each observation (must have same length as data)
#'     \item Integer vector: indices of observations to include (1-indexed)
#'   }
#'   This allows evaluating GAVS comparison on a subset of observations.
#' 
#' @return An object of class \code{gavs_compare_results} containing:
#' \itemize{
#'   \item unrestricted: \code{gavs_results} object for unrestricted strategy
#'   \item restricted: \code{gavs_results} object for restricted strategy
#'   \item difference: data.table with the difference (unrestricted - restricted)
#'     for each group, with properly computed standard errors
#'   \item top_bottom_diff: data.table with difference in top-bottom estimates
#'   \item strata_var: name of the stratification variable
#'   \item strata_levels: unique levels of the stratification variable
#'   \item n_groups: number of groups used
#'   \item outcome: the outcome variable used
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item n_used: number of observations used
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#' 
#' @examples
#' \dontrun{
#' # Fit ensemble
#' fit <- ensemble_hte(Y ~ X1 + X2, data = mydata, treatment = "D")
#' 
#' # Compare unrestricted vs restricted by income quintile
#' comparison <- gavs_compare(fit, strata = "income_quintile", n_groups = 5)
#' print(comparison)
#' 
#' # The difference shows how much GAVS estimates change when restricting
#' # ranking to within each income quintile
#' }
#' 
#' @export
gavs_compare <- function(ensemble_fit, strata, n_groups = 3, outcome = NULL, 
                         subset = NULL) {
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$ite[[m]])
    fit_type <- "hte"
    has_train_idx <- FALSE
  } else if (inherits(ensemble_fit, "ensemble_pred_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$predictions[[m]])
    fit_type <- "pred"
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
  
  # Determine which outcome to use
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
    stop("strata must be a character string (column name) or a numeric/factor vector")
  }
  
  # Convert strata to factor if not already
  if (!is.factor(strata_vec)) {
    strata_vec <- as.factor(strata_vec)
  }
  strata_levels <- levels(strata_vec)
  
  # Validate that strata is categorical (not too many unique values)
  # If strata has more than 20 unique values, it's likely a continuous variable
  max_strata_levels <- 20
  if (length(strata_levels) > max_strata_levels) {
    stop(paste0("strata variable '", strata_var, "' has ", length(strata_levels), 
                " unique values. This suggests it may be a continuous variable. ",
                "strata should be a categorical variable (factor, integer, or character) ",
                "with a small number of groups (e.g., low/medium/high, or discrete bins). ",
                "Consider using cut() to create discrete groups from a continuous variable."))
  }
  
  # Also warn if strata has very few observations per level
  min_obs_per_level <- min(table(strata_vec))
  if (min_obs_per_level < 10) {
    warning(paste0("Some strata levels have very few observations (min = ", 
                   min_obs_per_level, "). This may lead to unreliable estimates."))
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
  
  # Validate no NAs in strata for used observations
  if (any(is.na(strata_vec[use_idx]))) {
    stop("strata has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (fit_type == "pred" && has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("GAVS Comparison", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Compute GAVS for each repetition (both unrestricted and restricted)
  results_by_rep <- lapply(1:M, function(m) {
    # Unrestricted GAVS (no strata)
    gavs_unrest <- .gavs_single(
      Y = Y[use_idx],
      predicted_values = predictions_list[[m]][use_idx],
      fold = splits[[m]][use_idx],
      n_groups = n_groups,
      strata = NULL
    )
    
    # Restricted GAVS (with strata)
    gavs_rest <- .gavs_single(
      Y = Y[use_idx],
      predicted_values = predictions_list[[m]][use_idx],
      fold = splits[[m]][use_idx],
      n_groups = n_groups,
      strata = strata_vec[use_idx]
    )
    
    # Compute differences with proper SEs using helper
    diff_results <- .compute_reg_diff(gavs_unrest, gavs_rest, n_groups)
    
    list(
      unrestricted = gavs_unrest,
      restricted = gavs_rest,
      diff_results = diff_results
    )
  })
  
  # Aggregate unrestricted results
  all_unrest_group <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$group_estimates), 
                                 idcol = "repetition")
  unrest_combined <- all_unrest_group[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  unrest_combined[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_unrest_tb <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$top_bottom), 
                              idcol = "repetition")
  unrest_tb <- all_unrest_tb[, .(estimate = mean(estimate), se = mean(se))]
  unrest_tb[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_unrest_all <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$all), 
                               idcol = "repetition")
  unrest_all <- all_unrest_all[, .(estimate = mean(estimate), se = mean(se))]
  unrest_all[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_unrest_ta <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$top_all), 
                              idcol = "repetition")
  unrest_ta <- all_unrest_ta[, .(estimate = mean(estimate), se = mean(se))]
  unrest_ta[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  # Aggregate restricted results
  all_rest_group <- rbindlist(lapply(results_by_rep, function(x) x$restricted$group_estimates), 
                               idcol = "repetition")
  rest_combined <- all_rest_group[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  rest_combined[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_rest_tb <- rbindlist(lapply(results_by_rep, function(x) x$restricted$top_bottom), 
                            idcol = "repetition")
  rest_tb <- all_rest_tb[, .(estimate = mean(estimate), se = mean(se))]
  rest_tb[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_rest_all <- rbindlist(lapply(results_by_rep, function(x) x$restricted$all), 
                             idcol = "repetition")
  rest_all <- all_rest_all[, .(estimate = mean(estimate), se = mean(se))]
  rest_all[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_rest_ta <- rbindlist(lapply(results_by_rep, function(x) x$restricted$top_all), 
                            idcol = "repetition")
  rest_ta <- all_rest_ta[, .(estimate = mean(estimate), se = mean(se))]
  rest_ta[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  # Aggregate difference results
  all_diff_group <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$group_diff), 
                               idcol = "repetition")
  diff_combined <- all_diff_group[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  diff_combined[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_diff_tb <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$top_bottom_diff), 
                            idcol = "repetition")
  diff_tb <- all_diff_tb[, .(estimate = mean(estimate), se = mean(se))]
  diff_tb[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_diff_all <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$all_diff), 
                             idcol = "repetition")
  diff_all <- all_diff_all[, .(estimate = mean(estimate), se = mean(se))]
  diff_all[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_diff_ta <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$top_all_diff), 
                            idcol = "repetition")
  diff_ta <- all_diff_ta[, .(estimate = mean(estimate), se = mean(se))]
  diff_ta[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  # Construct unrestricted gavs_results-like object
  unrestricted <- structure(
    list(
      estimates = unrest_combined,
      top_bottom = unrest_tb,
      all = unrest_all,
      top_all = unrest_ta,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      n_used = n_used,
      M = M
    ),
    class = "gavs_results"
  )
  
  # Construct restricted gavs_results-like object
  restricted <- structure(
    list(
      estimates = rest_combined,
      top_bottom = rest_tb,
      all = rest_all,
      top_all = rest_ta,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      n_used = n_used,
      M = M
    ),
    class = "gavs_results"
  )
  
  # Construct gavs_compare_results object
  structure(
    list(
      unrestricted = unrestricted,
      restricted = restricted,
      difference = diff_combined,
      top_bottom_diff = diff_tb,
      all_diff = diff_all,
      top_all_diff = diff_ta,
      strata_var = strata_var,
      strata_levels = strata_levels,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      n_used = n_used,
      M = M,
      call = cl
    ),
    class = "gavs_compare_results"
  )
}


#' Print method for gavs_compare_results objects
#' @param x An object of class \code{gavs_compare_results} from \code{gavs_compare()}
#' @param ... Additional arguments (currently unused)
#' @export
print.gavs_compare_results <- function(x, ...) {
  cat("\n")
  cat("GAVS Comparison: Unrestricted vs Restricted Ranking\n")
  cat(paste0(rep("=", 52), collapse = ""), "\n\n")
  
  # Basic info
  cat("Strategy comparison:\n")
  cat("  - Unrestricted: Rank predictions across full sample within folds\n")
  cat(paste0("  - Restricted: Rank predictions within strata ('", x$strata_var, "')\n"))
  cat(paste0("  - Strata levels: ", paste(x$strata_levels, collapse = ", "), "\n\n"))
  
  pred_type <- if (x$fit_type == "hte") "predicted ITE" else "predicted Y"
  cat(paste0("Groups (", x$n_groups, ") defined by: ", pred_type, "\n"))
  cat(paste0("Outcome: ", x$outcome, "\n"))
  if (x$outcome != x$targeted_outcome) {
    cat(paste0("Note: outcome differs from targeted outcome (", x$targeted_outcome, ")\n"))
  }
  cat(paste0("Observations: ", x$n_used, "\n"))
  cat(paste0("Repetitions: ", x$M, "\n\n"))
  
  # Unrestricted estimates
  cat("Unrestricted GAVS Estimates:\n")
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  print_df <- as.data.frame(x$unrestricted$estimates[, .(group, estimate, se, t_value, p_value)])
  print_df$estimate <- sprintf("%.4f", print_df$estimate)
  print_df$se <- sprintf("%.4f", print_df$se)
  print_df$t_value <- sprintf("%.2f", print_df$t_value)
  print_df$p_value <- sprintf("%.4f", as.numeric(print_df$p_value))
  print_df$sig <- sapply(as.numeric(x$unrestricted$estimates$p_value), get_stars)
  names(print_df) <- c("Group", "Estimate", "SE", "t", "p-value", "")
  print(print_df, row.names = FALSE, right = FALSE)
  
  cat("\n")
  cat(paste0("Top-Bottom: ", sprintf("%.4f", x$unrestricted$top_bottom$estimate),
             " (SE: ", sprintf("%.4f", x$unrestricted$top_bottom$se), ", ",
             "p = ", sprintf("%.4f", x$unrestricted$top_bottom$p_value), ")\n"))
  if (!is.null(x$unrestricted$all)) {
    cat(paste0("All: ", sprintf("%.4f", x$unrestricted$all$estimate),
               " (SE: ", sprintf("%.4f", x$unrestricted$all$se), ", ",
               "p = ", sprintf("%.4f", x$unrestricted$all$p_value), ")\n"))
  }
  if (!is.null(x$unrestricted$top_all)) {
    cat(paste0("Top-All: ", sprintf("%.4f", x$unrestricted$top_all$estimate),
               " (SE: ", sprintf("%.4f", x$unrestricted$top_all$se), ", ",
               "p = ", sprintf("%.4f", x$unrestricted$top_all$p_value), ")\n"))
  }
  cat("\n")
  
  # Restricted estimates
  cat("Restricted GAVS Estimates:\n")
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  print_df <- as.data.frame(x$restricted$estimates[, .(group, estimate, se, t_value, p_value)])
  print_df$estimate <- sprintf("%.4f", print_df$estimate)
  print_df$se <- sprintf("%.4f", print_df$se)
  print_df$t_value <- sprintf("%.2f", print_df$t_value)
  print_df$p_value <- sprintf("%.4f", as.numeric(print_df$p_value))
  print_df$sig <- sapply(as.numeric(x$restricted$estimates$p_value), get_stars)
  names(print_df) <- c("Group", "Estimate", "SE", "t", "p-value", "")
  print(print_df, row.names = FALSE, right = FALSE)
  
  cat("\n")
  cat(paste0("Top-Bottom: ", sprintf("%.4f", x$restricted$top_bottom$estimate),
             " (SE: ", sprintf("%.4f", x$restricted$top_bottom$se), ", ",
             "p = ", sprintf("%.4f", x$restricted$top_bottom$p_value), ")\n"))
  if (!is.null(x$restricted$all)) {
    cat(paste0("All: ", sprintf("%.4f", x$restricted$all$estimate),
               " (SE: ", sprintf("%.4f", x$restricted$all$se), ", ",
               "p = ", sprintf("%.4f", x$restricted$all$p_value), ")\n"))
  }
  if (!is.null(x$restricted$top_all)) {
    cat(paste0("Top-All: ", sprintf("%.4f", x$restricted$top_all$estimate),
               " (SE: ", sprintf("%.4f", x$restricted$top_all$se), ", ",
               "p = ", sprintf("%.4f", x$restricted$top_all$p_value), ")\n"))
  }
  cat("\n")
  
  # Difference (Unrestricted - Restricted)
  cat("Difference (Unrestricted - Restricted):\n")
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  print_df <- as.data.frame(x$difference[, .(group, estimate, se, t_value, p_value)])
  print_df$estimate <- sprintf("%.4f", print_df$estimate)
  print_df$se <- sprintf("%.4f", print_df$se)
  print_df$t_value <- sprintf("%.2f", print_df$t_value)
  print_df$p_value <- sprintf("%.4f", as.numeric(print_df$p_value))
  print_df$sig <- sapply(as.numeric(x$difference$p_value), get_stars)
  names(print_df) <- c("Group", "Estimate", "SE", "t", "p-value", "")
  print(print_df, row.names = FALSE, right = FALSE)
  
  cat("\n")
  cat(paste0("Top-Bottom Diff: ", sprintf("%.4f", x$top_bottom_diff$estimate),
             " (SE: ", sprintf("%.4f", x$top_bottom_diff$se), ", ",
             "p = ", sprintf("%.4f", x$top_bottom_diff$p_value), ")\n"))
  if (!is.null(x$all_diff)) {
    cat(paste0("All Diff: ", sprintf("%.4f", x$all_diff$estimate),
               " (SE: ", sprintf("%.4f", x$all_diff$se), ", ",
               "p = ", sprintf("%.4f", x$all_diff$p_value), ")\n"))
  }
  if (!is.null(x$top_all_diff)) {
    cat(paste0("Top-All Diff: ", sprintf("%.4f", x$top_all_diff$estimate),
               " (SE: ", sprintf("%.4f", x$top_all_diff$se), ", ",
               "p = ", sprintf("%.4f", x$top_all_diff$p_value), ")\n"))
  }
  
  cat("\n")
  cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")
  
  invisible(x)
}


#' Plot method for gavs_compare_results objects
#'
#' @description
#' Creates a comparison plot showing GAVS estimates with confidence intervals
#' for both unrestricted and restricted strategies side by side.
#'
#' @param x An object of class \code{gavs_compare_results} from \code{gavs_compare()}
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   labs theme_minimal theme element_text scale_x_continuous position_dodge
#' @export
plot.gavs_compare_results <- function(x, alpha = 0.05, ...) {
  
  z_val <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  # Prepare unrestricted data
  unrest_data <- copy(x$unrestricted$estimates)
  unrest_data[, strategy := "Unrestricted"]
  unrest_data[, ci_lower := estimate - z_val * se]
  unrest_data[, ci_upper := estimate + z_val * se]
  
  # Prepare restricted data
  rest_data <- copy(x$restricted$estimates)
  rest_data[, strategy := "Restricted"]
  rest_data[, ci_lower := estimate - z_val * se]
  rest_data[, ci_upper := estimate + z_val * se]
  
  # Combine
  plot_data <- rbindlist(list(unrest_data, rest_data))
  plot_data[, strategy := factor(strategy, levels = c("Unrestricted", "Restricted"))]
  
  # Determine label for y-axis
  pred_type <- if (x$fit_type == "hte") "Predicted ITE" else "Predicted Y"
  
  # Calculate global mean for reference line
  global_mean <- mean(x$unrestricted$estimates$estimate)
  
  # Create comparison plot with dodge
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = estimate, 
                                                color = strategy, shape = strategy)) +
    ggplot2::geom_hline(yintercept = global_mean, linetype = "dashed", 
                        color = "blue", linewidth = 0.5, alpha = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3), size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15, linewidth = 0.5,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::annotate("text", x = x$n_groups + 0.6, y = global_mean, 
                      label = paste0("Mean = ", sprintf("%.3f", global_mean)),
                      hjust = 0, vjust = -0.5, size = 3, color = "blue") +
    ggplot2::scale_x_continuous(
      breaks = 1:x$n_groups,
      labels = paste0("G", 1:x$n_groups),
      limits = c(0.5, x$n_groups + 1.2)
    ) +
    ggplot2::scale_color_manual(values = c("Unrestricted" = "#2C3E50", "Restricted" = "#E74C3C")) +
    ggplot2::scale_shape_manual(values = c("Unrestricted" = 16, "Restricted" = 17)) +
    ggplot2::labs(
      title = "GAVS Comparison: Unrestricted vs Restricted",
      subtitle = paste0("Groups ranked by ", pred_type, 
                        " | Restricted by: ", x$strata_var),
      x = paste0("Group (", x$n_groups, " groups by ", pred_type, ")"),
      y = paste0("Average ", x$outcome),
      color = "Strategy",
      shape = "Strategy",
      caption = paste0(conf_level, "% confidence intervals. ", x$M, " repetitions. ",
                       x$n_used, " observations.")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  invisible(p)
}


#' Compare GATES: Unrestricted vs Restricted Ranking
#' 
#' @description
#' Compares Group Average Treatment Effects (GATES) between an unrestricted strategy 
#' (ranking predictions across the full sample) and a restricted strategy (ranking
#' predictions within strata of a stratification variable).
#' 
#' This function implements the analysis from Fava (2025) to test whether
#' ranking by predicted treatment effects within subgroups (e.g., income quintiles,
#' education levels) yields different treatment effect estimates than global ranking.
#' 
#' The key statistical challenge is computing correct standard errors for the
#' difference between the two strategies, since they use the same data. This
#' is handled by computing the cross-covariance between regression estimates.
#' 
#' @references
#' Fava, B. (2025). Training and Testing with Multiple Splits: A Central Limit
#' Theorem for Split-Sample Estimators. \emph{arXiv preprint arXiv:2511.04957}.
#' 
#' @param ensemble_fit An object of class \code{ensemble_hte_fit} from 
#'   \code{ensemble_hte()} or \code{ensemble_pred_fit} from \code{ensemble_pred()}.
#' @param strata The stratification variable defining groups for restricted ranking:
#'   \itemize{
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric/factor vector: strata indicator (must have same length as data)
#'   }
#' @param n_groups Number of groups to divide the sample into (default: 3)
#' @param outcome Either:
#'   \itemize{
#'     \item NULL (default): uses the same outcome as in the ensemble function
#'     \item Character string: column name in the \code{data} used in the ensemble function
#'     \item Numeric vector: custom outcome variable (must have appropriate length)
#'   }
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
#' @param controls Optional character vector of control variable names from \code{data}
#'   to include as covariates in the GATES regression.
#' @param subset Which observations to use for the GATES comparison. Options:
#'   \itemize{
#'     \item NULL (default): uses all observations (or training obs for \code{ensemble_pred_fit}
#'       when using the default outcome with subset training)
#'     \item \code{"train"}: uses only training observations (for \code{ensemble_pred_fit} only)
#'     \item \code{"all"}: explicitly uses all observations
#'     \item Logical vector: TRUE/FALSE for each observation (must have same length as data)
#'     \item Integer vector: indices of observations to include (1-indexed)
#'   }
#'   This allows evaluating GATES comparison on a subset of observations.
#' 
#' @return An object of class \code{gates_compare_results} containing:
#' \itemize{
#'   \item unrestricted: \code{gates_results} object for unrestricted strategy
#'   \item restricted: \code{gates_results} object for restricted strategy
#'   \item difference: data.table with the difference (unrestricted - restricted)
#'     for each group, with properly computed standard errors
#'   \item top_bottom_diff: data.table with difference in top-bottom estimates
#'   \item all_diff: data.table with difference in weighted average (all) estimates
#'   \item top_all_diff: data.table with difference in top-all estimates
#'   \item strata_var: name of the stratification variable
#'   \item strata_levels: unique levels of the stratification variable
#'   \item n_groups: number of groups used
#'   \item outcome: the outcome variable used
#'   \item targeted_outcome: the outcome used for prediction
#'   \item fit_type: "hte" or "pred" depending on input
#'   \item n_used: number of observations used
#'   \item M: number of repetitions
#'   \item call: the function call
#' }
#' 
#' @examples
#' \dontrun{
#' # Fit ensemble
#' fit <- ensemble_hte(Y ~ X1 + X2, data = mydata, treatment = "D")
#' 
#' # Compare unrestricted vs restricted by income quintile
#' comparison <- gates_compare(fit, strata = "income_quintile", n_groups = 5)
#' print(comparison)
#' }
#' 
#' @export
gates_compare <- function(ensemble_fit, strata, n_groups = 3, outcome = NULL, 
                          treatment = NULL, prop_score = NULL, controls = NULL, subset = NULL) {
  
  # Check input type and extract predictions accordingly
  if (inherits(ensemble_fit, "ensemble_hte_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$ite[[m]])
    fit_type <- "hte"
    has_train_idx <- FALSE
    # For HTE fits, treatment info comes from the fit
    D <- ensemble_fit$D
    ps <- ensemble_fit$prop_score
    W <- ensemble_fit$weights
  } else if (inherits(ensemble_fit, "ensemble_pred_fit")) {
    predictions_list <- lapply(1:ensemble_fit$M, function(m) ensemble_fit$predictions[[m]])
    fit_type <- "pred"
    has_train_idx <- !is.null(ensemble_fit$train_idx) && 
                     !is.null(ensemble_fit$n_train) && 
                     ensemble_fit$n_train < ensemble_fit$n
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
  
  # Determine which outcome to use
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
    W <- 1 / (ps * (1 - ps))
  }
  
  # Process strata argument
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
    stop("strata must be a character string (column name) or a numeric/factor vector")
  }
  
  # Convert strata to factor if not already
  if (!is.factor(strata_vec)) {
    strata_vec <- as.factor(strata_vec)
  }
  strata_levels <- levels(strata_vec)
  
  # Validate that strata is categorical (not too many unique values)
  # If strata has more than 20 unique values, it's likely a continuous variable
  max_strata_levels <- 20
  if (length(strata_levels) > max_strata_levels) {
    stop(paste0("strata variable '", strata_var, "' has ", length(strata_levels), 
                " unique values. This suggests it may be a continuous variable. ",
                "strata should be a categorical variable (factor, integer, or character) ",
                "with a small number of groups (e.g., low/medium/high, or discrete bins). ",
                "Consider using cut() to create discrete groups from a continuous variable."))
  }
  
  # Also warn if strata has very few observations per level
  min_obs_per_level <- min(table(strata_vec))
  if (min_obs_per_level < 10) {
    warning(paste0("Some strata levels have very few observations (min = ", 
                   min_obs_per_level, "). This may lead to unreliable estimates."))
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
  
  # Validate no NAs in treatment for used observations
  if (any(is.na(D[use_idx]))) {
    stop("treatment has NA values for the observations being used")
  }
  
  # Validate no NAs in strata for used observations
  if (any(is.na(strata_vec[use_idx]))) {
    stop("strata has NA values for the observations being used")
  }
  
  n_used <- sum(use_idx)
  n_fit <- if (fit_type == "pred" && has_train_idx) ensemble_fit$n_train else n
  
  # Print subset message
  .print_subset_message("GATES Comparison", n, n_fit, n_used, fit_type, has_train_idx)
  
  # Get controls data if specified (subset to used observations)
  if (!is.null(controls)) {
    if (!is.character(controls)) {
      stop("controls must be a character vector of column names")
    }
    missing_controls <- controls[!controls %in% names(ensemble_fit$data)]
    if (length(missing_controls) > 0) {
      stop(paste0("Control variable(s) not found in data: ", 
                  paste(missing_controls, collapse = ", ")))
    }
    control_data <- ensemble_fit$data[use_idx, ..controls]
  } else {
    control_data <- NULL
  }
  
  # Subset strata
  strata_used <- strata_vec[use_idx]
  
  # Extract components from ensemble_fit
  splits <- ensemble_fit$splits
  M <- ensemble_fit$M
  
  # Compute GATES for each repetition (both unrestricted and restricted)
  results_by_rep <- lapply(1:M, function(m) {
    # Unrestricted GATES (no strata)
    gates_unrest <- .gates_single(
      Y = Y[use_idx],
      D = D[use_idx],
      prop_score = ps[use_idx],
      weight = W[use_idx],
      predicted_values = predictions_list[[m]][use_idx],
      fold = splits[[m]][use_idx],
      n_groups = n_groups,
      controls = controls,
      control_data = control_data,
      strata = NULL
    )
    
    # Restricted GATES (with strata)
    gates_rest <- .gates_single(
      Y = Y[use_idx],
      D = D[use_idx],
      prop_score = ps[use_idx],
      weight = W[use_idx],
      predicted_values = predictions_list[[m]][use_idx],
      fold = splits[[m]][use_idx],
      n_groups = n_groups,
      controls = controls,
      control_data = control_data,
      strata = strata_used
    )
    
    # Compute differences with proper SEs using helper
    diff_results <- .compute_reg_diff(gates_unrest, gates_rest, n_groups)
    
    list(
      unrestricted = gates_unrest,
      restricted = gates_rest,
      diff_results = diff_results
    )
  })
  
  # Aggregate unrestricted results
  all_unrest_group <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$group_estimates), 
                                 idcol = "repetition")
  unrest_combined <- all_unrest_group[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  unrest_combined[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_unrest_tb <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$top_bottom), 
                              idcol = "repetition")
  unrest_tb <- all_unrest_tb[, .(estimate = mean(estimate), se = mean(se))]
  unrest_tb[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_unrest_all <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$all), 
                               idcol = "repetition")
  unrest_all <- all_unrest_all[, .(estimate = mean(estimate), se = mean(se))]
  unrest_all[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_unrest_ta <- rbindlist(lapply(results_by_rep, function(x) x$unrestricted$top_all), 
                              idcol = "repetition")
  unrest_ta <- all_unrest_ta[, .(estimate = mean(estimate), se = mean(se))]
  unrest_ta[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  # Aggregate restricted results
  all_rest_group <- rbindlist(lapply(results_by_rep, function(x) x$restricted$group_estimates), 
                               idcol = "repetition")
  rest_combined <- all_rest_group[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  rest_combined[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_rest_tb <- rbindlist(lapply(results_by_rep, function(x) x$restricted$top_bottom), 
                            idcol = "repetition")
  rest_tb <- all_rest_tb[, .(estimate = mean(estimate), se = mean(se))]
  rest_tb[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_rest_all <- rbindlist(lapply(results_by_rep, function(x) x$restricted$all), 
                             idcol = "repetition")
  rest_all <- all_rest_all[, .(estimate = mean(estimate), se = mean(se))]
  rest_all[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_rest_ta <- rbindlist(lapply(results_by_rep, function(x) x$restricted$top_all), 
                            idcol = "repetition")
  rest_ta <- all_rest_ta[, .(estimate = mean(estimate), se = mean(se))]
  rest_ta[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  # Aggregate difference results
  all_diff_group <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$group_diff), 
                               idcol = "repetition")
  diff_combined <- all_diff_group[, .(
    estimate = mean(estimate),
    se = mean(se),
    n_reps = .N
  ), by = group]
  diff_combined[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_diff_tb <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$top_bottom_diff), 
                            idcol = "repetition")
  diff_tb <- all_diff_tb[, .(estimate = mean(estimate), se = mean(se))]
  diff_tb[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_diff_all <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$all_diff), 
                             idcol = "repetition")
  diff_all <- all_diff_all[, .(estimate = mean(estimate), se = mean(se))]
  diff_all[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  all_diff_ta <- rbindlist(lapply(results_by_rep, function(x) x$diff_results$top_all_diff), 
                            idcol = "repetition")
  diff_ta <- all_diff_ta[, .(estimate = mean(estimate), se = mean(se))]
  diff_ta[, `:=`(t_value = estimate / se, p_value = 2 * pnorm(-abs(estimate / se)))]
  
  # Construct unrestricted gates_results-like object
  unrestricted <- structure(
    list(
      estimates = unrest_combined,
      top_bottom = unrest_tb,
      all = unrest_all,
      top_all = unrest_ta,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      controls = controls,
      n_used = n_used,
      M = M
    ),
    class = "gates_results"
  )
  
  # Construct restricted gates_results-like object
  restricted <- structure(
    list(
      estimates = rest_combined,
      top_bottom = rest_tb,
      all = rest_all,
      top_all = rest_ta,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      controls = controls,
      n_used = n_used,
      M = M
    ),
    class = "gates_results"
  )
  
  # Construct gates_compare_results object
  structure(
    list(
      unrestricted = unrestricted,
      restricted = restricted,
      difference = diff_combined,
      top_bottom_diff = diff_tb,
      all_diff = diff_all,
      top_all_diff = diff_ta,
      strata_var = strata_var,
      strata_levels = strata_levels,
      n_groups = n_groups,
      outcome = outcome_var,
      targeted_outcome = targeted_outcome,
      fit_type = fit_type,
      controls = controls,
      n_used = n_used,
      M = M,
      call = cl
    ),
    class = "gates_compare_results"
  )
}


#' Print method for gates_compare_results objects
#' @param x An object of class \code{gates_compare_results} from \code{gates_compare()}
#' @param ... Additional arguments (currently unused)
#' @export
print.gates_compare_results <- function(x, ...) {
  cat("\n")
  cat("GATES Comparison: Unrestricted vs Restricted Ranking\n")
  cat(paste0(rep("=", 53), collapse = ""), "\n\n")
  
  # Show fit type
  fit_label <- if (!is.null(x$fit_type) && x$fit_type == "hte") {
    "HTE (ensemble_hte)"
  } else if (!is.null(x$fit_type) && x$fit_type == "pred") {
    "Prediction (ensemble_pred)"
  } else {
    "HTE (ensemble_hte)"  # Default for backwards compatibility
  }
  cat("Fit type:", fit_label, "\n\n")
  
  # Basic info
  cat("Strategy comparison:\n")
  cat("  - Unrestricted: Rank predictions across full sample within folds\n")
  cat(paste0("  - Restricted: Rank predictions within strata ('", x$strata_var, "')\n"))
  cat(paste0("  - Strata levels: ", paste(x$strata_levels, collapse = ", "), "\n\n"))
  
  pred_label <- if (!is.null(x$fit_type) && x$fit_type == "pred") "predicted Y" else "predicted ITE"
  cat(paste0("Groups (", x$n_groups, ") defined by: ", pred_label, "\n"))
  cat(paste0("Outcome: ", x$outcome, "\n"))
  if (x$outcome != x$targeted_outcome) {
    cat(paste0("Note: outcome differs from targeted outcome (", x$targeted_outcome, ")\n"))
  }
  if (!is.null(x$controls)) {
    cat("Controls:", paste(x$controls, collapse = ", "), "\n")
  }
  cat(paste0("Observations: ", x$n_used, "\n"))
  cat(paste0("Repetitions: ", x$M, "\n\n"))
  
  # Unrestricted estimates
  cat("Unrestricted GATES Estimates:\n")
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  print_df <- as.data.frame(x$unrestricted$estimates[, .(group, estimate, se, t_value, p_value)])
  print_df$estimate <- sprintf("%.4f", print_df$estimate)
  print_df$se <- sprintf("%.4f", print_df$se)
  print_df$t_value <- sprintf("%.2f", print_df$t_value)
  print_df$p_value <- sprintf("%.4f", as.numeric(print_df$p_value))
  print_df$sig <- sapply(as.numeric(x$unrestricted$estimates$p_value), get_stars)
  names(print_df) <- c("Group", "Estimate", "SE", "t", "p-value", "")
  print(print_df, row.names = FALSE, right = FALSE)
  
  cat("\n")
  cat(paste0("Top-Bottom: ", sprintf("%.4f", x$unrestricted$top_bottom$estimate),
             " (SE: ", sprintf("%.4f", x$unrestricted$top_bottom$se), ", ",
             "p = ", sprintf("%.4f", x$unrestricted$top_bottom$p_value), ")\n"))
  if (!is.null(x$unrestricted$all)) {
    cat(paste0("All: ", sprintf("%.4f", x$unrestricted$all$estimate),
               " (SE: ", sprintf("%.4f", x$unrestricted$all$se), ", ",
               "p = ", sprintf("%.4f", x$unrestricted$all$p_value), ")\n"))
  }
  if (!is.null(x$unrestricted$top_all)) {
    cat(paste0("Top-All: ", sprintf("%.4f", x$unrestricted$top_all$estimate),
               " (SE: ", sprintf("%.4f", x$unrestricted$top_all$se), ", ",
               "p = ", sprintf("%.4f", x$unrestricted$top_all$p_value), ")\n"))
  }
  cat("\n")
  
  # Restricted estimates
  cat("Restricted GATES Estimates:\n")
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  print_df <- as.data.frame(x$restricted$estimates[, .(group, estimate, se, t_value, p_value)])
  print_df$estimate <- sprintf("%.4f", print_df$estimate)
  print_df$se <- sprintf("%.4f", print_df$se)
  print_df$t_value <- sprintf("%.2f", print_df$t_value)
  print_df$p_value <- sprintf("%.4f", as.numeric(print_df$p_value))
  print_df$sig <- sapply(as.numeric(x$restricted$estimates$p_value), get_stars)
  names(print_df) <- c("Group", "Estimate", "SE", "t", "p-value", "")
  print(print_df, row.names = FALSE, right = FALSE)
  
  cat("\n")
  cat(paste0("Top-Bottom: ", sprintf("%.4f", x$restricted$top_bottom$estimate),
             " (SE: ", sprintf("%.4f", x$restricted$top_bottom$se), ", ",
             "p = ", sprintf("%.4f", x$restricted$top_bottom$p_value), ")\n"))
  if (!is.null(x$restricted$all)) {
    cat(paste0("All: ", sprintf("%.4f", x$restricted$all$estimate),
               " (SE: ", sprintf("%.4f", x$restricted$all$se), ", ",
               "p = ", sprintf("%.4f", x$restricted$all$p_value), ")\n"))
  }
  if (!is.null(x$restricted$top_all)) {
    cat(paste0("Top-All: ", sprintf("%.4f", x$restricted$top_all$estimate),
               " (SE: ", sprintf("%.4f", x$restricted$top_all$se), ", ",
               "p = ", sprintf("%.4f", x$restricted$top_all$p_value), ")\n"))
  }
  cat("\n")
  
  # Difference (Unrestricted - Restricted)
  cat("Difference (Unrestricted - Restricted):\n")
  cat(paste0(rep("-", 40), collapse = ""), "\n")
  print_df <- as.data.frame(x$difference[, .(group, estimate, se, t_value, p_value)])
  print_df$estimate <- sprintf("%.4f", print_df$estimate)
  print_df$se <- sprintf("%.4f", print_df$se)
  print_df$t_value <- sprintf("%.2f", print_df$t_value)
  print_df$p_value <- sprintf("%.4f", as.numeric(print_df$p_value))
  print_df$sig <- sapply(as.numeric(x$difference$p_value), get_stars)
  names(print_df) <- c("Group", "Estimate", "SE", "t", "p-value", "")
  print(print_df, row.names = FALSE, right = FALSE)
  
  cat("\n")
  cat(paste0("Top-Bottom Diff: ", sprintf("%.4f", x$top_bottom_diff$estimate),
             " (SE: ", sprintf("%.4f", x$top_bottom_diff$se), ", ",
             "p = ", sprintf("%.4f", x$top_bottom_diff$p_value), ")\n"))
  if (!is.null(x$all_diff)) {
    cat(paste0("All Diff: ", sprintf("%.4f", x$all_diff$estimate),
               " (SE: ", sprintf("%.4f", x$all_diff$se), ", ",
               "p = ", sprintf("%.4f", x$all_diff$p_value), ")\n"))
  }
  if (!is.null(x$top_all_diff)) {
    cat(paste0("Top-All Diff: ", sprintf("%.4f", x$top_all_diff$estimate),
               " (SE: ", sprintf("%.4f", x$top_all_diff$se), ", ",
               "p = ", sprintf("%.4f", x$top_all_diff$p_value), ")\n"))
  }
  
  cat("\n")
  cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")
  
  invisible(x)
}


#' Plot method for gates_compare_results objects
#'
#' @description
#' Creates a comparison plot showing GATES estimates with confidence intervals
#' for both unrestricted and restricted strategies side by side.
#'
#' @param x An object of class \code{gates_compare_results} from \code{gates_compare()}
#' @param alpha Significance level for confidence intervals (default 0.05)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object (invisibly)
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar geom_hline
#'   labs theme_minimal theme element_text scale_x_continuous position_dodge
#' @export
plot.gates_compare_results <- function(x, alpha = 0.05, ...) {
  
  z_val <- qnorm(1 - alpha / 2)
  conf_level <- (1 - alpha) * 100
  
  # Prepare unrestricted data
  unrest_data <- copy(x$unrestricted$estimates)
  unrest_data[, strategy := "Unrestricted"]
  unrest_data[, ci_lower := estimate - z_val * se]
  unrest_data[, ci_upper := estimate + z_val * se]
  
  # Prepare restricted data
  rest_data <- copy(x$restricted$estimates)
  rest_data[, strategy := "Restricted"]
  rest_data[, ci_lower := estimate - z_val * se]
  rest_data[, ci_upper := estimate + z_val * se]
  
  # Combine
  plot_data <- rbindlist(list(unrest_data, rest_data))
  plot_data[, strategy := factor(strategy, levels = c("Unrestricted", "Restricted"))]
  
  # Calculate global ATE for reference line
  global_ate <- mean(x$unrestricted$estimates$estimate)
  
  # Create comparison plot with dodge
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = estimate, 
                                                color = strategy, shape = strategy)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = global_ate, linetype = "solid", 
                        color = "blue", linewidth = 0.5, alpha = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.3), size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15, linewidth = 0.5,
      position = ggplot2::position_dodge(width = 0.3)
    ) +
    ggplot2::annotate("text", x = x$n_groups + 0.6, y = global_ate, 
                      label = paste0("ATE = ", sprintf("%.3f", global_ate)),
                      hjust = 0, vjust = -0.5, size = 3, color = "blue") +
    ggplot2::scale_x_continuous(
      breaks = 1:x$n_groups,
      labels = paste0("G", 1:x$n_groups),
      limits = c(0.5, x$n_groups + 1.2)
    ) +
    ggplot2::scale_color_manual(values = c("Unrestricted" = "#2C3E50", "Restricted" = "#E74C3C")) +
    ggplot2::scale_shape_manual(values = c("Unrestricted" = 16, "Restricted" = 17)) +
    ggplot2::labs(
      title = "GATES Comparison: Unrestricted vs Restricted",
      subtitle = paste0("Groups ranked by Predicted ITE", 
                        " | Restricted by: ", x$strata_var),
      x = paste0("Group (", x$n_groups, " groups by Predicted ITE)"),
      y = "Treatment Effect",
      color = "Strategy",
      shape = "Strategy",
      caption = paste0(conf_level, "% confidence intervals. ", x$M, " repetitions. ",
                       x$n_used, " observations.")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  print(p)
  invisible(p)
}

