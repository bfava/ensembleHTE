# ==============================================================================
# Tests for ensemble_hte, ensemble_pred, and analysis functions (gates, gavs)
# ==============================================================================

# Create simulated data for tests
# Using n=200 by default to reduce zero-variance warnings with small samples
create_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)
  D <- rbinom(n, 1, 0.5)  # treatment
  # Heterogeneous treatment effect: 1 + 0.5 * X1
  tau <- 1 + 0.5 * X1
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  
  data.frame(Y = Y, D = D, X1 = X1, X2 = X2, X3 = X3)
}

# ==============================================================================
# Tests for ensemble_hte
# ==============================================================================

test_that("ensemble_hte creates proper object structure", {
  skip_on_cran()
  
  # Use n=200 and T-learner to avoid zero-variance issues
  data <- create_test_data(n = 200)
  
  # Fit basic model with T-learner which is more stable with small samples
  fit <- ensemble_hte(
    Y ~ X1 + X2 + X3,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Test class
  expect_s3_class(fit, "ensemble_hte_fit")
  
  # Test required components exist
  expect_true(!is.null(fit$formula))
  expect_true(!is.null(fit$data))
  expect_true(!is.null(fit$Y))
  expect_true(!is.null(fit$D))
  expect_true(!is.null(fit$X))
  expect_true(!is.null(fit$ite))
  expect_true(!is.null(fit$M))
  expect_true(!is.null(fit$K))
  expect_true(!is.null(fit$splits))
  expect_true(!is.null(fit$prop_score))
  expect_true(!is.null(fit$weights))
  
  # Test dimensions
  expect_equal(length(fit$Y), nrow(data))
  expect_equal(length(fit$D), nrow(data))
  expect_equal(nrow(fit$X), nrow(data))
  expect_equal(fit$n, nrow(data))
  expect_equal(fit$M, 2)
  expect_equal(fit$K, 2)
  
  # Test ITE predictions structure
  expect_equal(length(fit$ite), fit$M)
  for (m in 1:fit$M) {
    expect_equal(length(fit$ite[[m]]), nrow(data))
  }
  
  # Test splits structure
  expect_equal(length(fit$splits), fit$M)
  for (m in 1:fit$M) {
    expect_equal(length(fit$splits[[m]]), nrow(data))
    expect_true(all(fit$splits[[m]] %in% 1:fit$K))
  }
})


test_that("ensemble_hte print and summary work", {
  skip_on_cran()
  
  data <- create_test_data(n = 200)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Print should work without error
  expect_output(print(fit), "Ensemble HTE")
})


# ==============================================================================
# Tests for ensemble_pred
# ==============================================================================

test_that("ensemble_pred creates proper object structure", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  # Fit basic model
  fit <- ensemble_pred(
    Y ~ X1 + X2 + X3,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Test class
  expect_s3_class(fit, "ensemble_pred_fit")
  
  # Test required components exist
  expect_true(!is.null(fit$formula))
  expect_true(!is.null(fit$data))
  expect_true(!is.null(fit$Y))
  expect_true(!is.null(fit$X))
  expect_true(!is.null(fit$predictions))
  expect_true(!is.null(fit$M))
  expect_true(!is.null(fit$K))
  expect_true(!is.null(fit$splits))
  
  # Test dimensions
  expect_equal(length(fit$Y), nrow(data))
  expect_equal(nrow(fit$X), nrow(data))
  expect_equal(fit$n, nrow(data))
  
  # Test predictions structure
  expect_equal(length(fit$predictions), fit$M)
  for (m in 1:fit$M) {
    expect_equal(length(fit$predictions[[m]]), nrow(data))
  }
})


test_that("ensemble_pred works with train_idx", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  # Only use first 60 observations for training
  train_idx <- 1:60
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    train_idx = train_idx,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Test class
  expect_s3_class(fit, "ensemble_pred_fit")
  
  # Test train_idx is stored
  expect_true(!is.null(fit$train_idx))
  expect_equal(sum(fit$train_idx), 60)
  
  # Predictions should still be for all observations
  expect_equal(length(fit$predictions[[1]]), nrow(data))
})


test_that("ensemble_pred print and summary work", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Print should work without error
  expect_output(print(fit), "Ensemble Prediction")
  
  # Summary should work and return summary object
  summ <- summary(fit)
  expect_s3_class(summ, "summary.ensemble_pred_fit")
  expect_true(!is.null(summ$blp))
  expect_true(!is.null(summ$gavs))
})


# ==============================================================================
# Tests for gates() with both fit types
# ==============================================================================

test_that("gates works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 200)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Run GATES
  gates_result <- gates(fit, n_groups = 3)
  
  # Test class
  expect_s3_class(gates_result, "gates_results")
  
  # Test required components
  expect_true(!is.null(gates_result$estimates))
  expect_true(!is.null(gates_result$top_bottom))
  expect_equal(gates_result$n_groups, 3)
  expect_equal(gates_result$fit_type, "hte")
  expect_equal(gates_result$M, 2)
  
  # Test estimates structure
  expect_equal(nrow(gates_result$estimates), 3)
  expect_true("estimate" %in% names(gates_result$estimates))
  expect_true("se" %in% names(gates_result$estimates))
  expect_true("p_value" %in% names(gates_result$estimates))
  
  # Print should work
  expect_output(print(gates_result), "GATES")
  expect_output(print(gates_result), "HTE")
})


test_that("gates works with ensemble_pred_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Run GATES - must specify treatment
  gates_result <- gates(fit, treatment = "D", n_groups = 3)
  
  # Test class
  expect_s3_class(gates_result, "gates_results")
  
  # Test fit_type
  expect_equal(gates_result$fit_type, "pred")
  
  # Test required components
  expect_equal(nrow(gates_result$estimates), 3)
  expect_true(!is.null(gates_result$top_bottom))
  
  # Print should work
  expect_output(print(gates_result), "GATES")
  expect_output(print(gates_result), "Prediction")
})


test_that("gates with ensemble_pred_fit requires treatment", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Should error without treatment
  expect_error(
    gates(fit, n_groups = 3),
    "treatment must be specified"
  )
})


test_that("gates accepts treatment as vector", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Use treatment as vector
  gates_result <- gates(fit, treatment = data$D, n_groups = 3)
  
  expect_s3_class(gates_result, "gates_results")
  expect_equal(gates_result$fit_type, "pred")
})


test_that("gates accepts custom propensity scores", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Use custom propensity score
  gates_result <- gates(fit, treatment = "D", prop_score = 0.5, n_groups = 3)
  
  expect_s3_class(gates_result, "gates_results")
})


# ==============================================================================
# Tests for gavs() with both fit types
# ==============================================================================

test_that("gavs works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 200)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Run GAVS
  gavs_result <- gavs(fit, n_groups = 3)
  
  # Test class
  expect_s3_class(gavs_result, "gavs_results")
  
  # Test required components
  expect_true(!is.null(gavs_result$estimates))
  expect_true(!is.null(gavs_result$top_bottom))
  expect_equal(gavs_result$n_groups, 3)
  expect_equal(gavs_result$fit_type, "hte")
  
  # Print should work
  expect_output(print(gavs_result), "GAVS")
})


test_that("gavs works with ensemble_pred_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Run GAVS
  gavs_result <- gavs(fit, n_groups = 3)
  
  # Test class
  expect_s3_class(gavs_result, "gavs_results")
  
  # Test fit_type
  expect_equal(gavs_result$fit_type, "pred")
  
  # Print should work
  expect_output(print(gavs_result), "GAVS")
})


test_that("gavs works with train_idx subset", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  # Fit with train_idx
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    train_idx = 1:60,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # GAVS with subset=TRUE (default for train_idx fits)
  gavs_result <- gavs(fit, n_groups = 3)
  
  expect_s3_class(gavs_result, "gavs_results")
  # Should use only training observations
  expect_equal(gavs_result$n_used, 60)
  
  # GAVS with subset='all' (use all)
  gavs_all <- gavs(fit, n_groups = 3, subset = "all")
  expect_equal(gavs_all$n_used, 100)
})


# ==============================================================================
# Tests for plot methods
# ==============================================================================

test_that("plot methods work for gates and gavs results", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  data <- create_test_data(n = 200)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # GATES plot
  gates_result <- gates(fit, n_groups = 3)
  p_gates <- plot(gates_result)
  expect_s3_class(p_gates, "ggplot")
  
  # GAVS plot
  gavs_result <- gavs(fit, n_groups = 3)
  p_gavs <- plot(gavs_result)
  expect_s3_class(p_gavs, "ggplot")
})
