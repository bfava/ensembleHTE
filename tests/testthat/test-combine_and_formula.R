# ==============================================================================
# Tests for combine_ensembles() and reg_from_formula()
# ==============================================================================

# Shared test data creation helper
create_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)
  D <- rbinom(n, 1, 0.5)
  tau <- 1 + 0.5 * X1
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  
  data.frame(Y = Y, D = D, X1 = X1, X2 = X2, X3 = X3)
}

# ==============================================================================
# Tests for combine_ensembles()
# ==============================================================================

test_that("combine_ensembles works with ensemble_pred_fit objects", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  # Create two fits with same formula/data but different seeds
  fit1 <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  fit2 <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Combine
  combined <- combine_ensembles(fit1, fit2)
  
  # Test class preserved

  expect_s3_class(combined, "ensemble_pred_fit")
  
  # Test M is summed
  expect_equal(combined$M, fit1$M + fit2$M)
  
  # Test predictions combined
  expect_equal(ncol(combined$predictions), 4)  # 2 + 2
})


test_that("combine_ensembles works with ensemble_hte_fit objects", {
  skip_on_cran()
  
  data <- create_test_data(n = 200)
  
  # Use T-learner to avoid zero-variance warnings
  fit1 <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  fit2 <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  combined <- combine_ensembles(fit1, fit2)
  
  expect_s3_class(combined, "ensemble_hte_fit")
  expect_equal(combined$M, 4)
  expect_equal(length(combined$ite), 4)
})


test_that("combine_ensembles requires compatible fits", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  # Different formulas
  fit1 <- ensemble_pred(
    Y ~ X1,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  fit2 <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  expect_error(combine_ensembles(fit1, fit2))
})


test_that("combine_ensembles requires same class", {
  skip_on_cran()
  
  data <- create_test_data(n = 200)
  
  fit_pred <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Use T-learner to avoid zero-variance warnings
  fit_hte <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  expect_error(combine_ensembles(fit_pred, fit_hte))
})


test_that("combine_ensembles works with more than 2 fits", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit1 <- ensemble_pred(Y ~ X1 + X2, data = data, M = 2, K = 2, algorithms = c("lm"))
  fit2 <- ensemble_pred(Y ~ X1 + X2, data = data, M = 2, K = 2, algorithms = c("lm"))
  fit3 <- ensemble_pred(Y ~ X1 + X2, data = data, M = 2, K = 2, algorithms = c("lm"))
  
  combined <- combine_ensembles(fit1, fit2, fit3)
  
  expect_s3_class(combined, "ensemble_pred_fit")
  expect_equal(combined$M, 6)  # 2 + 2 + 2
})


# ==============================================================================
# Tests for reg_from_formula()
# ==============================================================================

test_that("reg_from_formula fits a regression model", {
  data <- data.frame(Y = rnorm(50), X1 = rnorm(50), X2 = rnorm(50))
  data$Y <- 0.5 * data$X1 + 0.3 * data$X2 + rnorm(50, sd = 0.5)
  
  result <- reg_from_formula(Y ~ X1 + X2, data = data)
  
  # Should return list with model components
  expect_true(is.list(result))
  expect_true(!is.null(result$coef))
  expect_true(!is.null(result$vcov))
  expect_true(!is.null(result$residuals))
  expect_true(!is.null(result$fitted))
})


test_that("reg_from_formula works with weights", {
  data <- data.frame(Y = rnorm(50), X1 = rnorm(50))
  weights <- runif(50, 0.5, 1.5)
  
  result <- reg_from_formula(Y ~ X1, data = data, weights = weights)
  
  expect_true(is.list(result))
  expect_true(!is.null(result$coef))
})


test_that("reg_from_formula handles character formula", {
  data <- data.frame(Y = rnorm(50), X1 = rnorm(50))
  
  result <- reg_from_formula("Y ~ X1", data = data)
  
  expect_true(is.list(result))
  expect_true(!is.null(result$coef))
})
