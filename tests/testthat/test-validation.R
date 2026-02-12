
library(ensembleHTE)
library(testthat)

test_that("Input validation works for ensemble_hte", {
  # Create dummy data
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  data <- data.frame(Y = Y, D = D, X)
  
  # Test M and K validation
  expect_error(ensemble_hte(Y = Y, X = X, D = D, M = 0), "M must be a positive integer")
  expect_error(ensemble_hte(Y = Y, X = X, D = D, K = 1), "K must be an integer >= 2")
  
  # Test algorithms validation
  expect_error(ensemble_hte(Y = Y, X = X, D = D, algorithms = character(0)), "algorithms must be a non-empty character vector")
  
  # Test prop_score validation
  expect_error(ensemble_hte(Y = Y, X = X, D = D, prop_score = rep(1.5, n)), "prop_score must be between 0 and 1")
  
  # Test ensemble_folds validation
  expect_error(ensemble_hte(Y = Y, X = X, D = D, ensemble_folds = 1), "ensemble_folds must be an integer >= 2")
  
  # Test n_cores validation
  expect_error(ensemble_hte(Y = Y, X = X, D = D, n_cores = 0), "n_cores must be a positive integer")
})

test_that("Treatment variable can be specified in multiple ways", {
  # Create dummy data
  set.seed(123)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- 1 + X1 + D * X2 + rnorm(n)
  data <- data.frame(Y = Y, treat = D, X1 = X1, X2 = X2)
  
  # Test unquoted treatment name
  fit1 <- ensemble_hte(Y ~ X1 + X2, treatment = treat, data = data, 
                       M = 2, K = 2, algorithms = "lm")
  expect_equal(fit1$treatment, "treat")
  
  # Test quoted treatment name
  fit2 <- ensemble_hte(Y ~ X1 + X2, treatment = "treat", data = data, 
                       M = 2, K = 2, algorithms = "lm")
  expect_equal(fit2$treatment, "treat")
  
  # Test treatment name from variable
  treat_col <- "treat"
  fit3 <- ensemble_hte(Y ~ X1 + X2, treatment = treat_col, data = data, 
                       M = 2, K = 2, algorithms = "lm")
  expect_equal(fit3$treatment, "treat")
  
  # Test error when treatment column doesn't exist
  expect_error(ensemble_hte(Y ~ X1 + X2, treatment = "nonexistent", data = data),
               "treatment 'nonexistent' not found in data")
  
  expect_error(ensemble_hte(Y ~ X1 + X2, treatment = nonexistent, data = data),
               "not found in data")
})

test_that("Input validation works for ensemble_pred", {
  # Create dummy data
  n <- 100
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)
  data <- data.frame(Y = Y, X)
  
  # Test M and K validation
  expect_error(ensemble_pred(Y = Y, X = X, M = 0), "M must be a positive integer")
  expect_error(ensemble_pred(Y = Y, X = X, K = 1), "K must be an integer >= 2")
  
  # Test algorithms validation
  expect_error(ensemble_pred(Y = Y, X = X, algorithms = character(0)), "algorithms must be a non-empty character vector")
  
  # Test ensemble_folds validation
  expect_error(ensemble_pred(Y = Y, X = X, ensemble_folds = 1), "ensemble_folds must be an integer >= 2")
  
  # Test n_cores validation
  expect_error(ensemble_pred(Y = Y, X = X, n_cores = 0), "n_cores must be a positive integer")
  
  # Test train_idx validation
  expect_error(ensemble_pred(Y = Y, X = X, train_idx = c(TRUE, FALSE)), "train_idx has length")
  expect_error(ensemble_pred(Y = Y, X = X, train_idx = c(0, n+1)), "train_idx contains indices outside")
  expect_error(ensemble_pred(Y = Y, X = X, train_idx = "invalid"), "train_idx must be NULL")
})
