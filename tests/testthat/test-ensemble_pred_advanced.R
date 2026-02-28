# ==============================================================================
# Tests for ensemble_pred advanced features
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
# Tests for ensemble_strategy parameter
# ==============================================================================

test_that("ensemble_pred works with ensemble_strategy = 'average'", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm", "ranger"),
    ensemble_strategy = "average"
  )
  
  expect_s3_class(fit, "ensemble_pred_fit")
  expect_equal(fit$ensemble_strategy, "average")
  
  # Print should show strategy
  expect_output(print(fit), "simple average")
})


test_that("ensemble_pred default ensemble_strategy is 'cv'", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  expect_equal(fit$ensemble_strategy, "cv")
  expect_output(print(fit), "cross-validated OLS")
})


test_that("ensemble_pred average strategy works with single algorithm", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("ranger"),
    ensemble_strategy = "average"
  )
  
  expect_s3_class(fit, "ensemble_pred_fit")
  expect_equal(fit$ensemble_strategy, "average")
})


# ==============================================================================
# Tests for matrix interface (Y, X)
# ==============================================================================

test_that("ensemble_pred works with matrix interface", {
  skip_on_cran()
  
  set.seed(123)
  n <- 100
  Y <- rnorm(n)
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Y <- 0.5 * X$X1 + 0.3 * X$X2 + rnorm(n, sd = 0.5)
  
  fit <- ensemble_pred(
    Y = Y,
    X = X,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  expect_s3_class(fit, "ensemble_pred_fit")
  expect_equal(length(fit$Y), n)
  expect_equal(nrow(fit$X), n)
})


test_that("ensemble_pred matrix interface creates proper formula", {
  skip_on_cran()
  
  set.seed(123)
  n <- 100
  Y <- rnorm(n)
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  
  fit <- ensemble_pred(
    Y = Y,
    X = X,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  expect_true(!is.null(fit$formula))
  expect_true(inherits(fit$formula, "formula"))
})


# ==============================================================================
# Tests for multiple algorithms
# ==============================================================================

test_that("ensemble_pred works with multiple algorithms", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm", "ranger", "grf")
  )
  
  expect_s3_class(fit, "ensemble_pred_fit")
  expect_equal(fit$algorithms, c("lm", "ranger", "grf"))
})


# ==============================================================================
# Tests for binary outcomes (classification)
# ==============================================================================

test_that("ensemble_pred auto-detects binary outcome", {
  skip_on_cran()
  
  set.seed(123)
  n <- 100
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Y <- rbinom(n, 1, plogis(0.5 * X$X1 + 0.3 * X$X2))
  
  fit <- ensemble_pred(
    Y = Y,
    X = X,
    M = 2,
    K = 2,
    algorithms = c("ranger")
  )
  
  expect_equal(fit$task_type, "classif")
})


test_that("ensemble_pred works with explicit task_type = 'classif'", {
  skip_on_cran()
  
  set.seed(123)
  n <- 100
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  Y <- rbinom(n, 1, 0.5)
  
  fit <- ensemble_pred(
    Y = Y,
    X = X,
    M = 2,
    K = 2,
    algorithms = c("ranger"),
    task_type = "classif"
  )
  
  expect_s3_class(fit, "ensemble_pred_fit")
  expect_equal(fit$task_type, "classif")
})


# ==============================================================================
# Tests for input validation
# ==============================================================================

test_that("ensemble_pred rejects factor Y", {
  data <- create_test_data(n = 100)
  data$Y <- factor(data$Y > 0)
  
  expect_error(
    ensemble_pred(Y ~ X1 + X2, data = data, M = 2, K = 2, algorithms = c("lm")),
    "factor"
  )
})


test_that("ensemble_pred validates M and K", {
  data <- create_test_data(n = 100)
  
  # K larger than n should error (or produce warning)
  expect_error(
    ensemble_pred(Y ~ X1 + X2, data = data, M = 2, K = 200, algorithms = c("lm"))
  )
})


test_that("ensemble_pred validates train_idx length", {
  data <- create_test_data(n = 100)
  
  expect_error(
    ensemble_pred(
      Y ~ X1 + X2,
      data = data,
      train_idx = c(TRUE, FALSE),  # Wrong length
      M = 2,
      K = 2,
      algorithms = c("lm")
    ),
    "train_idx"
  )
})


# ==============================================================================
# Tests for covariate scaling
# ==============================================================================

test_that("ensemble_pred respects scale_covariates = FALSE", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    scale_covariates = FALSE
  )
  
  expect_s3_class(fit, "ensemble_pred_fit")
  expect_false(fit$scale_covariates)
})


test_that("ensemble_pred stores original unscaled data", {
  skip_on_cran()
  
  set.seed(123)
  n <- 100
  X <- data.frame(X1 = rnorm(n, mean = 100, sd = 50), X2 = rnorm(n))
  Y <- 0.01 * X$X1 + X$X2 + rnorm(n)
  
  fit <- ensemble_pred(
    Y = Y,
    X = X,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    scale_covariates = TRUE
  )
  
  # Original X should be stored (unscaled)
  expect_equal(mean(fit$X$X1), mean(X$X1), tolerance = 0.01)
})
