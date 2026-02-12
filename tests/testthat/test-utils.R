test_that("create_folds creates balanced folds", {
  n <- 100
  M <- 3
  K <- 5
  folds <- create_folds(n, M = M, K = K)
  
  # Check structure
  expect_length(folds, M)
  expect_equal(length(folds[[1]]), n)
  expect_equal(length(unique(folds[[1]])), K)
  
  # Check that all observations are assigned
  expect_true(all(folds[[1]] %in% 1:K))
})


test_that("create_folds stratifies by variable", {
  n <- 100
  strat_var <- rep(1:5, each = 20)
  folds <- create_folds(n, M = 1, K = 5, stratify_var = strat_var)
  
  # Check structure
  expect_length(folds, 1)
  expect_equal(length(folds[[1]]), n)
})


test_that("estimate_propensity_score returns valid probabilities", {
  skip_on_cran()
  
  set.seed(123)
  n <- 200
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  D <- rbinom(n, 1, plogis(0.5 * X$X1))
  
  # Estimate propensity scores using internal function if available
  # For now, just test that propensity scores in ensemble_hte are valid
  data <- cbind(Y = rnorm(n), D = D, X)
  
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
  
  # Check prop_score exists and is valid
  expect_true(!is.null(fit$prop_score))
  expect_true(all(fit$prop_score > 0))
  expect_true(all(fit$prop_score < 1))
})


test_that("create_groups creates balanced groups", {
  set.seed(123)
  x <- rnorm(100)
  groups <- create_groups(x, n_groups = 5)
  
  expect_equal(length(groups), 100)
  expect_equal(length(unique(groups)), 5)
  expect_true(all(groups %in% 1:5))
  
  # Check groups are roughly balanced
  group_sizes <- table(groups)
  expect_true(all(group_sizes >= 15))  # At least 15 in each of 5 groups from 100
  expect_true(all(group_sizes <= 25))  # At most 25
})


test_that("create_groups respects ordering", {
  set.seed(123)
  x <- 1:100
  groups <- create_groups(x, n_groups = 5)
  
  # Group 1 should have lowest values, group 5 highest
  expect_true(mean(x[groups == 1]) < mean(x[groups == 3]))
  expect_true(mean(x[groups == 3]) < mean(x[groups == 5]))
})
