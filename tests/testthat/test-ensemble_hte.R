test_that("ensemble_hte creates proper object structure", {
  skip_on_cran()
  
  # Create simple simulated data with sufficient size
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  # Use T-learner to avoid zero-variance issues with small samples
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Test that object has correct class
  expect_s3_class(fit, "ensemble_hte_fit")
  
  # Test that required components exist
  expect_true(!is.null(fit$formula))
  expect_true(!is.null(fit$D))
  expect_true(!is.null(fit$ite))
  expect_true(!is.null(fit$M))
  expect_true(!is.null(fit$K))
})


test_that("ensemble_hte validates input correctly", {
  set.seed(123)
  n <- 200
  data <- data.frame(
    Y = rnorm(n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  
  # Should error with factor Y
  data_factor <- data
  data_factor$Y <- factor(data_factor$Y > 0)
  expect_error(
    ensemble_hte(Y ~ X1, treatment = D, data = data_factor, 
                 M = 2, K = 2, algorithms = c("lm")),
    "factor"
  )
  
  # Should error with missing treatment
  expect_error(
    ensemble_hte(Y ~ X1, data = data, M = 2, K = 2, algorithms = c("lm"))
  )
})


test_that("ensemble_hte handles different metalearners", {
  skip_on_cran()
  
  set.seed(123)
  n <- 250  # Larger sample for stability
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  tau <- 1 + 0.5 * X1  # True heterogeneous effect
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  # T-learner
  fit_t <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm"), metalearner = "t"
  )
  expect_s3_class(fit_t, "ensemble_hte_fit")
  expect_equal(fit_t$metalearner, "t")
  expect_true(var(fit_t$ite[[1]]) > 0)  # Non-zero variance in ITEs
  
  # S-learner
  fit_s <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm"), metalearner = "s"
  )
  expect_s3_class(fit_s, "ensemble_hte_fit")
  expect_equal(fit_s$metalearner, "s")
  expect_true(var(fit_s$ite[[1]]) > 0)  # Non-zero variance in ITEs
  
  # X-learner
  fit_x <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm"), metalearner = "x"
  )
  expect_s3_class(fit_x, "ensemble_hte_fit")
  expect_equal(fit_x$metalearner, "x")
  expect_true(var(fit_x$ite[[1]]) > 0)  # Non-zero variance in ITEs
  
  # R-learner with lm (not grf, for smaller sample compatibility)
  fit_r <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm"), metalearner = "r", r_learner = "lm"
  )
  expect_s3_class(fit_r, "ensemble_hte_fit")
  expect_equal(fit_r$metalearner, "r")
  expect_true(var(fit_r$ite[[1]]) > 0)  # Non-zero variance in ITEs
})


test_that("ensemble_hte R-learner with grf works with sufficient data", {
  skip_on_cran()
  
  set.seed(456)
  n <- 300  # Larger sample needed for grf
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  tau <- 1 + 0.5 * X1
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  # R-learner with grf (default)
  fit_r_grf <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm"), metalearner = "r"
  )
  expect_s3_class(fit_r_grf, "ensemble_hte_fit")
  expect_equal(fit_r_grf$r_learner, "grf")
  expect_true(var(fit_r_grf$ite[[1]]) > 0)
})

# Test ensemble_hte with grf algorithm ----
test_that("ensemble_hte works with grf algorithm", {
  skip_on_cran()
  set.seed(42)
  n <- 300  # grf needs larger samples
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  tau <- 1 + 0.5 * X1
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  # T-learner with grf
  fit_grf <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("grf"), metalearner = "t"
  )
  expect_s3_class(fit_grf, "ensemble_hte_fit")
  expect_true("grf" %in% fit_grf$algorithms)
  expect_true(var(fit_grf$ite[[1]]) > 0)
})

# Test ensemble_hte with multiple algorithms ----
test_that("ensemble_hte works with multiple algorithms", {
  skip_on_cran()
  set.seed(42)
  n <- 300
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  tau <- 1 + 0.5 * X1
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  # T-learner with both lm and grf
  fit_multi <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm", "grf"), metalearner = "t"
  )
  expect_s3_class(fit_multi, "ensemble_hte_fit")
  expect_true(all(c("lm", "grf") %in% fit_multi$algorithms))
  expect_true(var(fit_multi$ite[[1]]) > 0)
  
  # X-learner with both lm and grf
  fit_multi_x <- ensemble_hte(
    Y ~ X1 + X2, treatment = D, data = data,
    M = 2, K = 3, algorithms = c("lm", "grf"), metalearner = "x"
  )
  expect_s3_class(fit_multi_x, "ensemble_hte_fit")
  expect_true(var(fit_multi_x$ite[[1]]) > 0)
})

test_that("ensemble_hte respects scale_covariates = FALSE", {
  skip_on_cran()
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n, mean = 100, sd = 50)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    scale_covariates = FALSE,
    metalearner = "t"
  )
  
  expect_s3_class(fit, "ensemble_hte_fit")
  expect_false(fit$scale_covariates)
  # Check that stored X is unscaled (mean approx 100)
  expect_equal(mean(fit$X$X1), mean(X1), tolerance = 0.1)
})

test_that("ensemble_hte uses provided prop_score", {
  skip_on_cran()
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  data <- data.frame(Y = Y, D = D, X1 = X1)
  
  # Custom propensity scores
  ps <- runif(n, 0.4, 0.6)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    prop_score = ps,
    metalearner = "t"
  )
  
  expect_equal(fit$prop_score, ps)
})

# ---- store_baseline tests ----

# Shared fixture for store_baseline tests
make_baseline_data <- function(n = 250, seed = 42) {
  set.seed(seed)
  X1 <- rnorm(n); X2 <- rnorm(n)
  D  <- rbinom(n, 1, 0.5)
  Y  <- 0.5 * X1 + D * (1 + 0.5 * X1) + rnorm(n, sd = 0.5)
  list(data = data.frame(Y = Y, D = D, X1 = X1, X2 = X2), n = n)
}

test_that("store_baseline = 'none' gives NULL baseline", {
  skip_on_cran()
  d <- make_baseline_data()
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                      M = 2, K = 2, algorithms = c("lm"), metalearner = "t",
                      store_baseline = "none")
  expect_null(fit$baseline)
  expect_equal(fit$store_baseline, "none")
})

test_that("store_baseline = 'ensemble' (default) gives data.table n x M", {
  skip_on_cran()
  d <- make_baseline_data()
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                      M = 2, K = 2, algorithms = c("lm"), metalearner = "t")
  expect_equal(fit$store_baseline, "ensemble")
  expect_s3_class(fit$baseline, "data.table")
  expect_equal(dim(fit$baseline), c(d$n, 2L))  # n x M
  expect_false(anyNA(fit$baseline))
  expect_equal(names(fit$baseline), c("rep_1", "rep_2"))
})

test_that("store_baseline = 'all' gives 3D array n x A x M", {
  skip_on_cran()
  d  <- make_baseline_data()
  A  <- 2L  # two algorithms
  M  <- 2L
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                      M = M, K = 2, algorithms = c("lm", "grf"), metalearner = "t",
                      store_baseline = "all")
  expect_true(is.array(fit$baseline))
  expect_equal(dim(fit$baseline), c(d$n, A, M))
  expect_equal(dimnames(fit$baseline)[[2]], c("baseline_lm", "baseline_grf"))
  expect_equal(dimnames(fit$baseline)[[3]], c("rep_1", "rep_2"))
  expect_false(anyNA(fit$baseline))
  # Confirm no old $baseline_all field
  expect_null(fit$baseline_all)
})

test_that("is.array() distinguishes 'all' from 'ensemble' baseline", {
  skip_on_cran()
  d   <- make_baseline_data()
  f_ens <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                        M = 2, K = 2, algorithms = c("lm"), metalearner = "t",
                        store_baseline = "ensemble")
  f_all <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                        M = 2, K = 2, algorithms = c("lm"), metalearner = "t",
                        store_baseline = "all")
  expect_false(is.array(f_ens$baseline))
  expect_true(is.array(f_all$baseline))
})

test_that("ensemble_strategy = 'average' gives rowMeans baseline", {
  skip_on_cran()
  d <- make_baseline_data()
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                      M = 2, K = 2, algorithms = c("lm", "grf"), metalearner = "t",
                      ensemble_strategy = "average", store_baseline = "ensemble")
  expect_s3_class(fit$baseline, "data.table")
  expect_equal(dim(fit$baseline), c(d$n, 2L))
  expect_false(anyNA(fit$baseline))
})

test_that("ensemble_strategy = 'cv' baseline has no NAs and correct shape", {
  skip_on_cran()
  d <- make_baseline_data()
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                      M = 2, K = 2, algorithms = c("lm", "grf"), metalearner = "t",
                      ensemble_strategy = "cv", store_baseline = "ensemble")
  expect_s3_class(fit$baseline, "data.table")
  expect_equal(dim(fit$baseline), c(d$n, 2L))
  expect_false(anyNA(fit$baseline))
})

test_that("R-learner store_baseline uses all training obs (not control-only)", {
  skip_on_cran()
  # Use a large enough sample so the R-learner is stable
  set.seed(7)
  n <- 300
  X1 <- rnorm(n); X2 <- rnorm(n)
  D  <- rbinom(n, 1, 0.5)
  Y  <- X1 + D * (1 + X1) + rnorm(n)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = data,
                      M = 2, K = 2, algorithms = c("lm"), metalearner = "r",
                      r_learner = "lm", store_baseline = "ensemble")
  # Predictions should cover treated units too — non-NA for all rows
  expect_false(anyNA(fit$baseline))
  expect_equal(nrow(fit$baseline), n)
})

test_that("t-learner store_baseline fits on control units only", {
  skip_on_cran()
  d <- make_baseline_data()
  fit <- ensemble_hte(Y ~ X1 + X2, treatment = D, data = d$data,
                      M = 2, K = 2, algorithms = c("lm"), metalearner = "t",
                      store_baseline = "ensemble")
  # Baseline predictions should exist for all obs (including treated),
  # even though the ensemble was fitted on control units only
  expect_equal(nrow(fit$baseline), d$n)
  expect_false(anyNA(fit$baseline))
})
