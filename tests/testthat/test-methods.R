test_that("summary.ensemble_hte_fit works correctly", {
  skip_on_cran()
  
  set.seed(123)
  n <- 200
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  D <- rbinom(n, 1, 0.5)
  Y <- X1 + D * (1 + 0.5 * X1) + rnorm(n, sd = 0.5)
  data <- data.frame(Y = Y, D = D, X1 = X1, X2 = X2)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"  # Use T-learner for reliable ITEs
  )
  
  # Summary should produce output (it prints, doesn't return structured object)
  expect_output(summary(fit), "BLP|Ensemble HTE")
})


test_that("print.ensemble_hte_fit works correctly", {
  skip_on_cran()
  
  set.seed(123)
  n <- 200
  data <- data.frame(
    Y = rnorm(n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n)
  )
  data$Y <- data$X1 + data$D * 1 + rnorm(n, sd = 0.5)
  
  # Use T-learner to avoid zero-variance warnings
  fit <- ensemble_hte(
    Y ~ X1,
    treatment = D,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  expect_output(print(fit), "Ensemble HTE")
})


test_that("summary.ensemble_pred_fit provides complete information", {
  skip_on_cran()
  
  set.seed(123)
  n <- 100
  data <- data.frame(
    Y = rnorm(n),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  data$Y <- 0.5 * data$X1 + 0.3 * data$X2 + rnorm(n, sd = 0.5)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  summ <- summary(fit)
  
  expect_s3_class(summ, "summary.ensemble_pred_fit")
  expect_true(!is.null(summ$blp))
  expect_true(!is.null(summ$gavs))
  expect_true(!is.null(summ$metrics))  # It's 'metrics' not 'accuracy'
})
